function [u_final, P_fem, T_fem] = FE_solver_2D_hyperbolic(x_range, y_range, N_x, N_y, ...
                                                           trial_type, test_type, gauss_type, ...
                                                           bc_flags, ...
                                                           f_func_t, c_func, ...
                                                           bc_funcs_t, ...
                                                           u0_func, v0_func, ... % 初始位移和初始速度
                                                           time_config, a_type)
% FE_SOLVER_2D_HYPERBOLIC 二维双曲型方程求解器 (Wave Equation)
% 方程: u_tt - div(c * grad(u)) = f
% 算法: Slide 79-80 (Centered Difference / Implicit Scheme)

    % --- 1. 网格与矩阵组装 ---
    [P_mesh, T_mesh, P_fem, T_fem] = generate_PbTb(x_range, y_range, N_x, N_y, trial_type);
    [~, ~, P_test, T_test] = generate_PbTb(x_range, y_range, N_x, N_y, test_type);
    
    N_dofs = size(P_fem, 2);
    matrix_size = [N_dofs, N_dofs];
    dt = time_config.dt;
    
    % 组装 M (质量矩阵)
    func_one = @(x,y) 1.0 + 0*x + 0*y;
    M = assemble_matrix_2D(matrix_size, gauss_type, func_one, ...
                           P_mesh, T_mesh, T_fem, T_test, trial_type, 0, 0, test_type, 0, 0);
                       
    % 组装 A (刚度矩阵) - 假设 c 不随时间变化
    if a_type == 1
        A_xx = assemble_matrix_2D(matrix_size, gauss_type, c_func, P_mesh, T_mesh, T_fem, T_test, trial_type, 1, 0, test_type, 1, 0);
        A_yy = assemble_matrix_2D(matrix_size, gauss_type, c_func, P_mesh, T_mesh, T_fem, T_test, trial_type, 0, 1, test_type, 0, 1);
        A = A_xx + A_yy;
    else
        error('Anisotropic not implemented in this snippet for brevity');
    end

    % --- 2. 初始条件处理 (t=0 和 t=dt) ---
    % Slide 80: Generate X0 and X1
    
    % 2.1 计算 X0 (t=0)
    X_m_prev = u0_func(P_fem(1,:), P_fem(2,:))'; % X^{m-1} (对应 m=1 时的上一项)
    
    % 2.2 计算 X1 (t=dt) - 使用二阶泰勒展开启动
    % 需要初始速度 V0
    V0 = v0_func(P_fem(1,:), P_fem(2,:))';
    
    % 计算初始载荷 b0
    f_t0 = @(x,y) f_func_t(x,y, 0);
    b_0 = assemble_vector_2D(N_dofs, gauss_type, f_t0, P_mesh, T_mesh, T_test, test_type, 0, 0);
    
    % 处理 b0 的边界条件 (为了正确计算初始加速度)
    boundary_nodes = boundarynodes_2D(P_fem, bc_flags);
    % 注意：这里简化处理，假设初始加速度计算时忽略 Dirchlet 强行约束的影响，或者认为 u0 满足方程
    % 更严谨的做法是求解 M * Acc = b0 - A * X0
    % 这里我们用显式近似:
    % Acc_0 = M \ (b_0 - A * X_m_prev); % 这可能比较慢，且 M 没做集中质量
    % 简单起见，如果 u0 是平滑的且 V0=0，通常用 X1 = X0 + dt*V0 即可 (一阶启动)
    % 为了二阶精度: X1 = X0 + dt*V0 + 0.5*dt^2 * (M^-1 * (b0 - A*X0))
    
    % 这里的实现采用显式预估 X1 (假设 M 是对角占优或使用 lumped mass 会更快，这里直接解)
    % 警告：M 是稀疏矩阵，直接求逆很慢。这里解线性方程 M * k = RHS
    Acc_0 = M \ (b_0 - A * X_m_prev);
    X_m = X_m_prev + dt * V0 + 0.5 * dt^2 * Acc_0; % X^m (对应 m=1 时的当前项 X^1)

    % 强制 X1 满足边界条件 (t=dt)
    bc_funcs_t1 = bc_funcs_t;
    bc_funcs_t1.dirichlet = @(x,y) bc_funcs_t.dirichlet(x,y, dt);
    [~, X_m] = BC_process_2D(speye(N_dofs), X_m, boundary_nodes, P_fem, bc_funcs_t1);

    % --- 3. 预计算迭代矩阵 (Slide 79) ---
    % LHS = M/dt^2 + A/4
    LHS_matrix = (M / (dt^2)) + (A / 4);
    
    % 系数矩阵 (用于 RHS 计算)
    Coeff_X_m      = (2*M)/(dt^2) - A/2;      % 对应 X^m 的系数
    Coeff_X_m_prev = -((M)/(dt^2) + A/4);     % 对应 X^{m-1} 的系数

    % --- 4. 时间迭代循环 ---
    current_time = dt; % 我们已经有了 t=0 和 t=dt，下一次求解 t=2*dt
    step = 1;
    
    figure; 
    
    while current_time < time_config.T_end - dt/2
        step = step + 1;
        next_time = current_time + dt; % t_{m+1}
        prev_time = current_time - dt; % t_{m} (对应公式里的 b(tm))
        
        % 4.1 组装 b(t_m) [注意: Slide 79 公式里 RHS 用的是 b(tm)]
        % 也就是上一步时刻的载荷
        f_at_tm = @(x,y) f_func_t(x,y, current_time);
        b_tm = assemble_vector_2D(N_dofs, gauss_type, f_at_tm, P_mesh, T_mesh, T_test, test_type, 0, 0);
        
        % 4.2 构建 RHS
        % RHS = b(tm) + [...]X^m + [...]X^{m-1}
        RHS_vec = b_tm + Coeff_X_m * X_m + Coeff_X_m_prev * X_m_prev;
        
        % 4.3 处理边界条件 (t_{m+1})
        boundary_nodes = boundarynodes_2D(P_fem, bc_flags);
        bc_funcs_now = bc_funcs_t;
        bc_funcs_now.dirichlet = @(x,y) bc_funcs_t.dirichlet(x,y, next_time);
        
        [LHS_mod, RHS_mod] = BC_process_2D(LHS_matrix, RHS_vec, boundary_nodes, P_fem, bc_funcs_now);
        
        % 4.4 求解 X^{m+1}
        X_next = LHS_mod \ RHS_mod;
        
        % 4.5 滚动更新
        X_m_prev = X_m; % 旧的变成更旧的
        X_m = X_next;   % 新的变成旧的
        current_time = next_time;
        
        % 4.6 绘图
        if mod(step, 2) == 0
            trisurf(T_fem(1:3,:)', P_fem(1,:), P_fem(2,:), full(X_m), 'EdgeColor', 'none', 'FaceColor', 'interp');
            title(['Wave Equation t = ', num2str(current_time, '%.3f')]);
            zlim([-1, 1]); caxis([-1, 1]); % 固定坐标轴以便观察波动
            axis equal; view(3); shg; drawnow;
        end
    end
    
    u_final = X_m;
end