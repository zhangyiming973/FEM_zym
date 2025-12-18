function A = assemble_robin_matrix(A, P_mesh, T_mesh, P_trial, T_trial, P_test, T_test, ...
                                   trial_basis_type, test_basis_type, bc_flags, robin_q_funcs)
% ASSEMBLE_ROBIN_MATRIX 处理 Robin 边界条件的左端项 (刚度矩阵修正)
% 数学项: + int(q * u * v) ds 加到矩阵 A 中
% 输入:
%A: 当前的刚度矩阵
%bc_flags: 标记为 -3 的是 Robin 边界
%robin_q_funcs: 结构体, 包含各边界的 q(x,y) 系数函数 (.bottom, .right 等)

% --- 1. 准备一维高斯积分 ---
[gauss_w, gauss_pt] = gaussValues_1d(3); 

num_elements = size(T_mesh, 2);
num_local_trial = size(T_trial, 1);
num_local_test  = size(T_test, 1);

% 几何容差
tol = 1e-6;
x_min = min(P_mesh(1,:)); x_max = max(P_mesh(1,:));
y_min = min(P_mesh(2,:)); y_max = max(P_mesh(2,:));

% --- 2. 遍历所有单元 ---
for n = 1:num_elements
    geo_idx = T_mesh(:, n);
    P1 = P_mesh(:, geo_idx(1));
    P2 = P_mesh(:, geo_idx(2));
    P3 = P_mesh(:, geo_idx(3));
    
    % 定义三条边
    edges = [1, 2; 2, 3; 3, 1];
    
    % 获取自由度索引
    dofs_trial = T_trial(:, n);
    dofs_test  = T_test(:, n);
    
    % --- 3. 遍历单元的三条边 ---
    for e = 1:3
        idx_a = edges(e, 1);
        idx_b = edges(e, 2);
        Pa = P_mesh(:, geo_idx(idx_a));
        Pb = P_mesh(:, geo_idx(idx_b));
        
        % --- 4. 判断是否为 Robin 边界 (-3) ---
        is_robin = false;
        q_func = [];
        
        if abs(Pa(2)-y_min)<tol && abs(Pb(2)-y_min)<tol && bc_flags(1)==-3
            is_robin = true; q_func = robin_q_funcs.bottom;
        elseif abs(Pa(1)-x_max)<tol && abs(Pb(1)-x_max)<tol && bc_flags(2)==-3
            is_robin = true; q_func = robin_q_funcs.right;
        elseif abs(Pa(2)-y_max)<tol && abs(Pb(2)-y_max)<tol && bc_flags(3)==-3
            is_robin = true; q_func = robin_q_funcs.top;
        elseif abs(Pa(1)-x_min)<tol && abs(Pb(1)-x_min)<tol && bc_flags(4)==-3
            is_robin = true; q_func = robin_q_funcs.left;
        end
        
        % --- 5. 执行线积分并组装到 A ---
        if is_robin
            edge_len = norm(Pa - Pb);
            jac = edge_len / 2;
            
            for k = 1:length(gauss_w)
                xi = gauss_pt(k);
                w = gauss_w(k);
                
                % 物理坐标
                x_phys = Pa(1)*(1-xi)/2 + Pb(1)*(1+xi)/2;
                y_phys = Pa(2)*(1-xi)/2 + Pb(2)*(1+xi)/2;
                
                % 计算 q 系数
                q_val = q_func(x_phys, y_phys);
                
                % 几何顶点 (用于基函数计算)
                vertices = [P1, P2, P3];
                
                % 双重循环: Trial x Test
                for alpha = 1:num_local_trial
                    % 计算 Trial 基函数值
                    psi_trial = FE_basis_local_fun_2D(x_phys, y_phys, vertices, ...
                                                      trial_basis_type, alpha, 0, 0);
                    
                    for beta = 1:num_local_test
                        % 计算 Test 基函数值
                        psi_test = FE_basis_local_fun_2D(x_phys, y_phys, vertices, ...
                                                         test_basis_type, beta, 0, 0);
                        
                        % 积分: w * q * u * v * jac
                        val = w * q_val * psi_trial * psi_test * jac;
                        
                        % 组装到刚度矩阵
                        row = dofs_test(beta);
                        col = dofs_trial(alpha);
                        A(row, col) = A(row, col) + val;
                    end
                end
            end
        end
    end
end

end