function b = assemble_neumann_vector(b, P_mesh, T_mesh, P_test, T_test, ...
                                     test_basis_type, bc_flags, neumann_funcs)
% ASSEMBLE_NEUMANN_VECTOR 处理非零 Neumann 边界条件 (面外力)
% 原理: 计算边界线积分 int(g * psi) ds 并加到载荷向量 b 中
% 输入:
% b: 当前的载荷向量 (通常包含了体力的贡献)
% P_mesh, T_mesh: 几何网格 (用于判断边界位置)
% P_test, T_test: 测试空间节点和拓扑 (用于获取自由度)
% bc_flags: 1x4 向量 [下, 右, 上, 左], 标记边界类型 (-2)
% neumann_funcs: 结构体, 包含各边界的函数句柄
% .bottom, .right, .top, .left 
%
% 输出 b: 加上了面力贡献后的载荷向量
%∫_Γ g(x,y) · ψ(x,y) ds
% 准备一维高斯积分点 (用于线积分) 
[gauss_w, gauss_pt] = gaussValues_1d(3); 

num_elements = size(T_mesh, 2);
num_local_test = size(T_test, 1);

% 几何容差
tol = 1e-6;
x_min = min(P_mesh(1,:)); x_max = max(P_mesh(1,:));
y_min = min(P_mesh(2,:)); y_max = max(P_mesh(2,:));

% --- 2. 遍历所有单元 ---
for n = 1:num_elements
    % 获取几何顶点索引 (三角形的3个顶点)
    geo_idx = T_mesh(:, n);
    % 获取几何顶点坐标
    P1 = P_mesh(:, geo_idx(1));
    P2 = P_mesh(:, geo_idx(2));
    P3 = P_mesh(:, geo_idx(3));
    
    % 定义三角形的三条边 (顶点对)
    % Edge 1: P1 -> P2
    % Edge 2: P2 -> P3
    % Edge 3: P3 -> P1
    edges =[1, 2;
            2, 3; 
            3, 1];
    
    % 获取当前单元的所有测试基函数自由度序号
    dofs = T_test(:, n);
    
    %  遍历该单元的 3 条边
    for e = 1:3
        idx_a = edges(e, 1);
        idx_b = edges(e, 2);
        
        Pa = P_mesh(:, geo_idx(idx_a));
        Pb = P_mesh(:, geo_idx(idx_b));
        
        % 判断该边是否在 Neumann 边界上
        is_neumann = false;
        g_func = [];
        
        % 检查是否都在下边界 (y=y_min) 且标记为 Neumann (-2)
        if abs(Pa(2)-y_min)<tol && abs(Pb(2)-y_min)<tol && bc_flags(1)==-2
            is_neumann = true; g_func = neumann_funcs.bottom;
        % 检查是否都在右边界
        elseif abs(Pa(1)-x_max)<tol && abs(Pb(1)-x_max)<tol && bc_flags(2)==-2
            is_neumann = true; g_func = neumann_funcs.right;
        % 检查是否都在上边界
        elseif abs(Pa(2)-y_max)<tol && abs(Pb(2)-y_max)<tol && bc_flags(3)==-2
            is_neumann = true; g_func = neumann_funcs.top;
        % 检查是否都在左边界
        elseif abs(Pa(1)-x_min)<tol && abs(Pb(1)-x_min)<tol && bc_flags(4)==-2
            is_neumann = true; g_func = neumann_funcs.left;
        end
        
        % 如果是 Neumann 边，进行线积分
        if is_neumann
            % 边长 (Jacobian for line integral)
            edge_len = norm(Pa - Pb);
            jac = edge_len / 2; % 从 [-1,1] 映射到 [0, L] 的系数
            
            % 遍历高斯点
            for k = 1:length(gauss_w)
                xi = gauss_pt(k);
                w = gauss_w(k);
                
                % 物理坐标映射: Pa 和 Pb 之间的线性插值
                % P(xi) = Pa * (1-xi)/2 + Pb * (1+xi)/2
                x_phys = Pa(1)*(1-xi)/2 + Pb(1)*(1+xi)/2;
                y_phys = Pa(2)*(1-xi)/2 + Pb(2)*(1+xi)/2;
                
                % 计算 Neumann 边界函数值 g(x,y)
                g_val = g_func(x_phys, y_phys);
                
                % 计算该位置的所有测试基函数值
                % 直接调用 2D 基函数，输入物理坐标
                % 几何 vertices 必须是整个单元的顶点 (P1, P2, P3)
                vertices = [P1, P2, P3];
                
                for beta = 1:num_local_test
                    % 计算基函数值 (导数阶数为 0)
                    psi_val = FE_basis_local_fun_2D(x_phys, y_phys, vertices, ...
                                                    test_basis_type, beta, 0, 0);
                    
                    % 累加积分: result += w * g * psi * jac
                    a = w * g_val * psi_val * jac;
                    
                    % 加到全局向量
                    global_idx = dofs(beta);
                    b(global_idx) = b(global_idx) + a;
                end
            end
        end
    end
end

end