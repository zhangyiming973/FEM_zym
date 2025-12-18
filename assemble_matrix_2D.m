function A = assemble_matrix_2D(matrix_size, gauss_type, function_name, ...
                                P_mesh, T_mesh, T_trial, T_test, ...
                                trial_basis_type, trial_der_x, trial_der_y, ...
                                test_basis_type, test_der_x, test_der_y)
% ASSEMBLE_MATRIX_2D 组装有限元矩阵 A
% 
% 算法来源: Algorithm I-3
%
% 参数:
%   matrix_size: [N_test_dofs, N_trial_dofs]
%   P_mesh, T_mesh: 几何网格 (定义积分区域 En)
%   T_trial, T_test: 有限元空间拓扑 (定义 alpha, beta 和全局映射)
%   *_der_x, *_der_y: 定义导数阶数 (r, s, p, q)

% 1. 初始化矩阵
% Initialize A = sparse(...)
A = sparse(matrix_size(1), matrix_size(2));

% 2. 准备高斯点
[gauss_weight, gauss_xi, gauss_eta] = gaussValues_2d(gauss_type);

% 获取单元内基函数个数 (N_lb_trial, N_lb_test)
num_local_trial = size(T_trial, 1); 
num_local_test  = size(T_test, 1);
number_of_elements = size(T_mesh, 2);

% 3. 循环遍历所有单元
% FOR n = 1, ..., N
for n = 1:number_of_elements
    
    % 获取当前单元的几何顶点 (用于计算雅可比和坐标映射)
    % P_mesh 存储了网格节点坐标
    vertices = P_mesh(:, T_mesh(:, n)); 
    
    % 4. 循环遍历试探基函数 (alpha)
    % FOR alpha = 1, ..., N_lb_trial
    for alpha = 1:num_local_trial
        
        % 5. 循环遍历测试基函数 (beta)
        % FOR beta = 1, ..., N_lb_test
        for beta = 1:num_local_test
            
            % Compute r = integral(...)
            r = Gauss_quad_matrix_2D(gauss_xi, gauss_eta, gauss_weight, vertices, function_name, ...
                                     trial_basis_type, alpha, trial_der_x, trial_der_y, ...
                                     test_basis_type, beta, test_der_x, test_der_y);
            
            % 获取全局自由度索引
            global_i = T_test(beta, n);  % Row index (Test)
            global_j = T_trial(alpha, n); % Col index (Trial)
            
            % Add r to A(...)
            A(global_i, global_j) = A(global_i, global_j) + r;
            
        end
    end
end

end