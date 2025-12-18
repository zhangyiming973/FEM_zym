function b = assemble_vector_2D(vector_size, gauss_type, function_name, ...
                                P_mesh, T_mesh, T_test, ...
                                test_basis_type, test_der_x, test_der_y)
% ASSEMBLE_VECTOR_2D 组装有限元载荷向量 b
%
% 参数:
%   vector_size: 向量总长度 (通常是 N_test_dofs)
%   function_name: 载荷函数 f 的句柄
%   P_mesh, T_mesh: 几何网格 (用于定义积分区域 En)
%   T_test: 测试空间矩阵 (用于获取全局索引 global_i)
%   test_der_x, test_der_y: 导数阶数 (p, q)

% 1. 初始化向量 Initialize b = sparse(...)
b = sparse(vector_size, 1);

% 2. 准备高斯点
[gauss_weight, gauss_xi, gauss_eta] = gaussValues_2d(gauss_type);

% 获取基本信息
num_local_test = size(T_test, 1); % N_lb
number_of_elements = size(T_mesh, 2); % N

% 3. 循环遍历所有单元 FOR n = 1, ..., N
for n = 1:number_of_elements
    
    % 获取当前单元的几何顶点 (vertices)
    vertices = P_mesh(:, T_mesh(:, n)); 
    
    % 4. 循环遍历测试基函数 FOR beta = 1, ..., N_lb
    for beta = 1:num_local_test
        
        % Compute r = integral(...)
        r = Gauss_quad_vector_2D(gauss_xi, gauss_eta, gauss_weight, vertices, function_name, ...
                                 test_basis_type, beta, test_der_x, test_der_y);
        
        % 获取全局自由度索引 T_b(beta, n)
        global_idx = T_test(beta, n);
        
        % b(global_idx) = b(global_idx) + r
        b(global_idx) = b(global_idx) + r;
        
    end
end

end