function result = Gauss_quad_vector_2D(xi_list, eta_list, gauss_weight, vertices, function_name, ...
                                       test_basis_type, test_beta, test_der_x, test_der_y)
% GAUSS_QUAD_VECTOR_2D 计算二维单元载荷向量积分
%
% 依据: Compute r = integral(f * D^p D^q psi)
%
% 输入:
%   xi_list, eta_list, gauss_weight: 高斯点和权重
%   vertices: 单元顶点坐标 (2x3)
%   function_name: 载荷函数 f(x,y) 的句柄或名称
%   test_basis_type: 测试基函数类型 (如 201, 202)
%   test_beta: 测试基函数在单元内的局部索引 (beta)
%   test_der_x, test_der_y: 导数阶数 p, q

% 1. 获取几何信息 & 计算雅可比
x1 = vertices(1,1); y1 = vertices(2,1);
x2 = vertices(1,2); y2 = vertices(2,2);
x3 = vertices(1,3); y3 = vertices(2,3);

% 计算雅可比行列式 |J| (用于积分变量变换 dxdy = |J| dxi deta)
det_J = (x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1);
abs_det_J = abs(det_J); 

% 2. 高斯积分循环
integral_sum = 0;
num_points = length(gauss_weight);

for k = 1:num_points
    xi = xi_list(k);
    eta = eta_list(k);
    w = gauss_weight(k);
    
    % --- 坐标映射: 参考坐标 (xi, eta) -> 物理坐标 (x, y) ---
    % 使用线性形函数进行几何映射 (Isoparametric mapping)
    lambda1 = 1 - xi - eta;
    lambda2 = xi;
    lambda3 = eta;
    
    x_phys = x1*lambda1 + x2*lambda2 + x3*lambda3;
    y_phys = y1*lambda1 + y2*lambda2 + y3*lambda3;
    
    % --- 计算载荷函数 f(x,y) ---
    f_val = feval(function_name, x_phys, y_phys);
    
    % --- 计算测试基函数值/导数 ---
    % 调用之前写好的 FE_basis_local_fun_2D (支持 x,y 输入)
    test_val = FE_basis_local_fun_2D(x_phys, y_phys, vertices, ...
                                     test_basis_type, test_beta, test_der_x, test_der_y);
    
    % --- 累加积分 ---
    % r += weight * f * test_psi * |J|
    integral_sum = integral_sum + w * f_val * test_val * abs_det_J;
end

result = integral_sum;

end