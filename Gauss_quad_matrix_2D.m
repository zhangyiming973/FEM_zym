function element_int = Gauss_quad_matrix_2D(xi_list, eta_list, gauss_weight, vertices, function_name, ...
                                          trial_basis_type, trial_alpha, trial_der_x, trial_der_y, ...
                                          test_basis_type, test_beta, test_der_x, test_der_y)
% 计算单元刚度矩阵元素的积分值
% 依据: r = integral(c * D_trial * D_test)

% 1. 获取几何信息用于映射
x1 = vertices(1,1); y1 = vertices(2,1);
x2 = vertices(1,2); y2 = vertices(2,2);
x3 = vertices(1,3); y3 = vertices(2,3);

% 计算雅可比行列式 |J| (用于积分变量变换 dxdy = |J| dxi deta)
% |J| = (x2-x1)(y3-y1) - (x3-x1)(y2-y1)
det_J = (x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1);
abs_det_J = abs(det_J); 

% 2. 遍历高斯点进行数值积分
integral_sum = 0;
num_points = length(gauss_weight);

for k = 1:num_points
    xi = xi_list(k);
    eta = eta_list(k);
    w = gauss_weight(k);
    
    % 坐标映射: 参考坐标 (xi, eta) -> 物理坐标 (x, y) 
    % 必须先算出物理坐标，因为 coefficient function c(x,y) 和 
    % 我们修改后的 FE_basis_local_fun_2D 都需要物理坐标作为输入。
    % 使用线性形函数映射几何 (Isoparametric mapping for linear mesh):
    lambda1 = 1 - xi - eta;
    lambda2 = xi;
    lambda3 = eta;
    
    x_phys = x1*lambda1 + x2*lambda2 + x3*lambda3;
    y_phys = y1*lambda1 + y2*lambda2 + y3*lambda3;
    
    % 计算系数函数 c(x,y)
    c_val = feval(function_name, x_phys, y_phys);
    
    %  计算 Trial 基函数值/导数 
    % 注意：传入的是 x_phys, y_phys，内部会处理逆变换
    trial_val = FE_basis_local_fun_2D(x_phys, y_phys, vertices, trial_basis_type, trial_alpha, trial_der_x, trial_der_y);
    
    % 计算 Test 基函数值/导数 
    test_val = FE_basis_local_fun_2D(x_phys, y_phys, vertices, test_basis_type, test_beta, test_der_x, test_der_y);
    
    % --- 累加积分 ---
    % result += weight * c * trial * test * |J|
    % 注意：gaussValues_2d 返回的权重通常是针对面积 0.5 的三角形归一化的，或者直接是面积。
    % 严谨写法: integral over triangle = sum( w_k * f(x_k) ) * Area_Ref * |J|/Area_Ratio?
    % 简单记忆: int f dxdy = int f |J| dxi deta. 
    % gaussValues_2d 中的权重如果是针对标准单纯形(面积0.5)的，则不需要额外乘0.5，只需乘 |J|。
    
    integral_sum = integral_sum + w * c_val * trial_val * test_val * abs_det_J;
end

element_int = integral_sum;

end