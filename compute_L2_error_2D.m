function total_error = compute_L2_error_2D(solution, P_mesh, T_mesh, P_trial, T_trial, trial_basis_type, analytical_u_handle, der_x, der_y)
% 计算二维有限元解的 L2 误差或 H1 半范数误差

% 输入:
%solution:            有限元解向量 (系数向量)
%P_mesh, T_mesh:      几何网格信息 (用于定义积分区域 En)
%P_trial, T_trial:    有限元空间信息 (用于获取自由度索引和基函数)
%trial_basis_type:    基函数类型 (201, 202 等)
%analytical_u_handle: 解析解函数句柄 u(x,y) (若计算导数误差则是导数函数句柄)
%der_x, der_y:        导数阶数 alpha1, alpha2计算 L2 误差( 0, 0)
%
% 输出:
%   total_error: 计算得到的 L2 范数误差 (或 H1 半范数)

total_error_squared = 0;

% 准备高精度高斯积分
[gauss_weight, gauss_xi, gauss_eta] = gaussValues_2d(4); 

number_of_elements = size(T_mesh, 2);
num_local_basis = size(T_trial, 1);

% FOR n = 1, ..., N
for n = 1:number_of_elements
    
    % 获取当前单元的几何顶点 (用于坐标映射)
    geo_idx = T_mesh(:, n);
    x1 = P_mesh(1, geo_idx(1)); y1 = P_mesh(2, geo_idx(1));
    x2 = P_mesh(1, geo_idx(2)); y2 = P_mesh(2, geo_idx(2));
    x3 = P_mesh(1, geo_idx(3)); y3 = P_mesh(2, geo_idx(3));
    
    vertices = P_mesh(:, geo_idx);
    
    % 计算雅可比行列式 |J| (用于 dxdy = |J| dxi deta)
    det_J = (x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1);
    abs_det_J = abs(det_J);
    
    % 获取当前单元的数值解系数 u_Tb(k,n)
    % sum(u_Tb * psi)
    local_coeffs = solution(T_trial(:, n));
    
    %单元内高斯积分

    element_integral = 0;
    
    for k = 1:length(gauss_weight)
        xi = gauss_xi(k);
        eta = gauss_eta(k);
        w = gauss_weight(k);
        
        % 坐标映射: 参考 -》 物理
        lambda1 = 1 - xi - eta;
        lambda2 = xi;
        lambda3 = eta;
        
        x_phys = x1*lambda1 + x2*lambda2 + x3*lambda3;
        y_phys = y1*lambda1 + y2*lambda2 + y3*lambda3;
        
        % 计算解析解的值 (或其导数值)
        % part: partial^(a1+a2) u / ...
        u_exact_val = feval(analytical_u_handle,x_phys, y_phys);
        
        % B. 计算数值解 u_h 的值 (或其导数值)
        % part: sum( u_Tb * partial psi )
        u_fem_val = 0;
        for i = 1:num_local_basis
            % 调用基函数计算
            basis_val = FE_basis_local_fun_2D(x_phys, y_phys, vertices, ...
                                              trial_basis_type, i, der_x, der_y);
            u_fem_val = u_fem_val + local_coeffs(i) * basis_val;
        end
        
        % 计算被积函数 (差的平方)
        diff = u_exact_val - u_fem_val;
        
        % 累加积分: (u - u_h)^2 * |J| * w
        element_integral = element_integral + w * (diff^2) * abs_det_J;
    end
    
    % 将单元误差累加到总误差
    total_error_squared = total_error_squared + element_integral;
    
end

total_error = sqrt(total_error_squared);

end