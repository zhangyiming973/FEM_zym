function element_error_squared = Gauss_quad_error_1d(gauss_point, gauss_weight, vertices, local_solution_coeffs, trial_basis_type, s, analytical_deriv_handle)
%需要计算每个单元的积分值
%s是导数阶数，计算L2/H1范数误差
% 1. 获取高斯点和权重 (与 Gauss_quad_matrix_1d.m 相同)
xl = vertices(1);
xr = vertices(2);
gauss_point_local  = xl + (xr - xl)*(1+gauss_point)/2;
gauss_weight_local = gauss_weight*(xr - xl)/2;

% 在高斯点处计算 u_true^(s) (解析解)
u_true_deriv = feval(analytical_deriv_handle, gauss_point_local);


% 在高斯点处计算 u_fem^(s) (有限元数值解)

num_local_basis = length(local_solution_coeffs);
u_fem_deriv = zeros(size(gauss_point_local));

for k = 1:num_local_basis%所有基函数数值
    basis_deriv_val = FE_basis_local_fun_1D(gauss_point_local, vertices, trial_basis_type, k, s);%获取s阶导数
    u_fem_deriv = u_fem_deriv + local_solution_coeffs(k) * basis_deriv_val;
    %这是一个1*3or 1*4的行向量，由高斯点数量决定
end

integrand = (u_true_deriv - u_fem_deriv).^2;
element_error_squared = sum(gauss_weight_local .* integrand);


end