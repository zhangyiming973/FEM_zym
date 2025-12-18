function total_error=compute_L2_error_1D(solution, P_mesh, T_mesh, T_trial, trial_basis_type, s, analytical_deriv_handle)

[gauss_weight, gauss_point] = gaussValues_1d(4);
total_error_squared=0;
number_of_elements=size(T_mesh,2);
 for n=1:number_of_elements

    vertices=P_mesh(T_mesh(:,n));%去 P_mesh 里把第 5 个和第 6 个数取出来
    % 最重要的一环（遍历）：P括号中的内容存储在T矩阵里，第n个顶点的第k个坐标
    % 获取当前单元对应的 solution 系数 (u_Tb(k,n))
    local_solution_coeffs = solution(T_trial(:, n));% (u_Tb(k,n))提取出当前单元所对应的解


    element_integral = Gauss_quad_error_1d(gauss_point, gauss_weight, vertices, ...
        local_solution_coeffs, trial_basis_type, s, analytical_deriv_handle);
 
    total_error_squared = total_error_squared + element_integral;
end


total_error=sqrt(total_error_squared);

end