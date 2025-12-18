%保证程序通用性，试函数和
%1.input a b （左右界） N（网格数） generate mesh information martics P and T
%FE information matrices Pb Tb for trial and test func.s
%assemble A by algorithm 4
%assemble vector b  by algorithm 5
%assemble Drichlet BC by algorithm 6
%solve AX=b by using a direct method of matlab
%重要：网格节点与有限元节点分开
%程序开发流程：子程序搭结构，第二轮回填细节。
%basis_type_trial==101:1D linear
%basis_type_trial==102:二次元
function [solution, err]=FE_solver_1D_poisson(number_of_elements,trial_basis_type,test_basis_type,gauss_type)

[P_mesh, T_mesh, P_trial, T_trial] = generate_PbTb(number_of_elements, trial_basis_type);

% 同时生成几何网格(P_mesh, T_mesh)和试探空间(P_trial, T_trial)
[P_mesh, T_mesh, P_trial, T_trial] = generate_PbTb(number_of_elements, trial_basis_type);

% 同时生成几何网格(P_mesh, T_mesh)和测试空间(P_test, T_test)
[P_mesh, T_mesh, P_test, T_test] = generate_PbTb(number_of_elements, test_basis_type);


matrix_size = [length(P_test),length(P_trial)];
vector_size = length(P_test);

[gauss_weight,gauss_point] = gaussValues_1d(gauss_type); % both gauss_point and weight are row vector.


basis_der_x_test_A = 1;
basis_der_x_trial_A = 1;

% <-- 4. 将 P_mesh 和 T_mesh 作为 "几何网格" 传入
A = assemble_matrix_1D(matrix_size, gauss_point, gauss_weight, 'function_c', ...
                       P_mesh, T_mesh, T_trial, T_test, ...
                       trial_basis_type, basis_der_x_trial_A, ...
                       test_basis_type, basis_der_x_test_A);

test_basis_der_b = 0;
% <-- 5. 将 P_mesh 和 T_mesh 作为 "几何网格" 传入
b = assemble_vector_1D(gauss_point, gauss_weight, 'function_f', ...
                       test_basis_type, vector_size, ...
                       P_mesh, T_mesh, T_test, test_basis_der_b);

[A, b] = BC_process(A, b);

solution = A\b;

% <-- 6. 使用 P_trial (有限元节点坐标) 来绘图
plot(P_trial, solution, 'o', 'lineWidth', 1)
hold on 
plot(P_trial, P_trial.*cos(P_trial), '-*', 'lineWidth', 1)

% <-- 7. 使用 P_trial (有限元节点坐标) 来计算误差
% err = norm(solution'-P_trial.*cos(P_trial))/sqrt(num_of_element);  % L2 error
err = max(abs(solution' - P_trial.*cos(P_trial))); % L_inf error

legend("numerical solution", "real solution")
title(['L_inf error: ', num2str(err)]) % 在标题中显示误差
hold off

end