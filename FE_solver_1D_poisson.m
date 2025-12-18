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
function [solution, err]=FE_solver_1D_poisson(number_of_elements,trial_basis_type,test_basis_type,gauss_type,BC_type)

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

% <4. 将 P_mesh 和 T_mesh 作为 "几何网格" 传入
A = assemble_matrix_1D(matrix_size, gauss_point, gauss_weight, 'function_c', ...
                       P_mesh, T_mesh, T_trial, T_test, ...
                       trial_basis_type, basis_der_x_trial_A, ...
                       test_basis_type, basis_der_x_test_A);

test_basis_der_b = 0;
% <5. 将 P_mesh 和 T_mesh 作为 "几何网格" 传入
b = assemble_vector_1D(gauss_point, gauss_weight, 'function_f', ...
                       test_basis_type, vector_size, ...
                       P_mesh, T_mesh, T_test, test_basis_der_b);
f_handle_c = @function_c;
[A, b] = BC_process(A, b,BC_type,f_handle_c);

solution = A\b;

% %根据不同边界条件设计解析解用于误差分析
% 
% if BC_type==01
%     C=0;
%     D=0;
% elseif BC_type==02
%     C=0;
%     D=0;
% elseif BC_type==04
%     C=0;
%     D=0;
% elseif BC_type==03
%     C=(1 - cos(1) + sin(1)) * exp(1);
%     D=C;
% end
% 
% u_analytical = P_trial.*cos(P_trial) - C*exp(-P_trial) + D;
% 
% % <6. 使用 P_trial (有限元节点坐标) 来绘图
% plot(P_trial, solution, 'o', 'lineWidth', 1)
% hold on 
% plot(P_trial, u_analytical, '-*', 'lineWidth', 1)%使用新的解析解绘图！
% 
% hold off
% err = max(abs(solution' - u_analytical)); % <-- 使用新的解析解
% 
% % 更新图例和标题以反映当前情况
% legend("numerical solution", "real solution (BC_type=" + num2str(BC_type) + ")")
% title(['L_inf error: ', num2str(err)]) % 在标题中显示误差
%以上为无穷范误差

if BC_type == 1 || BC_type == 2 || BC_type == 4
    A = 0; B = 0;
elseif BC_type == 3
    A = (1 - cos(1) + sin(1)) * exp(1);
    B = A;
else
    error('未知的 BC_type。');
end

% 定义解析解 u(x) 和 u'(x) 的函数句柄
% 这是 s=0 时的解析解
u_analytical_s0 = @(x) x.*cos(x) - A*exp(-x) + B;
% 这是 s=1 时的解析解
u_analytical_s1 = @(x) cos(x) - x.*sin(x) + A*exp(-x);

% 在有限元节点处计算解析解 (用于绘图)
u_analytical_at_nodes = u_analytical_s0(P_trial);
plot(P_trial, solution, 'o', 'lineWidth', 1)
hold on 
plot(P_trial, u_analytical_at_nodes, '-*', 'lineWidth', 1) 
hold off

err_L2 = compute_L2_error_1D(solution, P_mesh, T_mesh, T_trial, trial_basis_type, 0, u_analytical_s0);
err_H1 = compute_L2_error_1D(solution, P_mesh, T_mesh, T_trial, trial_basis_type, 1, u_analytical_s1);


err = err_L2;

legend("numerical solution", "real solution (BC_type=" + num2str(BC_type) + ")")
title(['L2 error: ', num2str(err_L2), ' | H1 error: ', num2str(err_H1)]) % 在标题中显示新的误差

end

