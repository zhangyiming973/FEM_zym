%{
基本结构：
Algorithm 4 (upgraded version):

- Initialize the matrix: A = sparse(N_b^test, N_b^trial);
- Compute the integrals and assemble them into A:

  FOR n = 1, ..., N:网格总数
    FOR α = 1, ..., N_b^trial:每个单元上的试探函数基函数个数
      FOR j = 1, ..., N_b^test:每个单元上的测试函数基函数个数
        Compute r = ∫_{x_n}^{x_n+1} cφ_{nα}^{(r)}ψ_{nβ}^{(s)} dx;
        Add r to A(T_b(i, n), T_b(α, n));
    END
  END
END
T_b(i, n): 局部到全局自由度映射函数

φ_{nα}^{(r)}: 第n个单元上的第α个试探基函数

ψ_{nβ}^{(s)}: 第n个单元上的第β个测试基函数

%}
function A=assemble_matrix_1D(matrix_size,gauss_point,gauss_weight,function_name, P, T, T_trial, T_test, trial_basis_type, trial_basis_der, test_basis_type, test_basis_der)

A=sparse(matrix_size(1),matrix_size(2));

num_local_trial = size(T_trial,1); % how many basis within an element.
num_local_test  = size(T_test ,1);
number_of_elements=size(T,2);
%大规模计算spars更好,小规模用zero
for n=1:number_of_elements

    vertices=P(T(:,n));%最重要的一环（遍历）：P括号中的内容存储在T矩阵里，第n个顶点的第k个坐标
    
    for alpha=1:num_local_trial
        
        for beta=1:num_local_test
            
            int_value=Gauss_quad_matrix_1d(gauss_point,gauss_weight,vertices,function_name,trial_basis_type,alpha,trial_basis_der, test_basis_type,beta,test_basis_der);



            i=T_test(beta,n);
            j=T_trial(alpha,n);
            
            A(i,j)=A(i,j)+int_value;

        end

    end


end



