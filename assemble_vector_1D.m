%{
基本结构：
Algorithm 5 (upgraded version):

- Initialize the matrix: b = sparse(N_b^test, 1);
- Compute the integrals and assemble them into b:

  FOR n = 1, ..., N:网格总数
    FOR β = 1, ..., N_b^trial:每个单元上的试探函数基函数个数
      
        Compute r = ∫_{x_n}^{x_n+1} f*ψ_{nβ}^{(s)} dx;
        Add r to b(T_b(β, n), 1);
    
  END
END
T_b(i, n): 局部到全局自由度映射函数

ψ_{nβ}^{(s)}: 第n个单元上的第β个测试基函数

%}
function b=assemble_vector_1D(gauss_point, gauss_weight,function_name,test_basis_type,vector_size,P,T,T_test,test_basis_der)

b=sparse(vector_size(1),1);%1列的矩阵
%大规模计算spars更好,小规模用zero
num_local_test=size(T_test,1);
for n=1:size(T,2)

    vertices=P(T(:,n));%最重要的一环（遍历）：P括号中的内容存储在T矩阵里，第n个顶点的第k个坐标
    %与矩阵组装不同，只需遍历一遍β
   
        for beta=1:num_local_test
            
            val=Gauss_quad_vector_1d(gauss_point,gauss_weight,vertices,function_name, test_basis_type,beta,test_basis_der);
            
            b(T_test(beta,n))=b(T_test(beta,n))+val;

        end

end






