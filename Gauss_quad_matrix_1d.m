function result = Gauss_quad_matrix_1d(gauss_point,gauss_weight,vertices,function_name,trial_basis_type,trial_basis_index,trial_basis_der, test_basis_type,test_basis_index,test_basis_der)
%function_name就是c（x）
xl = vertices(1);
xr = vertices(2);
gauss_point_local  = xl + (xr - xl)*(1+gauss_point)/2;% 从[-1,1]到[xl,xr]
gauss_weight_local = gauss_weight*(xr - xl)/2;% 权重缩放

trial_val = FE_basis_local_fun_1D(gauss_point_local,vertices,trial_basis_type,trial_basis_index,trial_basis_der);
test_val  = FE_basis_local_fun_1D(gauss_point_local,vertices, test_basis_type, test_basis_index, test_basis_der);

func_val  = feval(function_name,gauss_point_local);
%根据 function_name 找到对应的函数c(x)

%将 gauss_point_local 作为参数传递给该函数

%返回函数在给定参数下的计算结果

%result =  gauss_weight_local*(trial_val.*test_val.*func_val)';
result = sum(gauss_weight_local.*trial_val.*test_val.*func_val) ;
end