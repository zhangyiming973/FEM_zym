%定义基函数
function result=FE_basis_local_fun_1D(x,vertices,basis_type,basis_index,basis_der_x) 
%vertices为顶点信息。--n
%basis_type用于区分test 和 trial，在积分时用一套调用逻辑，用test trial区分
%basis_der_x导数阶次
%basis_index  alpha 和 beta and theta，第几个基函数？--a b
%101 1D线性 102 1D二次……
%注意基函数顺序：左右中

h=vertices(2)-vertices(1);

if basis_type==101%1D线性

    if basis_index==1%对应α

        if basis_der_x==0
            
            result=(vertices(2)-x)/h;

        elseif basis_der_x==1

            result=-1/h;

        elseif basis_der_x>=2 
            result=0;
        else
        error('导数维数不符合要求')
        end
 
    elseif basis_index==2%对应β

        if basis_der_x==0
            
            result=(x-vertices(1))/h;

        elseif basis_der_x==1

            result=1/h;

        elseif basis_der_x>=2 
            result=0;
        else
        error('导数维数不符合要求')
        end


    else
        error('wrong input for basis_index,一维只有两个ab')
    
    end

elseif basis_type==102%二次元三个基函数
    h=vertices(2)-vertices(1);
    x_hat=(x-vertices(1))/h;
    d_x_hat=1/h;

if basis_index==1%对应α

        if basis_der_x==0
            
            result=2*x_hat.^2-3*x_hat+1;

        elseif basis_der_x==1

            result=(4*x_hat-3)*d_x_hat;

        elseif basis_der_x==2 
            result=4/(h*h);
        elseif basis_der_x>=3 
            result=0;  
        else
        error('导数维数不符合要求')
        end
 
    elseif basis_index==2%对应β

        if basis_der_x==0
            
            result=2*x_hat.^2-x_hat;

        elseif basis_der_x==1

            result=(4*x_hat-1)*d_x_hat;

        elseif basis_der_x==2 
            result=4/(h*h);
        elseif basis_der_x>=3 
            result=0;  
        else
        error('导数维数不符合要求')
        end
    elseif basis_index==3%对应gamma，二次元有三个基函数

        if basis_der_x==0
            
            result=-4*x_hat.^2+4*x_hat;

        elseif basis_der_x==1

            result=(-8*x_hat + 4) * d_x_hat;

        elseif basis_der_x==2 
            result=-8 * (d_x_hat)^2;
        elseif basis_der_x>=3 
            result=0;  
        else
        error('导数维数不符合要求')
        end

    else
        error('wrong input for basis_index,一维只有两个ab')
    

end

%elseif basis_type==103
end