function [A,b]=BC_process(A,b,BC_type,function_c)
    %三种边界条件，狄雷克雷-01，狄雷克雷+纽曼混合02，纽曼+狄雷克雷混合03，
    if BC_type==01
    g_a=0;
    g_b=cos(1);%位移限制边界条件
        A(1,:)=0;
        A(end,:)=0;
        A(1,1)=1;
        A(end,end)=1;
        b(1)=g_a;
        b(end)=g_b;
    elseif BC_type==02%u'a=ra ub=gb
        r_a=1;%纽曼边界
        g_b=cos(1);%位移限制边界条件
        c_a=function_c(0);
        A(end,:)=0;
        A(end,end)=1;
        b(end)=g_b;
        b(1) = b(1) - r_a * c_a;
    elseif BC_type==03%u'b=rb ua=ga
        r_b=1;%纽曼边界
        g_a=0;%位移限制边界条件
        A(1,:) = 0;
        A(1,1) = 1;
        b(1) = g_a;

        c_b=function_c(1);

        b(end)=b(end)+r_b*c_b;

        elseif BC_type==04%u'a+qaua=1 ub=gb混合边界条件
        % BC 1: u(1) = cos(1) (Dirichlet at x=b)
        g_b=cos(1);%v1=0

        c_a=function_c(0);
        A(end,:)=0;
        A(end,end)=1;
        A(1,1)=A(1,1)- c_a;
        b(1) =b(1) - c_a ;

        b(end)=g_b;

    end
end