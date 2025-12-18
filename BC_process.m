function [A,b]=BC_process(A,b)

g_a=0;
g_b=cos(1);%位移限制边界条件
    A(1,:)=0;
    A(end,:)=0;
    A(1,1)=1;
    A(end,end)=1;
    b(1)=g_a;
    b(end)=g_b;
end