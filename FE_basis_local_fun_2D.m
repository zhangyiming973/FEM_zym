%定义基函数
function result=FE_basis_local_fun_2D(x,y,vertices,basis_type,basis_index,basis_der_x,basis_der_y) 
% 输入:
%   x, y: 物理单元内的坐标点 物理坐标
%   vertices: 物理单元顶点坐标 [x1, x2, x3; y1, y2, y3] (2x3矩阵)
%   basis_type: 101/201(线性), 102/202(二次)
%   basis_index: 基函数编号
%   basis_der_x, basis_der_y: 对物理坐标 x, y 的导数阶数
%
% 依据:
%   - 几何变换公式 (x,y -> hat_x, hat_y)
%   - 导数变换链式法则
%   - 二次单元基函数定义
x1 = vertices(1,1); y1 = vertices(2,1);
x2 = vertices(1,2); y2 = vertices(2,2);
x3 = vertices(1,3); y3 = vertices(2,3);
% 计算差分项
x_21 = x2 - x1;
x_31 = x3 - x1;
x_13 = x1 - x3;
y_31 = y3 - y1;
y_21 = y2 - y1;
y_12 = y1 - y2;

% 计算雅可比行列式 |Jn|
det_J = x_21 * y_31 - x_31 * y_21;

%从物理坐标 (x,y) 计算参考坐标 (hat_x, hat_y)
hat_x = ( y_31 * (x - x1) - x_31 * (y - y1) ) / det_J;
hat_y = ( -y_21 * (x - x1) + x_21 * (y - y1) ) / det_J;


%3. 在参考单元上计算基函数值及其对 hat_x, hat_y 的导数

val_hat = 0;
d_val_d_hat_x = 0;   % d(psi_hat) / d(hat_x)
d_val_d_hat_y = 0;   % d(psi_hat) / d(hat_y)
d2_val_d_hat_x2 = 0; % d^2(psi_hat) / d(hat_x)^2
d2_val_d_hat_y2 = 0;
d2_val_d_hat_xy = 0; % 混合导数



if basis_type == 201  % 线性单元 (3节点)
    % N1 = 1 - hat_x - hat_y
    % N2 = hat_x
    % N3 = hat_y
    switch basis_index
        case 1
            val_hat = 1 - hat_x - hat_y;
            d_val_d_hat_x = -1;
            d_val_d_hat_y = -1;
        case 2
            val_hat = hat_x;
            d_val_d_hat_x = 1;
            d_val_d_hat_y = 0;
        case 3
            val_hat = hat_y;
            d_val_d_hat_x = 0;
            d_val_d_hat_y = 1;
        otherwise
            error('201类型单元只有3个基函数');
    end

elseif basis_type==202%二次元6个基函数
   switch basis_index
       case 1
            val_hat = 2*hat_x^2 + 2*hat_y^2 + 4*hat_x*hat_y - 3*hat_y - 3*hat_x + 1;
            d_val_d_hat_x = 4*hat_x + 4*hat_y - 3;
            d_val_d_hat_y = 4*hat_y + 4*hat_x - 3;
            d2_val_d_hat_x2 = 4; d2_val_d_hat_y2 = 4; d2_val_d_hat_xy = 4;
       case 2
            val_hat = 2*hat_x^2 - hat_x;
            d_val_d_hat_x = 4*hat_x - 1;
            d_val_d_hat_y = 0;
            d2_val_d_hat_x2 = 4; d2_val_d_hat_y2 = 0; d2_val_d_hat_xy = 0;
       case 3
            val_hat = 2*hat_y^2 - hat_y;
            d_val_d_hat_x = 0;
            d_val_d_hat_y = 4*hat_y - 1;
            d2_val_d_hat_x2 = 0; d2_val_d_hat_y2 = 4; d2_val_d_hat_xy = 0;
       case 4
            val_hat = -4*hat_x^2 - 4*hat_x*hat_y + 4*hat_x;
            d_val_d_hat_x = -8*hat_x - 4*hat_y + 4;
            d_val_d_hat_y = -4*hat_x;
            d2_val_d_hat_x2 = -8; d2_val_d_hat_y2 = 0; d2_val_d_hat_xy = -4;
       case 5
            val_hat = 4*hat_x*hat_y;
            d_val_d_hat_x = 4*hat_y;
            d_val_d_hat_y = 4*hat_x;
            d2_val_d_hat_x2 = 0; d2_val_d_hat_y2 = 0; d2_val_d_hat_xy = 4;
       case 6
            val_hat = -4*hat_y^2 - 4*hat_x*hat_y + 4*hat_y;
            d_val_d_hat_x = -4*hat_y;
            d_val_d_hat_y = -8*hat_y - 4*hat_x + 4;
            d2_val_d_hat_x2 = 0; d2_val_d_hat_y2 = -8; d2_val_d_hat_xy = -4;
       otherwise
           error('二次元六个基函数')
end
%以下分析导数，将参考坐标系映射回物理坐标系
if basis_der_x == 0 && basis_der_y == 0
    result = val_hat;
elseif basis_der_x == 1 && basis_der_y == 0
    result = (d_val_d_hat_x * y_31 + d_val_d_hat_y * y_12) / det_J;

elseif basis_der_x == 0 && basis_der_y == 1
    result = (d_val_d_hat_x * x_13 + d_val_d_hat_y * x_21) / det_J;

elseif basis_der_x == 2 && basis_der_y == 0
        a1 = d2_val_d_hat_x2 * (y_31^2);
        a2 = 2 * d2_val_d_hat_xy * (y_31 * y_12);
        a3 = d2_val_d_hat_y2 * (y_12^2);
        result = (a1+a2+a3) / (det_J^2);
elseif basis_der_x == 0 && basis_der_y == 2
        a1 = d2_val_d_hat_x2 * (x_13^2);
        a2 = 2 * d2_val_d_hat_xy * (x_13 * x_21);
        a3 = d2_val_d_hat_y2 * (x_21^2);
        result = (a1+a2+a3) / (det_J^2);
elseif basis_der_x == 1 && basis_der_y == 1
        a1 = d2_val_d_hat_x2 * (x_13 * y_31);
        a2 = d2_val_d_hat_xy * (x_13 * y_12); 
        a3 = d2_val_d_hat_xy * (x_21 * y_31);
        a4 = d2_val_d_hat_y2 * (x_21 * y_12);
        result = (a1+a2+a3+a4) / (det_J^2);
else
        error('不支持过高导数阶数.');
end
end

end