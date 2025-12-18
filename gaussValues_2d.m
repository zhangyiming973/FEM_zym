function [w, xi, eta] = gaussValues_2d(n)
% GAUSSVALUES_2D 三角形区域的高斯积分点和权重
% 积分区域: Reference Triangle {(0,0), (1,0), (0,1)}
% 输入 n: 精度等级
% 输出: w(权重), xi(坐标1), eta(坐标2)

switch n
    case 1 % 1点积分 (精度: Linear)
        w = 0.5; % 三角形面积 1/2
        xi = 1/3; 
        eta = 1/3;
        
    case 2 % 3点积分 (精度: Quadratic, 常用)
        w = [1/6, 1/6, 1/6];
        xi = [1/6, 2/3, 1/6];
        eta = [1/6, 1/6, 2/3];
        
    case 4 % 4点积分 (精度: Cubic) - Hammer-Stroud
        % 这是一个常用的高精度公式
        w = [-27/96, 25/96, 25/96, 25/96] * 0.5; % 注意总权重和应为0.5
        xi = [1/3, 0.6, 0.2, 0.2];
        eta = [1/3, 0.2, 0.6, 0.2];
        
    otherwise
        error('Gauss type not supported yet in 2D.');
end
end