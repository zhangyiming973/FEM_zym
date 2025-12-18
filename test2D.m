% % 清空环境
% clear; clc; close all;
% 
% % === 1. 准备函数句柄 (关键步骤) ===
% % 使用 @ 符号加上你的文件名
% % 这样 MATLAB 就会自动找到 function_c.m 等文件
% c_handle = @function_c;
% f_handle = @function_f;
% u_handle = @function_u;
% 
% % === 2. 定义边界条件 ===
% % 标记: [-1, -1, -1, -1] 全 Dirichlet
% bc_flags = [-1, -1, -1, -1]; 
% 
% % 定义边界函数结构体
% % 如果 Dirichlet 边界的值恰好就是解析解的值，可以直接复用 u_handle
% bc_funcs.dirichlet = u_handle; 
% 
% % 其他边界设为 0 (即便用不到也要定义，防止报错)
% bc_funcs.bottom = @(x,y) 0;
% bc_funcs.right  = @(x,y) 0;
% bc_funcs.top    = @(x,y) 0;
% bc_funcs.left   = @(x,y) 0;
% 
% % === 3. 设置参数 ===
% x_range = [0, 1];
% y_range = [0, 1];
% N = 20;
% trial_type = 202; % 二次单元
% test_type = 202;
% gauss_type = 4;
% a=1;%各向同性
% % === 4. 调用求解器 ===
% % 注意参数位置: f 在前, c 在后 (根据你的函数定义顺序)
% [sol, err, P, T] = FE_solver_2D_poisson(x_range, y_range, N, N, ...
%                                         trial_type, test_type, ...
%                                         gauss_type, bc_flags, ...
%                                         f_handle, c_handle, ...  % <--- 传入句柄
%                                         bc_funcs, u_handle,a);     % <--- 传入句柄
% 
% % === 5. 输出结果 ===
% fprintf('计算完成！\n');
% fprintf('L2 误差: %e\n', err);
%以上为各向同性，以下为各向异性
% 测试各向异性泊松方程: -div(C * grad(u)) = f

% test_anisotropy_case.m
% 目的: 测试 FE_solver_2D_poisson 的各向异性功能 (a=2)
% 物理方程: -div(C * grad(u)) = f









clear; clc; close all;

% === 1. 定义各项异性系数矩阵 C ===
% 我们定义一个 2x2 的张量 C = [2, 1; 1, 2]
% 这意味着 x, y 方向导通性不同，且存在耦合 (非对角项不为0)
% 必须构造结构体传入
c_struct.c11 = @(x,y) 2.0 + 0*x; 
c_struct.c12 = @(x,y) 1.0 + 0*x; % 交叉项 c12
c_struct.c21 = @(x,y) 1.0 + 0*x; % 交叉项 c21
c_struct.c22 = @(x,y) 2.0 + 0*x;





% === 2. 定义解析解 u (用于验证) ===
% u = sin(pi*x) * sin(pi*y)

u_handle = @function_u;
% === 3. 推导对应的源项 f ===
% 理论推导: f = -div(C * grad u)
% 对于常数矩阵 C，展开为: - ( c11*u_xx + (c12+c21)*u_xy + c22*u_yy )
% 
% 已知:
% u_xx = -pi^2 * sin(pi*x)sin(pi*y)
% u_yy = -pi^2 * sin(pi*x)sin(pi*y)
% u_xy =  pi^2 * cos(pi*x)cos(pi*y)
%
% 代入 C=[2 1; 1 2]:
% f = - [ 2*(-pi^2*u) + (1+1)*(pi^2*cos*cos) + 2*(-pi^2*u) ]
%   = 4*pi^2*u - 2*pi^2*cos(pi*x)cos(pi*y)

f_handle = @(x,y) (4 * pi^2) * (sin(pi*x).*sin(pi*y)) ...
                - (2 * pi^2) * (cos(pi*x).*cos(pi*y));

% === 4. 定义边界条件 (Dirichlet) ===
% 四周设为 Dirichlet，值取解析解的值 (即 0)
bc_flags = [-1, -1, -1, -1]; 

bc_funcs.dirichlet = u_handle; 
% 占位防止报错
bc_funcs.bottom = @(x,y) 0; bc_funcs.right = @(x,y) 0;
bc_funcs.top    = @(x,y) 0; bc_funcs.left  = @(x,y) 0;

% === 5. 设置求解参数 ===
x_range = [0, 1];
y_range = [0, 1];
N = 20;             % 网格密度
trial_type = 202;   % 二次单元 (精度高)
test_type  = 202;
gauss_type = 4;     % 高斯积分点数

% === 6. 调用求解器 (关键参数 a=2) ===
fprintf('开始各向异性测试 (a=2)...\n');

% 注意函数签名的最后一个参数是 2
[sol, err, P, T] = FE_solver_2D_poisson(x_range, y_range, N, N, ...
                                        trial_type, test_type, ...
                                        gauss_type, bc_flags, ...
                                        f_handle, c_struct, ...   % 传入结构体
                                        bc_funcs, u_handle, ...
                                        2);                       % <--- a=2 (各向异性)

% === 7. 结果可视化 ===
fprintf('计算完成。\n');
fprintf('L2 误差: %e\n', err);

