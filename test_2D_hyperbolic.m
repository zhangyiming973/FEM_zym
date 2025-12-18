% test_hyperbolic.m
% 测试二维双曲型方程 (波动方程) 求解器
% 物理场景: 四周固定的方形薄膜自由振动 (驻波)
% 方程: u_tt - div(1 * grad(u)) = 0

clear; clc; close all;

% === 1. 定义物理参数 ===
% 波速 c^2 = 1 -> c_val = 1
c_func = @(x,y) 1.0 + 0*x; 

% 自由振动，无外力 f = 0
f_func_t = @(x,y,t) 0 * x;

% === 2. 定义解析解 (用于验证和初始条件) ===
% 驻波解: u = sin(pi*x)sin(pi*y) * cos(omega*t)
% omega = c * sqrt(k^2 + l^2) = 1 * sqrt(pi^2 + pi^2) = sqrt(2)*pi
omega = sqrt(2) * pi;
u_exact_t = @(x,y,t) sin(pi*x) .* sin(pi*y) .* cos(omega * t);

% 初始位移 u0 (t=0)
u0_func = @(x,y) u_exact_t(x,y, 0);

% 初始速度 v0 (t=0) -> 对 u 求导 -> -omega * sin(...) * sin(...) * sin(0) = 0
v0_func = @(x,y) 0 * x; 

% === 3. 边界条件 ===
% 四周固定 (Dirichlet = 0)
bc_flags = [-1, -1, -1, -1];

% 边界函数随时间变化 (虽然这里始终是0，但为了通用性写作 t 的函数)
bc_funcs_t.dirichlet = @(x,y,t) u_exact_t(x,y,t); 
% 占位
bc_funcs_t.bottom = 0; bc_funcs_t.right = 0;
bc_funcs_t.top = 0;    bc_funcs_t.left = 0;

% === 4. 求解配置 ===
x_range = [0, 1];
y_range = [0, 1];
N = 32;              % 网格密度
trial_type = 201;    % 线性单元
test_type = 201;
gauss_type = 4;

% 时间步进配置
% 注意：双曲型方程显式/半隐式求解通常有 CFL 条件限制 (dt < h/c)
% 这里我们用隐式中心差分，理论上无条件稳定，但精度仍需 dt 较小
time_config.dt = 0.01; 
time_config.T_end = 3.0; % 模拟 2 秒

% === 5. 运行求解器 ===
fprintf('开始求解波动方程...\n');
[u_final, P, T] = FE_solver_2D_hyperbolic(x_range, y_range, N, N, ...
                                          trial_type, test_type, gauss_type, ...
                                          bc_flags, ...
                                          f_func_t, c_func, ...
                                          bc_funcs_t, ...
                                          u0_func, v0_func, ...
                                          time_config, 1); % 1=Isotropic

% === 6. 误差验证 (t = T_end) ===
u_exact_final = @(x,y) u_exact_t(x,y, time_config.T_end);
err = compute_L2_error_2D(u_final, P, T, P, T, trial_type, u_exact_final, 0, 0);

fprintf('L2 Error at t=%.2f: %e\n', time_config.T_end, err);