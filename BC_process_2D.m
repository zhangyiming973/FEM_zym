function [A, b] = BC_process_2D(A, b, boundary_nodes, P_fem, bc_value_functions)
% BC_PROCESS_2D 处理二维边界条件

% 输入:
%   A, b: 全局刚度矩阵和载荷向量
%   boundary_nodes: 由 boundarynodes_2D 生成的 (2 x k) 矩阵
%                   Row 1: 类型 (-1 Dirichlet, -2 Neumann, -3 Robin)
%                   Row 2: 全局节点索引
%   P_fem: 节点坐标 (用于计算边界函数值)
%   bc_value_functions: 一个结构体，包含不同类型边界的函数句柄
%       .dirichlet = @(x,y) ...
%       .neumann =暂不在此处理
%
% 输出:
%   A, b: 施加了 Dirichlet 约束后的系统矩阵

% 获取边界节点总数
num_boundary_nodes = size(boundary_nodes, 2);

% 遍历所有边界节点
for k = 1:num_boundary_nodes
    
    type = boundary_nodes(1, k);      % 边界类型
    global_idx = boundary_nodes(2, k); % 全局节点编号
    
    % 获取节点坐标
    x_node = P_fem(1, global_idx);
    y_node = P_fem(2, global_idx);
    
    % --- 处理 Dirichlet 边界 (-1) ---
    % 将 A 的第 i 行置为 0，A(i,i)=1，b(i)=g(x,y)
    if type == -1
        
        % 1. 计算 Dirichlet 边界上的已知值 u = g(x,y)
        g_val = bc_value_functions.dirichlet(x_node, y_node);
        
        % 2. 修改刚度矩阵 A
        % 将该节点对应的整行清零 
        
        A(global_idx, :) = 0; 
        A(global_idx, global_idx) = 1;
        
        % 3. 修改载荷向量 b
        b(global_idx) = g_val;
        
    %  处理 Neumann (-2) 和 Robin (-3)
    elseif type == -2 || type == -3

        continue; 
        % 保持矩阵原样 (即默认 satisfy du/dn = 0)
        
    end
end

end