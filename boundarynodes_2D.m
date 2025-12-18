function boundary_nodes = boundarynodes_2D(P_fem, bc_flags)

% 输入:
%   P_fem: 有限元节点坐标矩阵 (2 x N_nodes)
%   bc_flags: 1x4 向量，定义四条边的边界类型 [下, 右, 上, 左]
%             例如: [-1, -2, -1, -2] 表示上下为 Dirichlet(-1)，左右为 Neumann(-2)

% 输出:
%   boundary_nodes: (2 x N_boundary_nodes) 矩阵
%   第一行: 边界类型 (Type)
%   第二行: 全局节点索引 (Global Index)
%如果有冲突，优先保障狄雷克雷边界条件。

x = P_fem(1, :);
y = P_fem(2, :);

x_min = min(x);
x_max = max(x);
y_min = min(y);
y_max = max(y);

tol = 1e-8;

idx_bottom = find(abs(y - y_min) < tol);
idx_right = find(abs(x - x_max) < tol);
idx_top = find(abs(y - y_max) < tol);
idx_left = find(abs(x - x_min) < tol);
%顺序：下右上左
%识别边界点
node_types = zeros(1, size(P_fem, 2));
%初始化

type_bottom = bc_flags(1);
type_right  = bc_flags(2);
type_top    = bc_flags(3);
type_left   = bc_flags(4);
%想定义哪条边的边界条件就定义哪条边的


%赋予边界类型，注意狄雷克雷最优先，所以最后处理

if type_bottom ~= -1, node_types(idx_bottom) = type_bottom; end
if type_right  ~= -1, node_types(idx_right)  = type_right; end
if type_top    ~= -1, node_types(idx_top)    = type_top; end
if type_left   ~= -1, node_types(idx_left)   = type_left; end

% 最后处理所有 Dirichlet (-1) 的边界 (覆盖角点)
if type_bottom == -1, node_types(idx_bottom) = -1; end
if type_right  == -1, node_types(idx_right)  = -1; end
if type_top    == -1, node_types(idx_top)    = -1; end
if type_left   == -1, node_types(idx_left)   = -1; end


%组装输出矩阵


is_boundary = node_types ~= 0;
boundary_indices = find(is_boundary);
boundary_types = node_types(boundary_indices);


% 第一行: Type
% 第二行: Global Index
boundary_nodes = [boundary_types; boundary_indices];

end




