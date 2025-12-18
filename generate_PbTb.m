% function [Pb,Tb]=generate_PbTb(number_of_element,basis_type) %
function [P_mesh, T_mesh, P_fem, T_fem] = generate_PbTb(x_range, y_range, N_x, N_y, basis_type) 
% 输入:
%   x_range: [x_min, x_max]，例如 [0, 1]
%   y_range: [y_min, y_max]，例如 [0, 1]
%   N_x: x方向单元数量
%   N_y: y方向单元数量
%   basis_type: 201(线性), 202(二次)
x_min=x_range(1);
x_max=x_range(2);
y_min=y_range(1);
y_max=y_range(2);
%读出最大最小值
%划分网格

nx=N_x+1;
ny=N_y+1;
[X, Y] = meshgrid(linspace(x_min, x_max, nx), linspace(y_min, y_max, ny));

%转置后按列拉直,注意排列：按列
X = X'; Y = Y'; 
P_mesh = [X(:)'; Y(:)']; % 2行矩阵，先x后y

%201，202三角单元实现
% 单元总数 (每个矩形切成2个三角形)
num_elements = 2 * N_x * N_y;
T_mesh = zeros(3, num_elements);


idx = 0;
for j = 1:N_y
    for i = 1:N_x
        % 当前矩形单元的四个顶点索引 (基于 P_mesh 的线性索引)
        % 索引方式: (i-1)*ny + j
        n_bottom_left  = (i-1)*ny + j;
        n_bottom_right = i*ny + j;
        n_top_right    = i*ny + (j+1);
        n_top_left     = (i-1)*ny + (j+1);
        
        % 三角形 1 (下三角):下左-> 下右 -> 上左
        idx = idx + 1;
        T_mesh(:, idx) = [n_bottom_left; n_bottom_right; n_top_left];
        
        % 三角形 2 (上三角): 下右 -> 上右 -> 上左
        idx = idx + 1;
        T_mesh(:, idx) = [n_bottom_right; n_top_right; n_top_left];
    end
end
% --- 3. 生成有限元空间 (P_fem, T_fem) ---
if basis_type == 201 || basis_type == 101 % 线性单元
    P_fem = P_mesh;
    T_fem = T_mesh;
    
elseif basis_type == 202 || basis_type == 102 % 二次单元
    % 策略：生成一个 2倍密度的网格，包含了所有顶点和中点
    
    nx_fine = 2*N_x + 1;
    ny_fine = 2*N_y + 1;
    
    [X_fine, Y_fine] = meshgrid(linspace(x_min, x_max, nx_fine), linspace(y_min, y_max, ny_fine));
    X_fine = X_fine'; Y_fine = Y_fine';
    P_fem = [X_fine(:)'; Y_fine(:)'];
    
    % 二次单元有 6 个节点
    T_fem = zeros(6, num_elements);
    
    idx = 0;
    for j = 1:N_y
        for i = 1:N_x
            % 在细网格中的基准索引 (对应粗网格的 Bottom-Left)
            % 粗网格(i,j) -> 细网格(2i-1, 2j-1)
            % 计算细网格中的线性索引函数:
            get_node = @(r, c) (r-1)*ny_fine + c; 
            
            r_base = 2*i - 1; % Row index in fine grid (x direction logic in code above)
            c_base = 2*j - 1; % Col index in fine grid (y direction logic)
            %  meshgrid 转置后：行代表x变化，列代表y变化
            % X' 意味着第一维是 x 变化 (1..nx)，第二维是 y 变化 (1..ny)
            % 直接用 (i_fine, j_fine) 算线性索引：
            % node_idx = (i_fine - 1) * ny_fine + j_fine;
            
            idx_fine = @(I, J) (I-1)*ny_fine + J;
            
            % 获取当前矩形对应的 9 个关键点在 P_fem 中的索引
            n_BL = idx_fine(2*i-1, 2*j-1); % Bottom-Left Vertex
            n_BM = idx_fine(2*i,   2*j-1); % Bottom-Mid
            n_BR = idx_fine(2*i+1, 2*j-1); % Bottom-Right Vertex
            
            n_ML = idx_fine(2*i-1, 2*j);   % Mid-Left
            n_MM = idx_fine(2*i,   2*j);   % Center (Mid-Mid)
            n_MR = idx_fine(2*i+1, 2*j);   % Mid-Right
            
            n_TL = idx_fine(2*i-1, 2*j+1); % Top-Left Vertex
            n_TM = idx_fine(2*i,   2*j+1); % Top-Mid
            n_TR = idx_fine(2*i+1, 2*j+1); % Top-Right Vertex
            
            % 三角形 1 (下三角): BL -> BR -> TL
            % 6个节点顺序通常是：3个顶点 + 3个中点
            % 顺序约定：[V1, V2, V3, M12, M23, M31] (需与基函数定义对应)
            % 假设基函数顺序: 顶点1, 顶点2, 顶点3, 边1-2中点, 边2-3中点, 边3-1中点
            
            idx = idx + 1;
            % V1=BL, V2=BR, V3=TL
            % M12=BM, M23=Center(MM), M31=ML
            T_fem(:, idx) = [n_BL; n_BR; n_TL; n_BM; n_MM; n_ML];
            
            % 三角形 2 (上三角): BR -> TR -> TL
            % V1=BR, V2=TR, V3=TL
            % M12=MR, M23=TM, M31=Center(MM)
            idx = idx + 1;
            T_fem(:, idx) = [n_BR; n_TR; n_TL; n_MR; n_TM; n_MM];
        end
    end
    
else
    error('Unsupported basis type');
end

end






