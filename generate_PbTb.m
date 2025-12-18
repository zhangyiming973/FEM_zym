% function [Pb,Tb]=generate_PbTb(number_of_element,basis_type) % <-- 旧签名
function [P_mesh, T_mesh, P_fem, T_fem] = generate_PbTb(number_of_element, basis_type) % <-- 新签名

N = number_of_element;

% 1. 始终生成基础的 "线性几何网格" (P_mesh, T_mesh)
P_mesh = linspace(0, 1, N + 1); % [0, 1] 区间 N 个单元, N+1 个顶点
T_mesh = [1:N; 2:N+1];      % 存储 "顶点" 索引

% 2. 根据 basis_type 生成 "有限元空间" (P_fem, T_fem)
switch basis_type
    case 101 % 线性元
        % 有限元节点与网格节点重合
        P_fem = P_mesh;
        T_fem = T_mesh;
        
    case 102 % 二次元 (示例)
        % (这是您未来需要补充的)
        % 有限元节点 = 顶点 + 单元中点
        P_fem = linspace(0, 1, 2*N + 1); % N个单元, 2N+1 个节点
        
        % 自由度映射 (每个单元3个自由度)
        T_fem = zeros(3, N);
        for n = 1:N
            T_fem(:, n) = [2*n-1; 2*n+1; 2*n]; % [左顶点, 右顶点, 中点]
        end
        % 注意：FE_basis_local_fun_1D.m 中也要实现102的基函数
        
    case 103
        % ...
end
end