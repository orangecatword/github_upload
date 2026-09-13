% 综合IEEE33与IEEE123的网络规模,以及它们的网络结构信息 考虑到一阶节点度与二阶节点度
% 节点47为IEEE123网络最重要的节点
% 排行前10的节点为[47 48 49 44 86 76 72 77 65 66]
function [A, W] = node_importance(branch,Importance,pload) 
%BUILD_NODE_MATRIX 构造IEEE 123节点的一阶、二阶加权连接矩阵
%
% 输入：
%   branch - 2×M节点连接矩阵
%            branch(1,k)：第k条线路的输入节点
%            branch(2,k)：第k条线路的输出节点
%
% 输出：
%   W        - 节点重要性
%   A        - 33×33加权节点连接矩阵
%              一阶连接：1
%              二阶连接：0.5
%              其他情况：0
%
%   A_first  - 一阶邻接矩阵
%   A_second - 二阶邻接矩阵

    %% 1. 基本参数
    N = 123;

    %% 2. 构造一阶邻接矩阵
    A_first = zeros(N,N);

    for k = 1:size(branch,1)
        node_i = branch(k,1);
        node_j = branch(k,2);

        % 配电线路按照无向边处理
        A_first(node_i,node_j) = 1;
        A_first(node_j,node_i) = 1;
    end

    % 对角线元素保持为0
    A_first(1:N+1:end) = 0;

    %% 3. 利用邻接矩阵平方寻找二阶连接
    % (A_first^2)(i,j)>0说明i到j至少存在一条长度为2的路径
    twoStepPath = A_first^2;

    % 初步确定存在长度为2路径的节点对
    A_second = twoStepPath > 0;

    % 排除一阶直接连接的节点
    A_second(A_first == 1) = 0;

    % 排除节点自身
    A_second(1:N+1:end) = 0;

    % 转换为数值矩阵
    A_second = double(A_second);

    %% 4. 构造最终加权连接矩阵
    A = A_first + 0.5 * A_second;
    B = eye(123);

    %% 5. 计算最终的节点重要性
    % 改进想法-既要考虑连接的其他节点,也要考虑节点自身
    W = (A+B)*(pload.*Importance'); % 同时考虑负荷值和负荷权重,公式需验证
   
end