function [L_reach, componentSizes, componentID] = calculate_L_reach(Result, branch)
%CALCULATE_L_REACH 计算多条线路断开后的节点可达性损失指标
%
% 输入：
%   Result  - 断路边编号集合
%             支持单元格形式：Result = {a1, a2, ..., an}
%             支持数值形式：Result = [a1, a2, ..., an]
%
%   branch  - IEEE 33节点系统的线路连接矩阵
%             branch(k,1)和branch(k,2)是第k条线路的两个端点
%             大小为33×2，
%
% 输出：
%   L_reach        - 节点可达性损失指标
%   componentSizes - 故障后各连通分量的节点数
%   componentID    - 每个节点所属连通分量的编号
%
% 计算公式：
%   L_reach = 1 - sum_k[n_k(n_k-1)]/[N(N-1)]
%
% 其中：
%   N   为网络总节点数
%   n_k 为故障后第k个连通分量的节点数
    
    N = 33; % 节点数

    %% 1.删除故障线路
    branchStatus = true(32, 1);
    branchStatus(Result) = false;
    remainingBranches = branch(branchStatus, :);

    %% 2.建立故障后的无向图
    % 最后的参数N可以确保完全孤立的节点仍然保留在图中。
    Gf = graph(remainingBranches(:,1), remainingBranches(:,2), [], N);

    %% 3.确定每个节点所属的连通分量
    componentID = conncomp(Gf);
    numberOfComponents = max(componentID);
    % 统计各连通分量中的节点数n_k
    componentSizes = accumarray( ...
        componentID(:), 1, [numberOfComponents, 1]);

    %% 4.计算L_reach
    %
    % 故障后可达有序节点对数：
    %   sum_k n_k(n_k-1)
    reachablePairs = sum( ...
        componentSizes .* (componentSizes - 1));

    % 全部有序节点对数：
    %   N(N-1)
    totalPairs = N * (N - 1);

    % 节点可达性损失指标
    L_reach = 1 - reachablePairs / totalPairs;

    %% 5. 输出计算信息
    % fprintf('断路边编号：');
    % fprintf('%d ', Result);
    % fprintf('\n');
    % 
    % fprintf('网络节点总数 N = %d\n', N);
    % fprintf('故障后连通分量数 K = %d\n', numberOfComponents);
    % 
    % fprintf('各连通分量节点数：');
    % fprintf('%d ', componentSizes);
    % fprintf('\n');
    % 
    % fprintf('故障后可达有序节点对数 = %d\n', reachablePairs);
    % fprintf('全部有序节点对数 = %d\n', totalPairs);
    % fprintf('失去可达性的有序节点对数 = %d\n', ...
    %         totalPairs - reachablePairs);
    % 
    % fprintf('L_reach = %.6f\n', L_reach);
    % fprintf('节点可达性损失比例 = %.4f%%\n', 100 * L_reach);
end