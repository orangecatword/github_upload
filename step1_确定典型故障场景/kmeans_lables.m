function [C, k, labels] = kmeans_lables(data, K)  % 增加输出变量 labels
% 输入变量: data代表输入数据, K代表聚类数
% 输出变量: 
%   C - 聚类中心(K×n 矩阵)
%   k - 各簇样本占比(1×K 向量)
%   labels - 每个数据点对应的类别标签(m×1 向量)
    [m, n] = size(data);
    max_iter = 100;
    tol = 1e-6; % 更加通用的收敛阈值
    
    % 1. 随机初始化 (避免重复选择同一个点)
    rand_idx = randperm(m, K);
    C = data(rand_idx, :);
    labels = zeros(m, 1);
    
    for iter = 1:max_iter
        old_C = C;
        
        % 2. 分配簇 (向量化计算，提高速度)
        % 计算每个点到每个中心的欧氏距离平方
        dist_sq = sum(data.^2, 2) - 2 * data * C' + sum(C.^2, 2)';
        [~, labels] = min(dist_sq, [], 2);
        
        % 3. 更新聚类中心
        for y = 1:K
            relevant_points = data(labels == y, :);
            if ~isempty(relevant_points)
                C(y, :) = mean(relevant_points, 1);
            else
                % 处理空簇：如果某簇为空，重新随机选一个点作为中心
                C(y, :) = data(randi(m), :);
            end
        end
        
        % 4. 检查收敛 (判断总位移量)
        if norm(C - old_C, 'fro') < tol
            fprintf('在第 %d 次迭代时收敛。\n', iter);
            break;
        end
    end
    
    % 计算各类样本占比
    k = zeros(1, K);
    for i = 1:K
        k(i) = sum(labels == i) / m;
    end
end