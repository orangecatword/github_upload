function [C, k, labels] = kmeans_lables(data, K)  % 增加输出变量 labels
% 输入变量: data代表输入数据, K代表聚类数
% 输出变量: 
%   C - 聚类中心(K×n 矩阵)
%   k - 各簇样本占比(1×K 向量)
%   labels - 每个数据点对应的类别标签(m×1 向量)

[m, n] = size(data); % 求输入数据点的个数m 特征为n
labels = zeros(m, 1); % 初始化类别标签向量 (原变量名id改为更直观的labels)

% 随机初始化聚类中心
C = zeros(K, n);
for i = 1:K
    C(i, :) = data(randi(m, 1), :); % 随机选择数据点作为初始中心
end

for iter = 1:100 % 最大迭代次数
    % 分配簇 - 计算每个点到聚类中心的距离
    for x = 1:m
        d = zeros(1, K); % 预分配距离向量
        for y = 1:K
            d(y) = norm(data(x, :) - C(y, :)); 
        end
        [~, idx] = min(d); % 找到最近的聚类中心
        labels(x) = idx;   % 记录数据点所属类别
    end
    
    % 更新聚类中心
    new_C = zeros(K, n);
    num = zeros(K, 1); % 每类样本计数
    q = 0;             % 收敛计数器
    
    for y = 1:K
        % 计算新聚类中心
        for x = 1:m
            if labels(x) == y
                new_C(y, :) = new_C(y, :) + data(x, :);
                num(y) = num(y) + 1;
            end
        end
        new_C(y, :) = new_C(y, :) / num(y);
        
        % 检查是否收敛
        if norm(new_C(y, :) - C(y, :)) < 0.1
            q = q + 1;
        end
    end
    
    % 检查所有聚类中心是否收敛
    if q == K
        break;
    else
        C = new_C;
    end
end

% 计算各类样本占比
for i = 1:K
    k(i) = length(find(labels == i)) / m;
end