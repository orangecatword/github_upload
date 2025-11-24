% 采用kmeans聚类,将故障情况聚为三类
clc
clear
warning off

Result = zeros(1000,5); % 存储N-5断路情况
S = zeros(1000,1);      % 定义不同场景下的故障特征量

%% Step 1. 导入或定义矩阵 
for i = 1:1000
    % 从 1~32 中随机选择 5 个不重复数字
    randNums = randperm(32, 5);
    % 保存到矩阵中
    Result(i, :) = randNums;
    S(i) = function_fault_5(Result(i, :));
end
%% Step 2. 提取非零且 <=70 的元素及其位置-不考虑丢弃负荷的情况
[idx_row,idx_col,values] = find(S); % 找出非零元素

% 过滤掉大于70的值（不参与聚类,且这个值无解）
valid_idx = values <= 20000;
idx_row = idx_row(valid_idx);
values  = values(valid_idx);

fprintf('参与聚类的有效元素数: %d（已排除 >70 的元素）\n', length(values));

%% Step 3. 执行K-means聚类（3类）
k = 3;
[cluster_center,proportion,cluster_idx] = kmeans_lables(values, k);

%% Step 4. 输出结果
fprintf('\n===== K-means 聚类结果 =====\n');
for i = 1:k
    idx = find(cluster_idx == i);      % 该类样本索引
    vals = values(idx);                % 该类样本值

    % 找到最接近中心的样本
    [~, nearest] = min(abs(vals - cluster_center(i)));
    nearest_val = vals(nearest);
    nearest_pos = [idx_row(idx(nearest)), idx_col(idx(nearest))];

    % 输出信息
    fprintf('\n第 %d 类:\n', i);
    fprintf('  聚类中心值 (均值): %.4f\n', cluster_center(i));
    fprintf('  最接近中心的原始值: %.4f (行=%d, 列=%d)\n', ...
        nearest_val, nearest_pos(1), nearest_pos(2));
    fprintf('  具体故障情况:%d\n',Result(nearest_pos(1),:))
end

