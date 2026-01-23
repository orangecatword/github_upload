% 采用kmeans聚类,将故障情况聚为三类-确定好最佳聚类数是四类
%% 断路情况有问题 需要修改
% 2025.12.31版本
% 第一类：5 1 3 32 28 聚类中心值0.0087 比例：1029/1958
% 第二类：20 29 11 9 15 聚类中心值0.4491 比例：598/1958
% 第三类：5 23 17 4 24 聚类中心值1.9204 比例：331/1958

% 2025.1.15版本
% 第一类：10 24 1 2 19 聚类中心值0.7628 比例：0.243
% 第二类：3 32 6 28 23 聚类中心值6.0536 比例：0.055
% 第三类：5 4 29 32 30 聚类中心值2.4731 比例：0.217
% 第四类：10 23 13 25 27 聚类中心值0.0562 比例：0.485
clc
clear
warning off

% 抽样100个值和1000个值差距还是很大的
Result = zeros(1000,5); % 存储N-5断路情况
S = zeros(1000,1);      % 定义不同场景下的故障特征量
fault_line = zeros(1000,5); % 存放故障线路
%% Step 1. 导入或定义矩阵 
m = 1;
for i = 1:1000
    % 从 1~32 中随机选择 5 个不重复数字
    randNums = randperm(32, 5);
    % 保存到矩阵中
    Result(i, :) = randNums;
    S(i) = function_fault_5(Result(i, :));
    if S(i)>10
        fault_line(m,:) = Result(i, :);
        m = m+1;
    end
end
fault_line = value(fault_line);
%% Step 2. 提取非零且 <=10 的元素及其位置-不考虑丢弃负荷的情况
[idx_row,idx_col,values] = find(S); % 找出非零元素,idx_col始终为1-因为只有1列

% 过滤掉大于10的值（不参与聚类,故障特征值过大可能情况是无解的）
valid_idx = values <= 10;
idx_row = idx_row(valid_idx);
values  = values(valid_idx);

fprintf('参与聚类的有效元素数: %d（已排除 >10 的元素）\n', length(values));

x = zeros(8,1); 
%% Step 3. 执行K-means聚类（3类）
for k = 2:8
% k = 4;
[cluster_center,proportion,cluster_idx] = kmeans_lables(values, k);

clut = 0; % 用以确定最佳聚类数
%% Step 4. 输出结果
fprintf('\n===== K-means 聚类结果 =====\n');
for i = 1:k
    u = cluster_center(i);             % 每个聚类的均值
    idx = find(cluster_idx == i);      % 该类样本索引
    num = size(idx,1);                 % 记录每一类的样本数
    vals = values(idx);                % 该类样本值
    clut = sum(abs(vals - u*ones(num,1)))^2 + clut;

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
x(k) = clut;
end

iterations = 2:8;
figure; % 新建图形窗口
plot(iterations, x(2:8,:), 'b-', 'LineWidth', 1.5); 
% 添加标题和坐标轴标签
title('不同聚类数下的距离内和');
xlabel('聚类数');
ylabel('距离内和');