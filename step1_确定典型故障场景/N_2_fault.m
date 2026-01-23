% 采用kmeans聚类,将故障情况聚为三类
% 论文版:聚为3类 断开了支路5-28;15-30;14-31

% 2025.12.30更新
% 第一类：7 18 聚类中心值0.0042 比例：1-18/496-8/496
% 第二类：9 11 聚类中心值0.4225 比例：18/496
% 第三类：22 24 聚类中心值2.0395 比例：8/496

% 2025.1.16版本 0.0322580645161290	0.0282258064516129	0.0423387096774194	0.897177419354839
% 第一类：3 22 聚类中心值0.0081 比例：0.924
% 第二类：19 20 聚类中心值0.3589 比例：0.034
% 第三类：15 32 聚类中心值1.0041 比例：0.020
% 第四类：16 30 聚类中心值2.3329 比例：0.022
clc
clear
warning off

%% Step 1. 导入或定义矩阵 
S = zeros(32,32); % 定义不同场景下的故障特征量
tic
for i = 1:32  % 改为从1开始,可以通过分布式电源进行孤岛运行
    for j = i+1:32
    S(i,j) = function_fault_2(i,j);
    end
end
toc

%% Step 2. 提取非零且 <=10 的元素及其位置-不考虑丢弃负荷的情况
[idx_row, idx_col, values] = find(S); % 找出非零元素

% 过滤掉大于70的值（不参与聚类,故障特征值过大可能情况是无解的）
valid_idx = values <= 10;
idx_row = idx_row(valid_idx);
idx_col = idx_col(valid_idx);
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
