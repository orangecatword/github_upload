% 聚为3类 
% 第一类：71 92  聚类中心值174.1538 比例：0.6660
% 第二类：52 110 聚类中心值457.2970 比例：0.2220 
% 第三类：43 51  聚类中心值822.6367 比例：0.1120 

% 2026.3.11
% 第一类：18 33  聚类中心值591.4690 比例：0.8624
% 第二类：115 43 聚类中心值620.5553 比例：0.0275 
% 第三类：4 51   聚类中心值824.4931 比例：0.1101 

% clc
% clear
% warning off

Result = zeros(1000,2); % 存储N-5断路情况
S = zeros(1000,1);      % 定义不同场景下的故障特征量
fault_S = zeros(1000,1); % 存放故障线路对应的故障特征量
fault_line = zeros(1000,2); % 存放故障线路
%% Step 1. 导入或定义矩阵 
m = 1;
for i = 1:1000
    % 从 1~122 中随机选择 5 个不重复数字
    randNums = randperm(122, 2);
    % 保存到矩阵中
    Result(i, :) = randNums;
    S(i) = function_fault_2(Result(i, :));
    if S(i)>1000
        fault_line(m,:) = Result(i, :);
        fault_S(m) = S(i);
        m = m+1;
    end
end

%% Step 2. 提取非零且 <=1000 的元素及其位置-不考虑丢弃负荷的情况
[idx_row,idx_col,values] = find(S); % 找出非零元素

valid_idx = values <= 1000;
idx_row = idx_row(valid_idx);
idx_col = idx_col(valid_idx);
values  = values(valid_idx);

fprintf('参与聚类的有效元素数: %d（已排除 >1000 的元素）\n', length(values));

% x = zeros(8,1); 
%% Step 3. 执行K-means聚类（3类）
% for k = 2:8
k = 3;
[cluster_center,proportion,cluster_idx] = kmeans_lables(values, k);
% clut = 0; % 用以确定最佳聚类数
%% Step 4. 输出结果
fprintf('\n===== K-means 聚类结果 =====\n');
for i = 1:k
    u = cluster_center(i);             % 每个聚类的均值
    idx = find(cluster_idx == i);      % 该类样本索引
    num = size(idx,1);                 % 记录每一类的样本数
    vals = values(idx);                % 该类样本值
    % clut = sum(abs(vals - u*ones(num,1)))^2 + clut;
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
% x(k) = clut;
% end

% iterations = 2:8;
% figure; % 新建图形窗口
% plot(iterations, x(2:8,:)/(10*1000), 'b-', 'LineWidth', 1.5); 
% % 添加标题和坐标轴标签
% title('不同聚类数下的距离内和');
% xlabel('聚类数');
% ylabel('距离内和');

