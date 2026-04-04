% 采用kmeans聚类,将故障情况聚为三类-确定好最佳聚类数是四类

% r x 减5%时 
% 第一类：20 22 12 11 28 聚类中心值4.8519     比例：0.5763
% 第二类：21 11 15 32 25 聚类中心值108.3648   比例：0.1911
% 第三类：31 17 29 3 6   聚类中心值284.1119   比例：0.1709
% 第四类：21 32 22 26 28 聚类中心值643.6753   比例：0.0617

% r x 加5%时 
% 第一类：20 22 12 11 28 聚类中心值5.7846     比例：0.5657
% 第二类：21 11 15 32 25 聚类中心值108.4062   比例：0.2089
% 第三类：31 17 29 3 6   聚类中心值284.0483   比例：0.1634
% 第四类：21 32 22 26 28 聚类中心值643.6867   比例：0.0620

% r x 减10%时
% 第一类：20 22 12 11 28 聚类中心值4.5920     比例：0.5759
% 第二类：21 11 15 32 25 聚类中心值108.3454   比例：0.1937
% 第三类：31 17 29 3 6   聚类中心值284.1006   比例：0.1774
% 第四类：21 32 22 26 28 聚类中心值643.6698   比例：0.0530

% r x 加10%时 [0.200204290091931,0.586312563840654,0.053115423901941,0.160367722165475]
% 第一类：20 22 12 11 28 聚类中心值6.4274     比例：0.5863
% 第二类：21 11 15 32 25 聚类中心值108.4295   比例：0.2002
% 第三类：31 17 29 3 6   聚类中心值284.1624   比例：0.1604
% 第四类：21 32 22 26 28 聚类中心值643.6922   比例：0.0531

% m=1;
% for i = 1:10000
%     if values(i)<0.1
%         value(m) = values(i);
%         m=m+1;
%     end
% 
% end
% k = 4;
% [cluster_center,proportion,cluster_idx] = kmeans_lables(value', k);

% clc
% clear
% warning off
% 

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
    if S(i)>1000
        fault_line(m,:) = Result(i, :);
        m = m+1;
    end
end
fault_line = value(fault_line);
%% Step 2. 提取非零且 <=1000 的元素及其位置-不考虑丢弃负荷的情况
[idx_row,idx_col,values] = find(S); % 找出非零元素,idx_col始终为1-因为只有1列

% 过滤掉大于1000的值（不参与聚类,故障特征值过大可能情况是无解的）

valid_idx = values <= 1000;
idx_row = idx_row(valid_idx);
values  = values(valid_idx);

fprintf('参与聚类的有效元素数: %d（已排除 >1000 的元素）\n', length(values));

% x = zeros(8,1); 
%% Step 3. 执行K-means聚类（3类）
% for k = 2:8
k = 4;
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