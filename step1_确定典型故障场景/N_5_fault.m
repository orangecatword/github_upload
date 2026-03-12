% 采用kmeans聚类,将故障情况聚为三类-确定好最佳聚类数是四类

% 2025.2.7版本
% 第一类：1	15 27 4 13   聚类中心值11.2354  比例：0.5576
% 第二类：16 25 10 32 11 聚类中心值105.2897 比例：0.2130  改为：21	11 15 32 25
% 第三类：15 4 30 10 9   聚类中心值270.3424 比例：0.1733  改为：28 15 30 25 3
% 第四类：21 32	22 26 28 聚类中心值645.2953 比例：0.0561  改为：22 16 31 3 28

% 2026.3.4版本 0.5586	0.1713（3）	0.0561（4）	0.2141（2）
% 第一类：20 22	12 11 28 聚类中心值11.1837    比例：0.5586
% 第二类：21 11 15 32 25 聚类中心值105.7488   比例：0.2141
% 第三类：31 17	29 3 6   聚类中心值270.9940   比例：0.1713
% 第四类：21 32 22 26 28 聚类中心值645.2905   比例：0.0561


<<<<<<< HEAD

=======
% 2025.2.2版本
% 第一类：5	28 1 10 7    聚类中心值0.13 比例：0.5564
% 第二类：31 20	5 28 14  聚类中心值1.06 比例：0.2162
% 第三类：11 16	9 5	30   聚类中心值2.68 比例：0.1714  
% 第四类：29 32	14 23 27 聚类中心值6.28 比例：0.0560

% 2025.2.7版本
% 第一类：1	15 27 4 13   聚类中心值11.2354  比例：0.5576
% 第二类：16 25 10 32 11 聚类中心值105.2897 比例：0.2130  改为：21 11 15 32 25
% 第三类：15 4 30 10 9   聚类中心值270.3424 比例：0.1733  改为：28 15 30 25 3
% 第四类：21 32	22 26 28 聚类中心值645.2953 比例：0.0561  改为：22 16 31 3 28


>>>>>>> 73a3c4801881cb66c48104448e8046700252e1e6
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
<<<<<<< HEAD


=======
>>>>>>> 73a3c4801881cb66c48104448e8046700252e1e6
% clc
% clear
% warning off
% 
<<<<<<< HEAD
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
=======
% % 抽样100个值和1000个值差距还是很大的
% Result = zeros(1000,5); % 存储N-5断路情况
% S = zeros(1000,1);      % 定义不同场景下的故障特征量
% fault_line = zeros(1000,5); % 存放故障线路
% %% Step 1. 导入或定义矩阵 
% m = 1;
% for i = 1:1000
%     % 从 1~32 中随机选择 5 个不重复数字
%     randNums = randperm(32, 5);
%     % 保存到矩阵中
%     Result(i, :) = randNums;
%     S(i) = function_fault_5(Result(i, :));
%     if S(i)>10
%         fault_line(m,:) = Result(i, :);
%         m = m+1;
%     end
% end
% fault_line = value(fault_line);
%% Step 2. 提取非零且 <=10 的元素及其位置-不考虑丢弃负荷的情况
[idx_row,idx_col,values] = find(S); % 找出非零元素,idx_col始终为1-因为只有1列

% 过滤掉大于10的值（不参与聚类,故障特征值过大可能情况是无解的）
>>>>>>> 73a3c4801881cb66c48104448e8046700252e1e6
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