%% 故障边数选择
%% 故障概率也需要重新计算
% 选择故障概率p=0.05, 同时明确"极端事件造成了至少一条线路故障"-这一概率用来选择主要故障场景
% 概率前十的故障断路数为N-(6,5,7,4,8,3,9,10,2,11) 它们的概率和超过了95%,
% 以它们的概率和为分母,它们的所占概率分别为:17.03,16.59,14.85,13.36,11.24,8.53,7.49,4.46,4.05,2.39%

% 存储断路情况
Result1 = zeros(1000,6); 
Result2 = zeros(1000,5); 
Result3 = zeros(1000,7); 
Result4 = zeros(1000,4); 
Result5 = zeros(1000,8);
Result6 = zeros(1000,3); 
Result7 = zeros(1000,9); 
Result8 = zeros(1000,10); 
Result9 = zeros(1000,2); 
Result10 = zeros(1000,11);

% 定义故障特征值
S1 = zeros(1000,2); 
S2 = zeros(1000,2);
S3 = zeros(1000,2); 
S4 = zeros(1000,2); 
S5 = zeros(1000,2); 
S6 = zeros(1000,2); 
S7 = zeros(1000,2);
S8 = zeros(1000,2); 
S9 = zeros(1000,2); 
S10 = zeros(1000,2); 

% Step 1. 导入或定义矩阵 
m = 1;
for i = 1:200
    % 通过1000组数据的随机抽样,得到故障线路
    % 在选择故障线路时避开变电站故障情况
    randNums1 = randperm(121, 6) + 1;
    randNums2 = randperm(121, 5) + 1;
    randNums3 = randperm(121, 7) + 1;
    randNums4 = randperm(121, 4) + 1;
    randNums5 = randperm(121, 8) + 1;
    randNums6 = randperm(121, 3) + 1;
    randNums7 = randperm(121, 9) + 1;
    randNums8 = randperm(121,10) + 1;
    randNums9 = randperm(121, 2) + 1;
    randNums10 = randperm(121,11) + 1;

    Result1(i, :) = randNums1;
    Result2(i, :) = randNums2;
    Result3(i, :) = randNums3;
    Result4(i, :) = randNums4;
    Result5(i, :) = randNums5;
    Result6(i, :) = randNums6;
    Result7(i, :) = randNums7;
    Result8(i, :) = randNums8;
    Result9(i, :) = randNums9;
    Result10(i, :) = randNums10;

    [S1(i,1),S1(i,2)] = function_fault(Result1(i, :));
    [S2(i,1),S2(i,2)] = function_fault(Result2(i, :));
    [S3(i,1),S3(i,2)] = function_fault(Result3(i, :));
    [S4(i,1),S4(i,2)] = function_fault(Result4(i, :));
    [S5(i,1),S5(i,2)] = function_fault(Result5(i, :));
    [S6(i,1),S6(i,2)] = function_fault(Result6(i, :));
    [S7(i,1),S7(i,2)] = function_fault(Result7(i, :));
    [S8(i,1),S8(i,2)] = function_fault(Result8(i, :));
    [S9(i,1),S9(i,2)] = function_fault(Result9(i, :));
    [S10(i,1),S10(i,2)] = function_fault(Result10(i, :));
end


%% Step 2. 找到最佳聚类数

% x = zeros(20,1); % 数据存储空间
% 
% for k = 2:20
%     [cluster_center,proportion,cluster_idx] = kmeans_lables(S1, k);
%     clut = 0; % 用以确定最佳聚类数
% 
%     fprintf('\n===== K-means 聚类结果 =====\n');
%     for i = 1:k
%         u = cluster_center(i,:);             % 每个聚类的均值
%         idx = find(cluster_idx == i);      % 该类样本索引
%         num = size(idx,1);                 % 记录每一类的样本数
%         vals = S1(idx,:);                % 该类样本值
%         clut = sum(sum((vals-u).^2))+clut; % 计算类内平方和
%         distance = sqrt(sum((vals-u).^2,2)); % 计算每个样本到中心距离
%         [~,nearest]=min(distance);         % 找距离中心最近的样本
%         nearest_val = vals(nearest,:);
%         nearest_pos = idx(nearest);   
%     end
%     x(k) = clut;
% end
% 
% iterations = 2:20;
%     figure; % 新建图形窗口
%     plot(iterations, x(2:20,:), 'b-', 'LineWidth', 1.10); 
%     % 添加标题和坐标轴标签
%     title('不同聚类数下的距离内和');
%     xlabel('聚类数');
%     ylabel('距离内和');


%% 基于最佳聚类数进行聚类
% ChatGPT分析:从S1到S5 最佳聚类数分别为6 5 5 5 6
% 先统一为5 5 5 5 5 (已完成前两列的数据验证 13-7-17-4-3仍存在环路)
% N-1: 11; 23; 6; 3; 20                                                             对应概率: 0.298	0.258	0.224	0.128	0.092
% N-2: 29-2; 24-14; 22-7; 12-17; 18-26                                              对应概率: 0.089	0.126	0.231	0.280	0.274
% N-3: 30-28-25; 23-14-3; 28-5-20; 30-10-19; 19-30-23                               对应概率: 0.227	0.154	0.289	0.200	0.130
% N-4: 21-3-22-9; 27-16-19-3; 25-27-22-26; 5-4-21-2; 8-30-16-9                      对应概率: 0.181	0.318	0.150	0.055	0.296
% N-5: 13-7-17-4-3; 15-29-26-12-21; 20-21-4-8-23; 24-26-28-27-18; 13-23-19-14-26    对应概率: 0.267	0.200	0.338	0.088	0.107

[cluster_center,proportion,cluster_idx] = kmeans_lables(S1, 5);
    fprintf('\n===== K-means 聚类结果 =====\n');
    for i = 1:5
        u = cluster_center(i,:);             % 每个聚类的均值
        idx = find(cluster_idx == i);      % 该类样本索引
        num = size(idx,1);                 % 记录每一类的样本数
        vals = S1(idx,:);                % 该类样本值
        distance = sqrt(sum((vals-u).^2,2)); % 计算每个样本到中心距离
        [~,nearest]=min(distance);         % 找距离中心最近的样本
        nearest_val = vals(nearest,:);
        nearest_pos = idx(nearest);

        % 输出信息
        fprintf('\n第 %d 类:\n', i);
        fprintf('  聚类中心坐标: (%.4f, %.4f)\n',...
            cluster_center(i,1),...
            cluster_center(i,2));
        fprintf('  最接近中心的原始值: %.4f\n', nearest_val(:));
        fprintf('  最近中心的原始坐标:(%.4f, %.4f)\n',...
            nearest_val(1),...
            nearest_val(2));
        fprintf('  具体故障情况:%d\n',Result1(nearest_pos,:))
        fprintf('  故障情况所占概率:%.3f\n',proportion(i));
    end

[cluster_center,proportion,cluster_idx] = kmeans_lables(S2, 5);
    fprintf('\n===== K-means 聚类结果 =====\n');
    for i = 1:5
        u = cluster_center(i,:);             % 每个聚类的均值
        idx = find(cluster_idx == i);      % 该类样本索引
        num = size(idx,1);                 % 记录每一类的样本数
        vals = S2(idx,:);                % 该类样本值
        distance = sqrt(sum((vals-u).^2,2)); % 计算每个样本到中心距离
        [~,nearest]=min(distance);         % 找距离中心最近的样本
        nearest_val = vals(nearest,:);
        nearest_pos = idx(nearest);  
        
        % 输出信息
        fprintf('\n第 %d 类:\n', i);
        fprintf('  聚类中心坐标: (%.4f, %.4f)\n',...
            cluster_center(i,1),...
            cluster_center(i,2));
        fprintf('  最接近中心的原始值: %.4f\n', nearest_val(:));
        fprintf('  最近中心的原始坐标:(%.4f, %.4f)\n',...
            nearest_val(1),...
            nearest_val(2));
        fprintf('  具体故障情况:%d\n',Result2(nearest_pos,:))
        fprintf('  故障情况所占概率:%.3f\n',proportion(i));
    end

[cluster_center,proportion,cluster_idx] = kmeans_lables(S3, 5);
    fprintf('\n===== K-means 聚类结果 =====\n');
    for i = 1:5
        u = cluster_center(i,:);             % 每个聚类的均值
        idx = find(cluster_idx == i);      % 该类样本索引
        num = size(idx,1);                 % 记录每一类的样本数
        vals = S3(idx,:);                % 该类样本值
        distance = sqrt(sum((vals-u).^2,2)); % 计算每个样本到中心距离
        [~,nearest]=min(distance);         % 找距离中心最近的样本
        nearest_val = vals(nearest,:);
        nearest_pos = idx(nearest);   

        % 输出信息
        fprintf('\n第 %d 类:\n', i);
        fprintf('  聚类中心坐标: (%.4f, %.4f)\n',...
            cluster_center(i,1),...
            cluster_center(i,2));
        fprintf('  最接近中心的原始值: %.4f\n', nearest_val(:));
        fprintf('  最近中心的原始坐标:(%.4f, %.4f)\n',...
            nearest_val(1),...
            nearest_val(2));
        fprintf('  具体故障情况:%d\n',Result3(nearest_pos,:))
        fprintf('  故障情况所占概率:%.3f\n',proportion(i));
    end

[cluster_center,proportion,cluster_idx] = kmeans_lables(S4, 5);
    fprintf('\n===== K-means 聚类结果 =====\n');
    for i = 1:5
        u = cluster_center(i,:);             % 每个聚类的均值
        idx = find(cluster_idx == i);      % 该类样本索引
        num = size(idx,1);                 % 记录每一类的样本数
        vals = S4(idx,:);                % 该类样本值
        distance = sqrt(sum((vals-u).^2,2)); % 计算每个样本到中心距离
        [~,nearest]=min(distance);         % 找距离中心最近的样本
        nearest_val = vals(nearest,:);
        nearest_pos = idx(nearest);   

        % 输出信息
        fprintf('\n第 %d 类:\n', i);
        fprintf('  聚类中心坐标: (%.4f, %.4f)\n',...
            cluster_center(i,1),...
            cluster_center(i,2));
        fprintf('  最接近中心的原始值: %.4f\n', nearest_val(:));
        fprintf('  最近中心的原始坐标:(%.4f, %.4f)\n',...
            nearest_val(1),...
            nearest_val(2));
        fprintf('  具体故障情况:%d\n',Result4(nearest_pos,:))
        fprintf('  故障情况所占概率:%.3f\n',proportion(i));
    end

[cluster_center,proportion,cluster_idx] = kmeans_lables(S5, 5);
    fprintf('\n===== K-means 聚类结果 =====\n');
    for i = 1:5
        u = cluster_center(i,:);             % 每个聚类的均值
        idx = find(cluster_idx == i);      % 该类样本索引
        num = size(idx,1);                 % 记录每一类的样本数
        vals = S5(idx,:);                % 该类样本值
        distance = sqrt(sum((vals-u).^2,2)); % 计算每个样本到中心距离
        [~,nearest]=min(distance);         % 找距离中心最近的样本
        nearest_val = vals(nearest,:);
        nearest_pos = idx(nearest);   

        % 输出信息
        fprintf('\n第 %d 类:\n', i);
        fprintf('  聚类中心坐标: (%.4f, %.4f)\n',...
            cluster_center(i,1),...
            cluster_center(i,2));
        fprintf('  最接近中心的原始值: %.4f\n', nearest_val(:));
        fprintf('  最近中心的原始坐标:(%.4f, %.4f)\n',...
            nearest_val(1),...
            nearest_val(2));
        fprintf('  具体故障情况:%d\n',Result5(nearest_pos,:))
        fprintf('  故障情况所占概率:%.3f\n',proportion(i));
    end

