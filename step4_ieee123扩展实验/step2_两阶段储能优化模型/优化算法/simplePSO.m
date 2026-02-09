%% 33维PSO算法（每个个体仅3个非零变量）-simplePSO在x值的迭代变化中还存在一些问题，后续研究改进时需要注意

function[fy_max] = simplePSO(~,~)

%% 参数设置
N = 10;                 % 种群个数
d = 33;                 % 维度
ger = 10;               % 最大迭代次数
% limit = [0.1, 0.3];   % 位置边界
% vlimit = [-0.02, 0.02];   % 速度边界
limit = [0, 1];  % 位置边界
vlimit = [-0.1, 0.1];   % 速度边界
c1 = 0.5;               % 自我学习因子
c2 = 0.5;               % 群体学习因子
w  = 0.8;               % 惯性权重
forbidden_idx = [1 7 12 27];   % 禁止布置储能的位置
%% 新增自适应权重
% 定义权重范围
w_max = 1; 
%% 初始化种群
x = limit(1) + (limit(2) - limit(1)) .*rand(N, d);
v = 0.1*(rand(N, d)-0.5);
x_cun = 0; % x的中介值


% 强制每个个体去掉禁止位置外仅选择 3 个非零变量
for i = 1:N
    all_idx = 1:d;
    allowed_idx = setdiff(all_idx, forbidden_idx);
    idx = allowed_idx(randperm(length(allowed_idx), 3));
    mask = false(1, d);
    mask(idx) = true;
    x(i, ~mask) = 0;
end

x_max = x;                     % 个体历史最优位置
y_max = zeros(1, d);           % 全局最优位置
fx_max = ones(N, 1)*inf;       % 个体历史最优适应度
fy_max = inf;                  % 全局最优适应度
record = zeros(ger,1);      % 记录适应度变化

%% PSO迭代
tic
for iter = 1:ger
    % for i = 1:N
    %     % ===== 1️⃣ 找出非零元素位置与数值 =====
    %     idx_x(i,:) = find(x(i,:) ~= 0);      % 非零元素索引
    %     vals_x(i,:) = x(i,idx_x(i,:));  
    % end
    % 此时fx应该是外层目标函数-由于论文中没有提及,则定义w_F,ESS;w_F,load;w_F,net分别为1/3
    %% 外层目标函数
    fx = down(x);   
    % 更新个体最优
    % for i = 1:N
        % ===== 1️⃣ 找出非零元素位置与数值 =====
        % idx_x(i,:) = find(x(i,:) ~= 0);      % 非零元素索引
        % vals_x(i,:) = x(i,idx_x(i,:));  
    % end
    better = fx < fx_max;
    fx_max(better) = fx(better);
    x_max(better, :) = x(better, :);
    
    % for i = 1:N
        % ===== 1️⃣ 找出非零元素位置与数值 =====
        % idx_xmax(i,:) = find(x_max(i,:) ~= 0);      % 非零元素索引
        % vals_xmax(i,:) = x_max(i,idx_xmax(i,:));  
    % end
    % 更新全局最优
    [minval, nmin] = min(fx_max);
    if minval < fy_max
        fy_max = minval;
        y_max = x_max(nmin, :);
    end
    % idx_ymax(1,:) = find(y_max(1,:) ~= 0);      % 非零元素索引
    % vals_ymax(1,:) = y_max(1,idx_ymax(1,:));

    % 速度与位置更新
    v = w*v + c1*rand*(x_max - x) + c2*rand*(repmat(y_max, N, 1) - x);
    v(v > vlimit(2)) = vlimit(2);
    v(v < vlimit(1)) = vlimit(1);
    x = x_cun + v;
    x(x > limit(2)) = limit(2);
    x(x < limit(1)) = limit(1);
    
    % 保证每个个体只有3个非零元素-矩阵稀疏化,只保证矩阵中最大的三个元素非零-这不就导致每一行的储能位置固定住了;缺陷2:迭代后期是否应该逐渐收敛 
    % for i = 1:N
    %     [~, idx] = sort(abs(x(i,:)), 'descend');
    %     mask = false(1, d);
    %     mask(idx(1:3)) = true;
    %     x(i, ~mask) = 0;
    % end

    % 改进的矩阵稀疏化，使其可以避免选择分布式电源位置
    for i = 1:N
        % ---------- 1. 强制禁止位置为 0 ----------
        x(i, forbidden_idx) = 0;
        x_cun = x;
        % ---------- 2. 找可选位置 ----------
        allowed_idx = setdiff(1:d, forbidden_idx);
        % ---------- 3. 只在允许位置中选最大的 3 个 ----------
        [~, idx] = sort(abs(x(i, allowed_idx)), 'descend');
        mask = false(1, d);
        mask(allowed_idx(idx(1:3))) = true;
        % ---------- 4. 稀疏化 ----------
        x(i, ~mask) = 0;
    end


  % 记录与可视化
    figure(5)
    record(iter) = fy_max;
    subplot(1,2,1);
    bar(y_max);
    title(['迭代', num2str(iter), ' - 全局最优非零元素分布']);
    xlabel('变量索引'); ylabel('值');
    subplot(1,2,2);
    plot(record(1:iter), 'LineWidth', 1.5);
    title('最优适应度进化过程');
    xlabel('迭代次数'); ylabel('最优适应度');
    pause(0.05);
    %% --- 方案 C：基于成功率的自适应 (进阶) ---
    % 如果这一代全局最优 fy_max 没变，说明陷入停滞，适当增大 w 鼓励跳出
    if iter > 1 && record(iter) == record(iter-1)
        w = min(w + 0.05, w_max);
    end
end

toc
%% 输出结果
disp(['最优值：', num2str(fy_max)]);
disp(['最优变量（仅前10维显示）：', num2str(y_max(1:10))]);
disp(['非零变量索引：', num2str(find(y_max~=0))]);
disp(['非零变量取值：', num2str(y_max(y_max~=0))]);

figure(6);
subplot(1,2,1); bar(y_max); title('最终最优变量分布');
subplot(1,2,2); plot(record, 'LineWidth', 1.5); title('适应度收敛曲线');
end
