%% 离散粒子群算法（选址）-同时加入自适应权重与速度筛选
function [fy_max] = simpleBPSO(~,~)

%% 参数设置
N = 20;                 % 种群个数
d = 33;                 % 维度
ger = 20;               % 最大迭代次数
vlimit = [-2, 2];   % 速度边界
c1 = 0.5;               % 自我学习因子
c2 = 0.5;               % 群体学习因子
w  = 0.8;               % 惯性权重
%% 新增自适应权重
% 定义权重范围
w_max = 1; 

% 定义禁止位置
forbidden_idx = [1, 7, 12, 27]; 
%% ================ 初始化 =====================
%%% 【改动1】位置 x 变为 0/1 离散变量 x避开禁止位置
x = zeros(N,d);
for i = 1:N
    all_idx = 1:d;
    allowed_idx = setdiff(all_idx, forbidden_idx);
    idx = allowed_idx(randperm(length(allowed_idx), 3));
    x(i,idx) = 1;
end

%%% 【改动2】速度仍为连续变量
v = randn(N,d);

x_max  = x;                 % 个体历史最优
y_max  = zeros(1,d);        % 全局最优
fx_max = ones(N, 1)*inf;       % 个体历史最优适应度
fy_max = inf;
record = zeros(ger,1);

%% ================= BPSO 迭代 ==================
tic
for iter = 1:ger
    %% --------- 目标函数 ----------
    % 离散位置进入外层目标函数
    fx = down_position(x);
    % 更新个体最优
    better = fx < fx_max;
    fx_max(better) = fx(better);
    x_max(better,:) = x(better,:);

    % 更新全局最优
    [minval, idx] = min(fx_max);
    if minval < fy_max
        fy_max = minval;
        y_max  = x_max(idx,:);
    end

    %% --------- 更新速度 ----------
    v = w*v + c1*rand*(x_max - x) + c2*rand*(repmat(y_max, N, 1) - x);
    v(v > vlimit(2)) = vlimit(2);
    v(v < vlimit(1)) = vlimit(1);

    %% --------- 离散位置更新 ----------
    %%% 【改动3】Sigmoid 映射
    prob = 1 ./ (1 + exp(-v));

    %%% 【改动4】概率 → 0/1
    x = rand(N,d) < prob;

    %% --------- 强制 3 个 1 ----------
for i = 1:N
    % --- 新增逻辑：首先强制将禁止位置设为 0，确保万无一失 ---
    x(i, forbidden_idx) = 0;
    current_sum = sum(x(i,:));
    if current_sum > 3
        % 1. 在当前为 1 的位置中寻找
        idx = find(x(i,:) == 1);
        % 这里不需要额外排除，因为最开始已经把禁止位清零了，idx 里不会包含它们
        [~, order] = sort(v(i, idx), 'descend'); 
        keep = idx(order(1:3));
        x(i, :) = 0;
        x(i, keep) = 1;
    elseif current_sum < 3
        % 2. 在当前为 0 的位置中寻找可添加的位置
        % % --- 关键改动：依据速度 v 的大小进行排序补齐 ---
        zero_idx = find(x(i,:) == 0);
        % --- 关键改动：从可选的零位置中剔除禁止的位置 ---
        allowed_zero_idx = setdiff(zero_idx, forbidden_idx);
        % 我们提取这些可选 0 位置对应的速度值，按降序排列
        [~, order] = sort(v(i, allowed_zero_idx), 'descend');  
        % 计算需要补齐的数量
        add_count = 3 - current_sum; 
        % 挑选速度最大的前 add_count 个位置
        add = allowed_zero_idx(order(1:add_count));
        % 设为 1
        x(i, add) = 1;
    end
end

    %% --------- 记录与可视化 ----------
    figure(5)
    record(iter) = fy_max;
    subplot(1,2,1);
    bar(y_max);
    title(['迭代',num2str(iter),' - 选址分布']);
    xlabel('位置编号'); ylabel('是否选中');

    subplot(1,2,2);
    plot(record(1:iter),'LineWidth',1.5);
    title('适应度收敛');
    xlabel('迭代次数'); ylabel('最优值');
    pause(0.05);
    %% --- 方案 C：基于成功率的自适应 (进阶) ---
    % 如果这一代全局最优 fy_max 没变，说明陷入停滞，适当增大 w 鼓励跳出
    if iter > 1 && record(iter) == record(iter-1)
        w = min(w + 0.05, w_max);
    end
end
toc
%% ================= 输出结果 ==================
disp(['最优值：', num2str(fy_max)]);
disp(['选中位置索引：', num2str(find(y_max==1))]);

figure(6);
subplot(1,2,1); bar(y_max); title('simpleBPSO最终选址结果');
subplot(1,2,2); plot(record,'LineWidth',1.5); title('收敛曲线');

end
