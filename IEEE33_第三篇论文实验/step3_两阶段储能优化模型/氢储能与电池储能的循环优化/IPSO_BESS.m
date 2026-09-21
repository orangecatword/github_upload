%% IPSO算法
function[idx,vals,f1] = IPSO_BESS(idx_H,vals_H)
%% 参数设置
N = 10;                 % 种群个数
d = 33;                 % 维度
ger = 10;               % 最大迭代次数
limit = [0.05, 0.5];  % 位置边界
vlimit = [-0.05, 0.05];   % 速度边界
c1 = 0.5;               % 自我学习因子
c2 = 0.5;               % 群体学习因子
w  = 0.8;               % 惯性权重

candidate_nodes = [6 8 23 24 30]; % 5个节点重要性最高的点
allowed_idx = setdiff(candidate_nodes,idx_H); % 节点选择范围
%% 新增自适应权重
% 定义权重范围
w_max = 1; 
%% 初始化种群
% 初始化方法Halton法-具体解释见ChatGPT回答
p = haltonset(d);
p.Skip = 100;
p.Leap = 10;
p = scramble(p,'RR2');
H = net(p,N);
x1 = limit(1) + (limit(2) - limit(1)) .* H;

x = limit(1) + (limit(2) - limit(1)) .*rand(N, d);
v = 0.1*(rand(N, d)-0.5);

for i = 1:N
    num = randi(size(allowed_idx,2)); % 在 1、2、3、4、5 中等概率随机选择一个整数
    idx = allowed_idx(randperm(length(allowed_idx), num));
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
    %% 外层目标函数
    fx = function_up_BESS(x,idx_H,vals_H);   
    % 更新个体最优
    better = fx < fx_max;
    fx_max(better) = fx(better);
    x_max(better, :) = x(better, :);
    % 更新全局最优
    [minval, nmin] = min(fx_max);
    if minval < fy_max
        fy_max = minval;
        y_max = x_max(nmin, :);
    end

    % 速度与位置更新
    v = w*v + c1*rand*(x_max - x) + c2*rand*(repmat(y_max, N, 1) - x);
    v(v > vlimit(2)) = vlimit(2);
    v(v < vlimit(1)) = vlimit(1);
    x = x + v;
    x(x > limit(2)) = limit(2);
    x(x < limit(1)) = limit(1);
    

    % 改进的矩阵稀疏化，使其可以避免选择分布式电源位置
    for i = 1:N
        [~, idx] = sort(abs(x(i, allowed_idx)), 'descend');
        num = randi(size(allowed_idx,2));
        mask = false(1, d);
        mask(allowed_idx(idx(1:num))) = true;
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
    title('PSO最优适应度进化过程');
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
disp(['非零变量索引：', num2str(find(y_max~=0))]);
disp(['非零变量取值：', num2str(y_max(y_max~=0))]);

% figure(6);
% subplot(1,2,1); bar(y_max); title('最终最优变量分布');
% subplot(1,2,2); plot(record, 'LineWidth', 1.5); title('IPSO适应度收敛曲线');
idx  = find(y_max~=0);
vals = y_max(y_max~=0);
f1   = fy_max;
end
