%% 用灰狼优化算法替代 PSO-更改:加入了边界约束 (随机重置法)
%% 外层：组合位置搜索
%% 内层：down(x) 精确模型（不变）

% E_max设计时还是得大于0.05
function [Alpha_score, Alpha_pos, record] = EGWO(~,~)

%% ================== 参数设置 ==================
N   = 10;          % 种群规模
d   = 33;          % 维度（IEEE33 节点）
Max_iter = 10;     % 最大迭代次数
limit = [0.05, 0.5];   % 连续编码范围
% forbidden_idx = [1 7 12 27];   % 禁止布置储能的位置
allowed_idx = [6 8 23 24 30];

% 初始化方法Halton法-具体解释见ChatGPT回答
p = haltonset(d);
p.Skip = 100;
p.Leap = 10;
p = scramble(p,'RR2');
H = net(p,N);
x1 = limit(1) + (limit(2) - limit(1)) .* H;
%% ================== 初始化种群 ==================
x = limit(1) + (limit(2) - limit(1)) .* rand(N, d);
% 强制每个个体去掉禁止位置外仅选择 3 个非零变量
for i = 1:N
    % all_idx = 1:d;
    % allowed_idx = setdiff(all_idx, forbidden_idx);
    % idx = allowed_idx(randperm(length(allowed_idx), 3));
    num = randi(5); % 在 1、2、3、4、5 中等概率随机选择一个整数
    idx = allowed_idx(randperm(length(allowed_idx), num));
    mask = false(1, d);
    mask(idx) = true;
    x(i, ~mask) = 0;
end
%% ================== Alpha / Beta / Delta 初始化 ==================
Alpha_pos   = zeros(1, d);
Alpha_score = inf;

Beta_pos    = zeros(1, d);
Beta_score  = inf;

Delta_pos   = zeros(1, d);
Delta_score = inf;

record = zeros(Max_iter, 1);

%% ================== GWO 主循环 ==================
tic
for iter = 1:Max_iter
    %% ---------- 适应度计算 ----------
    tic
    fx = function_up(x);
    toc
    %% ---------- 更新 Alpha / Beta / Delta ----------
    for i = 1:N
        if fx(i) < Alpha_score
            Delta_score = Beta_score;
            Delta_pos   = Beta_pos;    
            Beta_score  = Alpha_score;
            Beta_pos    = Alpha_pos;  
            Alpha_score = fx(i);
            Alpha_pos   = x(i,:);
        elseif fx(i) < Beta_score
            Delta_score = Beta_score;
            Delta_pos   = Beta_pos;
            Beta_score  = fx(i);
            Beta_pos    = x(i,:);    
        elseif fx(i) < Delta_score
            Delta_score = fx(i);
            Delta_pos   = x(i,:);
        end
    end
    %% ---------- GWO 位置更新 ----------
    % a = 2 - iter * (2 / Max_iter);   % 收敛因子 a为从2到0线性递减的函数
    var = exp(-100*iter/Max_iter);
    % a = 2 * (1 - iter/Max_iter)^2; % 非线性收敛因子,使得a能更快的降为1

    % lamda = var * randn; % 模拟随机误差

    xp = 0.5*Alpha_pos + 0.3*Beta_pos + 0.2*Delta_pos;
    for i = 1:N
            % --- 更新位置 ---
        for j = 1:d
            x_t(i,j) = xp(j)-(4*rand-2)*(xp(j)-x(i,j)); % 模拟每只狼的绕圈行为;|r| < 1表示攻击猎物，|r| > 1表示寻找猎物
            if x_t(i,j) > limit(2)
                x(i,j) = x(i,j)+rand*(limit(2)-x(i,j));
            elseif x_t(i,j) < limit(1) % limit(1)改为0
                if x_t(i,j) > 0
                x(i,j) = x(i,j)+rand*(limit(1)-x(i,j));
                end
            else 
                x(i,j) = x_t(i,j);
            end
        end
    end

    %% 改进的矩阵稀疏化，使其可以避免选择分布式电源位置
    for i = 1:N
        [~, idx] = sort(abs(x(i, allowed_idx)), 'descend');
        num = randi(5);
        mask = false(1, d);
        mask(allowed_idx(idx(1:num))) = true;
        % ---------- 4. 稀疏化 ----------
        x(i, ~mask) = 0;
    end

    %% ---------- 记录与可视化 ----------
    record(iter) = Alpha_score;
    figure(5)
    subplot(1,2,1);
    bar(Alpha_pos);
    title(['迭代 ', num2str(iter), ' - 最优非零元素分布']);
    xlabel('节点编号'); ylabel('编码值');

    subplot(1,2,2);
    plot(record(1:iter), 'LineWidth', 1.5);
    title('EGWO 适应度收敛曲线');
    xlabel('迭代次数'); ylabel('最优适应度');

    pause(0.05);
end
toc
%% ================== 输出结果 ==================
disp(['最优目标值：', num2str(Alpha_score)]);
disp(['储能位置：', num2str(find(Alpha_pos ~= 0))]);
disp(['对应容量：', num2str(Alpha_pos(Alpha_pos ~= 0))]);

figure(6);
subplot(1,2,1);
bar(Alpha_pos);
title('最终最优变量分布');

subplot(1,2,2);
plot(record, 'LineWidth', 1.5);
title('EGWO 适应度收敛曲线');

end
