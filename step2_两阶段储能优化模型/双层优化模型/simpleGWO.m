%% 33维 GWO 算法（每个个体仅3个非零变量）
%% 用灰狼优化算法替代 PSO-更改:加入了边界约束 (随机重置法)
%% 外层：组合位置搜索
%% 内层：down(x) 精确模型（不变）

function [Alpha_score, Alpha_pos, record] = simpleGWO(~,~)

%% ================== 参数设置 ==================
N   = 10;          % 种群规模
d   = 33;          % 维度（IEEE33 节点）
Max_iter = 10;     % 最大迭代次数
limit = [0,1];   % 连续编码范围
forbidden_idx = [1 7 12 27];   % 禁止布置储能的位置
%% ================== 初始化种群 ==================
x = limit(1) + (limit(2) - limit(1)) .* rand(N, d);
% 强制每个个体去掉禁止位置外仅选择 3 个非零变量
for i = 1:N
    all_idx = 1:d;
    allowed_idx = setdiff(all_idx, forbidden_idx);
    idx = allowed_idx(randperm(length(allowed_idx), 3));
    
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
    % 上层模型目标函数新增了储能成本
    fx = down(x);   
    
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
    a = 2 - iter * (2 / Max_iter);   % 收敛因子 a为从2到0线性递减的函数
    % a = 2 * (1 - iter/Max_iter)^2; % 非线性收敛因子,使得a能更快的降为1
    for i = 1:N
        % --- Alpha ---
        r1 = rand; r2 = rand; % 对应论文中将r1赋给C1 r2赋给C2
        A1 = 2*a*r1 - a; 
        C1 = 2*r2;
        D_alpha = abs(C1*Alpha_pos - x(i,:));
        X1 = Alpha_pos - A1*D_alpha;
        
        % --- Beta ---
        r1 = rand; r2 = rand;
        A2 = 2*a*r1 - a;
        C2 = 2*r2;
        D_beta = abs(C2*Beta_pos - x(i,:));
        X2 = Beta_pos - A2*D_beta;

        % --- Delta ---
        r1 = rand; r2 = rand;
        A3 = 2*a*r1 - a;
        C3 = 2*r2;
        D_delta = abs(C3*Delta_pos - x(i,:));
        X3 = Delta_pos - A3*D_delta;
            
        % --- 更新位置 ---
        x(i,:) = (X1 + X2 + X3) / 3;
    end
    
    %% ---------- 边界约束 ----------
    % x(x > limit(2)) = limit(2);
    % x(x < limit(1)) = limit(1);
    %% ---------- 边界约束 (随机重置法) ----------
    for i = 1:N
        mask_upper = x(i,:) > limit(2);
        mask_lower = x(i,:) < limit(1);
        if any(mask_upper) || any(mask_lower)
            % 越限位置重新随机初始化
            x(i, mask_upper | mask_lower) = limit(1) + (limit(2)-limit(1)) * rand;
        end
    end
    %% ---------- 稀疏化：仅保留 3 个非零 ----------
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
        % ---------- 2. 找可选位置 ----------
        allowed_idx = setdiff(1:d, forbidden_idx);
        % ---------- 3. 只在允许位置中选最大的 3 个 ----------
        [~, idx] = sort(abs(x(i, allowed_idx)), 'descend');
        % % 引入扰动：50%概率选前3，50%概率选前6里的随机3个
            % if rand > 0.5
            %     keep_idx = idx(1:3); 
            % else
            % % 在前6名里随机挑3个，增加变数
            % pool = idx(1:min(3, length(idx)));
            % keep_idx = pool(randperm(length(pool), 3));
            % end
            % mask = false(1, d);
            % mask(allowed_idx(keep_idx)) = true;
            % x(i, ~mask) = 0;
        mask = false(1, d);
        mask(allowed_idx(idx(1:3))) = true;
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
    title('GWO 适应度收敛曲线');
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
title('GWO 适应度收敛曲线');

end
