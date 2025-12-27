%% 离散粒子群算法
function [fy_max] = simpleBPSO(~,~)

%% ================= 参数设置 =================
N = 20;              % 种群规模
d = 33;              % 维度（33个候选位置）
ger = 20;            % 迭代次数

w  = 0.8;            % 惯性权重
c1 = 0.5;            % 个体学习因子
c2 = 0.5;            % 群体学习因子

vlimit = [-4, 4];    % 速度边界（离散PSO中一般取较大）

%% ================ 初始化 =====================
%%% 【改动1】位置 x 变为 0/1 离散变量
x = zeros(N,d);
for i = 1:N
    idx = randperm(d,3);   % 每个粒子只选 3 个位置
    x(i,idx) = 1;
end

%%% 【改动2】速度仍为连续变量
v = randn(N,d);

x_max  = x;                 % 个体历史最优
fx_max = inf(N,1);          % 个体最优适应度

y_max  = zeros(1,d);        % 全局最优
fy_max = inf;

record = zeros(ger,1);

%% ================= PSO 迭代 ==================
for iter = 1:ger
    %% --------- 目标函数 ----------
    % 离散位置进入外层目标函数
    % fx = 2/3*down_position(x) + 1/3*100*sum(x,2); % 都是安装三个储能，所以无需考虑储能成本
    fx = down_position(x);

    %% --------- 更新个体最优 ----------
    better = fx < fx_max;
    fx_max(better) = fx(better);
    x_max(better,:) = x(better,:);

    %% --------- 更新全局最优 ----------
    [minval, idx] = min(fx_max);
    if minval < fy_max
        fy_max = minval;
        y_max  = x_max(idx,:);
    end

    %% --------- 更新速度 ----------
    v = w*v ...
        + c1*rand(N,d).*(x_max - x) ...
        + c2*rand(N,d).*(repmat(y_max,N,1) - x);

    v = max(min(v, vlimit(2)), vlimit(1)); % 不允许速度超过边界

    %% --------- 离散位置更新 ----------
    %%% 【改动3】Sigmoid 映射
    prob = 1 ./ (1 + exp(-v));

    %%% 【改动4】概率 → 0/1
    x = rand(N,d) < prob;

    %% --------- 强制 3 个 1 ----------
    %%% 【改动5】每个粒子只保留 3 个 1
    for i = 1:N
        if sum(x(i,:)) > 3
            idx = find(x(i,:)==1);
            % keep = idx(randperm(length(idx),3));
            [~,order] = sort(v(i,idx),'descend'); % 速度越大,被选中的概率越大
            keep = idx(order(1:3));

            x(i,:) = 0;
            x(i,keep) = 1;
        elseif sum(x(i,:)) < 3
            zero_idx = find(x(i,:)==0);
            add = zero_idx(randperm(length(zero_idx), 3-sum(x(i,:))));
            x(i,add) = 1;
        end
    end

    %% --------- 记录与可视化 ----------
    record(iter) = fy_max;
    subplot(1,2,1);
    bar(y_max);
    title(['迭代 ',num2str(iter),' - 选址分布']);
    xlabel('位置编号'); ylabel('是否选中');

    subplot(1,2,2);
    plot(record(1:iter),'LineWidth',1.5);
    title('适应度收敛');
    xlabel('迭代次数'); ylabel('最优值');
    pause(0.05);
end

%% ================= 输出结果 ==================
disp(['最优值：', num2str(fy_max)]);
disp(['选中位置索引：', num2str(find(y_max==1))]);

figure;
subplot(1,2,1); bar(y_max); title('最终选址结果');
subplot(1,2,2); plot(record,'LineWidth',1.5); title('收敛曲线');

end
