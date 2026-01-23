function [Alpha_score, Alpha_pos, record] = simpleGWO_PSO(~,~)
%% =========================================================
%%  GWO–PSO 混合算法（组合编码）
%%  外层：储能选址（33维，仅3个非零）
%%  内层：down(x) 不变
%% =========================================================

%% =============== 参数设置 ===============================
N = 20;                 % 种群规模
d = 33;                 % 维度（IEEE33节点）
Max_iter = 10;          % 最大迭代次数

limit  = [0.8, 1.2];    % 连续编码范围
vlimit = [-0.2, 0.2];   % 速度约束

forbidden_idx = [1 7 12 27];   % 禁止布置储能的位置

% PSO 参数
w  = 0.7;
c1 = 1.2;
c2 = 1.2;
c3 = 1.2;

%% =============== 初始化 ===============================
x = limit(1) + (limit(2)-limit(1))*rand(N,d);
v = zeros(N,d);

% 强制每个个体去掉禁止位置外仅选择 3 个非零变量
for i = 1:N
    all_idx = 1:d;
    allowed_idx = setdiff(all_idx, forbidden_idx);
    idx = allowed_idx(randperm(length(allowed_idx), 3));
    
    mask = false(1,d);
    mask(idx) = true;
    x(i,~mask) = 0;
end

%% =============== Alpha / Beta / Delta ===================
Alpha_pos = zeros(1,d); Alpha_score = inf;
Beta_pos  = zeros(1,d); Beta_score  = inf;
Delta_pos = zeros(1,d); Delta_score = inf;

record = zeros(Max_iter,1);
tic
%% =============== 主循环 ===============================
for iter = 1:Max_iter
    %% ---------- 适应度 ----------
    fx = down(x);

    %% ---------- 更新 α β δ ----------
    for i = 1:N
        if fx(i) < Alpha_score
            Delta_score = Beta_score;  Delta_pos = Beta_pos;
            Beta_score  = Alpha_score; Beta_pos = Alpha_pos;
            Alpha_score = fx(i);        Alpha_pos = x(i,:);
        elseif fx(i) < Beta_score
            Delta_score = Beta_score;  Delta_pos = Beta_pos;
            Beta_score  = fx(i);        Beta_pos = x(i,:);
        elseif fx(i) < Delta_score
            Delta_score = fx(i);
            Delta_pos   = x(i,:);
        end
    end

    %% ---------- GWO–PSO 混合更新 ----------
    a = 2 - 2*iter/Max_iter;   % GWO 收敛因子

    for i = 1:N
        r1 = rand(1,d); r2 = rand(1,d); r3 = rand(1,d);

        % Eq.(28)：速度更新
        v(i,:) = w*v(i,:) ...
            + c1*r1.*(Alpha_pos - x(i,:)) ...
            + c2*r2.*(Beta_pos  - x(i,:)) ...
            + c3*r3.*(Delta_pos - x(i,:));

        % 速度限幅
        v(i,v(i,:)>vlimit(2)) = vlimit(2);
        v(i,v(i,:)<vlimit(1)) = vlimit(1);

        % 位置更新
        x(i,:) = x(i,:) + v(i,:);
    end

    %% ---------- 边界 ----------
    x(x > limit(2)) = limit(2);
    x(x < limit(1)) = limit(1);

    %% ---------- 组合稀疏化（关键） ----------
    for i = 1:N
        % 禁止节点
        x(i,forbidden_idx) = 0;

        allowed = setdiff(1:d, forbidden_idx);
        [~,idx] = sort(abs(x(i,allowed)),'descend');

        mask = false(1,d);
        mask(allowed(idx(1:3))) = true;
        x(i,~mask) = 0;
    end

    %% ---------- 记录 ----------
    record(iter) = Alpha_score;
    subplot(1,2,1);
    bar(Alpha_pos);
    title(['迭代 ',num2str(iter),' 最优解']);
    subplot(1,2,2);
    plot(record(1:iter),'LineWidth',1.5);
    title('适应度收敛');
    pause(0.05)
end
toc
%% =============== 输出 ===============================
disp(['最优目标值: ',num2str(Alpha_score)])
disp(['储能位置: ',num2str(find(Alpha_pos~=0))])
disp(['对应容量: ',num2str(Alpha_pos(Alpha_pos~=0))])

figure;
subplot(1,2,1); bar(Alpha_pos); title('最终选址');
subplot(1,2,2); plot(record,'LineWidth',1.5); title('收敛曲线');

end
