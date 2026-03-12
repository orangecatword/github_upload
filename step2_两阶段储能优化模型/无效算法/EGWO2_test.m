%% 用灰狼优化算法替代 PSO-更改:加入了边界约束 (随机重置法)
%% 外层：组合位置搜索
%% 内层：down(x) 精确模型（不变）

function [Alpha_score, Alpha_pos, record] = EGWO2_test(~,~)

%% ================== 参数设置 ==================
N   = 10;          % 种群规模
d   = 3;          % 维度（IEEE33 节点）
Max_iter = 10;     % 最大迭代次数
limit = [-1.28, 1.28];   % 连续编码范围
% c1 = 0.5;                
% c2 = 0.3;                
% c3 = 0.2;               

%% ================== 初始化种群 ==================
x = limit(1) + (limit(2) - limit(1)) .* rand(N, d);
% v = 0.1*(rand(N, d)-0.5);
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
    % 单模态函数
    % fx = sum(x.^2,2);                    % F1 范围:[-100,100]
    % fx = sum(abs(x),2)+prod(abs(x),2);   % F2 范围:[-10,10]
    % fx=zeros(N,1);
    % for i = 1:N    
    %     for j=1:d
    %         fx(i,:)=fx(i,:)+sum(x(i,1:j))^2;% F3 范围:[-100,100]
    %     end
    % end
    % for i = 1:N
    %     fx(i,:) = max(abs(x(i,:)));      % F4 范围:[-100,100]
    % end
    % for i = 1:N
    %     fx(i,:)=sum(100*(x(i,2:d)-(x(i,1:d-1).^2)).^2+(x(i,1:d-1)-1).^2); % F5 范围:[-30,30]
    % end
    % for i = 1:N  
    %     fx(i,:)=sum(abs((x(i,:)+.5)).^2);  % F6 范围:[-100,100]
    % end
    for i = 1:N  
        fx(i,:)=sum([1:d].*(x(i,:).^4))+rand;% F7 范围:[-1.28, 1.28]
    end
    %% 多模态函数
    % for i = 1:N 
    %     fx(i,:)=sum(-x(i,:).*sin(sqrt(abs(x(i,:)))));% F8 范围:[-500, 500] 最小值：-418.98*（维度的n次方
    % end

    % for i = 1:N
    %     fx(i,:)=sum(x(i,:).^2-10*cos(2*pi.*x(i,:)))+10*d;% F9 范围:[-5.12, 5.12]
    % end

    % for i = 1:N
    %     fx(i,:)=-20*exp(-.2*sqrt(sum(x(i,:).^2)/d))-exp(sum(cos(2*pi.*x(i,:)))/d)+20+exp(1);% F10 范围:[-32,32]
    % end

    % for i = 1:N
    %     fx(i,:)=sum(x(i,:).^2)/4000-prod(cos(x(i,:)./sqrt([1:d])))+1;% F11 范围：[-600,600]
    % end

    % for i = 1:N
    %     fx(i,:)=0.1*((sin(3*pi*x(i,1)))^2+sum((x(i,1:d-1)-1).^2.*(1+(sin(3.*pi.*x(i,2:d))).^2))+...
    %             ((x(i,d)-1)^2)*(1+(sin(2*pi*x(i,d)))^2));% F13 范围：[-50,50]
    % end
    
    %% 复合基准测试函数-维度与设定的维度有冲突
    % 不使用F19函数进行性能测试
    % aH=[3 10 30;.1 10 35;3 10 30;.1 10 35];cH=[1 1.2 3 3.2];
    % pH=[.3689 .117 .2673;.4699 .4387 .747;.1091 .8732 .5547;.03815 .5743 .8828];
    % for i = 1:N
    %     o = 0;
    %     for j=1:4
    %         o=o-cH(j)*exp(-(sum(aH(j,:).*((x(i,:)-pH(j,:)).^2)))); % F19 范围:[1, 3] 最小值-3 与实际运行结果不符
    %     end
    %     fx(i,:) = o;
    % end

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
    %% ---------- GWO 位置更新 ----------
        % var = 1-(iter/Max_iter)^2;% 定义为标准差-对应公式（13）
        var = exp(-100*iter/Max_iter);
        % --- 更新位置 ---
        for j = 1:d
            lamda = var * randn; % 模拟随机误差
            x_t(i,j) = 0.5*Alpha_pos(j) + 0.3*Beta_pos(j) + 0.2*Delta_pos(j) + lamda;
            % x(i,j) = xp(i,j)-(4*rand-2)*(xp(i,j)-x(i,j));
            if x_t(i,j) > 1.28
                x(i,j) = x(i,j)+rand*(1.28-x(i,j));
            elseif x_t(i,j) < -1.28
                x(i,j) = x(i,j)+rand*(-1.28-x(i,j));
            else 
                x(i,j) = x_t(i,j);
            end
        end
    end
    %% ---------- 边界约束 (随机重置法) ----------
    % for i = 1:N
    %     mask_upper = x(i,:) > limit(2);
    %     mask_lower = x(i,:) < limit(1);
    %     if any(mask_upper) || any(mask_lower)
    %         % 越限位置重新随机初始化
    %         x(i, mask_upper | mask_lower) = limit(1) + (limit(2)-limit(1)) * rand;
    %     end
    % end
    %% ---------- 边界约束 ----------
    % x(x > limit(2)) = limit(2);
    % x(x < limit(1)) = limit(1);

    %% ---------- 记录与可视化 ----------
    record(iter) = Alpha_score;
    figure(5)
    subplot(1,2,1);
    bar(Alpha_pos);
    title(['迭代 ', num2str(iter), ' - 最优非零元素分布']);
    xlabel('节点编号'); ylabel('编码值');

    subplot(1,2,2);
    plot(record(1:iter), 'LineWidth', 1.5);
    title('EGWO2 适应度收敛曲线');
    xlabel('迭代次数'); ylabel('最优适应度');

    pause(0.05);
end
toc
%% ================== 输出结果 ==================
disp(['最优目标值：', num2str(Alpha_score)]);
disp(['储能位置：', num2str(find(Alpha_pos ~= 0))]);
disp(['对应容量：', num2str(Alpha_pos(Alpha_pos ~= 0))]);

% figure(6);
% subplot(1,2,1);
% bar(Alpha_pos);
% title('最终最优变量分布');
% 
% subplot(1,2,2);
% plot(record, 'LineWidth', 1.5);
% title('适应度收敛曲线');

end
