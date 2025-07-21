%% 初始化种群
% 个人想法:忽略主论文表8 表9的内容,按照论文算例部分解释 定死这几个位置, 光伏 4 7 27; 风力 13 17 25
function [LC_wt, LC_pv] = initialize_population(pop_size)

% 初始化参数
LC_wt = zeros(pop_size, 33);
LC_pv = zeros(pop_size, 33);

% 初始化容量
Cap_wt = zeros(pop_size, 33);
Cap_pv = zeros(pop_size, 33);

Loc_pv_initial = [0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0];
Loc_wt_initial = [0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0];

% 选址
LC_wt = repmat(Loc_wt_initial, pop_size, 1);
LC_pv = repmat(Loc_pv_initial, pop_size, 1);

% 定容 默认WT与PV装机容量比为5:3

    % 赋值容量结果是固定值了
    % for i = 1:pop_size
    %     for j = 1:33
    %         if LC_wt(i,j) == 1
    %         Cap_wt(i,j) = min(500, 1000 - sum(Cap_wt(i,:)));
    %         end
    %         if LC_pv(i,j) == 1
    %         Cap_pv(i,j) = min(300, 700 - sum(Cap_pv(i,:)));
    %         end
    %     end
    % end

for i = 1:pop_size
    for j = 1:33
        if LC_wt(i,j) == 1
            % 风电容量分配：在剩余容量范围内随机取值
            remaining_wt = 1000 - sum(Cap_wt(i,:));
            Cap_wt(i,j) = max(200,min(500, remaining_wt * rand()));
        end
        if LC_pv(i,j) == 1
            % 光伏容量分配：在剩余容量范围内随机取值
            remaining_pv = 700 - sum(Cap_pv(i,:));
            Cap_pv(i,j) = max(60,min(300, remaining_pv * rand()));
        end
    end
end

% 最终的初始容量
LC_wt = LC_wt .* Cap_wt;
LC_pv = LC_pv .* Cap_pv;



