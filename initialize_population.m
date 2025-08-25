%% 初始化种群
% 个人想法:忽略主论文表8 表9的内容,按照论文算例部分解释 定死这几个位置, 光伏 4 7 27; 风力 13 17 25
% 个人想法2:容量与风光出力的关系:100MVA的装机容量对应IEEE33BW中的mpc.pv与mpc.wt*100MW
% 暂时的取舍:100MVA为单机装机容量,而IEEE33BW中的mpc.pv mpc.wt是三个装机总的出力
% 2025.8.23新想法:是定死了这三个点 但这三个点中有容量为0也可以,现在问题是如何在有约束条件下随机赋值
function [LC_wt, LC_pv] = initialize_population(pop_size)

% 初始化参数
LC_wt = zeros(pop_size, 33);
LC_pv = zeros(pop_size, 33);

% 初始化容量
Cap_wt = zeros(pop_size, 33);
Cap_pv = zeros(pop_size, 33);

% 初始化容量和
Cap_wt_sum = zeros(pop_size, 1);
Cap_pv_sum = zeros(pop_size, 1);

Loc_pv_initial = [0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0];
Loc_wt_initial = [0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0];

% 选址
LC_wt = repmat(Loc_wt_initial, pop_size, 1);
LC_pv = repmat(Loc_pv_initial, pop_size, 1);

% 定容 假设WT与PV装机容量比为5:3

for i = 1:pop_size
    Cap_wt_sum = rand() * 1000; % 风电最大容量和
    Cap_pv_sum = rand() * (1700 - Cap_wt_sum); % 光伏最大容量和
    taken_wt = 0; % 前面风电占据的容量
    taken_pv = 0;% 前面光伏占据的容量
    for j = 1:33
        if LC_wt(i,j) == 1
            % 风电容量分配：在剩余容量范围内随机取值
            Cap_wt(i,j) = min(800,rand()*(Cap_wt_sum-taken_wt));
            taken_wt = taken_wt + Cap_wt(i,j);
        end
        if LC_pv(i,j) == 1
            % 光伏容量分配：在剩余容量范围内随机取值
            Cap_pv(i,j) = min(800,rand()*(Cap_pv_sum-taken_pv));
            taken_pv = taken_pv + Cap_pv(i,j);
        end
    end
end

% 最终的初始容量
LC_wt = LC_wt .* Cap_wt;
LC_pv = LC_pv .* Cap_pv;



