%% 初始化种群
% 个人想法:按照论文中来的话 固定好这几个位置, 光伏 4 7 27; 风力 13 17 25;有一定的概率不选其中的某个点
function [LC_wt, LC_pv] = initialize_population(pop_size)
% 初始化参数
LC_wt = zeros(pop_size, 33);
LC_pv = zeros(pop_size, 33);

% 初始化容量
Cap_wt = zeros(pop_size, 33);
Cap_pv = zeros(pop_size, 33);

Loc_pv_initial = [0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0]
Loc_wt_initial = [0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0]

    for i = 1:pop_size
        for j = 1:33
            if rand() > 0.1  % 90%概率执行此操作
            LC_wt(i,j) = Loc_wt_initial(j);
            LC_pv(i,j) = Loc_pv_initial(j);
            else
            LC_wt(i,j) = 0;
            LC_pv(i,j) = 0;
            end

            if LC_wt(i,j) == 1
            Cap_wt(i,j) = min(500, 1000 - sum(Cap_wt(i,:)));
            end
            if LC_pv(i,j) == 1
            Cap_pv(i,j) = min(300, 700 - sum(Cap_pv(i,:)));
            end
        end
    end
LC_wt = LC_wt .* Cap_wt;
LC_pv = LC_pv .* Cap_pv;



