% 参数设置
pop_size = 10;  % 种群大小
dim = 3;         % 变量维度（对应Re,W,Cs三助剂比例）
lb = [0.8, 0.8, 0.8];  % 下界（摩尔比最小值）
ub = [1.2, 1.2, 1.2];  % 上界（摩尔比最大值）

% 生成初始种群
population = tent_initialization(pop_size, dim, lb, ub);

% 可视化验证（三维空间分布）
scatter3(population(:,1), population(:,2), population(:,3), 'filled');
xlabel('Re'); ylabel('W'); zlabel('Cs');
title('Tent混沌初始化种群分布');