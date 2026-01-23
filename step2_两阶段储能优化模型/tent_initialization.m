function population = tent_initialization(pop_size, dim, lb, ub)
% 输入：pop_size（种群大小）, dim（变量维数）, lb/ub（变量下界/上界，标量或1×dim向量）
% 输出：初始种群矩阵（pop_size × dim）

% 处理边界输入（支持标量或向量）
if isscalar(lb)
    lb = lb * ones(1, dim);
end
if isscalar(ub)
    ub = ub * ones(1, dim);
end

% 初始化种群矩阵
population = zeros(pop_size, dim);

% 生成初始混沌值（避免0.5不动点）
x0 = rand();
while x0 == 0.5  % 排除0.5避免周期性
    x0 = rand();
end

% Tent混沌映射参数
epsilon = 1e-10;  % 扰动系数

% 遍历每个个体的所有维度
for i = 1:pop_size
    for j = 1:dim
        % Tent映射核心迭代
        if x0 < 0.5
            x0 = 2 * x0;
        else
            x0 = 2 * (1 - x0);
        end
        
        % 加入随机扰动避免周期性
        x0 = x0 + epsilon * rand();
        x0 = mod(x0, 1);  % 确保值域在[0,1)
        
        % 映射到解空间
        population(i, j) = lb(j) + x0 * (ub(j) - lb(j));
    end
end
end