%% 失败的方法_尝试上下层模型融合用求解器求解,未成功

function[objective] = up(~,~)

x = sdpvar(1,33);
xmin = 0.1;
xmax = 0.3;
position = binvar(1,33);

%% 设约束
Constraints = []; 
for j = 1:33
    Constraints = [Constraints,x(j)>=xmin*position(j)];
    Constraints = [Constraints,x(j)<=xmax*position(j)];
end
Constraints = [Constraints,sum(position) == 3]; 
objective = 2/3*(function2(x)) + 1/3*100*sum(x,2);


%% 设求解器
ops = sdpsettings('solver', 'cplex', 'verbose', 0);
ops.cplex.preprocessing.presolve = 1; % 启用预处理
% ops.cplex.workmem = 8192;  % 8GB
ops.cplex.workmem = 4096; 
% ops.cplex.mip.tolerances.mipgap = 0.02;  % 放宽最优间隙
% ops.cplex.parallel = 1;  % 启用并行计算
ops.cplex.mip.strategy.search = 1;  % 使用动态搜索
% ops.cplex.mip.strategy.heuristicfreq = 10; % 调整启发式频率

ops.cplex.mip.tolerances.mipgap = 0.05;  % 放宽最优间隙
ops.cplex.mip.tolerances.absmipgap = 1e-3; % 绝对间隙
ops.cplex.mip.strategy.heuristicfreq = 50; % 调整启发式频率

% 新加入
ops.cplex.threads = 1;        % ⭐ 关键 限制每个 CPLEX 求解器仅用 1 线程
ops.cplex.parallel = 0;       % ⭐ 关键 禁用 CPLEX 内部并行
ops.cplex.nodefileind = 2;    % 启用节点压缩磁盘文件
sol=optimize(Constraints,objective,ops);
