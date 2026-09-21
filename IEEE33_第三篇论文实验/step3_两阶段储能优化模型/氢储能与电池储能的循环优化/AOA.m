% 建立交替优化的大框架 交替优化算法

%% 1.设参数
idx_H = [6 23]; % 初始化氢储能配置情况
vals_H = [1 1];

% 循环优化后的最优结果
% [8,23,24,30]  [0.4823,0.4302,0.1023,0.4960] 
% 6 0.6947
% r_load 0.9929
% Rps 0.9611

for i =1:10
    [idx,vals,f1]     = EGWO_BESS(idx_H,vals_H);
    [idx_H,vals_H,f2] = EGWO_HESS(idx,vals);
    if abs(f1-f2) <= 1e-3
        break
    end
end

