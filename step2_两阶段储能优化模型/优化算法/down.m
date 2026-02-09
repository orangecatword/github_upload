function[objective] = down(x)

%% 功能：将1×33向量分解为3个独热向量 + 非零值列表
% 样本数N
N = 10;
objective = zeros(N,1);
% 建立10个样本
for i = 1:N
% ===== 1️⃣ 找出非零元素位置与数值 =====
idx = find(x(i,:) ~= 0);      % 非零元素索引
vals = x(i,idx);           % 非零元素取值

% ===== 2️⃣ 构造3×33的独热矩阵 =====
% Lc_pes = zeros(length(idx), length(x(i,:)));
% for j = 1:length(idx)
%     Lc_pes(j, idx(j)) = 1;
% end

% ===== 3️⃣ 非零值向量 =====
% Ees_max = vals(:);   % 转为列向量 (3×1)
% ===== 4️⃣ 拆分为单独变量 =====
% 独热矩阵逐行拆分
% Lc_pes1 = Lc_pes(1, :);
% Lc_pes2 = Lc_pes(2, :);
% Lc_pes3 = Lc_pes(3, :);
% 非零值逐个拆分
% Ees_max1 = Ees_max(1);
% Ees_max2 = Ees_max(2);
% Ees_max3 = Ees_max(3);

tic
% 未取消Fij约束前
%% GWO综合整体情况 ob1 ob2 X
% 允许环路
% idx = [18 22 24]; 
% vals = [0.4592 0.0841 0.2140];

% idx = [18 24 32];   0.0069+1.8399+2*1.2467=4.3402
% vals = [0.3438 0.2956 0.6073];

% idx = [18 24 32];  0.0069+1.8399+1.2467 = 3.0935
% vals = [0.3438 0.2956 0.6073];

% idx = [18 22 24];  0.0117+1.7633+1.6228 = 3.3979
% vals = [0.5506 0.3250 0.7472];
%% PSO综合整体情况
% 允许环路
% idx = [14 24 32]; 0.2359+1.5121+1.8470 = 3.5949
% vals = [0.7719 0.4226 0.6525];

% idx = [25 29 33]; 0.008+1.7911+1.2365 = 3.0356
% vals = [0.2 0.5365 0.5];

% 取消Fij约束后-无储能成本 13.9419
%% GWO综合整体情况 ob1 ob2 X
% idx = [25 29 33]; % 3.5352+0.3366+1.5234=5.3952
% vals = [0.3 0.628398 0.59505];
% idx = [2 5 30]; % 12.多 废弃数据
% vals = [0.644764537 0.422885689 0.146514911];

% idx = [21 22 33]; 
% vals = [0.7490 0.1335 0.5708];

%% PSO综合整体情况 ob1 ob2 X
% 同一次实验中 (两次实验中PSO结果代入后与目标函数不同
% idx = [25 29 33]; % 3.2120+0.5444+1.6997=5.4561
% vals = [0.2373 0.8055 0.6569];% (1-6次
% idx = [4 17 33]; % 1.1959+4.8173+1.1841=7.1973
% vals = [0.394707 0.394707 0.394707];(7-10次
tic
[f1,r_load1,V_bias1] = function1(idx,vals);
[f2,r_load2,V_bias2] = function2(idx,vals);
[f3,r_load3,V_bias3] = function3(idx,vals);
[f4,r_load4,V_bias4] = function4(idx,vals);
[f5,r_load5,V_bias5] = function5(idx,vals);
[f6,r_load6,V_bias6] = function6(idx,vals);
[f7,r_load7,V_bias7] = function7(idx,vals);
[f8,r_load8,V_bias8] = function8(idx,vals);
[f9,r_load9,V_bias9] = function9(idx,vals);
toc

% 归一化的方法
Ees_cost = 100; % 电池配置价格 单位:万美元/10MWh wESS = 0.1 % 储能成本权重

ob1 = 1/2*(f1+f2) + (0.9173*f3+0.0544*f4+0.0282*f5);
% objective(i,:) = 1/2*(f1+f2) + (0.940*f3+0.032*f4+0.028*f5) + 0.01*Ees_cost*sum(vals);
ob2 = 0.5576*f6+0.2130*f7+0.1733*f8+0.0561*f9;
% objective(i,:) = 0.5564*f6+0.2162*f7+0.1714*f8+0.0560*f9 + 0.01*Ees_cost*sum(vals);
X = 0.01*Ees_cost*sum(vals);
objective(i,:) = 1/2*(f1+f2) + (0.9173*f3+0.0544*f4+0.0282*f5) + 0.01*Ees_cost*sum(vals) +...
    0.5576*f6+0.2130*f7+0.1733*f8+0.0561*f9;

% r_load = (r_load1+r_load2+r_load3+r_load4+r_load5)/5;
% V_bias = (V_bias1+V_bias2+V_bias3+V_bias4+V_bias5)/5;
% r_load = (r_load6+r_load7+r_load8+r_load9)/4;
% V_bias = (V_bias6+V_bias7+V_bias8+V_bias9)/4;

r_load = (r_load1+r_load2+r_load3+r_load4+r_load5+r_load6+r_load7+r_load8+r_load9)/9;
V_bias = (V_bias1+V_bias2+V_bias3+V_bias4+V_bias5+V_bias6+V_bias7+V_bias8+V_bias9)/9;
end



