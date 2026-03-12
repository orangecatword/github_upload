function[objective] = down_capacity(x)

%% 功能：将1×33向量分解为3个独热向量 + 非零值列表
% 样本数N
N = 10;
objective = zeros(N,1);
vals = zeros(1,3); 
% 建立10个样本
for i = 1:N
% ===== 1️⃣ 找出非零元素位置与数值 =====
% 非零元素索引
vals = x(i,:);           % 非零元素取值

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

%% 同时考虑N-1 N-2 N-5场景时,一定要注意权重(N-5损失远大于N-1 N-2)
%% EGWO综合优化结果 ob1 ob2 X
vals = [0.2771 0.5357 0.6830]; % 0.0049+0.0027+1.4958=1.5034

%% PSO3综合优化结果 ob1 ob2 X
% vals = [0.3030 0.4149 0.8009]; % 0.0056+0.0067+1.5188=1.5311
%% GWO3综合优化结果 ob1 ob2 X
% vals = [0.2788 0.5478 0.6735]; % 0.0048+0.0027+1.5001=1.5076

tic
[f1,r_load1,V_bias1,R_load1] = function1(vals);
[f2,r_load2,V_bias2,R_load2] = function2(vals);
[f3,r_load3,V_bias3,R_load3] = function3(vals);
[f4,r_load4,V_bias4,R_load4] = function4(vals);
[f5,r_load5,V_bias5,R_load5] = function5(vals);
[f6,r_load6,V_bias6,R_load6] = function6(vals);
[f7,r_load7,V_bias7,R_load7] = function7(vals);
[f8,r_load8,V_bias8,R_load8] = function8(vals);
[f9,r_load9,V_bias9,R_load9] = function9(vals);
toc

% 归一化的方法
Ees_cost = 100; % 电池配置价格 单位:万美元/10MWh wESS = 0.1 % 储能成本权重

ob1 = 1/2*(f1+f2) + (0.9395*f3+0.0323*f4+0.0282*f5);
ob2 = 0.5586*f6+0.2141*f7+0.1713*f8+0.0561*f9;
X = 0.01*Ees_cost*sum(vals);
objective(i,:) = ob1 + ob2 + X;

% r_load = (r_load1+r_load2+r_load3+r_load4+r_load5)/5;
% V_bias = (V_bias1+V_bias2+V_bias3+V_bias4+V_bias5)/5;
% r_load = (r_load6+r_load7+r_load8+r_load9)/4;
% V_bias = (V_bias6+V_bias7+V_bias8+V_bias9)/4;
r_load = (r_load1+r_load2+r_load3+r_load4+r_load5+r_load6+r_load7+r_load8+r_load9)/9;
V_bias = (V_bias1+V_bias2+V_bias3+V_bias4+V_bias5+V_bias6+V_bias7+V_bias8+V_bias9)/9;
R_load = (R_load1+R_load2+R_load3+R_load4+R_load5+R_load6+R_load7+R_load8+R_load9)/9;
R_load = R_load';
V_bias = sum(V_bias)/4;
end




