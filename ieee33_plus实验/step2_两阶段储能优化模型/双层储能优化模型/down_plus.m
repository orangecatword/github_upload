function[objective] = down_plus(x)

%% 功能：将1×33向量分解为3个独热向量 + 非零值列表
% 样本数N
N = 20;
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

%% GWO综合整体情况 ob1 ob2 X
% 20-20情况
% idx = [19  23  26]; 
% vals = [0.19203     0.52222     0.70673]; 
%% PSOplus综合整体情况 ob1 ob2 X
% 20-20情况 1.4328
% idx = [18  30  33]; 
% vals = [0.14337     0.12719    0.026468]; 

tic
[f1,r_load1,V_bias1,R_load1] = function1_plus(idx,vals);
[f2,r_load2,V_bias2,R_load2] = function2_plus(idx,vals);
[f3,r_load3,V_bias3,R_load3] = function3_plus(idx,vals);
[f4,r_load4,V_bias4,R_load4] = function4_plus(idx,vals);
[f5,r_load5,V_bias5,R_load5] = function5_plus(idx,vals);
[f6,r_load6,V_bias6,R_load6] = function6_plus(idx,vals);
[f7,r_load7,V_bias7,R_load7] = function7_plus(idx,vals);
[f8,r_load8,V_bias8,R_load8] = function8_plus(idx,vals);
[f9,r_load9,V_bias9,R_load9] = function9_plus(idx,vals);
toc

% 归一化的方法
Ees_cost = 105; % 电池配置价格 单位:万美元/10MWh 
w_ess = 0.005; % 储能成本权重-每年可能发生的 1-2 次 电网断电事故;储能预期使用年限：10 – 15 年

ob1 = 1/2*(f1+f2) + (0.8085*f3+0.0766*f4+0.1149*f5);
ob2 = 0.4520*f6+0.3310*f7+0.1600*f8+0.0570*f9;
X = w_ess*Ees_cost*sum(vals);
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
r_load = value(r_load);
V_bias = value(V_bias);

end




