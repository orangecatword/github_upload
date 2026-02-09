function[objective] = down_capacity(x)

%% 功能：将1×33向量分解为3个独热向量 + 非零值列表
% 样本数N
N = 10;
objective = zeros(N,1);
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
% vals = vals';
%% 同时考虑N-1 N-2 N-5场景时,一定要注意权重(N-5损失远大于N-1 N-2)
%% EGWO综合优化结果 ob1 ob2 X
% vals = [0.3510;0.4813;0.4708]; % 0.2435+0.7182+1.3031 = 2.2648
% 取消Fij约束后
% vals = [0.5413 0.5326 0.5768]; % 0.0104+1.1537+1.6507 = 2.8148 (10-5
% vals = [0.5830 0.5830 0.5830]; %0.0376+1.1037+1.7490 = 2.8904 (10-5

%% PSO3综合优化结果 ob1 ob2 X
% 取消Fij约束后
vals = [0.5625 0.5837 0.5688]; % 0.0387+1.1044+1.7150 = 2.8581
tic
[f1,r_load1,V_bias1] = function1(vals);
[f2,r_load2,V_bias2] = function2(vals);
[f3,r_load3,V_bias3] = function3(vals);
[f4,r_load4,V_bias4] = function4(vals);
[f5,r_load5,V_bias5] = function5(vals);
[f6,r_load6,V_bias6] = function6(vals);
[f7,r_load7,V_bias7] = function7(vals);
[f8,r_load8,V_bias8] = function8(vals);
[f9,r_load9,V_bias9] = function9(vals);
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
r_load = value(r_load);
V_bias = value(V_bias);
end




