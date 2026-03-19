function[objective] = down(x)

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
% idx = [16 18 32]; 
% vals = [0.4269 0.5008 0.9967]; % 0.0136+2.6391+1.9244=4.5771
% 20-20情况
idx = [18 29 32]; 
vals = [0.2730 0.4987 0.7195]; % 0.0129+2.0668+1.4912=3.5709

%% PSO综合整体情况 ob1 ob2 X
% idx = [25 29 33];
% vals = [0.4 0.8423 0.4]; % 0.0066+2.0640+1.6423=3.7129
% 20-20情况
% idx = [17 29 30]; 
% vals = [0.2342 0.6089 0.6090]; % 0.6745+2.1246+1.4521=4.2512

tic
[f1,r_load1,V_bias1,R_load1] = function1(idx,vals);
[f2,r_load2,V_bias2,R_load2] = function2(idx,vals);
[f3,r_load3,V_bias3,R_load3] = function3(idx,vals);
[f4,r_load4,V_bias4,R_load4] = function4(idx,vals);
[f5,r_load5,V_bias5,R_load5] = function5(idx,vals);
[f6,r_load6,V_bias6,R_load6] = function6(idx,vals);
[f7,r_load7,V_bias7,R_load7] = function7(idx,vals);
[f8,r_load8,V_bias8,R_load8] = function8(idx,vals);
[f9,r_load9,V_bias9,R_load9] = function9(idx,vals);

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
r_load = value(r_load);
V_bias = value(V_bias);

end




