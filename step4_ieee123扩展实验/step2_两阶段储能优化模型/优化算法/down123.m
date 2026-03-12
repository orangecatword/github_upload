function[objective] = down123(x)

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

%% 同时考虑N-1 N-2 N-5场景时,一定要注意权重(N-5损失远大于N-1 N-2)
tic
%% 断开2条支路下-无储能成本0.8711
% 不允许环路
% idx =  [18 22 24] ;            % 比没有储能的更优解
% vals =  [0.7317 0.1707 0.3188];
% idx =  [2 16 26] ;
% vals = [0.1136 0.1144 0.0025];  % GWO的最优解 在储能成本上过于保守,使得对损失的优化效果并不理想
% idx =  [2 16 20];
% vals = [0.4249 0.5075 0.4712]; % GWO第一次优化后的最优解 反而在负荷损失上的效果更好,但相对来说并不注重储能成本
% 允许环路
% idx  =  [16 18 24];            % 允许环路情况下的解-此时能更快更好的得到优解
% vals =  [0.1344 0.2665 0.1278]; 储能成本0.6477-0.1190

%% 断开5条支路下-无储能成本7.0482
% 允许环路
% idx  =  [5 18 25]; 
% vals =  [0.6153     0.5442      0.6729];
%% GWO综合整体情况 ob1 ob2 X
% 允许环路
% idx = [18 22 24]; 
% vals = [0.4592 0.0841 0.2140];

% idx = [18 24 32];   0.0069+1.8399+2*1.2467=4.3402
% vals = [0.3438 0.2956 0.6073];

% idx = [18 24 32];  0.0069+1.8399+1.2467 = 3.0935
% vals = [0.3438 0.2956 0.6073];
%% PSO综合整体情况
% 允许环路
% idx = [14 24 32]; 0.2359+1.5121+1.8470 = 3.5949
% vals = [0.7719 0.4226 0.6525];

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

% X = 0.01*Ees_cost*sum(vals);

ob1 = 1/2*(f1+f2) + (0.940*f3+0.032*f4+0.028*f5);
% objective(i,:) = 1/2*(f1+f2) + (0.940*f3+0.032*f4+0.028*f5) + 0.01*Ees_cost*sum(vals);
ob2 = 0.5564*f6+0.2162*f7+0.1714*f8+0.0560*f9;
% objective(i,:) = 0.5564*f6+0.2162*f7+0.1714*f8+0.0560*f9 + 0.01*Ees_cost*sum(vals);
X = 0.01*Ees_cost*sum(vals);
objective(i,:) = 1/2*(f1+f2) + (0.940*f3+0.032*f4+0.028*f5) + 0.01*Ees_cost*sum(vals) +...
    0.5564*f6+0.2162*f7+0.1714*f8+0.0560*f9 ;

% r_load = (r_load1+r_load2+r_load3+r_load4+r_load5)/5;
% V_bias = (V_bias1+V_bias2+V_bias3+V_bias4+V_bias5)/5;
% r_load = (r_load6+r_load7+r_load8+r_load9)/4;
% V_bias = (V_bias6+V_bias7+V_bias8+V_bias9)/4;

r_load = (r_load1+r_load2+r_load3+r_load4+r_load5+r_load6+r_load7+r_load8+r_load9)/9;
V_bias = (V_bias1+V_bias2+V_bias3+V_bias4+V_bias5+V_bias6+V_bias7+V_bias8+V_bias9)/9;
r_load = value(r_load);
V_bias = value(V_bias);
end



