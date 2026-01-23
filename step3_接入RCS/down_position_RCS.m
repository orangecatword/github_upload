function[objective] = down_position_RCS(x)

%% 功能：将1×33向量分解为3个独热向量 + 非零值列表
% 样本数N
N = 10;
objective = zeros(N,1);
% 建立20个样本
for i = 1:N
% ===== 1️⃣ 找出非零元素位置与数值 =====
idx = find(x(i,:) ~= 0);      % 非零元素索引
% vals = x(i,idx);           % 非零元素取值

% ===== 2️⃣ 构造3×33的独热矩阵 =====
% Lc_pes = zeros(length(idx), length(x(i,:)));
% for j = 1:length(idx)
%     Lc_pes(j, idx(j)) = 1;
% end

% ===== 3️⃣ 非零值向量 =====
% Ees_max = vals(:);   % 转为列向量 (3×1)

% ===== 4️⃣ 拆分为单独变量 =====
% % 独热矩阵逐行拆分
% Lc_pes1 = Lc_pes(1, :);
% Lc_pes2 = Lc_pes(2, :);
% Lc_pes3 = Lc_pes(3, :);
% 
% % 非零值逐个拆分
% Ees_max1 = Ees_max(1);
% Ees_max2 = Ees_max(2);
% Ees_max3 = Ees_max(3);

%% 同时考虑N-1 N-2 N-5场景时,一定要注意权重(N-5损失远大于N-1 N-2)
tic
[f1,~,~] = function1_1(idx);
[f2,~,~] = function2_1(idx);
[f3,~,~] = function3_1(idx);
% [f4,~,~] = function4_1(idx);
[f5,~,~] = function5_1(idx);
[f6,~,~] = function6_1(idx);
% [f7,~,~] = function7_1(idx);
% [f8,~,~] = function8_1(idx);
% [f9,~,~] = function9_1(idx);
[f10,~,~] = function10_1(idx);
toc
 
% objective(i,:) = 1/2*(f1+f2)/(10^-1)+ (0.924*f3+0.034*f4+0.020*f5+0.022*f6)/(10^-1)+(0.243*f7+0.055*f8+0.217*f9+0.485*f10)/(10^2);
objective(i,:) = 1/2*(f1+f2)/+ (0.924*f3+0.020*f5+0.022*f6)+0.485*f10;
end



