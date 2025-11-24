function[objective] = down(x)

%% 功能：将1×33向量分解为3个独热向量 + 非零值列表
% 样本数N
N = 20;
objective = zeros(N,1);
% 建立50个样本
for i = 1:N
% ===== 1️⃣ 找出非零元素位置与数值 =====
idx = find(x(i,:) ~= 0);      % 非零元素索引
vals = x(i,idx);           % 非零元素取值

% ===== 2️⃣ 构造3×33的独热矩阵 =====
Lc_pes = zeros(length(idx), length(x(i,:)));
for j = 1:length(idx)
    Lc_pes(j, idx(j)) = 1;
end

% ===== 3️⃣ 非零值向量 =====
Ees_max = vals(:);   % 转为列向量 (3×1)

% ===== 4️⃣ 拆分为单独变量 =====
% 独热矩阵逐行拆分
Lc_pes1 = Lc_pes(1, :);
Lc_pes2 = Lc_pes(2, :);
Lc_pes3 = Lc_pes(3, :);

% 非零值逐个拆分
Ees_max1 = Ees_max(1);
Ees_max2 = Ees_max(2);
Ees_max3 = Ees_max(3);

%% 同时考虑N-1 N-2 N-5场景时,一定要注意权重(N-5损失远大于N-1 N-2)
tic
[f2,Psum_loss2,Psum_load2,U2] = function2(Lc_pes1,Lc_pes2,Lc_pes3, Ees_max1,Ees_max2,Ees_max3);
[f3,Psum_loss3,Psum_load3,U3] = function3(Lc_pes1,Lc_pes2,Lc_pes3, Ees_max1,Ees_max2,Ees_max3);
[f4,Psum_loss4,Psum_load4,U4] = function4(Lc_pes1,Lc_pes2,Lc_pes3, Ees_max1,Ees_max2,Ees_max3);
[f5,Psum_loss5,Psum_load5,U5] = function5(Lc_pes1,Lc_pes2,Lc_pes3, Ees_max1,Ees_max2,Ees_max3);
toc

% [U,Pload_total] = function_standard;
% % 故障恢复率
% rload = 1- (Psum_load2+Psum_load3+Psum_load4+Psum_load5)./(4*Pload_total);
% % 网损率
% rloss = 1- (Psum_loss2+Psum_loss3+Psum_loss4+Psum_loss5)./(4*Pload_total);
% % 电压偏差
% deltaU = sum(sum(abs(U-U2)+abs(U-U3)+abs(U-U4)+abs(U-U5)));
objective(i,:) = f2 + f3 + f4 + f5;
end



