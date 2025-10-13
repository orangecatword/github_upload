% 运用典型场景集法处理风光新能源的不确定性
% 改进的想法:赋予不同场景物理含义;将kmeans改为W距离
clc
clear

% 定义变量
% 多场景下的风光值
pv = zeros(4, 24);
wind = zeros(4, 24);

U = zeros(4, 33, 24); % 电压
P = zeros(4, 33, 24); % 节点功率
f = zeros(4,1);       % 适应度函数

[C_res,pk] = datap; % 获取在不同场景下的风光数据 pk:每种场景发生的概率
for i = 1:4
    pv(i,:) = C_res(i,1:24);
    wind(i,:) = C_res(i,25:48);
end


for i = 1:4
    [U(i,:,:),P(i,:,:),f(i)] = function1(pv(i,:),wind(i,:),pk(i));
end
f_sum = sum(f); % 最终的目标函数

