

clc
clear
warning off

%% 风光值不确定性的场景削减
pv = zeros(4, 24);
wind = zeros(4, 24);
[C_res,pk] = datap; % 获取在不同场景下的风光数据 pk:每种场景发生的概率
for i = 1:4
    pv(i,:) = C_res(i,1:24);
    wind(i,:) = C_res(i,25:48);
end


%% 考虑负荷需求的无故障配电网优化-function1
U = zeros(4, 33, 24); % 电压
P = zeros(4, 33, 24); % 节点功率
f = zeros(4,1);       % 适应度函数
%% function1
for i = 1:4
    [U(i,:,:),P(i,:,:),f(i)] = function1(pv(i,:),wind(i,:),pk(i));
end
f_sum = sum(f); % 最终的目标函数