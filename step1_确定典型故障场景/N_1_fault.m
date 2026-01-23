% 通过循环函数,得到N-1故障损失最大的边
% 论文版:断开了支路30
% 2025.1.15更新 损失最大的边应该是27   
clc
clear
warning off
S = zeros(32,1); % 定义不同场景下的故障特征量
tic
for i = 1:32
    S(i) = function_fault_1(i);
end
toc
[maxValue, maxIndex] = max(S);
