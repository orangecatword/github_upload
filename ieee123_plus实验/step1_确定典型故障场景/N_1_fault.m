% 上级电网和IEEE123断开的边为111

% 通过循环函数,得到N-1故障损失最大的边
% 2026.4.23 损失最大的边应该是110
clc
clear
warning off
S = zeros(122,1); % 定义不同场景下的故障特征量
tic
for i = 1:122
    S(i) = function_fault_1(i);
end
toc
[maxValue, maxIndex] = max(S);
