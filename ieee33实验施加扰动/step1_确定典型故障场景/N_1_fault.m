% r x 减5%时 损失最大的边是27   
% r x 加5%时 损失最大的边是28 
% r x 减10%时 损失最大的边是27 
% r x 加10%时 损失最大的边是28 
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
