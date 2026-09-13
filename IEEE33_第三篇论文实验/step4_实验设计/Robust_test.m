function[objective] = Robust_test(~)

%% 随机生成断路情况
n = 200;               % 分为3种情况,蒙特卡罗生成的场景数量分别为:200 500 1000
Zij = zeros(n,32); % 断路情况
for i = 1:n
    for j = 1:32
        if rand<0.05    % 分为三种情况,断路概率分别为:0.02 0.05 0.1
            Zij(i,j)=1;
        end
    end
end

f = zeros(n,1);
r_load = zeros(n,1);
V_bias = zeros(n,4);
R_load = zeros(n,33);
for i = 1:n
    [f(i),r_load(i),V_bias(i,:),R_load(i,:)] = optimized_model(Zij(i,:));
end

objective = sum(f)/n+ X;
r_load = sum(r_load)/n;