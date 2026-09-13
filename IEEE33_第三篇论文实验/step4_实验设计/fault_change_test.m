es_nodes = [8  23  30];
Ees_max  = [0.5225 0.51139 0.46168];

n = 10;   % 实验次数
f      = zeros(n,1); 
r_load = zeros(n,1); 
R_load = zeros(33,n);
V_bias = zeros(33,n);
Rps    = zeros(n,1);

for i=1:n
    [f(i),r_load(i),V_bias(:,i),R_load(:,i),Rps(i)] = fault_with_time(es_nodes,Ees_max);
end

f      = sum(f)/n;
r_load = sum(r_load)/n;
R_load = sum(R_load,2)/n;
V_bias = sum(V_bias,2)/n;
Rps    = sum(Rps)/n;