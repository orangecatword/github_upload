mpc = case33bw1;
branch = mpc.branch; 
upstream=zeros(33,37);%代表流入节点支路
dnstream=zeros(33,37);%代表流出节点支路

for i=1:37
    upstream(branch(i,2),i)=1;
end

for i=1:37
    dnstream(branch(i,1),i)=1;
end