mpc = case33bw1;
branch = mpc.branch; 
upstream=zeros(33,37);%代表流入节点支路
dnstream=zeros(33,37);%代表流出节点支路

for i=1:32
    upstream(i+1,i)=1;
end
upstream(21,33)=1;%支路33为21-8支路，流入节点21
upstream(15,34)=1;%支路34为15-9支路，流入节点15
upstream(22,35)=1;%支路35为22-12支路，流入节点22
upstream(33,36)=1;%支路36为33-18支路，流入节点33
upstream(29,37)=1;%支路37为29-25支路，流入节点29
for i=[1:17,19:21,23:24,26:32]
    dnstream(i,i)=1;
end
dnstream(2,18)=1;
dnstream(3,22)=1;
dnstream(6,25)=1;
% 5条流入，对应5条流出
dnstream(8,33)=1;
dnstream(9,34)=1;
dnstream(12,35)=1;
dnstream(18,36)=1;
dnstream(25,37)=1;