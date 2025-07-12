function y=wblLHS(ns,lamda,k)
% ns=1000;
for i=1:ns
    x(1,i)=wblinv((rand(1,1)/ns+(i-1)/ns),lamda,k);%是韦伯分布函数的反函数
end
[~,s1]=sort(rand(ns,1));
y=x(1,s1)';
%下为验证分布程序
% [x2,c]=hist(y,101);
% t=0:0.01:1;
% bar(t,x2)
