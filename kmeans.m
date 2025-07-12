function [C,k]=kmeans(data,K)
% 输入变量:data代表24小时内抽取不同样本的出力值,K代表聚类数
% 输出变量:C代表聚类中心(K×n 矩阵),k代表各簇样本占比(1×K 向量)
[m,n]=size(data); %求输入数据点的个数m 特征为n
id=zeros(m,1);
%K=22;
%随机初始化聚类中心
C=zeros(K,n);
for i=1:K
    C(i,:)=data(randi(m,1),:);     %随机初始化聚类中心%randi(n,1)产生一个1到n的伪随机整数
end
%C(i)表示第i类的聚类中心
 for i=1:100
%while 1
    %分配簇
    for x=1:m
        for y=1:K
            d(y)=norm(data(x,:)-C(y,:)); %计算数据点到每个聚类中心的距离
        end
        [~,idx]=min(d);
        id(x)=idx;
    end
    %更新聚类中心 
    new_C=zeros(K,n);
    num=zeros(K);
    q=0;
    for y=1:K
        for x=1:m
            if id(x,1)==y
            new_C(y,:)=new_C(y,:)+data(x,:);
            num(y)=num(y)+1;
            end
        end
        new_C(y,:)=new_C(y,:)/num(y);
        if norm(new_C(y,:)-C(y,:))<0.1  %判断是否收敛
            q=q+1;
        end
    end
% %更新聚类第二种写法
%     for y=1:K
%         datak=data(id==y,:);  %属于第k类的数据点
%         new_C(y,:)=mean(datak);
%         if norm(new_C(y,:)-C(y,:))<0.1
%             q=q+1;
%         end
%     end
%%%
    if q==K
        break;
    else
        C=new_C;
    end
end
center=C;
% figure
%[id,center]=k_means(data,K);
% for i=1:K
% plot(data(id==i,1),data(id==i,24),'*','Color',[rand(),rand(),rand()]);
% hold on;
% end

%for i=1:K
%plot(data(id==1,1),data(id==1,24),'*','Color',[rand(),rand(),rand()]);
%hold on;
%end
%plot(data(id==2,1),data(id==2,24),'*','Color',[rand(),rand(),rand()]);
%hold on;
%plot(data(id==3,1),data(id==3,24),'*','Color',[rand(),rand(),rand()]);
%hold on;
%plot(data(id==4,1),data(id==4,24),'*','Color',[rand(),rand(),rand()]);
%hold on;
figure;
x1=1:1:n;%y1=1:1:24;
for i=1:K
    y1=i*ones(1,n);
    plot3(x1,y1,center(i,:),'-*','Color',[rand(),rand(),rand()]);
    hold on
end
     xlabel('时段');

     ylabel('场景');

     zlabel('功率');
      grid on;
for i=1:K
k(i)=length(find(id==i))/m;
end