% 基础使用多个图像绘制，并添加图例-导入优化算法最后的结果数据，就能得到小提琴图数据
X1=[1];
Y1=randn(100,1)+sin(X1);
X2=[2];
Y2=randn(100,1)+cos(X2);

figure
Hdl1=violinChart(gca,X1,Y1,[0     0.447 0.741]);
Hdl2=violinChart(gca,X2,Y2,[0.850 0.325 0.098]);
legend([Hdl1.F_legend,Hdl2.F_legend],{'randn+sin(x)','randn+cos(x)'});
