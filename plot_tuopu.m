clear;
close all;
clc;
%%画一个含有N个节点的拓扑图，信号为USPS中成员函数mem_fn中的第一列
%即整个拓扑图节点归属于数字0的归属值，其中有m个采样点的归属值为1（标成红色），其余为0（标成黑色）
addpath(genpath('data'));
a = load(['set' num2str(1) '.mat']);%导入一组数据集，进行各种半监督学习
N = 70;%取其中的N维出来画图
x = a.mem_fn(:,1);%信号x
size(x(x==1),1); 
m = 7;
B = a.A(1:N,1:N);%取邻接矩阵A的N*N维
 ND_netplot(B,N,m);
