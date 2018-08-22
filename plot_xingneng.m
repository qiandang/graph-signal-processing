%%使用sp（EV）采样后的恢复误差
 addpath(genpath('error'));
  labelled_percentage = 0.04:0.01:0.1;
%labelled_percentage = [0.04,0.045,0.05,0.055,0.06,0.065,0.07];
 figure
%  load('e_sp_k13.mat');%当k为13时，使用EV采样后AMIA恢复方法的恢复精确度
%  plot(labelled_percentage,e_AMIA([1 2 3 4 5 6 7],:),'-<k');
%  hold on
 load('error_all.mat');%当k为13时，使用EV采样后AMIA恢复方法的恢复精确度
 plot(labelled_percentage, e_AMIA([1 3 5 7 9 11 13],:),'-p','LineWidth',1.5,'Color',[1,0,0])
 hold on 
 plot(labelled_percentage, e_MIA([1 3 5 7 9 11 13],:),'-s','LineWidth',1.5,'Color',[0,0,1])
 hold on 
 plot(labelled_percentage,e_pocs([1 3 5 7 9 11 13],:),'-d','LineWidth',1.5,'Color',[0,0,0]);
 hold on
 plot(labelled_percentage,e_RE([1 3 5 7 9 11 13],:),'->','LineWidth',1.5,'Color',[0,0.5,0]);

%%随机采样后恢复的误差
% addpath(genpath('error'));
% labelled_percentage = 0.04:0.01:0.1;
% load('e_rnd_all.mat');%导入一组数据集，进行各种半监督学习
% figure
% plot(labelled_percentage, e_MIA,'-<k')
% hold on 
% plot(labelled_percentage,e_pocs,'-^g');
% hold on
% plot(labelled_percentage,e_RE,'->b');
% hold on
% plot(labelled_percentage,e_AMIA,'-sr');
 legend('A-MIA','MIA','ILSR','LS');
 xlabel('Percentage of labeled data'); 
 ylabel('Accuracy');
 set(gca,'YLim',[0.88 0.93]);
%%不同k值的AMIA恢复误差进行比较
% addpath(genpath('error'));
%  labelled_percentage = 0.04:0.01:0.1;
% figure
% load('e_all.mat');
% plot(labelled_percentage,e_MIA,'-sk');
% hold on
% load('e_sp_k7.mat');%当k为7时，使用EV采样后AMIA恢复方法的恢复精确度
% plot(labelled_percentage,e_AMIA,'-dr');
% hold on
% load('e_sp_k12.mat');%当k为12时，使用EV采样后AMIA恢复方法的恢复精确度
% plot(labelled_percentage,e_AMIA,'-hb');
% hold on
% load('e_sp_k13.mat');%当k为13时，使用EV采样后AMIA恢复方法的恢复精确度
% plot(labelled_percentage,e_AMIA,'-dg');
% legend('矩阵逆近似(MIA)重建','加速的矩阵逆近似(A-MIA)重建,k = 7','k = 12','k = 13');
% xlabel('percentage of labeled data'); 
% ylabel('Accuracy');