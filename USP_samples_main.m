clear;
close all;
clc;

%% data and required toolboxes

addpath(genpath('data'));
addpath(genpath('sgwt_toolbox'));
%addpath(genpath('random_samples'));
addpath(genpath('EV_samples'));
%%
%1.先根据节点数目的要求，以最大化截止频率为目标选择节点集合，greedy的算法，并对应申请label
%2.根据选择的节点确定截止频率，在该截止频率之下的特征值（向量）作为采样集信号表示
%3.根据迭代的算法求解出原始的信号，也就是全部数据集的归属函数，归属函数也是10类的
%根据10类归属函数的恢复情况比较，选出它的label
%循环十个随机的数据集，求均值
%%

% Number of datasets (avg results reported) 数据集产生的随机的次数，在进行一种
num_datasets = 10;
%numan矩阵截止的数目
L = 10;

% Power of Laplacian
% k = 7; % higher k leads to better estimate of the cut-off frequency,尝试12看性能的变化
k = 13;

K = 60;%信号的带宽

%M = 300;%总的采样点数

% compare the classification accuracies
labelled_percentage = 0.04:0.005:0.1;
num_points = length(labelled_percentage);

% t_RE = zeros(num_points, num_datasets);%运行时间
% t_MIA = zeros(num_points, num_datasets);%运行时间
% t = zeros(num_points, num_datasets);

  error_list = zeros(num_points, num_datasets);
  error_list_MIA = zeros(num_points, num_datasets);
  error_list_RE = zeros(num_points, num_datasets);
 error_list_AMIA = zeros(num_points, num_datasets);

for iter = 1:num_datasets

    %% data
    
    fprintf(['\n\nloading set' num2str(iter) '...\n\n']);%\n空一行
    a = load(['set' num2str(iter) '.mat']);%导入一组数据集，进行各种半监督学习
    N = size(a.A,1);%节点的数目

    %% cells to store optimal sampling sets
    
    % We greedily select of batch of nodes to sample. Hence not necessary 
    % to start from scratch when a larger subset of nodes is to be sampled.
    % 当选择更大的采样集合的时候，其实前1%的采样点数就是之前的采样点，每次只需要
    %继续考虑增加的点数
    
     %% computation to be done only once 
    %拓扑结构已经固定对于每一组的数据来说

    
    % compute the symmetric normalized Laplacian matrix
    d = sum(a.A,2);
    d(d~=0) = d.^(-1/2);
    Dinv = spdiags(d,0,N,N);%spdiags 和简单的diag是什么区别呢？
    %A = spdiags(B,d,m,n) 
    %creates an m-by-n sparse matrix by taking the columns of B and placing them along the diagonals specified by d
    Ln = speye(N) - Dinv*a.A*Dinv;%都是稀疏的存储方式，1000000的数据以52万的空间
    clear Dinv;

    % make sure the Laplacian is symmetric
    Ln = 0.5*(Ln+Ln.');%之前存在非对称性质
    % make sure the Laplacian is symmetric
   
    Ln_k = Ln;
    for i = 1:(k-1)
        Ln_k = Ln_k*Ln;% s.^(3)矩阵的指数次是每个元素的单独指数次
    end
   Ln_k = 0.5*(Ln_k+Ln_k.');  
    %% Choosing random sampling sets of different sizes
    %以index的形式存储：queries = find(S_opt);%找出逻辑是1 的index
    
      Rt_RE =  zeros(num_points,1);
      Rt_MIA =  zeros(num_points,1);
      Rt =  zeros(num_points,1);
      Rt_AMIA =  zeros(num_points,1);

       error = zeros(num_points,1);
       error_MIA = zeros(num_points,1);
       error_RE = zeros(num_points,1);
       error_AMIA = zeros(num_points,1);
    
    for index_lp = 1:num_points%标记比例的标签
        fprintf('\n\n*** fraction of data labelled = %f ***\n\n', labelled_percentage(index_lp))
        %*** fraction of data labelled = 0.040000 ***
        s = load(['ev_samples' num2str(iter) '_' num2str(labelled_percentage(index_lp)) '.mat']);%导入一组数据集，进行各种半监督学习

        %save(['D:\matlab\仿真\shouxieshuzishibie\S_opt\rdm_samples' num2str(iter) num2str(index_lp) '.mat'],'S_opt')
        [Rt(index_lp), error(index_lp)] = proposed_active_ssl__samples_inc(a.mem_fn, Ln, k, Ln_k, s.S_opt);
          [error_MIA(index_lp),Rt_MIA(index_lp)] = compE_MIA_samples(a.mem_fn, Ln, L, s.S_opt);                                                                                           
           [error_RE(index_lp),Rt_RE] = compE_RE_samples(a.mem_fn, Ln, K, s.S_opt);
         [error_AMIA(index_lp),Rt_AMIA(index_lp)] = compE_AMIA_samples(a.mem_fn, Ln, Ln_k, K, L, k, s.S_opt);
         
         fprintf('classification error_MIA (proposed) = %f \n\n', error_MIA(index_lp));
         fprintf('classification error (proposed) = %f \n\n', error(index_lp));
         fprintf('classification error_RE (proposed) = %f \n\n', error_RE(index_lp));
         fprintf('classification error_AMIA (proposed) = %f \n\n', error_AMIA(index_lp));

    end
    
      error_list_RE(:,iter) = error_RE;
      error_list_MIA(:,iter) = error_MIA;
      error_list(:,iter) = error;
      error_list_AMIA(:,iter) = error_AMIA;

%     t_RE(:,iter) = Rt_RE;    
%     t_MIA(:,iter) = Rt_MIA;
%     t(:,iter) = Rt;
end
  e_RE = mean(1-error_list_RE,2);
  e_MIA = mean(1-error_list_MIA,2);
  e_pocs = mean(1-error_list,2);
  e_AMIA = mean(1-error_list_AMIA,2);

  save(['D:\matlab\仿真\USPS\USPS\error\error_all.mat'],'e_AMIA','e_RE','e_pocs','e_MIA')

