%%采样，产生S_opt
clear;
close all;
clc;

%% data and required toolboxes

addpath(genpath('data'));

%%
%1.先根据节点数目的要求，以最大化截止频率为目标选择节点集合，greedy的算法，并对应申请label
%2.根据选择的节点确定截止频率，在该截止频率之下的特征值（向量）作为采样集信号表示
%3.根据迭代的算法求解出原始的信号，也就是全部数据集的归属函数，归属函数也是10类的
%根据10类归属函数的恢复情况比较，选出它的label
%循环十个随机的数据集，求均值
%%

% Number of datasets (avg results reported) 数据集产生的随机的次数，在进行一种
num_datasets = 10;

% Power of Laplacian
k = 12; % higher k leads to better estimate of the cut-off frequency

% compare the classification accuracies
labelled_percentage = 0.045:0.01:0.095;
num_points = length(labelled_percentage);
for iter = 8:num_datasets

    %% data
    
    fprintf(['\n\nloading set' num2str(iter) '...\n\n']);%\n空一行
    a  = load(['set' num2str(iter) '.mat']);%导入一组数据集，进行各种半监督学习

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

%full函数功能：在MATLAB中，该函数用于把一个稀疏矩阵（sparse matrix）转换成一个全矩阵（full matrix）
    % higher power of Laplacian
    Ln_k = Ln;
    for i = 1:(k-1)
        Ln_k = Ln_k*Ln;% s.^(3)矩阵的指数次是每个元素的单独指数次
    end
   Ln_k = 0.5*(Ln_k+Ln_k.');  
    %% Choosing optimal sampling sets of different sizes
    
    prev_queries = []; % sampling set chosen in previous iteration
    %以index的形式存储：queries = find(S_opt);%找出逻辑是1 的index
    
    prev_nqueries = 0; % number of labels queried so far
    cur_nqueries = 0; % number of labels queried in current iteration
    
    for index_lp = 1:length(labelled_percentage)%标记比例的标签
        fprintf('\n\n*** fraction of data labelled = %f ***\n\n', labelled_percentage(index_lp))
        %*** fraction of data labelled = 0.010000 ***
        nqueries = round(labelled_percentage(index_lp) * N);
 
        cur_nqueries = nqueries - prev_nqueries;%每次只要求继续寻找新的节点，建立在之前的节点之上
        
        S_opt_prev = false(N,1);
        S_opt_prev(prev_queries) = true;%已经选择的为1逻辑
        [ S_opt, ~, ~ ] = compute_opt_set_inc( Ln_k, k, cur_nqueries, S_opt_prev);
        prev_queries = find(S_opt);                                                                                       
        prev_nqueries = nqueries;
        save(['D:\matlab\仿真\USPS\USPS\EV_samples\samples' num2str(iter) '_' num2str(labelled_percentage(index_lp)) '.mat'],'S_opt')
    end   
end
