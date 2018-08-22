clear;
close all;
clc;

%% data and required toolboxes
    load(['set' num2str(1) '.mat']);%导入一组数据集，进行各种半监督学习

    N = size(A,1);%节点的数目
    
    %拓扑结构已经固定对于每一组的数据来说

    % compute the symmetric normalized Laplacian matrix
    d = sum(A,2);
    d(d~=0) = d.^(-1/2);
    Dinv = spdiags(d,0,N,N);%spdiags 和简单的diag是什么区别呢？
    %A = spdiags(B,d,m,n) 
    %creates an m-by-n sparse matrix by taking the columns of B and placing them along the diagonals specified by d
    Ln = speye(N) - Dinv*A*Dinv;%都是稀疏的存储方式，1000000的数据以52万的空间
    clear Dinv;

    % make sure the Laplacian is symmetric
    Ln = 0.5*(Ln+Ln.');%之前存在非对称性质
    % make sure the Laplacian is symmetric
    [v,eigval] = eig(full(Ln));
    x = v'*mem_fn(:,1);
    e = eigval(eigval~=0);
    figure
%     bar(e',x',0.01);%1是width
    plot(e',x');
     set(gca,'XLim',[0 2]);%X轴的数据显示范围  
     set(gca,'YLim',[-2 6]);%X轴的数据显示范围  
     xlabel('λ_i');
     ylabel('f(λ_i)');
     eng_total = (norm(x))^2;
     eng = eng_total*0.9;
     eng_K = (norm(x(1:60,:)))^2;
%     hold on
%     scatter(e',x','r*');

