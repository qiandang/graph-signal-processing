%% 根据不同的采样信号点S计算误差
function [error_opt,Rt]=compE_MIA_samples(mem_fn, Ln, L, S_opt)
%compute error for the propsoed MIA sampling method whose reconstruction is
%MIA recovery (Section IV in paper)

%input:
%S_opt:已经采好的样本
%mem_fn: input signal(no noiseless and approximately bandlimited)
%L:诺依曼级数的截止值

%output
%error_opt: reconstruction error (原信号与恢复信号中不相等的数的个数和/未被采样个数和) 

%Tpoly: Chebychev matrix polynomial approximation of T, output for recovery
%gama: Neumann series of finally selected S and output for recovery

queries = find(S_opt);%找出逻辑是1 的index% S_opt = zeros(N,1);

tic;    
%%Tpoly: output for recovery
% approximate low pass filter using SGWT toolbox
lambda=eig(Ln);
lambda(30);
filterlen =10;
alpha = 8;
freq_range = [0 2];%由于归一化的矩阵决定
g = @(x)(1./(1+exp(alpha*(x-lambda(30)))));%这里带宽k设成了30
c = sgwt_cheby_coeff(g,filterlen,filterlen+1,freq_range);%1*11的系数
%如果需要计算矩阵的话
% rewrite sgwt_cheby_op.m to sgwt_cheby_matrix.m since we need polynomial
% matrix rather the result of matrix-vector product
Tpoly=sgwt_cheby_matrix(Ln,c,freq_range);%已经验证在各种矩阵稀疏与否的情况下的，都是double形式的最快。

% %gama:  output for recovery
eye_size=length(queries);
        B=eye(eye_size)-Tpoly(queries,queries);
        B_sum=B;
        for l=1:(L-1)
            B_sum=(eye(eye_size)+B_sum)*B;
        end
gama=B_sum+eye(eye_size);

% 近似带限+无噪声的采样信号;直接采样恢复，不涉及加噪以及噪声的均质化
% recovery

     x_S = mem_fn(queries,:);
     x_e=Tpoly(:,queries)*gama*x_S;% accompained elegent close-form recovery method 
     Rt=toc;
% predicted class labels  对于标签来说，需要考虑十种归属度的比例比较
[~,f_recon] = max(x_e,[],2);
%max(mem_fn_recon,[],2)直接返回每行的最大值，[a,b]=max(mem_fn_recon,[],2)则是a为具体值，b为index

% true class lables
[~,f] = max(mem_fn,[],2);

% reconstruction error 正确率的计算方法，估计的和原来的相等就是正确，只考虑未知标签的估计更加合理
error_opt = sum(f(~S_opt)~=f_recon(~S_opt))/sum(~S_opt); % error for unknown labels only

end