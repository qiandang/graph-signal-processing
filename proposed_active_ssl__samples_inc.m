function [Rt, error_opt] = proposed_active_ssl__samples_inc(mem_fn, Ln, k, Ln_k, S_opt)
%   AUTHOR: Akshay Gadde, USC
%   This function does the following: it finds the optimal set of nodes to 
%   add to the given sampling set. Using the labels on this net sampling
%   set, it predicts the unknown labels and reports the classification
%   error along with the net sampling set.

% % %
% PARAMETER DESCRIPTION
% 
% INPUT
%S_opt:已经采好的样本
% mem_fn:  ground truth membership functions of each
% Ln: normalized Laplacian
% k: Power of Laplacian while computing cutoff, higher the order,
% greater the accuracy, but the complexity is also higher.
% Ln_k: kth power of the Laplacian 
% 
% OUTPUT
% error_opt: classification error
% % %

%% options

num_iter = 100;%POCS的投影次数
%% compute its cutoff frequency 
[~,omega] = eigs(Ln_k(~S_opt,~S_opt),1,'sm');%只会在补集选择，每改变一个的变化
cutoff = abs(omega)^(1/k);

% S_opt = zeros(N,1);
% S_opt(Sample) = true;%已经选择的为1逻辑
% queries = Sample;
% [~,omega] = eigs(Ln_k(~S_opt,~S_opt),1,'sm');
% omega = abs(omega)^(1/k);%求出截止频率
%% reconstruction using POCS

tic;
norm_val = zeros(num_iter,1); % used for checking convergence

% reconstruction using POCS

% approximate low pass filter using SGWT toolbox
filterlen = 10;
alpha = 8;
freq_range = [0 2];%由于归一化的矩阵决定
g = @(x)(1./(1+exp(alpha*(x-cutoff))));
c = sgwt_cheby_coeff(g,filterlen,filterlen+1,freq_range);%1*11的系数


% initialization
mem_fn_du = mem_fn;
mem_fn_du(~S_opt,:) = 0;%down-up sampling
%只有选中的数据才保持原来的men-fun,也就是进行采样的样本，十个归属类一起考虑
mem_fn_recon = sgwt_cheby_op(mem_fn_du,Ln,c,freq_range);%向窄带信号空间投影的结果//f0

for iter = 1:num_iter % takes fewer iterations
    % projection on C1
    err_s = (mem_fn_du-mem_fn_recon); 
    err_s(~S_opt,:) = 0; % error on the known set
    
    % projection on C2
    mem_fn_temp = sgwt_cheby_op(mem_fn_recon + err_s,Ln,c,freq_range); % err on S approx LP
    
    norm_val(iter) = norm(mem_fn_temp-mem_fn_recon); % to check convergence，用上次的结果和这次的范数
    if (iter > 1 && norm_val(iter) > norm_val(iter-1) ), break; end % avoid divergence  发散，范数应该越来越小
    mem_fn_recon = mem_fn_temp;
end
Rt = toc;
% predicted class labels  对于标签来说，需要考虑十种归属度的比例比较
[~,f_recon] = max(mem_fn_recon,[],2);
%max(mem_fn_recon,[],2)直接返回每行的最大值，[a,b]=max(mem_fn_recon,[],2)则是a为具体值，b为index

% true class lables
[~,f] = max(mem_fn,[],2);

% reconstruction error 正确率的计算方法，估计的和原来的相等就是正确，只考虑未知标签的估计更加合理
error_opt = sum(f(~S_opt)~=f_recon(~S_opt))/sum(~S_opt); % error for unknown labels only