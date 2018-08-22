%% 根据不同的采样信号点S计算误差
function [error_opt,Rt,queries]=compE_RE(num_queries_to_add,mem_fn, Ln, K, k, Ln_k, prev_queries)

%compute error for random sampling or SP sampling
%they have the same reconstruction method (perfect recovery)

%input:
%num_queries_to_add: the number of selected samples(%*N-上次已经采过的样本数)  
%prev_queries: Samples have been taken before.
%mem_fn: input signal(no noiseless and approximately bandlimited)
%k:(L^k)S_c中的k次方，k越大越精确但是复杂度也越高
%K:近似带限信号的近似带宽

%output
%error_opt: reconstruction error (原信号与恢复信号中不相等的数的个数和/未被采样个数和) 
%queries:the index of the samples


    N = size(Ln,1);
    [v,~] = eig(full(Ln));
    %% compute optimal sampling set and store its cutoff frequency 
%     S_opt = zeros(N,1);
%     S_opt(Sample) = true;%已经选择的为1逻辑
%     queries = Sample;
S_opt_prev = false(N,1);
S_opt_prev(prev_queries) = true;%已经选择的为1逻辑
[S_opt, ~] = compute_opt_set_inc(Ln_k, k, num_queries_to_add, S_opt_prev);
queries = find(S_opt);%找出逻辑是1 的index% S_opt = zeros(N,1);

    % 考虑近似带限的信号    
    tic;
    x_S = mem_fn(queries,:);
    x_Ke=pinv(v(queries,1:K))*x_S;%pinv 计算伪逆 % x_Ke: estimate of x_K
    x_e=v(:,1:K)*x_Ke;
    % predicted class labels  对于标签来说，需要考虑十种归属度的比例比较
    [~,f_recon] = max(x_e,[],2);
    %max(mem_fn_recon,[],2)直接返回每行的最大值，[a,b]=max(mem_fn_recon,[],2)则是a为具体值，b为index
    % true class lables
    [~,f] = max(mem_fn,[],2);
    % reconstruction error 正确率的计算方法，估计的和原来的相等就是正确，只考虑未知标签的估计更加合理
    error_opt = sum(f(~S_opt)~=f_recon(~S_opt))/sum(~S_opt); % error for unknown labels only
    Rt = toc;
end