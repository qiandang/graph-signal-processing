function [ S_opt,omega] = rdm_sample( L_k, k, num_nodes_to_add )
%随机采样使用三种恢复方法进行比较
%

N = size(L_k,1);

s = randi(N,num_nodes_to_add,1);%产生num_nodes_to_add个随机的采样点
S_opt(s) = true;

[~,omega] = eigs(L_k(~S_opt,~S_opt),1,'sm');%只会在补集选择，每改变一个的变化
omega = abs(omega)^(1/k);


