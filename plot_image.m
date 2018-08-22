clear;
close all;
clc;

%% data and required toolboxes
    load('usps_all');
    b = zeros(16,16,12);%储存每个数字的像素图
    c = data(:,1:10);%构建3*4维的图
    for num = 1:10        %总共显示12个数字
        a = data(:,1,num);%data的第一维为256*1的像素列向量
            for n = 1:16
                b(n,:,num) = a((1+(n-1)*16):16*n,:);%将256*1转化为16*16的像素图，并存在b中
            end
    end     

    for num = 7:8        %总共显示12个数字
        a = data(:,1,num);%data的第一维为256*1的像素列向量
            for n = 1:16
                b(n,:,num+4) = a((1+(n-1)*16):16*n,:);%将256*1转化为16*16的像素图，并存在b中
            end
    end    
    for j = 1:12
        d = b(:,:,j);
        subplot(3,4,j);
        imshow(d');
    end

