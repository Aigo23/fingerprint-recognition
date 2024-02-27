clc;
clear;
close all;

%读取指纹库的图片并得到特征点
folder = 'data1\file\img\';
files = dir(fullfile(folder, '*.png')); %set中图片的顺序与files.name相对应

img_set = cell(1,numel(files));
for i=1:length(img_set)
    filename = strcat(folder,num2str(i),'.png');
    img_r = imread(filename);
    img_set{i} = img_r;
    % getIcn(img_r, i);
    % disp(i);
end


%检验原图上标记特征点
show = 1;
if show
for i=1:length(img_set)
    fig = figure(1);
    imshow(img_set{i});
    minu = LoadNeuFeature(strcat('IcnOut\Icn',num2str(i),'.txt'));
    minu(:,3) = -minu(:,3);
    DrawMinu(fig,minu);
end
end

% img1 =imread("data1\file\img\8.png");
% img2 = imread("data1\search\img\2.png");
% 
% minu1 = LoadNeuFeature('IcnOut\Icn8.txt');
% minu2 = LoadNeuFeature('IcnSearch\Icn.txt');
% 
% [match_pts1, match_pts2,PTS]  = match2(img1,minu1,img2,minu2);
% disp(PTS);
% 
% figure(3);
% showMatchedFeatures(img1, img2, match_pts1, match_pts2, 'montage');








