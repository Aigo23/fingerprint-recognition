% clc;
% clear;
% close all;

%用于给库中指纹提取特征点

function [] = getIcn(img, ID)

show = 0;
% img = imread("data1\file\img\2.png");
% ID = 2; 
[H,W]= size(img);

if show 
    figure(1),imshow(img),title("original"); end

[newimg,binimg,mask,orientim] = testfin(img);

if show
    figure(2),subplot(2,3,1),imshow(binimg),title("binary");end

se = strel('square', 2);
im = imopen(binimg, se);
im = imclose(im, se);

if show
    figure(2),subplot(2,3,2),imshow(im),title("ridge-smooth");end

im = bwmorph(im,'thin',inf);
if show
    figure(2),subplot(2,3,3),imshow(im),title("ridge-thin");end

%弥补断裂
se = strel('square',3);
im = imdilate(im,se);
if show
figure(2),subplot(2,3,4),imshow(im),title("ridge-dilated");end
im = bwmorph(im, 'thin', inf);
if show
figure(2),subplot(2,3,5),imshow(im),title("ridge-thin-again");end

%桥接目前没想到好的处理方法，决定在特征点提取的过程中去减少桥接带来的影响，让距离的较近的特征点只保留一个
%去除短线
im = bwareaopen(im, 30);
if show
    figure(2),subplot(2,3,6),imshow(im),title("ridge-delsmall");end

%去除毛刺
im = Pruning(im, 5);
if show
    figure(3),subplot(2,3,1),imshow(im),title("ridge-prunned");end

%细节点检测,Icn为细节点标记阵
I = im;
Icn=zeros(size(I));
C = I;
for i=2:H-1
    for j=2:W-1
        if C(i,j) == 1
            cn=abs(C(i-1,j)-C(i-1,j-1))+abs(C(i-1,j+1)-C(i-1,j))+abs(C(i,j+1)-C(i-1,j+1))...
            +abs(C(i+1,j+1)-C(i,j+1))+abs(C(i+1,j)-C(i+1,j+1))+abs(C(i+1,j-1)-C(i+1,j))...
            +abs(C(i,j-1)-C(i+1,j-1))+abs(C(i-1,j-1)-C(i,j-1));
            Icn(i,j)=cn/2;
        end
    end
end

%标注
if show
[row,col]=find(Icn==1);
hold on, plot(col,row,'gs','MarkerSize',5)
[row,col]=find(Icn==3);
hold on, plot(col,row,'rs','MarkerSize',5)
end

%细节点筛选，筛除距离较近的细节点
if show
    figure(3),subplot(2,3,2),imshow(im),title("cn-2");end

ind = find(Icn == 1 | Icn == 3);

[rowind, colind] = ind2sub(size(Icn), ind);

rowind_new = [];
colind_new = [];
for i = 1:length(colind)
    flag = 0;
    for j = 1:length(colind_new)
        dis = [rowind(i) - rowind_new(j), colind(i)-colind_new(j)];
        distance = norm(dis);
        if(distance < 16)
            flag = 1;
            break;
        end
    end
    if(flag == 1)
        continue;
    end
    rowind_new(end+1) = rowind(i);
    colind_new(end+1) = colind(i);
end

Icn_new = zeros(size(Icn));%不区分两种细节点
for i = 1:length(rowind_new)
    Icn_new(rowind_new(i),colind_new(i)) = 1;
end
%标注
if show
[row,col]=find(Icn_new==1);
hold on, plot(col,row,'gs','MarkerSize',5)
end

%去除mask边缘的干扰细节点
se = strel('disk',10);
maskforedge = ~imdilate(~mask,se);
if show
    figure(3),subplot(2,3,3),imshow(~maskforedge),title("edge");end
Icn_new(maskforedge == 0) = 0;
if show
    figure(3),subplot(2,3,3),imshow(im),title("final-Icn");
    [row,col]=find(Icn_new==1);hold on, plot(col,row,'gs','MarkerSize',5)
end

%orient场的方向定义为0-pi，参考坐标系，x右，y下
filename = strcat('IcnOut\Icn',num2str(ID),'.txt');
fileID = fopen(filename,'w');
fprintf(fileID,'');
ind = find(Icn_new > 0);
[rowind, colind] = ind2sub(size(Icn_new),ind);
for i = 1:length(rowind)
    pt = [rowind(i),colind(i)];
    or = pi - orientim(pt(1),pt(2));
move = 3;bk_c1 = [fix(pt(1)- move*sin(or)),fix(pt(2)+move*cos(or))];
bk_c2 = [fix(pt(1)+move*sin(or)),fix(pt(2)-move*cos(or))];
bk1 = im(bk_c1(1)-move:bk_c1(1)+move,bk_c1(2)-move:bk_c1(2)+move);
bk2 = im(bk_c2(1)-move:bk_c2(1)+move,bk_c2(2)-move:bk_c2(2)+move);
if(sum(bk1(:))<sum(bk2(:)))
    or = or + pi;
end
fprintf(fileID, '%d %d %d\n', colind(i), rowind(i), or*180/pi);
end
fclose(fileID);

show1 = 0;
if show1
    fig = figure;imshow(im);
minu = LoadNeuFeature('IcnOut\Icn2.txt');
minu(:,3) = -minu(:,3);
DrawMinu(fig,minu);
end





%剪枝算法
function C = Pruning(A, len)

B = CreateEndpointSE();
X1 = A;
for k = 1:len
    endpoints = false(size(A));
    for m = 1:size(B,1)
        endpoints = endpoints | bwhitmiss(X1, B{m,1}, B{m,2});
    end
    X1(endpoints) = 0;
end

X2 = false(size(A));
for m = 1:size(B,1)
    endpoints = bwhitmiss(X1, B{m,1}, B{m,2});
    X2(endpoints) = 1;
end

se = strel(ones(3,3));
X3 = X2;
for k = 1:len
    X3 = imdilate(X3, se) & A;
end

C = X3 | X1; 
end
%构造8种毛刺se
function B = CreateEndpointSE()
B{1,1} = [0 0 0; 1 1 0; 0 0 0];
B{1,2} = [0 1 1; 0 0 1; 0 1 1];
for k = 2:4
    B{k,1} = rot90(B{k-1,1});
    B{k,2} = rot90(B{k-1,2});
end
B{5,1} = [1 0 0; 0 1 0; 0 0 0];
B{5,2} = [0 1 1; 1 0 1; 1 1 1];
for k = 6:8
    B{k,1} = rot90(B{k-1,1});
    B{k,2} = rot90(B{k-1,2});
end
end

end