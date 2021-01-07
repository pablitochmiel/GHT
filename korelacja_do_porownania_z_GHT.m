close all; clear; clc;

image=imread('shape_moj.png');
image=rgb2gray(image);

imageEdg=edge(image);         %nie wywala
%imshow(imageEdg);
shape=imread('tr.png');

shape=(shape(:,:,1)==0);

res=filter2(shape,imageEdg);
res=res./max(max(res));
[y,x]=find(res==1);

imtool(res)
 
