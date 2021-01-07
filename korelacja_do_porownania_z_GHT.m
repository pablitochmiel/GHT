close all; clear; clc;

image_rgb=imread('shape_moj.png');
image=rgb2gray(image_rgb);
%test
imageEdg=edge(image);         %nie wywala
%imshow(imageEdg);
shape=imread('eksperyment.png');

shape=(shape(:,:,1)==0);

tic;
res=filter2(shape,imageEdg);
time=toc;
res=res./max(max(res));
[y,x]=find(res==1);
s=size(x,1);
[a,b,~]=size(image);
[c,d]=size(shape);
temp=zeros([a,b,3],'uint8');
temp(:,:,:)=image_rgb(:,:,:);
for k=1:s
    for i=1:c
        for j=1:d
            if (shape(i,j))
                temp(i+y(k)-round(0.5*c),j+x(k)-round(0.5*d),1)=255;
                temp(i+y(k)-round(0.5*c),j+x(k)-round(0.5*d),2)=0;
                temp(i+y(k)-round(0.5*c),j+x(k)-round(0.5*d),3)=0;
            end
        end
    end
end
Text="found "+string(s)+' shapes, time: '+string(time)
imshow(temp);
%imtool(res)
 
