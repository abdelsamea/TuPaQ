function [Epitheluim_bw, Stroma_bw]= FCMSegm(filename)
tic
II=imread(filename);

%% create a te
% LAB=rgb2lab(II);
% L=LAB(:,:,1);
% A=LAB(:,:,2);
% B=LAB(:,:,3);
%  IA=(L+A+B)/2;
% %%
% Img = double(rgb2gray(II));
% [aa1 aa2] = size(Img);
II = imresize(II,.4);
%Img1=Img;

% XR=a1;
% YR=a2;
%P=zeros(a1,a2);

%%
%%
im_red = II;
im_green = II;
im_blue = II;

% Red channel only
im_red(:,:,2) = 0; 
im_red(:,:,3) = 0; 
I1 = double(rgb2gray(im_red));
% Green channel only
im_green(:,:,1) = 0; 
im_green(:,:,3) = 0; 
I2 = double(rgb2gray(im_green));
% Blue channel only
im_blue(:,:,1) = 0; 
im_blue(:,:,2) = 0; 
I3 = double(rgb2gray(im_blue));
%%
IA=(I1+I2+I3)/3;
[a1, a2 ] = size(IA);
%tic
for maxClassnum=3:3
            [centers,U,of] = fcm(IA(:),maxClassnum);
            error(maxClassnum)=min(of);
            [q , indexx]=min(centers);
            
            EP1(:,:,maxClassnum)=reshape(U(indexx,:),a1,a2,1);
            A2 = sort(centers);
            out = A2(2);
            maxIndex = ismember(centers,out);
            ST1(:,:,maxClassnum)=reshape(U(maxIndex,:),a1,a2,1);
          [q , indexx]=max(centers);
            Back(:,:,maxClassnum)=reshape(U(indexx,:),a1,a2,1);
end
 [v in]=min(error(3:3));
 NumOfClasses=3
 Epitheluim(:,:)=EP1(:,:,NumOfClasses);
 Stroma(:,:)=ST1(:,:,NumOfClasses);
 %figure(1), imshow(Epitheluim,[]);
 %figure(2), imshow(Stroma,[]);
 Background(:,:)=Back(:,:,NumOfClasses);
 %figure(10), imshow(Background,[]);
     threshold = graythresh(Epitheluim);
     Epitheluim_bw = im2bw(Epitheluim,threshold);
 %figure(3), imshow(Epitheluim_bw,[]);
 threshold = graythresh(Stroma);
 Stroma_bw = im2bw(Stroma,threshold);
 %figure(4), imshow(Stroma_bw,[]);
 toc
end