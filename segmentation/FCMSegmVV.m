function [Epitheluim_bw, Stroma_bw]= FCMSegmVV(II)
%II=imread(filename);

%% create a te


%%
Img = double(rgb2gray(II));
rangefilt_channel=rangefilt(Img);
% contrastfilt_channel=Contrastfilter (Img,9);
[aa1 aa2] = size(Img);
%Img1 = imresize(Img,.1);
Img1=Img;
[a1 a2] = size(Img1);
XR=a1;
YR=a2;
%P=zeros(a1,a2);
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
%I1=(I1+I2+I3)/3;
%I2=rangefilt_channel;
%I3=Img;
%%

for maxClassnum=2:3
            [centers,U,of] = fcm([I1(:) I2(:) I3(:)],maxClassnum);
            error(maxClassnum)=min(of);
            %
            centers=mean(centers,2);
            %
            [q , indexx]=min(centers);
            
            EP1(:,:,maxClassnum)=reshape(U(indexx,:),a1,a2,1);
            %
            centers=mean(centers,2);
            %
            A2 = sort(centers);
            out = A2(2);
            maxIndex = ismember(centers,out);
            ST1(:,:,maxClassnum)=reshape(U(maxIndex,:),a1,a2,1);
           %
           centers=mean(centers,2);
           %
           [q , indexx]=max(centers);
            Back(:,:,maxClassnum)=reshape(U(indexx,:),a1,a2,1);
end
 [v in]=min(error(3:3));
 NumOfClasses=in+1
 Epitheluim(:,:)=EP1(:,:,NumOfClasses);
 Stroma(:,:)=ST1(:,:,NumOfClasses);
  Background(:,:)=Back(:,:,NumOfClasses);
 %figure(10), imshow(Background,[]);
 %figure(1), imshow(Epitheluim,[]);
 %%
 Ep=Epitheluim-Stroma-Background;
 Epitheluim_bw=(Ep>0);
 %imtool(Epitheluim_bw);
 %figure(2), imshow(Epitheluim_bw,[]);
 
 %%
  St=Stroma-Epitheluim-Background;
 Stroma_bw=(St>0);
 %imtool(Stroma_bw);
 %figure(3), imshow(Stroma_bw,[]);
 %%
  Ba=Background-Stroma-Epitheluim;
 Background_bw=(Ba>0);
 %imtool(Background_bw);
 %figure(4), imshow(BA_bw,[]);
 %%
 %figure(2), imshow(Stroma,[]);

   %  threshold = graythresh(Epitheluim);
   %  Epitheluim_bw = im2bw(Epitheluim,threshold);
 %figure(3), imshow(Epitheluim_bw,[]);
 %threshold = graythresh(Stroma);
 %Stroma_bw = im2bw(Stroma,threshold);
 %figure(4), imshow(Stroma_bw,[]);

end