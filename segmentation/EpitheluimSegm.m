function [Epitheluim_bw, Stroma]= EpitheluimSegm(II,NumOfClasses,SIGMA)
%tic
%II=imread(filename);
%II = imresize(II,.4);
%imtool(rgb2gray(II))
Epitheluim_RGB=II;
%%
im_red = II;
im_green = II;
im_blue = II;
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
% %  sigma=.1;%3
% %  KernelSize = 2*round(2*sigma)+1;
% %  
% % K = fspecial('gaussian', KernelSize, sigma); 
% %    
% %  %      ProcessedHighResolutionImage1 = conv2(ProcessedHighResolution1,K,'same');     
% %         IA = conv2(double(rgb2gray(Epitheluim_RGB)),K,'same');  
% IA=double(rgb2gray(Epitheluim_RGB));

[a1, a2 ] = size(IA);
%tic
NumOfClasses1=3;
if NumOfClasses == 3
    for maxClassnum=NumOfClasses:NumOfClasses%
            [centers,U,of] = fcm(IA(:),maxClassnum,[NaN 100 0.00001 0]);
            error(maxClassnum)=min(of);
            [q , indexx]=min(centers);
            
            EP1(:,:,maxClassnum)=reshape(U(indexx,:),a1,a2,1);
            A2 = sort(centers);
            out = A2(2);
            maxIndex = ismember(centers,out);
            ST1(:,:,maxClassnum)=reshape(U(maxIndex,:),a1,a2,1);
            Stroma(:,:)=ST1(:,:,NumOfClasses);
          [q , indexx]=max(centers);
            Back(:,:,maxClassnum)=reshape(U(indexx,:),a1,a2,1);
    end
else
        for maxClassnum=NumOfClasses1:NumOfClasses1%
            [centers,U,of] = fcm(IA(:),maxClassnum,[NaN 100 0.00001 0]);
            error(maxClassnum)=min(of);
            [q , indexx]=min(centers);
            
            EP1(:,:,maxClassnum)=reshape(U(indexx,:),a1,a2,1);
            A2 = sort(centers);
            out = A2(2);
            maxIndex = ismember(centers,out);
            ST1(:,:,maxClassnum)=reshape(U(maxIndex,:),a1,a2,1);
            Stroma(:,:)=ST1(:,:,NumOfClasses1);
            [q , indexx]=max(centers);
            Back(:,:,maxClassnum)=reshape(U(indexx,:),a1,a2,1);
        end
        for maxClassnum=NumOfClasses:NumOfClasses%
            [centers,U,of] = fcm(IA(:),maxClassnum,[NaN 100 0.00001 0]);
            error(maxClassnum)=min(of);
            [q , indexx]=min(centers);
            
            EP1(:,:,maxClassnum)=reshape(U(indexx,:),a1,a2,1);
           % A2 = sort(centers);
           % out = A2(2);
           % maxIndex = ismember(centers,out);
           % ST1(:,:,maxClassnum)=reshape(U(maxIndex,:),a1,a2,1);
           % Stroma(:,:)=ST1(:,:,NumOfClasses);
            %[q , indexx]=max(centers);
            %Back(:,:,maxClassnum)=reshape(U(indexx,:),a1,a2,1);
        end
    
end
 %[v in]=min(error(3:3));
 
 Epitheluim(:,:)=EP1(:,:,NumOfClasses);
 %if NumOfClasses == 3
 %Stroma(:,:)=ST1(:,:,NumOfClasses);
%  Background(:,:)=Back(:,:,NumOfClasses);
     threshold = graythresh(Epitheluim);
     Epitheluim_bw = im2bw(Epitheluim,threshold);
   threshold = graythresh(Stroma);
  Stroma = im2bw(Stroma,threshold);
if(SIGMA == 0)
else
    Epitheluim_bw= SPFLFCM(Epitheluim_bw,Epitheluim, Epitheluim_RGB,5000,SIGMA);
end
                 %%
                 [labeled,N] = bwlabel(Epitheluim_bw,4);
                    tempor = regionprops(labeled,'Area');
                   idx = find([tempor.Area] > 20);
                    bw = ismember(labeled,idx);
                     Epitheluim_bw=bw>0; 
                     
% imtool(Epitheluim_RGB)
% 
% imtool(Epitheluim_bw);
% imtool(Stroma);
% toc
end