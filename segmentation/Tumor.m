function [] = TumorTrainingandTesting(Normal,Tumor, HighResolutionImage,TESTHighResolutionImage,TESTMask)
%tic
%[Mask, Stroma]=FCMSegm('MMR-5-MSH2_2015-05-08_21.12.05_x2.5_z0.tif');
%[LowResolutionImage, Stroma]=FCMSegmVV('MMR-3-MSH2_2015-05-08_19.25.43_x0.625_z0.tif');
%MMR-3-MSH2_2015-05-08_19.25.43_x10_z0.tif
%MMR-3-MSH2_2015-05-08_19.25.43_x0.625_z0.tif
%MMR-5-MSH2_2015-05-08_21.12.05_x0.625_z0
%MMR-6-MLH1_2015-05-08_22.09.09_x0.625_z0
%MMR-9-MLH1_2015-05-09_01.31.18_x0.625_z0
%MMR-11-MLH1_2015-05-09_03.53.17_x0.625_z0
%LowResolutionImage=imread('Mask.bmp');
%[nrowLR, ncolLR]=size(LowResolutionImage);
%MMR-3-MSH2_2015-05-08_19.25.43_x2.5_z0.tif
%MMR-5-MSH2_2015-05-08_21.12.05_x2.5_z0
%MMR-6-MLH1_2015-05-08_22.09.09_x2.5_z0
%MMR-9-MLH1_2015-05-09_01.31.18_x2.5_z0
%MMR-11-MLH1_2015-05-09_03.53.17_x2.5_z0
%[Mask]=SPFLFCM('MMR-11-MLH1_2015-05-09_03.53.17_x2.5_z0.tif',100,0.5);
HighResolutionImage=rgb2gray(HighResolutionImage);
[nrowHR, ncolHR]=size(HighResolutionImage);
TESTHighResolutionImage=rgb2gray(TESTHighResolutionImage);
[nrowTHR, ncolTHR]=size(TESTHighResolutionImage);
%                  Mask=imresize(LowResolutionImage,[nrowHR, ncolHR]);
%                  Mask=imbinarize(double(Mask));
                 %%
%                  [labeled,N] = bwlabel(Mask,4);
%                     tempor = regionprops(labeled,'Area');
%                    idx = find([tempor.Area] > 30);
%                     bw = ismember(labeled,idx);
%                      Mask=bw>0; 
% % %                  
%                    imtool(Mask)
%Mask=imread('M3.bmp');
%imtool(Mask)
%       Unkown=zeros(nrowHR,ncolHR);
  %     Normal=zeros(nrowHR,ncolHR);
  %     Tumor=zeros(nrowHR,ncolHR);
        NORMAL=zeros(nrowTHR,ncolTHR);
       TUMOR=zeros(nrowTHR,ncolTHR);
%        Indicator2=zeros(nrowHR,ncolHR);
%        Indicator22=zeros(nrowHR,ncolHR);
%        Indicator222=zeros(nrowHR,ncolHR);
%        Indicator000=zeros(nrowHR,ncolHR);
ProcessedHighResolution=HighResolutionImage;
ProcessedHighResolutiontest=TESTHighResolutionImage;
% 
sigma=3.;%3
KernelSize = 2*round(2*sigma)+1;
K = fspecial('gaussian', KernelSize, sigma); 

     ProcessedHighResolutionImage = conv2(ProcessedHighResolution,K,'same');      
        ProcessedHighResolutionIma = conv2(ProcessedHighResolutiontest,K,'same');      
clc;    % Clear the command window.
close all;  % Close all figures (except those of imtool.)
workspace;  % Make sure the workspace panel is showing.
                [bw, numberOfitems] = bwlabel(Normal);
                measures = regionprops(bw, 'Area','Image');
                        FVnormal=[];
                        index=[];
                        for k = 1 : numberOfitems % Loop through  all items.
                                 M1=measures(k).Image();
                                 f1  =mean(ProcessedHighResolutionImage((ismember(bw, k))==1));
                                 AV=mmomALI(M1,2);%  AVabs(255-f2)abs(255-f2) momALI(M1)
                                 LEN=ALID(M1);
                                            FVnormal=[FVnormal  ; abs(255-f1)  AV LEN];
                                            index=[index; 1];
                        end
                   %%
                        FS=[FVnormal index];
                        [row c]=size(FS);
                         A=FS(:,1:end-1);
                         Anorm = A;%(A - min(min(A)))/(max(max(A)) - min(min(A)));
                         %                   [idx,C] = kmeans(Anorm,2);
                         [neuronsNormal]=SOM(Anorm,30,10000,0.1);
%                          idx=[];
%                              for i = 1:row
%                                     vector = A(i,1:c-1);
%                                      newmatrix = [neuronsNormal;vector];
%                                      distance = pdist(newmatrix);
%                                      s = squareform(distance);
%                                      [mindist, bmu] = getMinimum(s, size(s,1));
%                                      idx = [idx;bmu];
%                              end
                             
                [bbw, numberOfitems] = bwlabel(Tumor);
                measures = regionprops(bbw, 'Area','Image');
                        FVtumor=[];
                        index=[];
                        for k = 1 : numberOfitems % Loop through  all items.
                            M1=measures(k).Image();
                            Are=measures(k).Area;
                            f1  =mean(ProcessedHighResolutionImage((ismember(bbw, k))==1));
                                  AV=mmomALI(M1,2);%  AVabs(255-f2)abs(255-f2) momALI(M1)
                                    LEN=ALID(M1);
                                            FVtumor=[FVtumor  ; abs(255-f1)  AV LEN];
                                            index=[index; 2];
                        end
                   %%
                        FS=[FVtumor index];
                        [row c]=size(FS);
                         A=FS(:,1:end-1);
                         Anorm = A;%(A - min(min(A)))/(max(max(A)) - min(min(A)));
                         %                   [idx,C] = kmeans(Anorm,2);
                         [neuronsTumor]=SOM(Anorm,30,10000,0.1);
%                          idx=[];
%                              for i = 1:row
%                                     vector = A(i,1:c-1);
%                                      newmatrix = [neuronsTumor;vector];
%                                      distance = pdist(newmatrix);
%                                      s = squareform(distance);
%                                      [mindist, bmu] = getMinimum(s, size(s,1));
%                                      idx = [idx;bmu];
%                              end        
                             
                        %% TEST PHASE
               [bww, numberOfitems] = bwlabel(TESTMask);
                measures = regionprops(bww, 'Area','Image');
                        FVunkown=[];
                        index=[];
                        for k = 1 : numberOfitems % Loop through  all items.
                            M1=measures(k).Image();
                            Are=measures(k).Area;
                            f1  =mean(ProcessedHighResolutionIma((ismember(bww, k))==1));
                                  AV=mmomALI(M1,2);%  AVabs(255-f2)abs(255-f2) momALI(M1)
                                    LEN=ALID(M1);
                                            FVunkown=[FVunkown  ; abs(255-f1)  AV LEN];
                                            index=[index; k];
                        end
                   %%
                        FS=[FVunkown index];
                        [row c]=size(FS);
                         A=FS(:,1:end-1);
                         Anorm = A;%(A - min(min(A)))/(max(max(A)) - min(min(A)));
                         %                   [idx,C] = kmeans(Anorm,2);
                       %  [neuronsTumor]=SOM(Anorm,30,10000,0.1);
                         idx=[];
                             for i = 1:row
                                    vector = A(i,1:c-1);
                                     newmatrix1 = [neuronsNormal;vector];
                                     distance1 = pdist(newmatrix1);
                                     s1 = squareform(distance1);
                                     [mindist1, bmu1] = getMinimum(s1, size(s1,1));
                                      newmatrix2 = [neuronsTumor;vector];
                                     distance2 = pdist(newmatrix2);
                                     s2 = squareform(distance2);
                                     [mindist2, bmu2] = getMinimum(s2, size(s2,1));
                                     if(mindist1 < mindist2)
                                         idx = [idx;1];
                                     else
                                         idx = [idx;2];
                                     end
                             end        
                             [q1 w1]=find(idx==1);
                            q11=index(q1);
                        [q2 w2]=find(idx==2);
                            q22=index(q2);
                              NORMAL((ismember(bww, q11))==1)=1; %merge
                              TUMOR((ismember(bww, q22))==1)=1; %merge
                              imtool(NORMAL)
                              imtool(TUMOR)
                              