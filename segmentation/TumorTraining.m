function [] = TumorTraining(Normal,Tumor,Tumor1, HighResolutionImage, HighResolutionImage1,TESTHighResolutionImage,TESTMask)
%function [] = TumorTrainingandTesting(Normal,Tumor, HighResolutionImage, Normal1,Tumor1, HighResolutionImage1,TESTHighResolutionImage,TESTMask)
tic

HighResolutionImage=rgb2gray(HighResolutionImage);
[nrowHR, ncolHR]=size(HighResolutionImage);
 HighResolutionImage1=rgb2gray(HighResolutionImage1);
 [nrowHR1, ncolHR1]=size(HighResolutionImage1);
TESTHighResolutionImage=rgb2gray(TESTHighResolutionImage);
[nrowTHR, ncolTHR]=size(TESTHighResolutionImage);

        NORMAL=zeros(nrowTHR,ncolTHR);
       TUMOR=zeros(nrowTHR,ncolTHR);

ProcessedHighResolution=HighResolutionImage;
 ProcessedHighResolution1=HighResolutionImage1;
ProcessedHighResolutiontest=TESTHighResolutionImage;
% 
sigma=2.;%3
KernelSize = 2*round(2*sigma)+1;
K = fspecial('gaussian', KernelSize, sigma); 

     ProcessedHighResolutionImage = conv2(ProcessedHighResolution,K,'same');   
      ProcessedHighResolutionImage1 = conv2(ProcessedHighResolution1,K,'same');     
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
                                 Are=measures(k).Area;
                                 if(Are > 5)
                                    f1  =mean(ProcessedHighResolutionImage((ismember(bw, k))==1));
                                    % AV=mmomALI(M1,2);%  AVabs(255-f2)abs(255-f2) AV 
                                     %LEN=momALI(M1);%AV LENLEN    Are./momALI(M1)
                                                FVnormal=[FVnormal  ;f1 Are./mmomALI(M1,2) ];
                                                index=[index; 1];
                                 end
                        end
                   %%
                          FS=[FVnormal index];
                          [row c]=size(FS);
                           A=FS(:,1:end-1);%normc(A);%
                           Anorm = A;%(A - min(min(A)))/(max(max(A)) - min(min(A)));
% %                          %                   [idx,C] = kmeans(Anorm,2);
                           [neuronsNormal]=SOMM(Anorm,3,3,10000,0.1)
 %                         neuronsNormal=mean(neuronsNormal,1)
%                                          [bw1, numberOfitems] = bwlabel(Normal1);
%                 measures = regionprops(bw1, 'Area','Image');
%                         FVnormal1=[];
%                         index1=[];
%                         for k = 1 : numberOfitems % Loop through  all items.
%                                  M1=measures(k).Image();
%                                  Are=measures(k).Area;
%                                  if(Are > 5)
%                                  %    f1  =mean(ProcessedHighResolutionImage1((ismember(bw1, k))==1));
%                                      AV=mmomALI(M1,1);%  AVabs(255-f2)abs(255-f2) AV 
%                                     LEN=momALI(M1);%AV LENLEN   abs(255-f1)
%                                                 FVnormal1=[FVnormal1  ;  AV mmomALI(M1,2) ];
%                                                 index1=[index1; 1];
%                                  end
%                         end
%                    %%
%                         FS=[FVnormal index;FVnormal1 index1];
%                         [row c]=size(FS);
%                          A=FS(:,1:end-1);%normc(A);%
%                          Anorm = A;%(A - min(min(A)))/(max(max(A)) - min(min(A)));
%                          %                   [idx,C] = kmeans(Anorm,2);
%                          [neuronsNormal]=SOMM(Anorm,3,3,10000,0.1)
                        % neuronsNormal=mean(neuronsNormal,1)
                         %neuronsNormal=[neuronsNormal;neuronsNormal1]
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
                            if(Are > 5)
                                   f1  =mean(ProcessedHighResolutionImage((ismember(bbw, k))==1));
                                         % AV=mmomALI(M1,2);%  AVabs(255-f2)abs(255-f2)AV  
                                          % LEN=momALI(M1);%AV LEN    mmomALI(M1,2) Are./momALI(M1)
                                                    FVtumor=[FVtumor  ;f1 Are./mmomALI(M1,2)];
                                                    index=[index; 2];
                            end
                        end
                   %%
                         FS=[FVtumor index];
                         [row c]=size(FS);
                          A=FS(:,1:end-1);%normc(A);%
                          Anorm = A;%(A - min(min(A)))/(max(max(A)) - min(min(A)));
%                          %                   [idx,C] = kmeans(Anorm,2);
                          [neuronsTumor]=SOMM(Anorm,3,3,10000,0.1);
%                            neuronsTumor=mean(neuronsTumor,1)
                                        [bbw1, numberOfitems] = bwlabel(Tumor1);
                 measures = regionprops(bbw1, 'Area','Image');
                         FVtumor1=[];
                         index1=[];
                         for k = 1 : numberOfitems % Loop through  all items.
                             M1=measures(k).Image();
                             Are=measures(k).Area;
                             if(Are > 5)
                                     f1  =mean(ProcessedHighResolutionImage1((ismember(bbw1, k))==1));
%                                           AV=mmomALI(M1,1);%  AVabs(255-f2)abs(255-f2)AV  
%                                           LEN=momALI(M1); %AV LEN LEN   abs(255-f1)
                                                     FVtumor1=[FVtumor1  ;  f1 Are./mmomALI(M1,2)];
                                                     index1=[index1; 2];
                             end
                         end
%                    %%
                         FS=[FVtumor1 index1];
                         [row c]=size(FS);
                          A=FS(:,1:end-1);%normc(A);%
                          Anorm = A;%(A - min(min(A)))/(max(max(A)) - min(min(A)));
%                          %                   [idx,C] = kmeans(Anorm,2);
                          [neuronsTumor1]=SOMM(Anorm,3,3,10000,0.1)
                        % neuronsTumor=mean(neuronsTumor,1)
                        neuronsTumor=[neuronsTumor;neuronsTumor1]
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
                            if(Are > 5)
                                    f1  =mean(ProcessedHighResolutionIma((ismember(bww, k))==1));
                                       %   AV=mmomALI(M1,2);%  AVabs(255-f2)abs(255-f2)  AV 
                                         %  LEN=momALI(M1);%AV LEN   mmomALI(M1,2) Are./momALI(M1)
                                                    FVunkown=[FVunkown  ;f1 Are./mmomALI(M1,2) ];
                                                    index=[index; k];
                            end
                        end
                   %%
                        FS=[FVunkown index];
                        [row c]=size(FS);
                         A=FS(:,1:end-1);%%normc(A);%
                         Anorm = A;%(A - min(min(A)))/(max(max(A)) - min(min(A)));
                         %                   [idx,C] = kmeans(Anorm,2);
                       %  [neurons]=SOMM(Anorm,3,3,10000,0.1);
                       size(neuronsNormal);
                         idx=[];
                             for i = 1:row
                                    vector = Anorm(i,1:c-1);
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
toc                              