function [Tumor] = TumorDetection(Mask, HighResolutionImage,Postprocessing,Postprocessingmiddle)
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
       Unkown=zeros(nrowHR,ncolHR);
       Normal=zeros(nrowHR,ncolHR);
       Tumor=zeros(nrowHR,ncolHR);
       Indicator2=zeros(nrowHR,ncolHR);
       Indicator22=zeros(nrowHR,ncolHR);
       Indicator222=zeros(nrowHR,ncolHR);
       Indicator000=zeros(nrowHR,ncolHR);
ProcessedHighResolution=HighResolutionImage;
% 
sigma=3.;%3
KernelSize = 2*round(2*sigma)+1;
K = fspecial('gaussian', KernelSize, sigma); 
   sigma=21;%11 3
   KernelSize = 2*round(2*sigma)+1;
  KK = fspecial('gaussian', KernelSize, sigma); 
     ProcessedHighResolutionImage = conv2(ProcessedHighResolution,K,'same');      
        ProcessedHighResolutionIma = conv2(ProcessedHighResolution,KK,'same');      
clc;    % Clear the command window.
close all;  % Close all figures (except those of imtool.)
workspace;  % Make sure the workspace panel is showing.
                [bw, numberOfitems] = bwlabel(Mask);
                measures = regionprops(bw, 'Area','Image');
                        FVnormal=[];
                        index=[];
                        for k = 1 : numberOfitems % Loop through  all items.
                            M1=measures(k).Image();
                            Are=measures(k).Area;
                         %   M=mmomALI(M1,2);
                            
                           if (Are > 25) % remove noise 19.8	7.0373	1.6	0.3425	1.4 %29.9255	15.6	2.3413	0.5	2.0173
%61.0812	35.1751	4.7265	0.7141	3.9006% 146.723625	24.288	10.62765	1.94785	0.3995	1.6935

                                    f1  =mean(ProcessedHighResolutionImage((ismember(bw, k))==1));
                              %      f2  =mean(ProcessedHighResolutionIma((ismember(bw, k))==1));
                              % if(M(3)<1.94785 && M(4)<0.3995 && M(5)<1.6935)     abs(255-f2) ALID(M1) ALIR(M1,1)ALIR(M1,2)
                                  AV=mmomALI(M1,2);%  AVabs(255-f2)abs(255-f2) momALI(M1)
                              %    if(any(AV))
                                    LEN=ALID(M1);
                                            FVnormal=[FVnormal  ; abs(255-f1)  AV LEN];
                                            index=[index; k];
                               %   end
                                            

                           end
                        end
                   %%
                        FS=[FVnormal index];
                        [row c]=size(FS);
                         A=FS(:,1:end-1);
                         Anorm = A;%(A - min(min(A)))/(max(max(A)) - min(min(A)));
                         %                   [idx,C] = kmeans(Anorm,2);
                         [neurons]=SOM(Anorm,3,10000,0.1)
                         idx=[];
                             for i = 1:row
                                    vector = A(i,1:c-1);
                                     newmatrix = [neurons;vector];
                                     distance = pdist(newmatrix);
                                     s = squareform(distance);
                                     [mindist, bmu] = getMinimum(s, size(s,1));
                                     idx = [idx;bmu];
                             end
                      if(length(unique(idx))==3)
                          NEWDATA=neurons(idx,:);
                             size(NEWDATA);
                             [iidx,C] = kmeans(NEWDATA,3);
                          [q3 w3]=find(iidx==3);
                            q33=index(q3);
                           Test(3)= sum((neurons(3,2:4)));%5:7
                            [q2 w2]=find(iidx==2);
                            q22=index(q2);
                           Test(2)=sum((neurons(2,2:4)));
                            [q1 w1]=find(iidx==1);
                            q11=index(q1);
                            Test(1)=sum((neurons(1,2:4)));
                            [v inMin]=min(Test(:));
                            [vvv inMax]=max(Test(:));
                            ATest=[1;2;3];
                            ATest(ATest==inMin)=[];
                            ATest(ATest==inMax)=[];
                            BET=ATest(1);
                            ALDImin=mean(neurons(:,5))
                       v=mean([neurons(inMin,2:4)./ALDImin])
                            vv=mean([neurons(BET,2:4)./ALDImin])
                       
                        vvv=mean([neurons(inMax,2:4)./ALDImin])
                           [q00 w3]=find(idx==inMin);
                            q000=index(q00);
                            [q11 w2]=find(idx==inMax);
                            q111=index(q11);
                            [q22 w2]=find(idx==BET);
                            q222=index(q22);
                              %Indicator000((ismember(bw, q000))==1)=1; %merge  Unkown
                          %    neurons(inMin,2)
                          %    vv=mean(neurons(inMin,5:7))
                            %  vvv=mean(neurons(inMax,5:7))
                            %  vv/(vv+vvv)
                          %    vv=neurons(inMin,5)/neurons(inMin,6)
                          %    if(vv<100)
                                Normal((ismember(bw, q000))==1)=1; %merge
                                %%
                             %%
                                 [bww, numberOfitems] = bwlabel(Normal);
                                    measures = regionprops(bww, 'Image','Area');
                                            FVnormal=[];
                                            index=[];
                                            for k = 1 : numberOfitems % Loop through  all items.
                                                M1=measures(k).Image();
                                              %  Are=measures(k).Area;
                                              %  if(Are>70)
                                                                FVnormal=[FVnormal  ;   ALIDD(M1)];
                                              %  end               %index=[index; k];
                                            end
                                if(max(mean(FVnormal(:,1)), mean(FVnormal(:,2)))/min(mean(FVnormal(:,1)), mean(FVnormal(:,2))) > 3.5)
                                    Postprocessingmiddle=1;
                                end

   %                           else
   %                                Tumor((ismember(bw, q000))==1)=1; %merge
   %                           end
                              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                              Unkown((ismember(bw, q222))==1)=1; %merge
                              Indicator2((ismember(bw, q111))==1)=1; %merge
                             % Tumor((ismember(bw, q222))==1)=1; %merge
                             % Tumor((ismember(bw, q111))==1)=1; %merge
                            Normal=Normal+Unkown;
                             Tumor=Indicator2;
                               imtool((Normal))
%                            imtool((Indicator000))
  %                             imtool((Unkown))
  %                           imtool((Indicator2))
                             %Unkown=Tumor;
                            imtool((Tumor))

                   %%
                      elseif(length(unique(idx))>1)
                             NEWDATA=neurons(idx,:);
                             size(NEWDATA);
                             [iidx,C] = kmeans(NEWDATA,2);
                            [q2 w]=find(iidx==2);
                            q22=index(q2);
                            Test(2)=sum(C(2,:));
                            [q1 w]=find(iidx==1);
                            q11=index(q1);
                            Test(1)= sum(C(1,:));
                            [v inMin]=min(Test(:));
                            [vv inMax]=max(Test(:));
                            v
                            vv
%                            Indicator1((ismember(bw, q11))==1)=1; %merge
%                            Indicator1((ismember(bw, q22))==1)=2; %merge
                            [q00 w3]=find(iidx==inMin);
                            q000=index(q00);
                            [q11 w2]=find(iidx==inMax);
                            q111=index(q11);
                            Normal((ismember(bw, q000))==1)=1; %merge
                            Tumor((ismember(bw, q111))==1)=1; %merge
                            imtool((Normal))
                          imtool((Tumor))
                     else
                          
                             Indicator2((ismember(bw, index))==1)=1; %merge
                             imtool(Indicator2)
                      end
                      
                   if(Postprocessingmiddle && ~Postprocessing) 
                       %Normal=Normal+Unkown;
                                     [bbw, numberOfitems] = bwlabel(Normal);
                                    measures = regionprops(bbw, 'Image');
                                            %FVnormal=[];
                                            index=[];
                                            for k = 1 : numberOfitems % Loop through  all items.
                                                M1=measures(k).Image();
                                               % Are=measures(k).Area;
                                               %IND=abs(ALIR(M1,2)- ALIR(M1,1));IND <vvv+1
                                               IND=mean([mmomALI(M1,2)./ALID(M1)]);
                                         %      IND=ALIC(M1,2);
                                           %    IND=ALID(M1);
                                                if ( IND < v ) % remove noiseALID(M1) nmomALI(M1,2)  ALI(M1,2) mmomALI(M1,2)momALI(M1)
                                                                %FVnormal=[FVnormal  ;   Are];
                                                                index=[index; k];


                                               end
                                            end
                                %%
%                          FS=[FVnormal index];
%                         [row c]=size(FS);
%                          A=FS(:,1:end-1);
%                          Anorm = A;%(A - min(min(A)))/(max(max(A)) - min(min(A)));
%                         [idx,C] = kmeans(Anorm,2);
%                         Test=zeros(2,1);
%                             [q2 w]=find(idx==2);
%                             q22=index(q2);
%                             Test(2)=sum(C(2,:));
%                             [q1 w]=find(idx==1);
%                             q11=index(q1);
%                             Test(1)= sum(C(1,:));
%                             [v inMin]=min(Test(:));
%                            
%                             [vv inMax]=max(Test(:));
%                             
%                             [q00 w3]=find(idx==inMin);
%                             q000=index(q00);
%                             [q11 w2]=find(idx==inMax);
%                             q111=index(q11);
                            Indicator000((ismember(bbw, index))==1)=1; %merge
                          %  Indicator000((ismember(bbw, q111))==1)=1; %merge
                            %Indicator222((ismember(bbw, q111))==1)=1; %merge
                          %  Indicator222((ismember(bw, index))==1)=1; %merge
                          %  Normal=Normal+Indicator222;
                          %  Tumor=Tumor-Indicator222;
                          Indicator222=Normal-Indicator000;
                          Normal=Normal-Indicator222;
                          Tumor=Indicator2+Indicator222;
                            imtool(Normal)
                            imtool(Tumor)
                    end     
                    if(Postprocessing && ~Postprocessingmiddle) 
                                     [bbw, numberOfitems] = bwlabel(Normal);
                                    measures = regionprops(bbw, 'Image');
                      %                      FVnormal=[];
                                            index=[];
                                            for k = 1 : numberOfitems % Loop through  all items.
                                                M1=measures(k).Image();
                     %                           Are=measures(k).Area;
                                                IND=ALIR(M1,2);%;abs(ALIR(M1,2)- ALIR(M1,1));
                                                if (IND <.2 ) % remove noiseALID(M1) nmomALI(M1,2)  ALI(M1,2) mmomALI(M1,2)momALI(M1)
                     %                                           FVnormal=[FVnormal  ;   Are];
                                                                index=[index; k];


                                               end
                                            end
                                %%
%                          FS=[FVnormal index];
%                         [row c]=size(FS);
%                          A=FS(:,1:end-1);
%                          Anorm = A;%(A - min(min(A)))/(max(max(A)) - min(min(A)));
%                         [idx,C] = kmeans(Anorm,2);
%                         Test=zeros(2,1);
%                             [q2 w]=find(idx==2);
%                             q22=index(q2);
%                             Test(2)=sum(C(2,:));
%                             [q1 w]=find(idx==1);
%                             q11=index(q1);
%                             Test(1)= sum(C(1,:));
%                             [v inMin]=min(Test(:));
%                            
%                             [vv inMax]=max(Test(:));
%                             
%                             [q00 w3]=find(idx==inMin);
%                             q000=index(q00);
%                             [q11 w2]=find(idx==inMax);
%                             q111=index(q11);
                            Indicator000((ismember(bbw, index))==1)=1; %merge
                          %  Indicator000((ismember(bbw, q111))==1)=1; %merge
                            %Indicator222((ismember(bbw, q111))==1)=1; %merge
                          %  Indicator222((ismember(bw, index))==1)=1; %merge
                          %  Normal=Normal+Indicator222;
                          %  Tumor=Tumor-Indicator222;
                          Indicator222=Normal-Indicator000;
                          Normal=Indicator000;
                          Tumor=Tumor+Indicator222;
                            imtool(Normal)
                            imtool(Tumor)
                    end     
                     %%
                                
                                
%                     
%       
%                         
%                 [bw, numberOfitems] = bwlabel(Normal);
%                 measures = regionprops(bw, 'Image','Area');
%                         FVnormal=[];
%                         index=[];
%                         for k = 1 : numberOfitems % Loop through  all items.
%                             M1=measures(k).Image();
%                             Are=measures(k).Area;
%                            % if (Are >60 ) % remove noiseALID(M1) nmomALI(M1,2)  ALI(M1,2) mmomALI(M1,2)momALI(M1)
%                                             FVnormal=[FVnormal  ;   momALI(M1) ];
%                                             index=[index; k];
%                                       
% 
%                            %end
%                         end
%                        
%                         FS=[FVnormal index];
%                         [row c]=size(FS);
%                          A=FS(:,1:end-1);
%                          Anorm = A;%(A - min(min(A)))/(max(max(A)) - min(min(A)));
%                          [neurons]=SOM(Anorm,3,10000,0.1)
%                          idx=[];
%                              for i = 1:row
%                                     vector = A(i,1:c-1);
%                                      newmatrix = [neurons;vector];
%                                      distance = pdist(newmatrix);
%                                      s = squareform(distance);
%                                      [mindist, bmu] = getMinimum(s, size(s,1));
%                                      idx = [idx;bmu];
%                              end
%                       if(length(unique(idx))>1)
%                              NEWDATA=neurons(idx,:);
%                              size(NEWDATA);
%                              [iidx,C] = kmeans(NEWDATA,2);
%                              Test=zeros(2,1);
%                              [q2 w]=find(iidx==2);
%                             q22=index(q2);
%                             Test(2)=sum(C(2,:));
%                             [q1 w]=find(iidx==1);
%                             q11=index(q1);
%                             Test(1)= sum(C(1,:));
%                             [v, inMin]=min(Test(:));
%                             [v, inMax]=max(Test(:));
%                         
% %                            Indicator1((ismember(bw, q11))==1)=1; %merge
% %                            Indicator1((ismember(bw, q22))==1)=2; %merge
%                             [q00 w3]=find(iidx==inMin);
%                             q000=index(q00);
%                             [q11 w2]=find(iidx==inMax);
%                             q111=index(q11);
%                             Normal((ismember(bw, q000))==1)=1; %merge
%                             Indicator22((ismember(bw, q111))==1)=1; %merge
%                             Normal=Normal+Indicator000;
%                             imtool(Normal)
%                            % imtool(Indicator22)
%                             Tumor=Indicator22+Indicator2;
%                            imtool(Tumor)
% 
%                             
%                       else
%                           
%                              Normal((ismember(bw, index))==1)=1; %merge
%                              imtool(Normal)
%                       end
%                      %%PreProcessing Area Filtering
%                      
%                      [bw, numberOfitems] = bwlabel(Normal);
%                       measures = regionprops(bw, 'Area');
%                       FVnormal=[];
%                       index=[];
%                        for k = 1 : numberOfitems % Loop through  all items.
%                             Are=measures(k).Area;
%                            % M1=measures(k).Image();
%                             FVnormal=[FVnormal  ;ALID(M1)];
%                             index=[index; k];
%                         end
%                        LowerSize= min((FVnormal))
%                        UpperSize= max((FVnormal))
%                                 
% %                      [bw, numberOfitems] = bwlabel(Tumor);
% %                       measures = regionprops(bw, 'Area');
% %                       FVnormal=[];
% %                       index=[];
% %                        for k = 1 : numberOfitems % Loop through  all items.
% %                            Are=measures(k).Area;
% %                            %Are=ALID(measures(k).Image());
% %                               if(Are<UpperSize && Are>LowerSize)
% %                                             FVnormal=[FVnormal  ;Are];
% %                                             index=[index; k];
% %                               end
% %                        end
%                         FS=[FVnormal index];
%                         [row c]=size(FS);
%                          A=FS(:,1:end-1);
%                          Anorm = A;%(A - min(min(A)))/(max(max(A)) - min(min(A)));
%                         [idx,C] = kmeans(Anorm,2);
%                         Test=zeros(2,1);
%                             [q2 w]=find(idx==2);
%                             q22=index(q2);
%                             Test(2)=sum(C(2,:));
%                             [q1 w]=find(idx==1);
%                             q11=index(q1);
%                             Test(1)= sum(C(1,:));
%                             [v inMin]=min(Test(:));
%                             v
%                             [vv inMax]=max(Test(:));
%                             vv
%                             [q00 w3]=find(iidx==inMin);
%                             q000=index(q00);
%                             [q11 w2]=find(iidx==inMax);
%                             q111=index(q11);
%                             Indicator000((ismember(bw, q000))==1)=1; %merge
%                             Indicator222((ismember(bw, q111))==1)=1; %merge
%                           %  Indicator222((ismember(bw, index))==1)=1; %merge
%                           %  Normal=Normal+Indicator222;
%                           %  Tumor=Tumor-Indicator222;
%                             imtool(Indicator000)
%                             imtool(Indicator222)
%                             
%                      %%

 %toc
% end