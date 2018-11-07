function [ProcessedHighResolutionn, TUMOR, TESTMask, Stroma] = TumorTesting(neuronsNormal,neuronsTumour,TESTHighResolutionImage,NumOfClasses,SIGMA,flag)

TESTHighResolutionImageImage=TESTHighResolutionImage;
if(flag ==1)

  load neuronsNormal.mat
  load neuronsTumour.mat 

end

[TESTMask, Stroma]= EpitheluimSegm(TESTHighResolutionImage,NumOfClasses,SIGMA);


TESTHighResolutionImage=rgb2gray(TESTHighResolutionImage);
[nrowTHR, ncolTHR]=size(TESTHighResolutionImage);

        NORMAL=zeros(nrowTHR,ncolTHR);
       TUMOR=zeros(nrowTHR,ncolTHR);
       MASK=zeros(nrowTHR,ncolTHR);

% ProcessedHighResolution1=HighResolutionImage1;
ProcessedHighResolutiontest=TESTHighResolutionImage;
Ref=mean(mean(TESTHighResolutionImage(TESTMask >0)));
sigma=2.;%2
KernelSize = 2*round(2*sigma)+1;
K = fspecial('gaussian', KernelSize, sigma); 
       
        ProcessedHighResolutionIma = conv2(ProcessedHighResolutiontest,K,'same');      


                        %% TEST PHASE
               [bww, numberOfitems] = bwlabel(TESTMask);
                measures = regionprops(bww, 'Area','Image');
                        FVunkown=[];
                        index=[];
                        for k = 1 : numberOfitems % Loop through  all items.
                            M1=measures(k).Image();
                            Are=measures(k).Area;
                            if(Are > 50)
                                    f1  =mean(ProcessedHighResolutionIma((ismember(bww, k))==1));
                                     % f2  =median(ProcessedHighResolutionIma((ismember(bww, k))==1));
                                     %f3  =var(ProcessedHighResolutionIma((ismember(bww, k))==1));
                                     % f4  =skewness(ProcessedHighResolutionIma((ismember(bww, k))==1));
                                     % f5  =kurtosis(ProcessedHighResolutionIma((ismember(bww, k))==1));
% % %                                     %f6  =moment(ProcessedHighResolutionIma((ismember(bww, k))==1),2);
                                       %f=[abs(255-f1) abs(255-f2) abs(255-f3) abs(255-f4) abs(255-f5)];
                                       f=[f1];% f2 f3 f4 f5];
                                         %FVunkown=[FVunkown  ;(100*f)./Ref  Are Are./mmomALI(M1,1) Are./mmomALI(M1,2) mmomALI(M1,1) mmomALI(M1,2)];
                                        % FVunkown=[FVunkown  ;f Are mmomALI(M1,2)];
                                         FVunkown=[FVunkown  ;f mmomALI(M1,2)];
                                                    index=[index; k];
                            end
                        end
                   %%
                        FS=[FVunkown index];
                        [row c]=size(FS);
                         A=FS(:,1:end-1);%%zscore(A);%normc(A);%
                         Anorm = A;%normc(A);%(A - min(min(A)))/(max(max(A)) - min(min(A)));
                         %                   [idx,C] = kmeans(Anorm,2);
                       %  [neurons]=SOMM(Anorm,10,10,10000,0.1);
                       size(neuronsNormal);
                         idx=[];
                             for i = 1:row
                                    vector = Anorm(i,1:c-1);
                                     newmatrix1 = [neuronsNormal;vector];
                                     distance1 = pdist(newmatrix1);
                                     s1 = squareform(distance1);
                                     [mindist1, bmu1] = getMinimum(s1, size(s1,1));
                                      newmatrix2 = [neuronsTumour;vector];
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
 
ProcessedHighResolutionn=TESTHighResolutionImageImage;%imread(TESTHighResolutionImagefilename);%TESTHighResolutionImage;%
%imtool(TUMOR)
    ProcessedHighResolutionn(TUMOR>0)=0;
 
                             