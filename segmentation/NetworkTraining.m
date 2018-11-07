function [Anorm] = NetworkTraining(Normal, HighResolutionImage, mapWidth,mapHeight,startLearningRate)
%imtool(Normal)
%imtool(HighResolutionImage)
%HighResolutionImage=rgb2gray(HighResolutionImage);
[nrowHR, ncolHR]=size(HighResolutionImage);
ProcessedHighResolution=HighResolutionImage;
Ref=mean(mean(HighResolutionImage(Normal >0)));
% 
sigma=2.;%3
KernelSize = 2*round(2*sigma)+1;
K = fspecial('gaussian', KernelSize, sigma); 
     ProcessedHighResolutionImage = conv2(ProcessedHighResolution,K,'same');   

                [bw, numberOfitems] = bwlabel(Normal);
                measures = regionprops(bw, 'Area','Image');
                        FVnormal=[];
                      %  index=[];
                        for k = 1 : numberOfitems % Loop through  all items.
                                 M1=measures(k).Image();
                                 Are=measures(k).Area;
                                 if(Are > 50)
                                    f1  =mean(ProcessedHighResolutionImage((ismember(bw, k))==1));
                                   %    f2  =median(ProcessedHighResolutionImage((ismember(bw, k))==1));
                                   %    f3  =var(ProcessedHighResolutionImage((ismember(bw, k))==1));
                                   %  f4  =skewness(ProcessedHighResolutionImage((ismember(bw, k))==1));
                                   %   f5  =kurtosis(ProcessedHighResolutionImage((ismember(bw, k))==1));
% % %                                     %f6  =moment(ProcessedHighResolutionIma((ismember(bww, k))==1),2);
                                     %  f=[abs(255-f1) abs(255-f2) abs(255-f3) abs(255-f4) abs(255-f5)];
                                       f=[f1 ];%f2 f3 f4 f5];
                                         %FVnormal=[FVnormal  ;(100*f)./Ref  Are Are./mmomALI(M1,1) Are./mmomALI(M1,2) mmomALI(M1,1) mmomALI(M1,2)];
                                         %FVnormal=[FVnormal  ;f  Are mmomALI(M1,2)];
                                         FVnormal=[FVnormal  ;f  mmomALI(M1,2)];
                                              %  index=[index; 1];
                                 end
                        end
                   %%
                          FS=[FVnormal];
                          [row c]=size(FS);
                           A=FS(:,1:end);%zscore(A);%normc(A);%
                           Anorm = A;%normc(A);%(A - min(min(A)))/(max(max(A)) - min(min(A)));
% %                          %                   [idx,C] = kmeans(Anorm,2);
                          
                               