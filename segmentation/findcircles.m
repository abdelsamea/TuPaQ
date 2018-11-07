function [] = DemoClassification()
tic
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
HighResolutionImage=rgb2gray(imread('MMR-3-MSH2_2015-05-08_19.25.43_x2.5_z0.tif'));
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
Mask=imread('UN.bmp');
% Maskk=imresize(Stroma,[nrowHR, ncolHR]);
% Maskk=imbinarize(double(Maskk));
 imtool(Mask)
[x,y] = find (Mask) ;
CenterOfMassXY = [mean(x) mean(y)] ;
       Indicator0=zeros(nrowHR,ncolHR);
       Indicator1=zeros(nrowHR,ncolHR);
       Indicator2=zeros(nrowHR,ncolHR);
%%  
%ProcessedHighResolution(Mask==0)=0; 
%ProcessedHighResolutionImage=medfilt2(HighResolutionImage,[5 5]);
  ProcessedHighResolution=HighResolutionImage;
% 
sigma=1;
KernelSize = 2*round(2*sigma)+1;
K = fspecial('gaussian', KernelSize, sigma); 
     ProcessedHighResolutionImage = conv2(ProcessedHighResolution,K,'same');
% ProcessedHighResolutionImage = imgaussfilt(ProcessedHighResolution,10);  
%    sigma=15;
% KernelSize = 2*round(2*sigma)+1;
% K = fspecial('gaussian', KernelSize, sigma); 
%    ProcessedHighResolutionImageh = conv2(ProcessedHighResolution,K,'same');
%ProcessedHighResolutionImage=HighResolutionImage;
%%
 %
 

rgbImage=Mask;

%%
% Demo to divide a color image up into blocks.
clc;    % Clear the command window.
close all;  % Close all figures (except those of imtool.)
workspace;  % Make sure the workspace panel is showing.
fontSize = 30;
%%

%                         x = ones(nrowHR,1)*[1:ncolHR];    % Matrix with each pixel set to its x coordinate
%                         y = [1:nrowHR]'*ones(1,ncolHR);   % Matrix with each pixel set to its y coordinate
%                     area = sum(sum(Mask));
%                     meanx = sum(sum(double(Mask).*x))/area;
%                     meany = sum(sum(double(Mask).*y))/area;
%                     centroid = [meanx, meany];
                [bw, numberOfitems] = bwlabel(Mask);
                measures = regionprops(bw, 'all');
                  c1=0;
                  c2=0;
                        FVnormal=[];
                        index=[];
                        for k = 1 : numberOfitems % Loop through  all items.
                            M1=measures(k).Image();
                            %  figure(k*10), imshow(M1,[]);
                            Are=measures(k).Area;
                            if (Are >100 ) % remove noise
                                % M1= imfill(M1,'holes');
                                       % [FV] = momentALI(M1,2)
                                       %  FVnormal=[FVnormal  ;momentALI(M1,2)];
                                       %%
                                        c=measures(k).Centroid;
                                        %Are=measures(k).Area;
                                       % distC = pdist([c;CenterOfMassXY],'euclidean');
                                      %  distC=sqrt(sum((CenterOfMassXY-c).^2,2));
%                                            A=(centroid-c).^2;
%                                             dis= sqrt(sum(A,2));
                                       %%
                                      % averageIntensity=mean((ProcessedHighResolutionImage((ismember(bw, k))==1)));
                                  
                          %      f1  =mean(ProcessedHighResolutionImage((ismember(bw, k))==1));
                                %   f0   =mean(hampel(double(ProcessedHighResolution((ismember(bw, k))==1))));
%momALI(M1)
                                %      [f0,x0]=ksdensity([momALI(M1);m0(:)]);
                                %segm=5 4*f1 mmomALI(M1,2) DALI'/2 mmomALI(M1,2) mmomALI(M1,1) nmomALI(M1,2) nmomALI(M1,1)
                               
                                
                                        c1=c1+1;
                                      
                                  %    NNM=[nnmomALI(M1,2)];
                                  %    NNM(NNM<1)=0;
                                  %     if(all(NNM))
%                                            NM1=[ALI(M1,2) ];
%                                            NM2=[ALI(M1,1)];
%                                            NM=mean([NM1;NM2],1);
%                                    NM(NM<2)=0;
%                                    NM(NM>12)=0;
                                   %NM(abs(mean(NM1)-mean(NM2))>)
                                 %  if(mean(NM1)<mean(NM2)/2 ||  mean(NM2)<mean(NM1)/2)
                                 %      NM=0;
                                 %  end
%                                  NM2=nmomALI(M1,2);
%                                  NM1=nmomALI(M1,1);
%                                  NM=mean([NM1; NM2],1);
%                                  NM(NM>20)=0;
%                                   if(NM>5)
%                                       NM=0;
%                                   end
                                [centers, radiiBright] = imfindcircles(M1,[4 18],'ObjectPolarity', ...
                                'bright','Sensitivity',0.92,'Method','TwoStage');
                              if(~isempty(centers))    %  NM2  NM1   momALI(M1)
%                                 R1=mean(ALIR(M1,1));
%                                 R2=mean(ALIR(M1,2));
%                                      mi=min(R1,R2);
%                                      ma=max(R1,R2);
%                                    FR=(ma-mi)/ma ;% Are  nmomALI(M1,2)  ALI(M1,2) mmomALI(M1,2)
                                            FVnormal=[FVnormal  ;Are];
                                            index=[index; k];
                                            c2=c2+1;
                               end
                                       
                                       
                                       %FVnormal=[FVnormal abs(BB(2,1)- BB(2,2))];
                                      %  FVnormal=[FVnormal max([momentALI(M1,1) momentALI(M1,2)])];
                                  %  M=measures(k).Image();
                                 %   figure(k), imshow(M,[]);



                            end
                        end
                        %c1
                        %c2
                        
                        FS=[FVnormal index];
                        [row c]=size(FS);
                        %%
                       % min((FVnormal))
                     %   max((FVnormal))
                         A=FS(:,1:end-1);
                       %  Anorm= normc(A);
                           
                         Anorm = A;%(A - min(min(A)))/(max(max(A)) - min(min(A)));
%                   [idx,C] = kmeans(Anorm,2);
                        %%

                         [neurons]=SOM(Anorm,3,10000,0.1);
                         idx=[];
                             for i = 1:row
                                    vector = A(i,1:c-1);
                                     newmatrix = [neurons;vector];
                                     distance = pdist(newmatrix);
                                     s = squareform(distance);
                                     [mindist, bmu] = getMinimum(s, size(s,1));
                                     idx = [idx;bmu];
                             end
                      if(length(unique(idx))==33)
                          NEWDATA=neurons(idx,:);
                             size(NEWDATA);
                             [iidx,C] = kmeans(NEWDATA,3);
                             
                          %%
                        %[C, idx]=Cluster;
                            [q3 w3]=find(iidx==3);
                            q33=index(q3);
                           Test(3)= sum((neurons(3,:)));
                            [q2 w2]=find(iidx==2);
                            q22=index(q2);
                           Test(2)=sum((neurons(2,:)));
                            [q1 w1]=find(iidx==1);
                            q11=index(q1);
                            Test(1)=sum((neurons(1,:)));
                            [v inMin]=min(Test(:));
                            [v inMax]=max(Test(:));
                            ATest=[1;2;3];
                            ATest(ATest==inMin)=[];
                            ATest(ATest==inMax)=[];
                            BET=ATest(1);
                          
                           % Indicator1((ismember(bw, q11))==1)=1; %merge
                            Indicator1((ismember(bw, q33))==1)=1; %merge
                            Indicator2((ismember(bw, q22))==1)=1; %merge
                            Indicator0((ismember(bw, q11))==1)=1; %merge
                            
                            imtool((Indicator1))
                            imtool((Indicator2))
                            imtool((Indicator0))
                           [q00 w3]=find(idx==inMin);
                            q000=index(q00);
                            [q11 w2]=find(idx==inMax);
                            q111=index(q11);
                            [q22 w2]=find(idx==BET);
                            q222=index(q22);
                                   Indicator0=zeros(nrowHR,ncolHR);
                           
                               Indicator2=zeros(nrowHR,ncolHR);
                              Indicator0((ismember(bw, q000))==1)=1; %merge
                            Indicator2((ismember(bw, q111))==1)=1; %merge
                            imtool((Indicator0))
                            imtool((Indicator2))
%                             Indicator1((ismember(bw, q000))==1)=1; %merge
%                             Indicator1((ismember(bw, q111))==1)=2; %merge
%                             Indicator2=zeros(nrowHR,ncolHR);
%                             Indicator2((ismember(bw, q222))==1)=1; %merge
%                             %%
%                           %Indicator1= imfill(Indicator1,'holes');
%                                                      [bw, numberOfitems] = bwlabel(Indicator1);
%                          measures = regionprops(bw, 'Centroid','Area','Image');
%                         Group=[];
%                         Training=[];
%                         Tr=[];
%                         index=[];
%                         for k = 1 : numberOfitems % Loop through  all items.
%                             c=measures(k).Centroid;
%                           %    M1=measures(k).Image();
%                            %  index=[index; k];
% %                              M2= imfill(M1,'holes');
%                          %   Tr=[Tr;c(1) c(2)];
%                            
%                           %if (measures(k).Area > 200  ) % remove noise
%                                     Training=[Training; c(1) c(2)  measures(k).Area ];
%                                     Group=[Group;mode(Indicator1((ismember(bw, k))==1))];    
%                            %end
%                         end
%                        [bw, numberOfitems] = bwlabel(Indicator2);
%                         measures = regionprops(bw, 'Centroid','Area','Image');
%                        Tr=[];
%                        index=[];
%                        for k = 1 : numberOfitems % Loop through  all items.
%                          %  M1=measures(k).Image();
%                             index=[index; k];
%                             c=measures(k).Centroid;
% %                            M1=measures(k).Image();
%                             Tr=[Tr;c(1) c(2)  measures(k).Area ];
%                         end
%                         Sample=Tr;
%                          Class = knnclassify(Sample, Training, Group, 9);
%                          [q2 w]=find(Class==2);
%                             q22=index(q2);
%                             [q1 w]=find(Class==1);
%                             q11=index(q1);
%                                Indicator1=zeros(nrowHR,ncolHR);
%                                Indicator2=zeros(nrowHR,ncolHR);
%                             Indicator1((ismember(bw, q11))==1)=1; %merge
%                             Indicator2((ismember(bw, q22))==1)=1; %mergetumor
%                             imtool(Indicator1)
%                             imtool(Indicator2)
%                             %%
                      elseif(length(unique(idx))>1)
                             NEWDATA=neurons(idx,:);
                             size(NEWDATA);
                             [iidx,C] = kmeans(NEWDATA,2);
                             
                          %%
                        %[C, idx]=Cluster;
                            [q2 w]=find(iidx==2);
                            q22=index(q2);
                            Test(2)=sum(C(2,:));
                            [q1 w]=find(iidx==1);
                            q11=index(q1);
                            Test(1)= sum(C(1,:));
                            [v inMin]=min(Test(:));
                            [v inMax]=max(Test(:));
                            
                            Indicator1((ismember(bw, q11))==1)=1; %merge
                            Indicator1((ismember(bw, q22))==1)=2; %merge

                            [q00 w3]=find(iidx==inMin);
                            q000=index(q00);
                            [q11 w2]=find(iidx==inMax);
                            q111=index(q11);
                            Indicator0((ismember(bw, q000))==1)=1; %merge
                            Indicator2((ismember(bw, q111))==1)=1; %merge
                            imtool((Indicator0+Indicator2))
                        %   imtool((Indicator2))
%                         %%
%                             [bw, numberOfitems] = bwlabel(Indicator0);
%                          measures = regionprops(bw, 'Area');
%                         Areaa=[];
%                         for k = 1 : numberOfitems % Loop through  all items.
%                            % Are=Are+measures(k).Area;
%                             Areaa=[Areaa;measures(k).Area];
%                         end
%                         Are=mean(Areaa(:))-min(Areaa(:));%Are/numberOfitems;
%                         %%
%                             %%
%                         %    Indicator1= imfill(Indicator1,'holes');
%                            %% Indicator1 = imcomplement(imfill(imcomplement(Indicator1),'holes'));
%                          [bw, numberOfitems] = bwlabel(Indicator1);
%                          measures = regionprops(bw, 'Area','Centroid','Image');
%                         Group=[];
%                         Training=[];
%                         Tr=[];
%                         index=[];
%                         for k = 1 : numberOfitems % Loop through  all items.
%                             c=measures(k).Centroid;
%                             
%                              M1=measures(k).Image();
% %                             M2= imfill(M1,'holes');
%                          %   Tr=[Tr;c(1) c(2)];
%                            
%                           if (measures(k).Area > 250  ) % remove noise
%                              %  f1  =mean(ProcessedHighResolutionImage((ismember(bw, k))==1));
%                                     Training=[Training; c(1) c(2) measures(k).Area];
%                                     Group=[Group;mode(Indicator1((ismember(bw, k))==1))];    
%                           end
%                         end
%                         [bw, numberOfitems] = bwlabel(Mask);
%                          measures = regionprops(bw, 'Area','Centroid','Image');
%                         Tr=[];
%                        for k = 1 : numberOfitems % Loop through  all items.
%                             index=[index; k];
%                             c=measures(k).Centroid;
%                             M1=measures(k).Image();
%                          %  f1  =mean(ProcessedHighResolutionImage((ismember(bw, k))==1));
%                             Tr=[Tr; c(1) c(2) measures(k).Area];
%                         end
%                         Sample=Tr;
%                          Class = knnclassify(Sample, Training, Group, 13);
%                          [q2 w]=find(Class==2);
%                             q22=index(q2);
%                             [q1 w]=find(Class==1);
%                             q11=index(q1);
%                                Indicator1=zeros(nrowHR,ncolHR);
%                                Indicator2=zeros(nrowHR,ncolHR);
%                             Indicator1((ismember(bw, q11))==1)=1; %merge
%                             Indicator2((ismember(bw, q22))==1)=1; %merge
%                             imtool(Indicator1)
%                             imtool(Indicator2)
                            
                      else
                          
                             Indicator0((ismember(bw, index))==1)=1; %merge
                             imtool(Indicator0)
                     end
% %%
% % Read in a standard MATLAB color demo image.
% % folder = fullfile('D:\Segm-counting code');
% % baseFileName = 'MMR-3-MSH2_2015-05-08_19.25.43_x0.625_z0.tif';
% % % Get the full filename, with path prepended.
% % fullFileName = fullfile(folder, baseFileName);
% % if ~exist(fullFileName, 'file')
% % 	% Didn't find it there.  Check the search path for it.
% % 	fullFileName = baseFileName; % No path this time.
% % 	if ~exist(fullFileName, 'file')
% % 		% Still didn't find it.  Alert user.
% % 		errorMessage = sprintf('Error: %s does not exist.', fullFileName);
% % 		uiwait(warndlg(errorMessage));
% % 		return;
% % 	end
% % end
% % % Read the image from disk.
% % rgbImage = imread(fullFileName);
% % % Test code if you want to try it with a gray scale image.
% % % Uncomment line below if you want to see how it works with a gray scale image.
% % % rgbImage = rgb2gray(rgbImage);
% % % Display image full screen.
% imshow(rgbImage);
% % Enlarge figure to full screen.
% set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
% drawnow;
% % Get the dimensions of the image.  numberOfColorBands should be = 3.
% [rows columns numberOfColorBands] = size(rgbImage)
% %==========================================================================
% % The first way to divide an image up into blocks is by using mat2cell().
% blockSizeR = 700; % Rows in block.
% blockSizeC = 700; % Columns in block.
% % Figure out the size of each block in rows. 
% % Most will be blockSizeR but there may be a remainder amount of less than that.
% wholeBlockRows = floor(rows / blockSizeR);
% blockVectorR = [blockSizeR * ones(1, wholeBlockRows), rem(rows, blockSizeR)];
% % Figure out the size of each block in columns. 
% wholeBlockCols = floor(columns / blockSizeC);
% blockVectorC = [blockSizeC * ones(1, wholeBlockCols), rem(columns, blockSizeC)];
% % Create the cell array, ca.  
% % Each cell (except for the remainder cells at the end of the image)
% % in the array contains a blockSizeR by blockSizeC by 3 color array.
% % This line is where the image is actually divided up into blocks.
% if numberOfColorBands > 1
% 	% It's a color image.
% 	ca = mat2cell(rgbImage, blockVectorR, blockVectorC, numberOfColorBands);
% else
% 	ca = mat2cell(rgbImage, blockVectorR, blockVectorC);
% end
% % Now display all the blocks.
% plotIndex = 1;
% numPlotsR = size(ca, 1);
% numPlotsC = size(ca, 2);
%  counter=1;
% NumberOfCells=[];
% Position=[];
% Sizes=[];
% StepR=blockSizeR; 
% StepC=blockSizeC;
% 
% for r = 1 : numPlotsR
%     
% 	for c = 1 : numPlotsC
%         
% % 		fprintf('plotindex = %d', plotIndex, c, r);
% % 		% Specify the location for display of the image.
% % 		subplot(numPlotsR, numPlotsC, plotIndex);
% 		% Extract the numerical array out of the cell
% 		% just for tutorial purposes.
% 		rgbBlock = ca{r,c};
%         [pr,pc]=size(rgbBlock);
% %         imshow(rgbBlock);
%         center=size(rgbBlock)/2+.5;
%         Cx=center(1)+(r-1)*StepR;
%         Cy=center(2)+(c-1)*StepC;
%         Position = [Position; [Cy Cx]];
%         Sizes = [Sizes; [(c-1)*StepC (r-1)*StepR pc pr]];
%         
%         %%
%          Img = rgbBlock;
% %         if mean(mean(Img))==0
% %             number_of_objects=0;
% %    
% %              NumberOfCells=[NumberOfCells number_of_objects];
% %         else
% % 
% %             [number_of_objects bwBlock]=automated_cell_counting(Img);
% % 
% %             NumberOfCells=[NumberOfCells number_of_objects];
% %         end
%         %%
% %         length(Img>0)<200
% %         bw = logical(Img);
% %        measures = regionprops(bw, 'all');
% %     numberOfitems = size(measures, 1);
% %    
% %        for k = 1 : numberOfitems % Loop through all items.
% %             M1=measures(k).Image();
% %             %  figure(k*10), imshow(M1,[]);
% %             if (measures(k).Area > 1000 ) % remove noise
% %         
%             format bank;    
%         if (mean(mean(Img))==0 || length(Img>0)<500)
%             FVnormal=0;
%         else
%             %%
% %                 bw = logical(Img);
% %                 measures = regionprops(bw, 'all');
% %                     numberOfitems = size(measures, 1);
%  
% %                         FVnormal=[];
% %                         for k = 1 : numberOfitems % Loop through all items.
% %                             M1=measures(k).Image();
% %                             %  figure(k*10), imshow(M1,[]);
% %                             if (measures(k).Area > 500 ) % remove noise
% % 
% %                                         a=[momentALI(M1,1); momentALI(M1,2)];
% %                                         [BB,III] = sortrows(a,1);
% %                                            FVnormal=[FVnormal mean([momALI(M1)])];
% %                                         FVnormal=[FVnormal abs(BB(2,1)- BB(2,2))];
% %                                       %  FVnormal=[FVnormal max([momentALI(M1,1) momentALI(M1,2)])];
% %                                   %  M=measures(k).Image();
% %                                  %   figure(k), imshow(M,[]);
% % 
% % 
% % 
% %                             end
% %                         end
% % 
% %             %%
% %              FVnormal=mean(FVnormal);
% %              a=[momentALI(Img,1); momentALI(Img,2)];
% %              [BB,III] = sortrows(a,2);
% %              FVnormal=BB(1,2);%abs(BB(2,1)- BB(2,2));
%             %FVnormal=max([momentALI(Img,1) momentALI(Img,2)]);
%             FVnormal=momALI(Img)
%             fprintf('plotindex = %d', plotIndex);
%             % Specify the location for display of the image.
%             subplot(numPlotsR, numPlotsC, plotIndex);
%             % Extract the numerical array out of the cell
%             % just for tutorial purposes.
%             rgbBlock = ca{r,c};
%             [pr,pc]=size(rgbBlock);
%             imshow(rgbBlock);
%             		caption = sprintf(' %d',FVnormal);
%                 title(caption);
%                 drawnow;
%         end
% 
% % 		[rowsB columnsB numberOfColorBandsB] = size(rgbBlock);
% 		% Make the caption the block number.
% 
%                 %% save image rgbBlock
%        
% %         scounter = num2str(counter);
% %         FileName11=strcat(scounter,'.tif');
% %         saveDataName = fullfile(FileName11);
% %         imwrite(rgbBlock,saveDataName,'compression','none') ;
% %         counter=counter+1;
%         %%
% 		% Increment the subplot to the next location.
% 		plotIndex = plotIndex + 1;
% 	end
% end
% % Display the original image in the upper left.
% subplot(4, 6, 1);
% imshow(rgbImage);
% title('Original Image');
% %HighResolutionImageWithText = insertText(ProcessedHighResolutionImage,Position,NumberOfCells,'FontSize',30);
% %HighResolutionImageWithText = insertShape(HighResolutionImageWithText,'Rectangle',Sizes,'LineWidth',5);
% 
% %figure(9999999),imshow(HighResolutionImageWithText,[]);
% %sum(NumberOfCells)
% 
% %%
% 
% %%
% 
% %%
% % imagefiles = dir('*.tif');      
% % nfiles = length(imagefiles);    % Number of files found
% % info=[];
% % info1=[];
% % for ii=1:nfiles
% %    currentfilename = imagefiles(ii).name;
% %    [CPUtime, Noiteration, x, y] = SPFLFCM(currentfilename, 100, 1.5);
% %    info1=[ii x y CPUtime Noiteration];
% %    info=[info; info1];
% % end
% %   
% % save information.txt info -ASCII
 toc
% end