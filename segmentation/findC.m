tic
Unkown=imread('test.bmp');
%imtool(Unkown)
%Unkown= imfill(Unkown,'holes');
%imtool(Unkown)
%[centers, radiiBright] = imfindcircles(Unkown,[18/2 56/2],'Sensitivity',0.92);
[centers, radiiBright] = imfindcircles(Unkown,[5 10],'ObjectPolarity','bright','Sensitivity',0.92);
length(centers)
% Group=repmat(1, length(centers), 1);
% %%
%                    [bw, numberOfitems] = bwlabel(Unkown);
%                          measures = regionprops(bw, 'Centroid','Area','Image','Extrema');
%                         Test=[];
%                         index=[];
%                         for k = 1 : numberOfitems % Loop through  all items.
%                        %     M1=measures(k).Image();
%                              index=[index; k];
%                              c=measures(k).Centroid;
%                              %EX=measures(k).Extrema;
%                            % M1=measures(k).Image();
%                                % CO=corner(M1,2);
%                                 %f1  =mean(ProcessedHighResolutionImage((ismember(bw, k))==1));
%                              Test=[Test;  c(1) c(2)];
%                         end
%                          Class = knnclassify(Test, centers, Group, 1);
%                           [q2 w]=find(Class==2);
%                              q22=index(q2);
%%

figure, imshow(Unkown)
%      hold;
%plot(centers(:,1), centers(:,2),'r*')

hBright = viscircles(centers, radiiBright,'Color','b');
toc