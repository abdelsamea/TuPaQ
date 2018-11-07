function [bw]=merge(Bimage, Indicator)

unmerge=Bimage;
[UnRegion, numberOfUnregions]=bwlabel(unmerge);
uns = regionprops(UnRegion, 'Eccentricity');
unmerge(ismember(UnRegion, find(0.95> [uns(1:numberOfUnregions).Eccentricity]  ))==1)=0; % remove the region
%figure(333), imshow(unmerge,[]);
[UnRegion2, numberOfUnregions2]=bwlabel(unmerge);
unss = regionprops(UnRegion2,'MinorAxisLength','MajorAxisLength');


AcceptedMinLength=max([unss(1:numberOfUnregions2).MinorAxisLength]);%-min([m(1:numberOfEl).MinorAxisLength])
AcceptedMaxLength=max([unss(1:numberOfUnregions2).MajorAxisLength]);%+min([m(1:numberOfEl).MinorAxisLength])

   % Indicator=double(Indicator);
  %  whos Indicator
    %figure(99), imshow(Indicator,[]);
    %%
  %   [Regio, numberOfEl]=bwlabel(Bimage);
  %   m = regionprops(Regio, 'MinorAxisLength','MajorAxisLength');
%     [Reg, numberOfElj]=bwlabel(Bimage);
 %    mj = regionprops(Reg, 'MajorAxisLength');
%     AcceptedMinLength= median([m(1:numberOfEl).MinorAxisLength]);
     
     [Ind, numberOfIndicatoritems]=bwlabel(Indicator);
     measures = regionprops(Ind, 'Area');
    


              %  AcceptedMinLength=min([m(1:numberOfEl).MinorAxisLength]);%-min([m(1:numberOfEl).MinorAxisLength])
              %  AcceptedMaxLength=mean([m(1:numberOfEl).MajorAxisLength]);%+min([m(1:numberOfEl).MinorAxisLength])
                Indicator(ismember(Ind, find(AcceptedMinLength/2 > [measures(1:numberOfIndicatoritems).Area]  ))==1)=0; %split
                Indicator(ismember(Ind, find(AcceptedMaxLength/2 < [measures(1:numberOfIndicatoritems).Area]  ))==1)=0; %split
 
%
%      [index]=find([measures(1:numberOfIndicatoritems).Area] > 5);
% index0=ismember(Ind, index);
% Indicator(index0==1)=1;
% for k = 1:numberOfIndicatoritems
%       %find the two objects in L which are most close
%     % to the set of locations [r c]
%     if(measures(k).Area>12)
%         
%         Indicator(Ind==k)=1; % do not seperate
%     else
%         Indicator(Ind==k)=0;% seperate
%     end
% end




UndividedImage = Bimage+Indicator;
%figure(11), imshow(Bimage,[]);
bw=UndividedImage;
end