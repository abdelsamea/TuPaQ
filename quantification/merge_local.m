function [bw]=merge_local(Bimage, Indicator)
%%
     %[Ind, numberOfIndicatoritems]=bwlabel(Indicator);
     %measures = regionprops(Ind, 'Area');
lines=Indicator;
unmerge=Bimage;
merge=unmerge+lines;
[r c]=size(merge);
[Region, numberOfregions]=bwlabel(merge);
%s = regionprops(Region, 'Eccentricity','Solidity','MinorAxisLength','MajorAxisLength');
% phi = linspace(0,2*pi,50);
% cosphi = cos(phi);
% sinphi = sin(phi);
NewIndicator=zeros(r,c);
%%
for i=1:numberOfregions
    temp=zeros(r,c);
    
    temp(Region==i)=1;
    %figure(1), imshow(temp);
    subregions=double(temp) & double(unmerge);
    sublines=double(temp) & double(lines);
    %figure(2), imshow(subregions);
    %figure(3), imshow(sublines);
    
    %xbar = s(i).Centroid(1);
    %ybar = s(i).Centroid(2);
    % a = s(i).MajorAxisLength/2;
    % b = s(i).MinorAxisLength/2;
%     theta = pi*s(i).Orientation/180;
%     R = [ cos(theta)   sin(theta)
%          -sin(theta)   cos(theta)];
% 
%     xy = [a*cosphi; b*sinphi];
%     xy = R*xy;
% 
%     x = xy(1,:) + xbar;
%     y = xy(2,:) + ybar;
    
    %err = immse([xbar ybar],[a b]);
  
    if (nnz(sublines)>0 )
      %  if(s(i).Eccentricity < .9)
         [Regio, numberOfEl]=bwlabel(subregions);
         m = regionprops(Regio, 'MinorAxisLength','MajorAxisLength');
             [Ind, numberOfIndicatoritems]=bwlabel(sublines);
             measures = regionprops(Ind, 'Area');
             %if (s(i).Eccentricity < mean([m(1:numberOfEl).Eccentricity]))
                AcceptedMinLength=max([m(1:numberOfEl).MinorAxisLength]);%-min([m(1:numberOfEl).MinorAxisLength])
                AcceptedMaxLength=max([m(1:numberOfEl).MajorAxisLength]);%+min([m(1:numberOfEl).MinorAxisLength])
                sublines(ismember(Ind, find(AcceptedMinLength/2 > [measures(1:numberOfIndicatoritems).Area]  ))==1)=0; %split
                sublines(ismember(Ind, find(AcceptedMaxLength/2 < [measures(1:numberOfIndicatoritems).Area]  ))==1)=0; %split
                % figure(4), imshow(sublines);
                NewIndicator=double(NewIndicator) + double(sublines);
               % figure(5), imshow(NewIndicator);
       %  end
        %end
    end

    
end
%%





%imtool(Indicator)
%[r c]=size(Bimage);
   % Indicator=double(Indicator);
  %  whos Indicator
    %figure(99), imshow(Indicator,[]);
    %%
     %[Regio, numberOfEl]=bwlabel(Bimage);
     %m = regionprops(Regio, 'MinorAxisLength','MajorAxisLength');
%     [Reg, numberOfElj]=bwlabel(Bimage);
 %    mj = regionprops(Reg, 'MajorAxisLength');
%     AcceptedMinLength= median([m(1:numberOfEl).MinorAxisLength]);
     
     %[Ind, numberOfIndicatoritems]=bwlabel(Indicator);
     %measures = regionprops(Ind, 'Area');
    


       %         AcceptedMinLength=mean([m(1:numberOfEl).MinorAxisLength]);%-min([m(1:numberOfEl).MinorAxisLength])
       %         AcceptedMaxLength=mean([m(1:numberOfEl).MajorAxisLength]);%+min([m(1:numberOfEl).MinorAxisLength])
       %         Indicator(ismember(Ind, find(AcceptedMinLength/2 > [measures(1:numberOfIndicatoritems).Area]  ))==1)=0; %split
       %         Indicator(ismember(Ind, find(AcceptedMaxLength/2 < [measures(1:numberOfIndicatoritems).Area]  ))==1)=0; %split
 
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




UndividedImage = double(Bimage) | double(NewIndicator);
%figure(11), imshow(UndividedImage,[]);
bw=UndividedImage;
end