% MOMENTS
%
% Function calculates the moments of a binary image and returns 
% the centroid, the angle of axis of minimum inertia, and a measure 
% of 'roundness'.  The function assumes that there is only one object
% in the binary image.
%
% function [centroid, thetamin, roundness] = moments(im)
%
% Argument:  im   - a binary image containing values of 0 or 1
%
% Returns:   centroid  - a 2 element vector
%            thetamin  - the angle of axis of minimum inertia (radians)
%            roundness - ratio of minimum inertia/maximum inertia.
%
% Note that positive x is to the right and positive y is downwards
% thus angles are positive clockwise.
%
% The function also displays the image and overlays the position of
% the centroid and the axis of minimum inertia.

 function [FV] = nmomALI(im,flage)
%im=imread(filename); 
im=double(im);
%scale = .5;
%im = imresize(im,scale);
%sigma=1.5;
%KernelSize = 2*round(2*sigma)+1;
%G = fspecial('gaussian', KernelSize, sigma);
%  im = conv2(im, G, 'same'); % regularization
  
s=strel('disk',1,0);%Structuring element
F=imerode(im,s);%Erode the image by structuring element
%Bound = bwmorph(im,'skel',Inf);
% Bound = bwmorph(im,'remove');
% im=Bound;
Bound=im;%-F;

%Bound = bwperim(im,8);

%figure(100), imshow(Bound,[]);
    [rows,cols] = size(im);
    x = ones(rows,1)*[1:cols];    % Matrix with each pixel set to its x coordinate
    y = [1:rows]'*ones(1,cols);   % Matrix with each pixel set to its y coordinate

    area = sum(sum(im));
    meanx = sum(sum(double(im).*x))/area;
    meany = sum(sum(double(im).*y))/area;
    centroid = [meanx, meany];
    
    % coordinates changed with respect to centre of mass
    x = x - meanx;
    y = y - meany;
    
    a = sum(sum(double(im).*x.^2));
    b = sum(sum(double(im).*x.*y))*2;
    c = sum(sum(double(im).*y.^2));
    
    denom = b^2 + (a-c)^2;
    
    if denom == 0
        % let thetas equal arbitrary angles
        thetamin = 2*pi*rand;
        thetamax = 2*pi*rand;
        roundness = 1;
    else
        sin2thetamin = b/sqrt(denom);       %positive solution
        sin2thetamax = -sin2thetamin;
        cos2thetamin = (a-c)/sqrt(denom);   %positive solution
        cos2thetamax = -cos2thetamin;
    
        thetamin = atan2(sin2thetamin, cos2thetamin)/2;
        thetamax = atan2(sin2thetamax, cos2thetamax)/2;
        Imin = 0.5*(c+a) - 0.5*(a-c)*cos2thetamin - 0.5*b*sin2thetamin;
        Imax = 0.5*(c+a) - 0.5*(a-c)*cos2thetamax - 0.5*b*sin2thetamax;
        roundness = Imin/Imax;
    end
    
    % draw an axis proportional to object size
    % 0.5 takes into acount lines with roundness = 0
    % 5   takes into acount small objects, so axis is still visible.
     %rho = sqrt(area)/(roundness + 0.5) ; 
   rho = area;%/5;%(roundness + 0.5) ;
    [X1,Y1] = pol2cart(thetamin,      rho);
    [X2,Y2] = pol2cart(thetamin + pi, rho);
    %[X1 + meanx, X2 + meanx],[Y1 + meany, Y2 + meany]
    
 %   figure(222),imshow(Bound);title('Boundary extracted Image');
%    FileName11=strcat('Proposed', filename);
%saveDataName = fullfile(FileName11);
%imwrite(Bound,saveDataName);
[rrr, ccc]=find(Bound);
%rc = bwtraceboundary(im,[y1 x1],'N');
L=[ccc(:) rrr(:)];

%[Bound]=Binary_trans(im);
%entropy(Bound)
[r1, c1]=find(Bound);
%rc = bwtraceboundary(im,[y1 x1],'N');
%Location=[rc(:,2) rc(:,1)];


[r, c] = poly2cw(r1, c1);
Location=[c ,r];

%%
%figure(11),

l=1:length(c);
C=num2cell(l);
%plot(c,r)
%text(c,r,C)

%%
%figure(10),
%plot(c,r,'k.-')
%s = getcontourlines(Bound);
%Locatio=[c r];
Location1 = sortrows(Location,1);
Location2 = sortrows(Location,2);

x1 = Location2(end,1);
x2 = Location2(1,1);
y1 = Location2(end,2);
y2 = Location2(1,2);

%rc = bwtraceboundary(im,[y1 x1],'N');
%Location=[rc(:,2) rc(:,1)];
x11 = X1 + meanx;
x22 = X2 + meanx;
y11 = Y1 + meany;
y22 = Y2 + meany;

%    aD=([x1 y1]-[x2 y2]).^2;
%    D = sqrt(sum(aD,2));
n =1;  %10 segments
deltax = (x2 - x1)/n;
deltay = (y2 - y1)/n;
halfdx = deltax/2;
halfdy = deltay/2;
%radius = sqrt(halfdx.^2 + halfdy.^2);
xcents = x1 + halfdx + (0:n-1)*deltax;
ycents = y1 + halfdy + (0:n-1)*deltay;

     %line1 = [x11 y11 x22-x11 y22-y11];
     %x = projPointOnLine([xcents(:) ycents(:)], line1); %x's are point on ALI
%%
       m = (diff([Y1 + meany, Y2 + meany ])/diff([X1 + meanx, X2 + meanx]));
    % Slope of new line
    slope = -1/m;
    x = meanx;
    y = meany;
    lengthLine = rho;
    xLine = x-lengthLine:x+lengthLine;
    yLine = slope*(xLine-x) + y;
if(flage==1)
    line1 = [floor(xLine(1,1)) floor(yLine(1,1)) floor(xLine(1,floor(rho)))-floor(xLine(1,1)) floor(yLine(1,floor(rho)))-floor(yLine(1,1))];
elseif(flage==2)
    line1 = [x11 y11 x22-x11 y22-y11];
end
     [nrow, ncol]=size(Location);
     
     xPro = projPointOnLine([Location(:,1) Location(:,2)], line1); %x's are poit on ALI from the shape
     [nsmaples po] = size(xPro);

    DALI=sqrt((Location2(end,1)-Location2(1,1))^2+(Location2(end,2)-Location2(1,2))^2);
    



    xPro = projPointOnLine([Location(:,1) Location(:,2)], line1); %x's are poit on ALI from the shape
     Q=xPro;
    D=Q-Location;  

  c1=1;
  c2=1;

  Data1=[];
  Data2=[];
 D1=1;
 D2=1;
  for kk=1:length(Location)
      if(D(kk,1)>0 && D(kk,2)<0 || D(kk,1)>0 && D(kk,2)>0)
          Data1(c1,:)=Location(kk,:);
          c1=c1+1;
      elseif(D(kk,1)<0 && D(kk,2)>0 || D(kk,1)<0 && D(kk,2)<0)
          Data2(c2,:)=Location(kk,:);
          c2=c2+1;
 
      end
  end
  if(numel(Data1>6))
      [qqq, www]=poly2cw(Data1(:,1), Data1(:,2));
      Data1=[qqq, www];
  end
  if(numel(Data2>6))
      [aaa, sss]=poly2cw(Data2(:,1), Data2(:,2));
      Data2=[aaa, sss];
  end
   %%
LevelOfCloseness=.5;
%%
[l1, m1]=size(Data1);
if(l1<=1)
     E_distanceData1=[];
     E_distanceDataC1=[];
          E_distanceDataMin1=0;
     E_distanceDataax1=0;
     D1=0;
     M1=[];
else
xPro1 = projPointOnLine([Data1(:,1) Data1(:,2)], line1); %x's are poit on ALI from the shape
projected=[xPro1(:,1),xPro1(:,2)];
 
     Q=projected;
     A=(Q-Data1).^2;
     A= sqrt(sum(A,2));
     E_distanceData1=A;
     
     
     
   QQ=centroid;
      QQ= repmat(QQ,length(Data1),1);
     A=(QQ-Q).^2;
     A= sqrt(sum(A,2));
     E_distanceDataC1=A;
     [M1 I]=find(E_distanceDataC1<LevelOfCloseness);
     IndexOfClosePoints=M1;
     
     Count=sum(I);
     
     if(~isempty(M1))
 %           E_distanceDataMin1=1000;
 %           E_distanceDataMax1=1000;
 %    else
           [M ,Imin1]=min(E_distanceData1(IndexOfClosePoints));
             QQ=centroid;
             Q=[Data1(Imin1,1) Data1(Imin1,2)];
             A=(QQ-Q).^2;
             A= sqrt(sum(A,2));
             E_distanceDataMin1=A;
             [M, Imax1]=max(E_distanceData1(IndexOfClosePoints));
             QQ=centroid;
             Q=[Data1(Imax1,1) Data1(Imax1,2)];
             A=(QQ-Q).^2;
             A= sqrt(sum(A,2));
           E_distanceDataMax1=A;
     end
         
%      [M I]=min(E_distanceDataC1)
%           QQ=centroid;
%      Q=[Data1(I,1) Data1(I,2)];
%      A=(QQ-Q).^2;
%      A= sqrt(sum(A,2));
%      E_distanceData1=A;
end

%%
[l2, m1]=size(Data2);
if(l2<=1)
     E_distanceData2=[];
     E_distanceDataC2=[];
     E_distanceDataMin2=0;
     E_distanceDataax2=0;
     D2=0;
     M2=[];
else
xPro1 = projPointOnLine([Data2(:,1) Data2(:,2)], line1); %x's are poit on ALI from the shape
projected=[xPro1(:,1),xPro1(:,2)];
     Q=projected;
     A=(Q-Data2).^2;
     A= sqrt(sum(A,2));
      E_distanceData2=A;
     QQ=centroid;
       QQ= repmat(QQ,length(Data2),1);
      A=(QQ-Q).^2;
      A= sqrt(sum(A,2));
      E_distanceDataC2=A;
%     if(M==[])
%            E_distanceDataMin2=1000;
%            E_distanceDataMax2=1000;
        
     
 %    [M2 II]=min(E_distanceDataC2);
 %    QQ=centroid;
 %    Q=[Data2(II,1) Data2(II,2)];
 %    A=(QQ-Q).^2;
 %    A= sqrt(sum(A,2));
 %    E_distanceData2=A;
     %     [MM IIM]=find(M==E_distanceDataC2);
     %sum(IIM)
     [M2 I]=find(E_distanceDataC2<LevelOfCloseness);
     IndexOfClosePoints=M2;
     Count=sum(I);
    
     if(~isempty(M2))
           [M, Imin2]=min(E_distanceData2(IndexOfClosePoints));
             QQ=centroid;
             Q=[Data2(Imin2,1) Data2(Imin2,2)];
             A=(QQ-Q).^2;
             A= sqrt(sum(A,2));
             E_distanceDataMin2=A;
             [M, Imax2]=max(E_distanceData2(IndexOfClosePoints));
             QQ=centroid;
             Q=[Data2(Imax2,1) Data2(Imax2,2)];
             A=(QQ-Q).^2;
             A= sqrt(sum(A,2));
           E_distanceDataMax2=A;
     end 
end
  
     if(isempty(M1) && ~isempty(M2))
                     E_distanceDataMin1=E_distanceDataMin2;
         E_distanceDataMax1=E_distanceDataMax2;
         Imin1=Imin2;
         Imax1=Imax2;
     elseif(~isempty(M1) && isempty(M2))

                  E_distanceDataMin2=E_distanceDataMin1;
         E_distanceDataMax2=E_distanceDataMax1;
         Imin2=Imin1;
         Imax2=Imax1;
     elseif(isempty(M1) && isempty(M2))
         E_distanceDataMin2=0;
         E_distanceDataMax2=0;
         Imin2=0;
         Imax2=0;
         E_distanceDataMin1=0;
         E_distanceDataMax1=0;
         Imin1=0;
         Imax1=0;
     end
%  line2 = [floor(xLine(1,1)) floor(yLine(1,1)) floor(xLine(1,floor(rho)))-floor(xLine(1,1)) floor(yLine(1,floor(rho)))-floor(yLine(1,1))];
%  line1 = [x11 y11 x22-x11 y22-y11];
% xPro = projPointOnLine([L(:,1) L(:,2)], line1); %x's are poit on ALI from the shape
% xProo = projPointOnLine([L(:,1) L(:,2)], line2); %x's are poit on ALI from the shape
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    figure, imshow(im)
%      hold;
%      line([X1 + meanx, X2 + meanx],[Y1 + meany, Y2 + meany])
%     
%      plot(Data1(:,1), Data1(:,2),'g.')
%           plot(Data2(:,1), Data2(:,2),'r.')
%      plot(Data1(Imin1,1), Data1(Imin1,2),'g*')
%           plot(Data1(Imax1,1), Data1(Imax1,2),'go')
%     plot(Data2(Imin2,1), Data2(Imin2,2),'r*')
%         plot(Data2(Imax2,1), Data2(Imax2,2),'ro')
% %      plot(xProo(:,1), xProo(:,2),'g.','markers',1)
%       plot(xPro(:,1), xPro(:,2),'r.','markers',1)
%     plot(meanx, meany,'go')%centroid
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %                                           plot(c, r,'r.')
%     text(rr1,cc1,C)
%     plot(ref(1), ref(2),'go')
%     plot(xcents, ycents, 'r*');   % unprojects points on ALI
%     plot(x(:,1), x(:,2), 'r*');  
%    plot(Location2(end,1), Location2(end,2),'go')%  exteremepoint 1
%    plot(Location2(1,1), Location2(1,2),'go')% exteremepoint 2
%     %%%%%%%%%%%%
% %    %plot(c(1), r(1),'ro')
% %                   %  plot(xPro(1,1), xPro(1,2),'ro')
% %     %%%%%%%%%%%%%%%%%%
    %plot(x11, y11,'go')% Last point on ALI
    %plot(x22, y22,'go')% first point on ALI
  
  %  plot([X1 + meanx, X2 + meanx],[Y1 + meany, Y2 + meany],'r+')
%     hold;


%[f1,x0]=ksdensity([E_distanceData1 ;E_distanceData2]);
%[fC1,xC0]=ksdensity(E_distanceDataC1);
%figure(100),
%plot(f1)
%[f2,x1]=ksdensity(E_distanceData2);
%[fC2,xC1]=ksdensity(E_distanceDataC2);
%[aqa indf1]=max(f1);
%[aqa indfC1]=max(fC1);
%[aqa indf2]=max(f2);
%[aqa indfC2]=max(fC2);
%sqrt(sum((x0(indf1)-x1(indf2))^2));
%[x0(indf1) x1(indf2)]
%[mean(E_distanceData1) mean(E_distanceData2)]
%sqrt(sum((indf1-indf2)^2))
%[min(x0) max(x0)]
%[min(x1) max(x1)]
 %          figure(6)
 %     plot(f1,'k.-')
 %     hold on;
 %     plot(f2,'r.-')
% figure(10),
% plot(x0,f0,'k-')
%   figure(10),
%     plot((E_distanceData1(:)),'b.-')
%  
%   figure(11),
%      plot((E_distanceData2(:)),'r.-')
%      
%      figure(12),
%      plot((E_distanceDataC1(:)),'g.-')
%      
%      figure(13),
%      plot((E_distanceDataC2(:)),'y.-')
   % A=(f1-f2).^2;
   % error= sqrt(sum(A,2));
  %FV=[mean(x0(indf1))  1-std(E_distanceDataC1) 1-std(E_distanceDataC2)];
   % FV=[mean(x0(indf1)) mean(x1(indf2)) 1-std(E_distanceDataC1) 1-std(E_distanceDataC2)];
 %FV=[f1 ; f2];

 if(D1 == 0)
     E_distanceData1=0;
     E_distanceDataC1=0;
 end
 if(D2 == 0)
     E_distanceData2=0;
     E_distanceDataC2=0;
 end
%     [f0,x0]=ksdensity([E_distanceData1(:);E_distanceData2(:)]);
%    [aqa indf0]=max(f0);
  %  [f1,x1]=ksdensity(E_distanceData2(:));
%    [aqa indf1]=max(f1);
%    [f2,x2]=ksdensity(E_distanceDataC1(:));
%   [aqa indf2]=max(f2);
%     [f3,x3]=ksdensity(E_distanceDataC2(:));
%   [aqa indf3]=max(f3);
% FeatureVector =[x0(indf0) x1(indf1) x2(indf2)];
%  [f1,x1]=ksdensity(E_distanceData1);
%  [f2,x1]=ksdensity(E_distanceData2);
%  [f3,x1]=ksdensity(E_distanceDataC1);
%  [f4,x1]=ksdensity(E_distanceDataC2);mean(E_distanceDataC1(:))  mean(E_distanceDataC2(:))
     %FV=[mean(E_distanceData1(:)) mean(E_distanceData2(:))  ];   
  %    FV=[ abs(E_distanceDataMin1 -E_distanceDataMax1) abs(E_distanceDataMin2 -E_distanceDataMax2) numel(E_distanceData1(:))/(numel(E_distanceData1(:))+numel(E_distanceData2(:))) numel(E_distanceData2(:))/(numel(E_distanceData1(:))+numel(E_distanceData2(:)))]; 
        FV=[ abs(numel(E_distanceData1(:))-numel(E_distanceData2(:))) ]; 
   %   FV=[ E_distanceDataMin1 E_distanceDataMax1 E_distanceDataMin2 E_distanceDataMax2]; 
      %mean(E_distanceDataC1(:)) mean(E_distanceDataC2(:))
  %FV=[x0 ];   
   %   FV=[x0(indf0)  x1(indf1) x2(indf2) x3(indf3)];   
%     figure(12),
%     [r lag]=xcorr(f1,f2);
%     plot(lag,r,'g.-')
%%


     
%     if(length(E_distanceData1)>length(E_distanceData2))
%        E_distanceData2(length(E_distanceData2):length(E_distanceData1),1)=0;
%    elseif(length(E_distanceData1)<length(E_distanceData2))
%        E_distanceData1(length(E_distanceData1):length(E_distanceData2),1)=0;
%     end

 %FV=[xcorr((E_distanceData1),(E_distanceData2),0,'coeff') ];
  %  a=[abs(mean(E_distanceData1)-mean(E_distanceData2)) abs(mean(E_distanceData1)-mean(flipud(E_distanceData2)))];
 % aa=[mean([std(E_distanceData1),std(E_distanceData2)]) mean([std(E_distanceData1),std(flipud(E_distanceData2))])];
 
  %a=[(mean(E_distanceData1+E_distanceData2)) (mean(E_distanceData1+flipud(E_distanceData2)))];
 % aa=[mean([std(E_distanceData1),std(E_distanceData2)]) mean([std(E_distanceData1),std(flipud(E_distanceData2))])];
 
 
 
 %aa=[max([std(E_distanceData1),std(E_distanceData2)]) max([std(E_distanceData1),std(flipud(E_distanceData2))])];
% 
% if a(1)>a(2)
%    % FS=[E_distanceData1(:);E_distanceData2(:)];
%     figure(flage),
%     plot((E_distanceData1(:)),'b.-')
%     hold on;
%     plot((E_distanceData2(:)),'r.-')
%     
%    %  aa=[max(abs(E_distanceData1(2:end)-E_distanceData1(1:end-1))),max(abs(E_distanceData2(2:end)-E_distanceData2(1:end-1)))];
% else
%   %   FS=[E_distanceData1(:);flipud(E_distanceData2(:))];
%         figure(flage),
%     plot((E_distanceData1(:)),'b.-')
%     hold on;
%     plot(flipud(E_distanceData2(:)),'r.-')
%   %  aa=[max(abs(E_distanceData1(2:end)-E_distanceData1(1:end-1))),max(abs(flipud(E_distanceData1(2:end))-flipud(E_distanceData1(1:end-1))))];
% end
 % [f0,x0]=ksdensity(FS);
  %   sttd= max(aa);
  %   coor= min(a);
 %    FV=f0;
 %[r l] = xcorr(E_distanceData1,flipud(E_distanceData2));
 %FeatureVector=r(find(l==0))/(rows*cols);
%FeatureVector =mean(r)/max(r);
 %FeatureVector=abs(std(E_distanceData1)-std(E_distanceData2));%sqrt((out1-out2).^2+(out3-out4).^2);

%FeatureVector=[entropy(E_distance(:,1)) mean(E_distance(:,1)) skewness(E_distance(:,1))];
%for i=2:n
%FeatureVector=[FeatureVector entropy(E_distance(:,i)) mean(E_distance(:,i)) skewness(E_distance(:,i))];
%end
%FeatureVector;
