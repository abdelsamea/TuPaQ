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

 function [FeatureVector] = momALI(im)
im=double(im);
  
s=strel('disk',1,0);%Structuring element
F=imerode(im,s);%Erode the image by structuring element
Bound=im;%-F;


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
    rho = sqrt(area)/(roundness + 0.5) ; 
    %rho = sqrt(area)/5;%(roundness + 0.5) ;
    [X1,Y1] = pol2cart(thetamin,      rho);
    [X2,Y2] = pol2cart(thetamin + pi, rho);
    
    [Xo1,Yo1] = pol2cart(1/thetamin,      rho);
    [Xo2,Yo2] = pol2cart(1/(thetamin - pi), rho);


%[Bound]=Binary_trans(im);
%entropy(Bound)
[r, c]=find(Bound);
%rc = bwtraceboundary(im,[y1 x1],'N');
%Location=[rc(:,2) rc(:,1)];


%[r, c] = poly2cw(r1, c1);
Location=[c r];


l=1:length(c);
C=num2cell(l);

Location1 = sortrows(Location,1);
Location2 = sortrows(Location,2);

x1 = Location2(end,1);
x2 = Location2(1,1);
y1 = Location2(end,2);
y2 = Location2(1,2);

x11 = X1 + meanx;
x22 = X2 + meanx;
y11 = Y1 + meany;
y22 = Y2 + meany;

xo11 = Xo1 + meanx;
xo22 = Xo2 + meanx;
yo11 = Yo1 + meany;
yo22 = Yo2 + meany;
n =1;  %10 segments
deltax = (x2 - x1)/n;
deltay = (y2 - y1)/n;
halfdx = deltax/2;
halfdy = deltay/2;
%radius = sqrt(halfdx.^2 + halfdy.^2);
     [nrow, ncol]=size(Location);
     line1 = [x11 y11 x22-x11 y22-y11];
    %%
       m = (diff([Y1 + meany, Y2 + meany ])/diff([X1 + meanx, X2 + meanx]));
    % Slope of new line
    slope = -1/m;
    x = meanx;
    y = meany;
    lengthLine = rho;
    xLine = x-lengthLine:x+lengthLine;
    yLine = slope*(xLine-x) + y;
    line2 = [floor(xLine(1,1)) floor(yLine(1,1)) floor(xLine(1,floor(rho)))-floor(xLine(1,1)) floor(yLine(1,floor(rho)))-floor(yLine(1,1))];
   DALI=sqrt((Location2(end,1)-Location2(1,1))^2+(Location2(end,2)-Location2(1,2))^2);
    
  
%%
xPro = projPointOnLine([Location(:,1) Location(:,2)], line1); %x's are poit on ALI from the shape
projected1=[xPro(:,1),xPro(:,2)];
 Q1=projected1;
    
    A=(Q1-Location).^2;
    E_distance1= sqrt(sum(A,2));
%%
xProo = projPointOnLine([Location(:,1) Location(:,2)], line2); %x's are poit on ALI from the shape
projected2=[xProo(:,1),xProo(:,2)];
 Q2=projected2;
     
    A=(Q2-Location).^2;
    E_distance2= sqrt(sum(A,2));
%%
   Q3=centroid;
      Q3= repmat(Q3,length(Location),1);
     A=(Q3-Location).^2;
     E_distance3= sqrt(sum(A,2));
%%
  %%
% % 
%     A=(Q3-projected1).^2;
%     E_distance4= sqrt(sum(A,2));
% % %%
% % %%
%     A=(Q3-projected2).^2;
%     E_distance5= sqrt(sum(A,2));
%%
%      if(length(E_distance1)>length(E_distance2))
%         E_distance2(length(E_distance2):length(E_distance1),1)=0;
%     elseif(length(E_distance1)<length(E_distance2))
%         E_distance1(length(E_distance1):length(E_distance2),1)=0;
%      end
%  FeatureSpace=[E_distance1(:);E_distance2(:)];% ;E_distance5(:) ;E_distance4(:)];
  %   FeatureSpace=[E_distance3(:);E_distance4(:);E_distance5(:)];%E_distance5(:)];% ;E_distance3(:)];   


%  figure, imshow(im)
%     hold;
%     line([X1+ meanx , X2 + meanx],[Y1 + meany, Y2+ meany ])
%  % plot(xLine,yLine,'k.')
%     plot(meanx, meany,'y*','markers',10)%centroid
%     %plot(c(:), r(:),'y.','markers',1)
%     plot(c(5), r(5),'r*','markers',10)
%      plot(c(900), r(900),'g*','markers',10)
%     plot(Location2(end,1), Location2(end,2),'b+')%  exteremepoint 1
%     plot(Location2(1,1), Location2(1,2),'b+')% exteremepoint 2
%      plot(xProo(:,1), xProo(:,2),'g.','markers',1)
%      plot(xPro(:,1), xPro(:,2),'r.','markers',1)
%    plot(xPro(5,1), xPro(5,2),'r*','markers',10)
%     plot(xProo(900,1), xProo(900,2),'g*','markers',10)
    
    %%
%      figure, imshow(im)
%     hold;
%     line([X1+ meanx , X2 + meanx],[Y1 + meany, Y2+ meany ])
%   plot(xLine,yLine,'k.')
%     plot(meanx, meany,'y*','markers',20)%centroid
%     %plot(c(:), r(:),'y.','markers',1)
%     plot(c(20000), r(20000),'r*','markers',20)
%      plot(c(10000), r(10000),'g*','markers',20)
%     plot(Location2(end,1), Location2(end,2),'b+')%  exteremepoint 1
%     plot(Location2(1,1), Location2(1,2),'b+')% exteremepoint 2
%      plot(xProo(:,1), xProo(:,2),'g.','markers',1)
%      plot(xPro(:,1), xPro(:,2),'r.','markers',1)
%    plot(xPro(20000,1), xPro(20000,2),'r*','markers',20)
%     plot(xProo(10000,1), xProo(10000,2),'g*','markers',20)
    %%
%%          figure(3)
%      plot(E_distance1(:),'r.-')
%     figure(4)
%      plot(E_distance2(:),'b.-')
%       figure(5)
%      plot(FeatureSpace,'g.-')
%      FeatureSpace = [E_distance1(:) ; E_distance2(:); E_distance3(:)];
%    [f0,x0]=ksdensity(FeatureSpace);
%    [aqa indf0]=max(f0);
%    [f1,x1] = ksdensity(E_distance1);
%     [f2,x2] = ksdensity(E_distance2);
%      [f3,x3] = ksdensity(E_distance3);
%    figure(1),
%    plot(x0,f0,'k-')
%    hold on,
%    plot(x1,f1,'r-')
%    plot(x2,f2,'b-')
%   plot(x3,f3,'g-')
%FeatureVector=f0;
%[bandwidth,FeatureVector,xmesh,cdf]=kde(FeatureSpace);
%    [f0,x0]=ksdensity(E_distance1(:));
%    [aqa indf0]=max(f0);
%    [f1,x1]=ksdensity(E_distance2(:));
%    [aqa indf1]=max(f1);
%    [f2,x2]=ksdensity(E_distance3(:));
%   [aqa indf2]=max(f2);
%  % FeatureVector =[f0 f1 f2];
%  FeatureVector =[x0(indf0) x1(indf1) x2(indf2)];
 % FeatureVector =[mean(FeatureSpace(:))];
% FeatureVector =[mean(FeatureSpace) ]; mean(E_distance2(:)) 
 FeatureVector =[  mean(E_distance1(:))  mean(E_distance3(:))];
  %FeatureVector=(FeatureVector);
%FeatureVector =[x0 ];
 %   FeatureVector = SOM(FeatureSpace,10,10000,0.9);
%FeatureVector=  sqrt(sum((f1-f2).^2,2))