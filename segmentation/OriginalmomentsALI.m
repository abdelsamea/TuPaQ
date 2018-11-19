function [] = OriginalmomentsALI()
imagefiles = dir('*.png');      
nfiles = length(imagefiles);    % Number of files found
FeatureVector=[];
for ii=1:nfiles
    currentfilename = imagefiles(ii).name;
%    currentimage = imread(currentfilename);
%    images{ii} = currentimage;
% end
im=imread(currentfilename); 
im=double(im);
%scale = .5;
%im = imresize(im,scale);
%sigma=1.5;
%KernelSize = 2*round(2*sigma)+1;
%G = fspecial('gaussian', KernelSize, sigma);
%  im = conv2(im, G, 'same'); % regularization
  
s=strel('disk',1,0);%Structuring element
F=imerode(im,s);%Erode the image by structuring element
Bound=im-F;


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
    rho = sqrt(area)/5;%(roundness + 0.5)+ 5 ; 
    [X1,Y1] = pol2cart(thetamin,      rho);
    [X2,Y2] = pol2cart(thetamin + pi, rho);
    
    
    %figure(2),imshow(Bound);title('Boundary extracted Image');
%     FileName11=strcat('Proposed', filename);
% saveDataName = fullfile(FileName11);
% imwrite(Bound,saveDataName) ;
% [r, c]=find(Bound);
% Location=[c r];
% Locatio=[c r];
% Location1 = sortrows(Location,1);
% Location2 = sortrows(Location,2);
% 
% x1 = Location2(end,1);
% x2 = Location2(1,1);
% y1 = Location2(end,2);
% y2 = Location2(1,2);
% 
% x11 = X1 + meanx;
% x22 = X2 + meanx;
% y11 = Y1 + meany;
% y22 = Y2 + meany;

%    aD=([x1 y1]-[x2 y2]).^2;
%    D = sqrt(sum(aD,2));
% n = 5;  %10 segments
% deltax = (x2 - x1)/n;
% deltay = (y2 - y1)/n;
% halfdx = deltax/2;
% halfdy = deltay/2;
% %radius = sqrt(halfdx.^2 + halfdy.^2);
% xcents = x1 + halfdx + (0:n-1)*deltax;
% ycents = y1 + halfdy + (0:n-1)*deltay;
% 
%      line1 = [x11 y11 x22-x11 y22-y11];
%      x = projPointOnLine([xcents(:) ycents(:)], line1); %x's are point on ALI
% 
%      line1 = [x11 y11 x22-x11 y22-y11];
%      x=zeros(n,2);
%      x = projPointOnLine([xcents(:) ycents(:)], line1); %x's are point on ALI
% 
% %%
% x1=x;
% x1(:,1) = x(:,1) + 50*cos(-1/thetamin);
% x1(:,2) = x(:,2) + 50*sin(-1/thetamin);
%%
%imshow(im);
%    hold;
%    line([X1 + meanx, X2 + meanx],[Y1 + meany, Y2 + meany])
%    line([x1(1,1) , x(1,1)],[x1(1,2), x(1,2)])
%    plot(meanx, meany,'r+')%centroid
    %plot(c, r,'g.')
 %   plot(xcents, ycents, 'r*');   % unprojects points on ALI
 %   plot(Location2(end,1), Location2(end,2),'go')% exteremepoint 1
 %   plot(Location2(1,1), Location2(1,2),'go')% exteremepoint 2
    %%%%%%%%%%%%%
 %   plot(x(:,1), x(:,2),'go')
    
    %%%%%%%%%%%%%%%%%%
    %plot(x11, y11,'go')% Last point on ALI
    %plot(x22, y22,'go')% first point on ALI
    
   % plot([X1 + meanx, X2 + meanx],[Y1 + meany, Y2 + meany],'r+')
   % hold;
    %plot(r(1), c(1),'go')
    %hold;
    %Q=[meanx meany];
    
   % Q= repmat(Q,length(Location),1);
   % A=(Q-Location).^2;
   % E_distance = sqrt(sum(A,2));
    %E_distance = sort(E_distance);
    
    %A1=(Q-Location1).^2;
    %E_distance1 = sqrt(sum(A1,2));
    
   % A2=(Q-Location2).^2;
   % E_distance2 = sqrt(sum(A2,2));
    
    
    %E_distance3 =E_distance1./E_distance2;
    
     %figure(2)
     %plot(E_distance)
%     
%         figure(3)
%     plot(E_distance1)
%     
%         figure(4)
%     plot(E_distance2)
%     
%             figure(5)
%     plot(E_distance3)
%     
    %FeatureVector=[mean(E_distance) median(E_distance) std(E_distance) skewness(E_distance) kurtosis(E_distance)]
%    FeatureVector=[mean(E_distance1) median(E_distance1) std(E_distance1) skewness(E_distance1) kurtosis(E_distance1)]
%    FeatureVector=[mean(E_distance2) median(E_distance2) std(E_distance2) skewness(E_distance2) kurtosis(E_distance2)]
%    FeatureVector=[mean(E_distance3) median(E_distance3) std(E_distance3) skewness(E_distance3) kurtosis(E_distance3)]

%     Q1=[X1 + meanx, Y1 + meany];
%     
%     Q1= repmat(Q1,length(Locatio),1);
%     A1=(Q1-Locatio).^2;
%     E_distanceExtreme1 = sqrt(sum(A1,2));
%     [i index1]=min(E_distanceExtreme1)
%     
%     
%     Q2=[X2 + meanx, Y2 + meany];
%     
%     Q2= repmat(Q2,length(Locatio),1);
%     A2=(Q2-Locatio).^2;
%     E_distanceExtreme2 = sqrt(sum(A2,2));
%     [i index2]=min(E_distanceExtreme2)
%     
%     
% rho1=(X1+meanx)*cos(thetamin)+(Y1+meany)*sin(thetamin);
% [X11,Y11] = pol2cart(1/thetamin,      rho1);
% [X22,Y22] = pol2cart(1/(thetamin + pi), rho1);
% 
% thetamin
% rho2=Location(index2,2)*cos(thetamin)+Location(index2,1)*sin(thetamin);
% [X111,Y111] = pol2cart(1/thetamin,      rho2);
% [X222,Y222] = pol2cart(1/(thetamin + pi), rho2);
FeatureVector1=[centroid thetamin roundness];
%for i=2:n
%FeatureVector1=[FeatureVector1 mean(E_distance(:,i)) median(E_distance(:,i)) std(E_distance(:,i)) skewness(E_distance(:,i)) kurtosis(E_distance(:,i))];
%end
%FeatureVector1=[mean(E_distance(:,:)) var(E_distance(:,:))];
FeatureVector=[FeatureVector;FeatureVector1];
%moment(E_distance(:,:),2)
%curvature_central(E_distance(:,1))
%     
%         figure(4)
%     plot(E_distance2)
%     
%             figure(5)
%     plot(E_distance3)
%     Distributions
%    FeatureVector=[mean(E_distance) median(E_distance) std(E_distance) skewness(E_distance) kurtosis(E_distance)]
%    FeatureVector=[mean(E_distance1) median(E_distance1) std(E_distance1) skewness(E_distance1) kurtosis(E_distance1)]
%    FeatureVector=[mean(E_distance2) median(E_distance2) std(E_distance2) skewness(E_distance2) kurtosis(E_distance2)]
%    FeatureVector=[mean(E_distance3) median(E_distance3) std(E_distance3) skewness(E_distance3) kurtosis(E_distance3)]

%     Q1=[X1 + meanx, Y1 + meany];
%     
%     Q1= repmat(Q1,length(Locatio),1);
%     A1=(Q1-Locatio).^2;
%     E_distanceExtreme1 = sqrt(sum(A1,2));
%     [i index1]=min(E_distanceExtreme1)
%     
%     
%     Q2=[X2 + meanx, Y2 + meany];
%     
%     Q2= repmat(Q2,length(Locatio),1);
%     A2=(Q2-Locatio).^2;
%     E_distanceExtreme2 = sqrt(sum(A2,2));
%     [i index2]=min(E_distanceExtreme2)
%     
%     
% rho1=(X1+meanx)*cos(thetamin)+(Y1+meany)*sin(thetamin);
% [X11,Y11] = pol2cart(1/thetamin,      rho1);
% [X22,Y22] = pol2cart(1/(thetamin + pi), rho1);
% 
% thetamin
% rho2=Location(index2,2)*cos(thetamin)+Location(index2,1)*sin(thetamin);
% [X111,Y111] = pol2cart(1/thetamin,      rho2);
% [X222,Y222] = pol2cart(1/(thetamin + pi), rho2);
%sqrt(sum((a-b).^2)/10)
end
save aALI.txt FeatureVector -ASCII
