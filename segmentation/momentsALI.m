function [FeatureVector] = momentsALI(im)
%im=imread(filename); 
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
    rho = sqrt(area)/5;%(roundness + 0.5) ;
    [X1,Y1] = pol2cart(thetamin,      rho);
    [X2,Y2] = pol2cart(thetamin + pi, rho);
%    FileName11=strcat('Proposed', filename);
%saveDataName = fullfile(FileName11);
%imwrite(Bound,saveDataName);
[r1, c1]=find(Bound);
[r, c] = poly2cw(r1, c1);
Location=[c r];
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

n =15;  %10 segments
deltax = (x2 - x1)/n;
deltay = (y2 - y1)/n;
halfdx = deltax/2;
halfdy = deltay/2;
xcents = x1 + halfdx + (0:n-1)*deltax;
ycents = y1 + halfdy + (0:n-1)*deltay;
    line1 = [x11 y11 x22-x11 y22-y11];
    x=zeros(n,2);
    x = projPointOnLine([xcents(:) ycents(:)], line1); %x's are point on ALI
    [nrow ncol]=size(Location);
    line1 = [x11 y11 x22-x11 y22-y11];
    xPro=zeros(nrow,2);
    xPro = projPointOnLine([Location(:,1) Location(:,2)], line1); %x's are poit on ALI from the shape
    [nsmaples po] = size(xPro);
    I = round(linspace(1,nsmaples,n));
    PoALI= xPro(I,:);
    size(PoALI);
    E_distance = zeros(nsmaples,n) ;
    Distributions = zeros(nsmaples-1,n) ;
    DALI=sqrt((Location2(end,1)-Location2(1,1))^2+(Location2(end,2)-Location2(1,2))^2);
    
for i=1:n
    Q=PoALI(i,:);
    Q= repmat(Q,length(Location),1);
    A=(Q-Location).^2;
    E_distance(:,i) = sqrt(sum(A,2));
 end
 for i=1:n
    E_distance(:,i) = E_distance(:,i);
 end
 
% imshow(im);
%     hold;
%     line([X1 + meanx, X2 + meanx],[Y1 + meany, Y2 + meany])
%     plot(meanx, meany,'ro')%centroid
%     plot(Location2(end,1), Location2(end,2),'go')% exteremepoint 1
%     plot(Location2(1,1), Location2(1,2),'go')% exteremepoint 2
%     hold;
%     figure(3)
%     plot(E_distance(:))
FeatureVector=[mean(E_distance(:,1)) median(E_distance(:,1)) std(E_distance(:,1)) skewness(E_distance(:,1)) kurtosis(E_distance(:,1))];
for i=2:n
FeatureVector=[FeatureVector mean(E_distance(:,i)) median(E_distance(:,i)) std(E_distance(:,i)) skewness(E_distance(:,i)) kurtosis(E_distance(:,i))];
end
%FeatureVector