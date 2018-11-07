
 function [FV] = mmomALI(im,flage)
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
   
   rho = area;%/5;%(roundness + 0.5) ;
    [X1,Y1] = pol2cart(thetamin,      rho);
    [X2,Y2] = pol2cart(thetamin + pi, rho);
  
[rrr, ccc]=find(Bound);

%entropy(Bound)
[r, c]=find(Bound);
%rc = bwtraceboundary(im,[y1 x1],'N');
%Location=[rc(:,2) rc(:,1)];


%[r, c] = poly2cw(r1, c1);
Location=[c r];
if(r~=2)
    Location=[0 0];
end
%%

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

n =1;  %10 segments
deltax = (x2 - x1)/n;
deltay = (y2 - y1)/n;
halfdx = deltax/2;
halfdy = deltay/2;

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
     xPro = projPointOnLine([Location(:,1) Location(:,2)], line1); %x's are poit on ALI from the shape
    DALI=sqrt((xPro(end,1)-xPro(1,1))^2+(xPro(end,2)-xPro(1,2))^2);
elseif(flage==2)
    line1 = [x11 y11 x22-x11 y22-y11];
    DALI=sqrt((Location2(end,1)-Location2(1,1))^2+(Location2(end,2)-Location2(1,2))^2);
end




    xPro = projPointOnLine([Location(:,1) Location(:,2)], line1); %x's are poit on ALI from the shape
     Q=xPro;

    D=Q-Location;  

  c1=1;
  c2=1;

  Data1=[];
  Data2=[];
 D1=1;
 D2=1;
 [le,col]=size(Location);
  for kk=1:le
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
      [aaa, sss]=poly2ccw(Data2(:,1), Data2(:,2));
      Data2=[aaa, sss];
  end
   %%

%%
[l1, m1]=size(Data1);
if(l1<=1)
     E_distanceData1=[];
     E_distanceDataC1=[];
     D1=0;
else
xPro1 = projPointOnLine([Data1(:,1) Data1(:,2)], line1); %x's are poit on ALI from the shape
projected=[xPro1(:,1),xPro1(:,2)];
 
     Q=projected;
     A=(Q-Data1).^2;
     A= sqrt(sum(A,2));
     E_distanceData1=A;
   QQ=centroid;
      QQ= repmat(QQ,length(Data1),1);
     A=(QQ-Data1).^2;
     A= sqrt(sum(A,2));
     E_distanceDataC1=A;
             QQ=centroid;
      QQ= repmat(QQ,length(Data1),1);
     A=(QQ-Q).^2;
     A= sqrt(sum(A,2));
     E_distanceC1=A;
end

%%
[l2, m1]=size(Data2);
if(l2<=1)
     E_distanceData2=[];
     E_distanceDataC2=[];
     D2=0;
else
xPro1 = projPointOnLine([Data2(:,1) Data2(:,2)], line1); %x's are poit on ALI from the shape
projected=[xPro1(:,1),xPro1(:,2)];
     Q=projected;
     A=(Q-Data2).^2;
     A= sqrt(sum(A,2));
     E_distanceData2=A;
    QQ=centroid;
      QQ= repmat(QQ,length(Data2),1);
     A=(QQ-Data2).^2;
     A= sqrt(sum(A,2));
     E_distanceDataC2=A;
      QQ=centroid;
      QQ= repmat(QQ,length(Data2),1);
     A=(QQ-Q).^2;
     A= sqrt(sum(A,2));
     E_distanceC2=A;
end
  



 if(D1 == 0)
     E_distanceData1=0;
     E_distanceDataC1=0;
     E_distanceC1=0;
     FV=[100 100 100];
 elseif(D2 == 0)
     E_distanceData2=0;
     E_distanceDataC2=0;
      E_distanceC2=0;
      FV=[100 100 100];
 else
 
   if(length(E_distanceC1)>=length(E_distanceC2))
       E_distanceC2(length(E_distanceC2):length(E_distanceC1),1)=0;
   else
       E_distanceC1(length(E_distanceC1):length(E_distanceC2),1)=0;
    end
    if(length(E_distanceData1)>=length(E_distanceData2))
       E_distanceData2(length(E_distanceData2):length(E_distanceData1),1)=0;
   else
       E_distanceData1(length(E_distanceData1):length(E_distanceData2),1)=0;
    end
    if(length(E_distanceDataC1)>=length(E_distanceDataC2))
       E_distanceDataC2(length(E_distanceDataC2):length(E_distanceDataC1),1)=0;
   else
       E_distanceDataC1(length(E_distanceDataC1):length(E_distanceDataC2),1)=0;
    end
     A=(E_distanceC1-E_distanceC2).^2;
     d1= mean(sqrt(sum(A,2)));
          A=(E_distanceData1-E_distanceData2).^2;
     d2= mean(sqrt(sum(A,2)));
          A=(E_distanceDataC1-E_distanceDataC2).^2;
     d3= mean(sqrt(sum(A,2)));
    if(d1==0)
        d1=.01;
    end
    if(d2==0)
        d2=.01;
    end
    if(d3==0)
        d3=.01;
    end
    FV=[d1 d2 d3];

 end
