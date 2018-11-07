% This Matlab file demomstrates the proposed the Sign Pressure Force Level Set Function based Local Fuzzy C-means (SPF-LFCM) 
%Regio Detection and Segmentation model
%
%
% the main function to be called is SPFLFCM(filename, Iter, sigma);
%
% filename is the file name of the original image to be segmented
% MaxIter is the maximum number of iteration
% segma is to control the smothness of the level set function in each
% iteration 
% 
% By
% Mohammed Abdelsamea (Email: mohammed.abdelsamea@nottingham.ac.uk)
%% the SPFLFCM parameters
%e.g., SPFLFCM('11.tif', 100, 1.5);
function [LL]= SPFLFCM(Mas, P, II, MaxIter, sigma)
%tic
KernelSize = 2*round(2*sigma)+1;
%%
%II=imread(filename);
%%
im_red = II;
im_green = II;
im_blue = II;

% Red channel only
im_red(:,:,2) = 0; 
im_red(:,:,3) = 0; 
I1 = double(rgb2gray(im_red));
% Green channel only
im_green(:,:,1) = 0; 
im_green(:,:,3) = 0; 
I2 = double(rgb2gray(im_green));
% Blue channel only
im_blue(:,:,1) = 0; 
im_blue(:,:,2) = 0; 
I3 = double(rgb2gray(im_blue));
%%
Img=(I1+I2+I3)/3;
%%
%Img = double(rgb2gray(II));
%Img=rangefilt(Img);
%Img = imresize(Img,.1);
[aa1 aa2] = size(Img);

Img1=Img;
[a1 a2] = size(Img1);
XR=a1;
YR=a2;
%P=zeros(a1,a2);



%Img1 = Img1.*(1-P);
    %Img1 = Img1.*P;

phi = ones(a1,a2);
%Mask= imfill(Mask,'holes');
%Mas=imread('Mas.bmp');
phi(Mas>0)=-1;
%%%% for particular initialization using mask image "FileN" %%%%
%    B(:,:)= imread('mask.bmp');% mask image for contour initialization
%    for i=1:a1
%        for j1=1:a2
%            if(B(i,j1) >0 )
%                phi(i,j1) = -1;
%            end
%        end
%    end
%%%% classical rectangular initialization  %%%%
 %r=150;
 %phi(r:a1-r,r:a2-r) = -1;
u =  -phi;
%%
%figure(1);
%title('Initial contour');
G = fspecial('gaussian', KernelSize, sigma);
u1=zeros(a1, a2);
t = cputime;
for n = 1:MaxIter
    %% convergence condition
         if (n > 1 && immse(u, u1)<.0001)%isequal(u1,u))
             break;
         end
        
    u1 = u;
   
    %%
    [ux, uy] = gradient(u);

     c2 =  mean(Img1(u<0));
     c1 = mean(Img1(u>=0));

    A=Img1-((c1^2-c2^2)/(2*c1-2*c2));
    spf = sign(c1-c2)*sign(A).*abs(A.^2);%.*A.^2;

%% level set equarion
u = u + (spf.*sqrt(ux.^2 + uy.^2));
%%
  if mod(n,100)==0
    imagesc(Img1,[0 255]); colormap(gray);hold on;
    [c, h] = contour(u, [0 0], 'r','LineWidth',2);
    iterNum = [num2str(n), 'iterations'];
    title(iterNum);
    pause(.0001);
  end
%%
  u = (u >= 0) - ( u< 0); % reintialization  
  u = conv2(u, G, 'same'); % regularization
end

% figure(1)
% imshow(Img1,[],'initialmagnification',200,'displayrange',[0 255]); hold on;
%   contour(u, [0 0], 'y','LineWidth',2);
%   contour(phi, [0 0], 'r','LineWidth',2);
%   hold off; title([num2str(n) ' Iterations']); drawnow;
%   FileName00=strcat('color', filename);
%   saveas(1,FileName00)
%% segmented result
LL = zeros(a1,a2);
LL(u>=0)=1;

%toc
end