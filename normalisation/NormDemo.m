% function [] = readwrite(filename)
% I=imread(filename);
% saveDataName = fullfile(filename);
% imwrite(I,saveDataName,'compression','none') ;
% end

function [] = NormDemo()
%verbose = 0;
TargetImage = imread('Ref.bmp');
imagefiles = dir('*.bmp');      
nfiles = length(imagefiles);    % Number of files found
%FeatureVector=[];
for ii=1:nfiles
    currentfilename = imagefiles(ii).name;
%    currentimage = imread(currentfilename);
%    images{ii} = currentimage;

SourceImage=imread(currentfilename);
%I = imresize(I,.5);


%%
[ NormHS ] = Norm( SourceImage, TargetImage, 'RGBHist' );

FileName1=strcat(currentfilename,'NormHS.png');
saveDataName = fullfile(FileName1);
imwrite(NormHS,saveDataName,'compression','none') ;
%%
[ NormRH ] = Norm( SourceImage, TargetImage, 'Reinhard' );
FileName2=strcat(currentfilename,'NormRH.png');
saveDataName = fullfile(FileName2);
imwrite(NormRH,saveDataName,'compression','none') ;
%%
[ NormMM ] = Norm(SourceImage, TargetImage, 'Macenko', 255, 0.15, 1);
FileName3=strcat(currentfilename,'NormMM.png');
saveDataName = fullfile(FileName3);
imwrite(NormMM,saveDataName,'compression','none') 

%%
%saveDataName = fullfile(currentfilename);
%imwrite(I,saveDataName,'compression','none') ;
end
