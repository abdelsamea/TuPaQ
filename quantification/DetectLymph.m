
function [Inew, numberOfPositiveStainPix] = DetectLymph(OI,Istogram, step_size, theta) 
Result = OI;
SIZEOI = size(OI);

% go through every single pixel of the image
for i=1:1:SIZEOI(1)
	for j=1:1:SIZEOI(2)
        % Get the location of the histogram
        block_red = ceil(double(OI(i,j,1)) / step_size);
        block_green = ceil(double(OI(i,j,2)) / step_size);
        block_blue = ceil(double(OI(i,j,3)) / step_size);
        if (block_red == 0) block_red = 1; end;
        if (block_blue == 0) block_blue = 1; end;
        if (block_green == 0) block_green = 1; end;  
        skin = Istogram(block_red, block_green, block_blue, 1);
        not_skin = Istogram(block_red, block_green, block_blue, 2);
        psgivenx = (skin /(skin + not_skin));
        if (psgivenx > theta)
            Result(i,j,3)= 255;
            Result(i,j,2)= 255; 
            Result(i,j,2)= max(255, skin);          
        else
            Result(i,j,3)= 0; 
            Result(i,j,2)= 0;  
            Result(i,j,1)= 0;            
        end;
        
    end;
end;
Result=rgb2gray(Result);
 Result=imbinarize(Result);
 numberOfPositiveStainPix = nnz(Result);
 Inew = OI.*repmat(uint8(Result),[1,1,3]);
  Inew(Inew==0)=255;
%imtool(Inew)


