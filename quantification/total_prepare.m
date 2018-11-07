
function [bw,threshold,graynuclei]=total_prepare(nuclei,thrmethod)

if nargin < 4
    backsub=1;
end

 

    threshold = graythresh(nuclei);
    if strcmp('otsu',thrmethod)
        bw = im2bw(nuclei,threshold);
    elseif strcmp('localvar',thrmethod)
        bw = adthr(nuclei);
    end

   bw=~bw;
end