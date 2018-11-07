function [labimg, MergIma,number_of_nuclei]=total_cutter(x,small_nuclei,hmin,do_cut)

if ~do_cut
    [labimg,number_of_nuclei] = bwlabel(x);
    return;
end


    if (size(small_nuclei,3)==3)
        small_nuclei=rgb2gray(small_nuclei);
    end

    small_nuclei=double(small_nuclei);

    
    
   
    small_nuclei=small_nuclei/max(max(small_nuclei));

     if(hmin==9999)
           Foreground=small_nuclei.*x;
           hmin=mode(mode(Foreground(Foreground>0)))/2;
         %  hmin=0.01;
     else
         
     end
    lines=ones(size(small_nuclei));
   
    markers=imhmin(small_nuclei,hmin);

    lines(watershed(markers)==0)=0;
    y=lines&x;
    [labimg,number_of_nuclei] = bwlabel(y);

    MergIma=abs(y-x);
    


