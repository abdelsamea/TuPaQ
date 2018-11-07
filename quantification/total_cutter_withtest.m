function [labimg, MergIma,number_of_nuclei]=total_cutter_withtest(x,small_nuclei,hmin,do_cut)

if ~do_cut
    [labimg,number_of_nuclei] = bwlabel(x);
    return;
end


    if (size(small_nuclei,3)==3)
        small_nuclei=rgb2gray(small_nuclei);
    end

  
 %%%%%%%%%%%%%%%%%%%%%%%%   
 
   Foreground=double(small_nuclei);
   BW = edge(Foreground,'Canny');
   
   
 %%%%%%%%%%%%%%%%%%%%%%%
  % imtool(small_nuclei)
     if(hmin==9999)
           Foreground=double(small_nuclei).*double(x);
           [hmin, number]=mode(Foreground(Foreground>0));
           hmin=hmin/2;
                lines=ones(size(small_nuclei));
                markers=imhmin(small_nuclei,hmin);
                lines(watershed(markers)==0)=0;
                yStraight=lines&x;
                MergImaStraight=abs(yStraight-x);
                Straight=sum(sum( imabsdiff(logical(BW),logical(MergImaStraight))==1) );
                
           
            small_nuclei=double(small_nuclei).*double(x);
            small_nuclei=small_nuclei/max(max(small_nuclei));%%%
            small_nuclei=imcomplement(small_nuclei);
            Foreground=double(small_nuclei).*double(x);
            hmin=mode(Foreground(Foreground>0))/2;
                    lines=ones(size(small_nuclei));
                markers=imhmin(small_nuclei,hmin);
                lines(watershed(markers)==0)=0;
                yComplement=lines&x;
                MergImaComplement=abs(yComplement-x);
                Complement=sum(sum( imabsdiff(logical(BW),logical(MergImaComplement))==1) );
            
     end
     if(Straight>Complement)
         [labimg,number_of_nuclei] = bwlabel(yStraight);
         MergIma=MergImaStraight;
     else
         [labimg,number_of_nuclei] = bwlabel(yComplement);
         MergIma=MergImaComplement;
     end
         
%      if (hmin > aver)
%           
%             small_nuclei=double(small_nuclei).*double(x);
%             small_nuclei=small_nuclei/max(max(small_nuclei));%%%
%             small_nuclei=imcomplement(small_nuclei);
%             Foreground=double(small_nuclei).*double(x);
%             hmin=mode(Foreground(Foreground>0))/2;
%             clc
%           hmin 
%           aver 
%      end
%     lines=ones(size(small_nuclei));
%     markers=imhmin(small_nuclei,hmin);
%     lines(watershed(markers)==0)=0;
%     y=lines&x;
%     [labimg,number_of_nuclei] = bwlabel(y);
% 
%     MergIma=abs(y-x);
%     imtool(MergIma)
%    sum(sum( imabsdiff(logical(BW),logical(MergIma))==1) )


