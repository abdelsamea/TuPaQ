function [FVnormal]= FExtraction(filenamenormal)
Normal=imread(filenamenormal);
 bw = logical(Normal);
measures = regionprops(bw, 'all');
    numberOfitems = size(measures, 1);
   
        FVnormal=[];
        for k = 1 : numberOfitems % Loop through all items.
            M1=measures(k).Image();
            %  figure(k*10), imshow(M1,[]);
          %  if (measures(k).Area > 1000 ) % remove noise
                    
                        
                        FVnormal=[FVnormal max([momentALI(M1,1) momentALI(M1,2)])];
                  %  M=measures(k).Image();
                 %   figure(k), imshow(M,[]);
                    

                   
           % end
        end

% a=[mmomALI(Normal,1); mmomALI(Normal,2)];
 a1=[mmomALI(Normal,1)];
 a2=[mmomALI(Normal,2)];
     
 % [f0,x0]=ksdensity(a1);
 %[f1,x1]=ksdensity(a2);
 %A=(a1-a2).^2;
     FVnormal= [a1 a2];%sqrt(sum(A,2));
% figure(10),
% plot(a1(1,:),'r-')
% hold on
% plot(a1(2,:),'g-')
% 
% plot(a2(1,:),'b-')
% 
% plot(a2(2,:),'y-')
% [BB,III] = sortrows(a,1);
%FVnormal=abs(BB(2,1)- BB(2,2));
%a1=BB(2,1);% corr
%a2=BB(2,2);%std
%FVnormal=momentALI(Normal,2);
%FVnormal=momALI(Normal);


end