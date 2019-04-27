% by soleimani h.soleimani@ec.iut.ac.ir
%input---> im1 and im2... they should be in gray scale,[0 255], and have the same size
function MI=mutual_information_images(im1, im2)
im1=double(im1)+1;
im2=double(im2)+1;

% find joint histogram
joint_histogram=zeros(256,256);

for i=1:min(size(im1,1),size(im2,1))
    for j=1:min(size(im1,2),size(im2,2))
       joint_histogram(im1(i,j),im2(i,j))= joint_histogram(im1(i,j),im2(i,j))+1;
    end
end


 JPDF=joint_histogram/sum(joint_histogram(:)); % joint pdf of two images
 pdf_im1=sum(JPDF,1); % pdf of im1
 pdf_im2=sum(JPDF,2); % pdf of im2
 
 % find MI
 MI=0;
 for i=1:256
     for j=1:256
         if JPDF(i,j)>0
             MI=MI+JPDF(i,j)*log2(JPDF(i,j)/(pdf_im1(i)*pdf_im2(j)));
         end
     end
 end
end