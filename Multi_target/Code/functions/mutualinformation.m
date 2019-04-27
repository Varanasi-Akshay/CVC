function MI=mutualinformation(im1, im2)
% by soleimani h.soleimani@ec.iut.ac.ir
%input---> im1 and im2... they should be in gray scale,[0 255], and have the same size


im1=double(im1)+1;
im2=double(im2)+1;

% MI =entropy(im1) + entropy(im2) - jointentropy(im1,im2)

%	find joint histogram
joint_histogram=zeros(256,256);

for i=1:min(size(im1,1),size(im2,1))
    for j=1:min(size(im1,2),size(im2,2))
       joint_histogram(im1(i,j),im2(i,j))= joint_histogram(im1(i,j),im2(i,j))+1;
    end
end


 JPDF=joint_histogram/sum(joint_histogram(:)); % joint pdf of two images
%  pdf_im1=sum(JPDF,1); % pdf of im1
%  pdf_im2=sum(JPDF,2); % pdf of im2
 
%  % find MI
%  MI=0;
%  for i=1:256
%      for j=1:256
%          if JPDF(i,j)>0
%              MI=MI+JPDF(i,j)*log2(JPDF(i,j)/(pdf_im1(i)*pdf_im2(j)));
%          end
%      end
%  end
% % end


indNoZero = joint_histogram ~= 0;
jointProb1DNoZero = JPDF(indNoZero);

jointEntropy = -sum(jointProb1DNoZero.*log2(jointProb1DNoZero));
MI=entropy(im1) + entropy(im2) - jointEntropy;

% indrow = double(im1(:)) + 1;
% indcol = double(im2(:)) + 1; %// Should be the same size as indrow


% [~,~,indrow] = unique(im1(:)); %// Change here
% [~,~,indcol] = unique(im2(:)); %// Change here
% jointHistogram = accumarray([indrow indcol], 1);
% jointProb = jointHistogram / numel(indrow);

% indNoZero = jointHistogram ~= 0;
% jointProb1DNoZero = jointProb(indNoZero);

% jointEntropy = -sum(jointProb1DNoZero.*log2(jointProb1DNoZero));
% MI=entropy(im1) + entropy(im2) - jointEntropy;