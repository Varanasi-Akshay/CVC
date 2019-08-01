function [C,S,level_zoom] = zoom_level(patch_img,basis)
img_ycbcr = rgb2ycbcr(patch_img);
s = size(img_ycbcr(:,:,1));
l = wmaxlev(s,basis);
[C,S] = wavedec2(img_ycbcr(:,:,1),l,basis);
details = ['h','v','d'];
threshold = 0.99;
level_zoom = 0; 
% disp(l)
Ct=C;
for i=1:l
%     disp('Removing')
    N = i;
    NC = wthcoef2(details(1),Ct,S,N);
    %dec=plotwavelet2(NC,S,l,basis,255,'square');
    NC = wthcoef2(details(2),NC,S,N);
%     dec=plotwavelet2(NC,S,l,basis,255,'square');
    NC = wthcoef2(details(3),NC,S,N);
%     dec=plotwavelet2(NC,S,l,basis,255,'square');

    low_res = waverec2(NC,S,basis);
    %low_res = reshape(low_res,s);

%     imtool(low_res,[]);
%     value = sum(NC.^2)/sum(C.^2);
    value = sumsqr(NC)/sumsqr(C);

    Ct=NC;
    if value < threshold
        level_zoom = l-i;
        break
    end    
    
end

end
