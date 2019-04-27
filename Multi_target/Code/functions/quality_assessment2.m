function [stand_dev,He,ce,sf,fmi,fqi,fsm,wfqi,efqi] = quality_assessment2(MS_image,HS_image,Fused_image,ignore_edges)
%MS_Image is for Multi spectral image
%HS_Image is for Hyper spectral image
%Fused_Image is for Fused Hyper spectral image

% This is for referenceless metrics

% std = standard deviation
% He = entropy
% ce = cross entropy
% sf=spatial frequency

% Ignore borders
z = MS_image(ignore_edges+1:end-ignore_edges, ignore_edges+1:end-ignore_edges, :);
y = HS_image(ignore_edges+1:end-ignore_edges, ignore_edges+1:end-ignore_edges, :);
x = Fused_image(ignore_edges+1:end-ignore_edges, ignore_edges+1:end-ignore_edges, :);



% Size, bands, samples 
sz_x = size(x);
n_bands = sz_x(3);
n_samples = sz_x(1)*sz_x(2);

% mb = (mean(z(:)) - mean(y(:)))/mean(z(:))

% Standard deviation of 
stand_dev=0;
for i=1:n_bands
	stand_dev=stand_dev+ std2(x(:,:,i));
end	
stand_dev=stand_dev/n_bands;

% Entropy 
He=0;
for i=1:n_bands
	He=He+entropy(x(:,:,i));
end

He=He/n_bands;	

% %cross entropy - https://www.mathworks.com/help/nnet/ref/crossentropy.html
% ce_fm=0;
% ce_fh=0;
% for i=1:n_bands
% 	ce_fh=ce_fh+crossentropy(x(:,:,i),y(:,:,i));
	
% end
% ce_fh=ce_fh/n_bands;	

% size_z = size(z)

% for i=1:size_z(3)
% 	ce_fm=ce_fm+crossentropy(x(:,:,i),z(:,:,i));
	
% end
% ce_fm=ce_fm/size_z(3);	


% !ce_fm=crossentropy(x(:,:,i),z(:,:,i));
% ce=(ce_fm+ce_fh)/2;



% FMI- https://ieeexplore.ieee.org/document/7036000/metrics 
% and https://www.mathworks.com/matlabcentral/fileexchange/45926-feature-mutual-information--fmi--metric-for-non-reference-image-fusion

% Spatial frequency
sf=0;
mean_sf=0;
for i=1:n_bands
	sum_rf=0;
	sum_cf=0;
	for j=1:sz_x(1)
		for k=2:sz_x(2)
			sum_rf=sum_rf+norm(x(j,k,i)-x(j,k-1,i))^2;
			%sum_cf=sum_cf+norm(x(j,k,i)-x(j-1,k-1,i))^2;
		end
	end
	for j=2:sz_x(1)
		for k=1:sz_x(2)
			%sum_rf=sum_rf+norm(x(j,k,i)-x(j,k-1,i))^2;
			sum_cf=sum_cf+norm(x(j,k,i)-x(j-1,k,i))^2;
		end
	end
	cf=sqrt(sum_cf/n_samples);
	rf=sqrt(sum_rf/n_samples);
	sf=rf^2+cf^2;
	mean_sf=mean_sf+sf;		

end
mean_sf=mean_sf/n_bands;


% For other metrics-https://www.mathworks.com/help/images/image-quality-metrics.html



% % RMSE
% aux = sum(sum((x - y).^2, 1), 2)/n_samples;
% rmse_per_band = sqrt(aux);
% rmse = sqrt(sum(aux, 3)/n_bands);

% %sum(A,n) sum the matrix along that dimension n


% % ERGAS
% mean_y = sum(sum(y, 1), 2)/n_samples;
% ergas = 100*ratio_ergas*sqrt(sum((rmse_per_band ./ mean_y).^2)/n_bands);




% % SAM
% sam= SpectAngMapper( ground_truth, estimated );
% sam=sam*180/pi;
% % num = sum(x .* y, 3);
% % den = sqrt(sum(x.^2, 3) .* sum(y.^2, 3));
% % sam = sum(sum(acosd(num ./ den)))/(n_samples);

% % UIQI - calls the method described in "A Universal Image Quality Index"
% % by Zhou Wang and Alan C. Bovik
% q_band = zeros(1, n_bands);
% for idx1=1:n_bands
%     q_band(idx1)=img_qi(ground_truth(:,:,idx1), estimated(:,:,idx1), 32);
% end
% uiqi = mean(q_band);
% ssim=cal_ssim(ground_truth, estimated,0,0);
% DD=norm(ground_truth(:)-estimated(:),1)/numel(ground_truth);
% CCS = CC(ground_truth,estimated);
% CCS=mean(CCS);
% psnr=csnr(ground_truth, estimated,0,0);




%Mutual information
a=0;

for i=1:n_bands
	mi=mutualinformation(y(:,:,i),x(:,:,i));
	a=a+mi;
	if i==1
		entropy(y(:,:,i))
		entropy(x(:,:,i))
		mi
	end	


end
mi=a/n_bands;