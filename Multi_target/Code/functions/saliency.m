function s=saliency(im1,w)

% Definition taken from wiki https://en.wikipedia.org/wiki/Saliency_map


[m,n]=size(im1);
sal=zeros(m,n);
for i=1:m
	for j=1:n
		for k=i:i+w
			for l=j:j+w	
				sal(i,j)=sal(i,j)+norm(im1(i,j)-im1(k,l),1)
			end
		end
	end
end				

