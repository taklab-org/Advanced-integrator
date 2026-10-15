function [Wn]= getting_Wn(Phi,r1,m,nu,j,Un_endpt)
arguments
  Phi,r1,m,nu,j
  Un_endpt = eye
end
cheblength = size(Phi,1);
Un = intval(zeros(size(Phi)));
Un_sup = intval(zeros(size(Phi,[2,3])));
Phi = intval(Phi);

if j==1
  for i=1:m
    Un_sup(:,i) = chebmag_oneside(Phi(:,:,i));
  end
  Wn = wopnorm(Un_sup+r1,nu);
elseif j>1
  Phi(1,:,:) = Phi(1,:,:) + r1*infsup(-1,1);
  for i=1:m
    for j=1:m
      Un(:,i,j) = sum(permute(Phi(:,i,:),[1,3,2]).*intval(repmat(Un_endpt(:,j).',[cheblength,1])),2);
    end
  end
  for i=1:m
    Un_sup(:,i) = chebmag_oneside(Un(:,:,i));
  end
  Wn = wopnorm(Un_sup,nu);
end
%
