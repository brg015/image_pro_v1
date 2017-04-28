function R2=square_image(R,p)

% Create blank image with a little extra space/padding around the edges
R2=uint8(ones([max(size(R))+p,max(size(R))+p,3])*255);
S=size(R);

% Compute X assignments
X=round(size(R2,1)/2-S(1)/2)+1:round(size(R2,1)/2-S(1)/2)+S(1);
% Compute Y assignments
Y=round(size(R2,2)/2-S(2)/2)+1:round(size(R2,2)/2-S(2)/2)+S(2);
R2(X,Y,:)=R;
