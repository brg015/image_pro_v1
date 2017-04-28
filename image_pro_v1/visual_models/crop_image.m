function CI=crop_image(I)
% Assumptions
% 1) Image has white background (or very nearly so!)
% 2) Image has rgb basis
% 3) White is 255,255,255

% White limit, if average color is this, then it is considered that all
% pixels included in the mean are white
Wlim=250;
STDmin=0.5;
PAD=5;

% Convert Image to BW
BW=rgb2gray(I);
H=round(mean(BW,1)); 
Hx=std(double(BW),[],1); 
W=round(mean(BW,2));
Wx=std(double(BW),[],2); 

H_cropped=(H<Wlim | Hx>STDmin);
W_cropped=(W<Wlim | Wx>STDmin);

yv=pad_vector(H_cropped,PAD);
xv=pad_vector(W_cropped,PAD);
try
CI=I(xv,yv,:);
catch
    keyboard
end
end

