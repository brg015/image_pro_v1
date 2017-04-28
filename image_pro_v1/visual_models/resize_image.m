function R=resize_image(I,Resize_to)

R=uint8(ones(Resize_to)*255);
S=size(I);

if S(1)<=Resize_to(1)
    X=round(Resize_to(1)/2-S(1)/2)+1:round(Resize_to(1)/2-S(1)/2)+S(1);
    Y=round(Resize_to(2)/2-S(2)/2)+1:round(Resize_to(2)/2-S(2)/2)+S(2);
    R(X,Y,:)=I;
else
    X=round(S(1)/2-Resize_to(1)/2)+1: round(S(1)/2-Resize_to(1)/2)+Resize_to(1);
    Y=round(S(2)/2-Resize_to(2)/2)+1: round(S(2)/2-Resize_to(2)/2)+Resize_to(2);
    
    R=I(X,Y,:);
end