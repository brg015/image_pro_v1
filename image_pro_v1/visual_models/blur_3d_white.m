function bI=blur_3d_white(I,M,addpad)

[x,y,~]=size(I);
I=resize_image(I,[x+addpad,y+addpad,3]);

Ia=squeeze(I(:,:,1));
Ib=squeeze(I(:,:,2));
Ic=squeeze(I(:,:,3));
f1=fspecial('gaussian',size(Ia),M);
Iaf=imfilter(Ia,f1);
Ibf=imfilter(Ib,f1);
Icf=imfilter(Ic,f1);

Ifshell=nan([size(Iaf),3]);
Ifshell(:,:,1)=Iaf;
Ifshell(:,:,2)=Ibf;
Ifshell(:,:,3)=Icf;
Ifshell=uint8(Ifshell);

bI=resize_image(Ifshell,[x,y,3]);

end