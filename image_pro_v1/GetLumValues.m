%% Compute luminance
% see: 
% https://stackoverflow.com/questions/596216/formula-to-determine-brightness-of-rgb-color
image_path='';      % Add image path
image_name={'',''}; % List of image names

for iStim=1:length(im_name)
    X=double(imread(fullfile(image_path,image_name{ii})));
    % Step 1: Normalize
    R=reshape(X(:,:,1)/255,1,[]);
    G=reshape(X(:,:,2)/255,1,[]);
    B=reshape(X(:,:,3)/255,1,[]);
    % Step 2: Convert
    Rnew=NaN(1,length(R));
    Rnew(R<=0.04045)=R(R<=0.04054)/12.92;
    Rnew(R>0.04045)=((R(R>0.04045)+0.055)/1.055).^2.4;
    
    Gnew=NaN(1,length(G));
    Gnew(G<=0.04045)=G(G<=0.04054)/12.92;
    Gnew(G>0.04045)=((G(G>0.04045)+0.055)/1.055).^2.4;
    
    Bnew=NaN(1,length(B));
    Bnew(B<=0.04045)=B(B<=0.04054)/12.92;
    Bnew(B>0.04045)=((B(B>0.04045)+0.055)/1.055).^2.4;
    
    Y1=Rnew*0.299 + Gnew*0.7152 + Bnew*0.0722;
    % Perception
    Y2=NaN(1,length(Y1));
    Y2(Y1<=0.008856)=Y1(Y1<=0.008856)*903.3;
    Y2(Y1>0.008856)=(Y1(Y1>0.008856).^(1/3))*116-16;
    
    ImInfo(iStim).name=data1{3}.col{iStim};
    ImInfo(iStim).lum=mean2(Y1);
    ImInfo(iStim).perlum=mean2(Y2); % L*
    clear Y1 X R G B;
    clear Rnew Gnew Bnew Y2;
end
