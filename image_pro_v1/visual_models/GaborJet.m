function GaborJet(stim)
% Implements GaborJet by Biederman
% Notes on this
% 256*256 images
% 10x10 grid w/ 20 pixels between nodes
% 40 Gabor filters at 5 scales and 8 orientations
% Requires squared images...

% Used magnitude output

% im    : image file
% cs    : complex (defualt) or simple cells    (0 1) 
% grid  : 10x10, 12x12 or image size (0 1 2)
% sigma : Gaussian envelope (2pi) is defualt

% Similarity measure
% dot(A,B)/(norm
bad_stim=zeros(1,length(stim.exp));
for ii=1:length(stim.exp)
    try
        im=fancy_read(stim,stim.exp{ii}); 
        % [sJetsMagnitude{ii},sJetsPhase{ii},~]=GWTWgrid_Simple(im,1,1,2*pi);   
        [cJetsMagnitude{ii},~,~]=GWTWgrid_Simple_bg(im,0,2,2*pi);    
        data{1}.col{ii}=stim.exp{ii};
    catch err
        bad_stim(ii)=1;
        continue;
    end
    if rem(ii,100)==0, display(ii); end
end

% Index 37 -> low contrast horozontal engergy
% Index 33 -> low contrast vertical engergy
% Index 1-8 -> high-contrast engergy

% Consider only the engery in the mask
I=imread('D:\Data\SAM\PresentationScripts\circle.jpg');
I=imresize(I,[256,256]); 
I(I<130)=0; I=im2bw(I); I=reshape(~I,1,[]);
keyboard;
for ii=1:length(stim.exp)
    Ijets=cJetsMagnitude{ii};
    He(ii)=mean(Ijets(I,37));
    Ve(ii)=mean(Ijets(I,33));
    Hc(ii)=mean(mean(Ijets(I,1:8)));
    Lc(ii)=mean(mean(Ijets(I,32:40)));
end

keyboard;
csv='D:\Data\SAM\PresentationScripts\stimuli00_sqr\image_list_io.csv';
X=excel_reader(csv); 
Indoor=cell2num(X{1}.col);

X=[Hc;Lc;He;Ve];
mdl=fitglm(X',Indoor,'distribution','binomial');

b=mean(mdl.coefCI');
Y=exp(X'*b(2:end)'+b(1))./(1+exp(X'*b(2:end)'+b(1)));
 
[~,a1]=min(Y)
imagesc(imread(fullfile(stim.dir,stim.exp{a1})))
[~,a2]=max(Y)
imagesc(imread(fullfile(stim.dir,stim.exp{a2})))
keyboard;
Yd=Y; Yd(Y>=0.5)=1; Yd(Y<0.5)=0;

[In,Ii]=sort(Y);

keyboard;

M=40; L=0;
for ii=1:96
Ijets=reshape(cJetsMagnitude{Ii(ii)},[256 256 40]);
figure(1);
subplot(3,4,1); imagesc(imread(fullfile(stim.dir,stim.exp{Ii(ii)})));
subplot(3,4,2); imagesc(Ijets(:,:,33)',[0 M]); title('Ve');
subplot(3,4,3); imagesc(Ijets(:,:,37)',[0 M]); title('He');
subplot(3,4,5); imagesc(Ijets(:,:,1+L)',[0 M]);
subplot(3,4,6); imagesc(Ijets(:,:,2+L)',[0 M]);
subplot(3,4,7); imagesc(Ijets(:,:,3+L)',[0 M]);
subplot(3,4,8); imagesc(Ijets(:,:,4+L)',[0 M]);
subplot(3,4,9); imagesc(Ijets(:,:,5+L)',[0 M]);
subplot(3,4,10); imagesc(Ijets(:,:,6+L)',[0 M]);
subplot(3,4,11); imagesc(Ijets(:,:,7+L)',[0 M]);
subplot(3,4,12); imagesc(Ijets(:,:,8+L)',[0 M]);
pause; close(gcf);
end

keyboard;
% Similarity
for ii=1:length(stim.exp)
    for jj=1:length(stim.exp)
        v1=reshape(cJetsMagnitude{ii},numel(cJetsMagnitude{ii}),1);
        v2=reshape(cJetsMagnitude{jj},numel(cJetsMagnitude{jj}),1);
        dot_val=sum(v1.*v2);
        norm1=norm(v1);
        norm2=norm(v2);
        R(ii,jj)=dot_val/(norm1*norm2);
    end
end

% Save models
if ~exist(fullfile(stim.save,'gbjet'),'dir'),
    mkdir(fullfile(stim.save,'gbjet'));
end
if ~exist(fullfile(stim.save,'gbjet_nan'),'dir'),
    mkdir(fullfile(stim.save,'gbjet_nan'));
end
save(fullfile(stim.save,'gbjet','model.mat'),'R','stim_ID_num');
x=eye(size(R));
R(x==1)=nan;
save(fullfile(stim.save,'gbjet_nan','model.mat'),'R','stim_ID_num');
