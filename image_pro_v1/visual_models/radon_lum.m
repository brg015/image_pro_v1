function radon_lum(stim)

stim_ID_num=stim.stim_ID_num;
for ii=1:length(stim.exp) % Loop through images
    im=fancy_read(stim,stim.exp{ii});
    hsv=rgb2hsv(im);         % Capture lum
    lum=hsv(:,:,3);
    lum_array(:,ii)=reshape(lum,numel(lum),1);
    rad=radon(lum);
    rad_array(:,ii)=reshape(rad,numel(rad),1);
end

% Save data
R=corr(lum_array);
x=eye(size(R));
save(fullfile(stim.save,'lum','model.mat'),'R','stim_ID_num');
R(x==1)=nan;
save(fullfile(stim.save,'lum_nan','model.mat'),'R','stim_ID_num');
R=corr(rad_array);
save(fullfile(stim.save,'radon','model.mat'),'R','stim_ID_num');
R(x==1)=nan;
save(fullfile(stim.save,'radon_nan','model.mat'),'R','stim_ID_num');