function im=fancy_read(stim,im)
    im=imread(fullfile(stim.dir,im));
    if stim.grey==1, im = rgb2gray(im); end
    if stim.downsample~=0, im=imresize(im,stim.downsample); end
    if stim.square==1,
        % set up cent
        im_size=size(im);
        [~,I]=max(im_size);
        resize_to=im_size(setdiff([1,2],I));
        center_to=round(im_size(I)/2);
        S=round(center_to-resize_to/2);
        sample_from=S+1:S+resize_to;
        if I==1
            im=im(sample_from,:);
        else
            im=im(:,sample_from);        
        end
    end
end