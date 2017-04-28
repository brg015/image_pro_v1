function image_modifier(image_array,params)

% Make output directories if images are bing saved
if (~exist(params.crop.dir,'dir') && params.crop.wrte==1), mkdir(params.crop.dir); end
if (~exist(params.resz.dir,'dir') && params.resz.wrte==1), mkdir(params.resz.dir); end
if (~exist(params.sqre.dir,'dir') && params.sqre.wrte==1), mkdir(params.sqre.dir); end

for ii=1:length(image_array)

    Image=fullfile(params.raw_dir,image_array{ii});

    if ~exist(Image,'file')
        display(['Missing: ' Image]); continue;
    end
    display([n2sp(ii,3) ': ' Image]);

    % Load in Image
    I=imread(Image);
    
    % GOOD: Most of background removed
    CI=crop_image(I);
    if params.crop.wrte==1, imwrite(CI,fullfile(params.crop.dir,image_array{ii})); end
    
    % GOOD: Fixed to properly square image
    % * second parameter allows for addition of white space around the
    % image too
    CIs=square_image(CI,0);
    if params.sqre.wrte==1, imwrite(CIs,fullfile(params.sqre.dir,image_array{ii})); end
    
    % Very basic resize with built in scripts
    if params.resz.wrte==1, 
        IS = imresize(CIs, params.resz.size);
        imwrite(IS,fullfile(params.resz.dir,image_array{ii})); 
    end

end