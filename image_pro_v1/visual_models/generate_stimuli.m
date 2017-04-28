function generate_stimuli(design)
%=========================================================================%
% Basic preperation
%=========================================================================%

% Setup the out directoriesjh
if ~exist(design.out_dir1,'dir'), mkdir(design.out_dir1); end
if ~exist(design.crp_dir,'dir'), mkdir(design.crp_dir); end
data=excel_reader(design.file,design.col,design.head);
% This is just the first 430 items as the database has been trimmed

file_idx=strcmp(design.head,design.set);

fimg=dir([design.raw_dir,'*.jpg']);
%=========================================================================%
% Setup output images
%=========================================================================%
setup_images=1; Resize_to=[500 500 3];
if setup_images==1
    % Crop them images
    for ii=1:460
        
        switch design.set
            case 'practice', Image=fullfile(design.raw_dir,fimg(ii).name);
            otherwise, Image=fullfile(design.raw_dir,data{file_idx}.col{ii});
        end

        if ~exist(Image,'file')
            display(['Missing: ' Image]); continue;
        end
        display([n2sp(ii,3) ': ' Image]);
        
        I=imread(Image);
        
        % Blurred image
        % bI=blur_3d_white(I,2,40);
        % imwrite(bI,fullfile(design.out_dir3,data{file_idx}.col{ii}));    
        
        CI=crop_image(I);
        % Find size of image and save image with square dimensions based
        % upon the largest dimension that the image has
        S(ii,:)=size(CI); x=S(ii,1); y=S(ii,2);

        if x>y
           CIs=resize_image(CI,[x,x,3]); 
        else
           CIs=resize_image(CI,[y,y,3]); 
        end
        switch design.set
            case 'practice', imwrite(CIs,fullfile(design.crp_dir,fimg(ii).name));
            otherwise, imwrite(CIs,fullfile(design.crp_dir,data{file_idx}.col{ii}));
        end
        
%         
        % Now take that cropped image and resize it to a small and large
        % dimension
        % IS = imresize(CIs, [195 195]);
        % IL = imresize(CIs, [465 465]);
        % RS=resize_image(IS,Resize_to);
        % RL=resize_image(IL,Resize_to);
        % imwrite(RS,fullfile(design.out_dir1,data{file_idx}.col{ii}));
        % imwrite(RL,fullfile(design.out_dir2,data{file_idx}.col{ii}));
        Ir = imresize(CIs, [450 450]);
        R=resize_image(Ir,Resize_to);        
         switch design.set
            case 'practice', imwrite(R,fullfile(design.out_dir1,fimg(ii).name));
            otherwise, imwrite(R,fullfile(design.out_dir1,data{file_idx}.col{ii}));
        end    
        
%         
    end

%     for ii=1:length(item_idx)
%        I=uint8(ones(Resize_to)*255); % White image 
%        save_name=fullfile(design.wrd_dir,data{file_idx}.col{ii});
%        Letter=data{wrd_idx}.col{ii};
%        [~,~,~,~]=AddTextToImage(I,Letter,[Resize_to(1)/2 Resize_to(2)/2],...
%             [1,1,1],'Serif',54,save_name,1);
%     end
end
return;
%=========================================================================%
% Write out code 
%=========================================================================%


    
% Now write out html code
fid_img_name=fullfile(design.txt_dir,'Mturk_arrays_img.txt');
fid_txt_name=fullfile(design.txt_dir,'Mturk_arrays_txt.txt');
fid_let1_name=fullfile(design.txt_dir,'Mturk_arrays_let1.txt');
fid_let2_name=fullfile(design.txt_dir,'Mturk_arrays_let2.txt');
fid_img=fopen(fid_img_name,'w+');
fid_txt=fopen(fid_txt_name,'w+');
fid_let1=fopen(fid_let1_name,'w+');
fid_let2=fopen(fid_let2_name,'w+');

for ii=1:4
    switch ii
        case 1
            f=fid_img; pull=file_idx;
            root_dir='http://cabezalab.org/wp-content/gallery/sefer_img/';
            La{1}='<div class="enc" id="enc_trials">';
            Lb{1}='\t<a class="enc_displays" id="e';
            Lb{2}='" href="';
            Lb{3}='">1</a>';
            Lc{1}='\t<a href="javascript:StartTrial()" id="startTrial">Start Trial</a>';
            Lc{2}='</div>';
        case 2
            f=fid_txt; pull=file_idx;
            root_dir='http://cabezalab.org/wp-content/gallery/sefer_wrd/';
            La{1}='<div class="ret" id="ret_trials">';
            Lb{1}='\t<a class="ret_displays" id="r';
            Lb{2}='" href="';
            Lb{3}='">1</a>';
            Lc{1}='</div>';
        case 3
            f=fid_let1;
            root_dir=''; pull=let_idx;
            La{1}='<div class="enc" id="enc_let1">';
            Lb{1}='\t<p class="enc_let_displays" id="la';
            Lb{2}='">';
            Lb{3}='</p>';
            Lc{1}='</div>';
        case 4
            f=fid_let2;
            root_dir=''; pull=let2_idx;
            La{1}='<div class="enc" id="enc_let2">';
            Lb{1}='\t<p class="enc_let_displays" id="lb';
            Lb{2}='">';
            Lb{3}='</p>';
            Lc{1}='</div>';
    end
    
    for jj=1:length(La), fprintf(f,[La{jj} '\n']); end
    for kk=1:length(item_idx)
        obj=[root_dir data{pull}.col{kk}];
        for jj=1:length(Lb)
            if jj==1, 
                fprintf(f,[Lb{jj} num2str(kk)]); 
            elseif jj==2
                fprintf(f,[Lb{jj} obj]);
            else
                fprintf(f,[Lb{jj} '\n']);
            end
        end 
    end
    for jj=1:length(Lc), fprintf(f,[Lc{jj} '\n']); end
    clear La Lb Lc
end
fclose('all');

end
%=========================================================================%
% Subfunctions (will be moved later)
%=========================================================================%







