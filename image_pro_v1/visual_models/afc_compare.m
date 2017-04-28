function afc_compare(afc)

write_to='D:\Data\SEFER\Presentation\Stimuli06192014\afc\';
data=excel_reader(afc.file,afc.col,afc.head);
I1_dir='D:\Data\SEFER\Presentation\Stimuli06192014\stimuli\pro\';
stim_idx=strcmp(afc.head,'stim');
I2_dir='D:\Data\SEFER\Presentation\Stimuli06192014\exemplar\pro\';
xmpl_idx=strcmp(afc.head,'exemplar');
name_idx=strcmp(afc.head,'name');
catg_idx=strcmp(afc.head,'category');

category_all=data{catg_idx}.col;
category_u=unique(category_all);
for ii=1:length(category_u)
    if ~exist(fullfile(write_to,category_u{ii}),'dir')
        mkdir(fullfile(write_to,category_u{ii}));
    end
end

for ii=1:430,
    I1=imread(fullfile(I1_dir,data{stim_idx}.col{ii}));
    I2=imread(fullfile(I2_dir,data{xmpl_idx}.col{ii}));
    I=[I1,I2];
    cat=data{catg_idx}.col{ii};
    imwrite(I,fullfile(write_to,cat,[data{name_idx}.col{ii} '.jpg']));
end
% OI=length(data)+1;
% data{OI}.header='Old';
% data{OI}.col=nan(1,430);
% for ii=1:length(category_u)
%     lookin1=fullfile(write_to,category_u{ii},'Old');
%     old_images=dir([lookin1 filesep '*.jpg']);
%     for jj=1:length(old_images)
%         [~,name,~]=fileparts(old_images(jj).name);
%         data{OI}.col(find(strcmp(name,data{name_idx}.col)))=1;
%     end
%     lookin2=fullfile(write_to,category_u{ii},'New');
%     new_images=dir([lookin2 filesep '*.jpg']);
%     for jj=1:length(new_images)
%         [~,name,~]=fileparts(new_images(jj).name);
%         data{OI}.col(find(strcmp(name,data{name_idx}.col)))=0;
%     end
% end
% 
% write_struct(data,afc.file);