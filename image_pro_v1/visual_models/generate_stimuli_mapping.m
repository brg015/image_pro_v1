% BRG2014(Winter)
function generate_stimuli_mapping(stim)

Alphabet={'a' 'b' 'c' 'd' 'e' 'f' 'g' 'h' 'i' 'j' 'k' 'l' 'm' 'n' 'o' ...
    'p' 'q' 'r' 's' 't' 'u' 'v' 'w' 'x' 'y' 'z'};

stim_data=excel_reader(stim.file,stim.col,stim.head);

stim_data{length(stim.head)+1}.header='catch';
stim_data{length(stim.head)+2}.header='code1';
stim_data{length(stim.head)+3}.header='code2';

name_idx=strcmp(stim.head,'name');
file_idx=strcmp(stim.head,'stim');
xmpl_idx=strcmp(stim.head,'exemplar');
lett_idx=strcmp(stim.head,'letter');

%===================
% Generate catch letters
%===================
let=stim_data{lett_idx}.col;
let=let(~cellfun('isempty',let)); % Remove empites

% Double spaces?
for ii=1:length(let)
    if isspace(let{ii}(end)); let{ii}=let{ii}(1:end-1); end
end
for ii=1:length(let)
    if isspace(let{ii}(end)); let{ii}=let{ii}(1:end-1); end
end

% Now generate catch letters
lett2_idx=length(stim.head)+1;
for ii=1:length(let)
    Letter=stim_data{lett_idx}.col{ii};
    same_let=1;
    while same_let
        catch_letter=Alphabet(randi(26,1));
        if ~strcmp(catch_letter,Letter)
            same_let=0;
        end
    end
    stim_data{lett2_idx}.col{ii}=char(catch_letter);
end
%==================

% Open file to save array mappings
fid_img_name1=fullfile(stim.out_dir,'presentation_arrays_img1.txt');
name1dir='';
fid_img_name2=fullfile(stim.out_dir,'presentation_arrays_img2.txt');
name2dir='';
fid_txt_name=fullfile(stim.out_dir,'presentation_arrays_txt.txt');
fid_let1_name=fullfile(stim.out_dir,'presentation_arrays_let1.txt');
fid_let2_name=fullfile(stim.out_dir,'presentation_arrays_let2.txt');

fid_img1=fopen(fid_img_name1,'w+');
fid_img2=fopen(fid_img_name2,'w+');
fid_txt=fopen(fid_txt_name,'w+');
fid_let1=fopen(fid_let1_name,'w+');
fid_let2=fopen(fid_let2_name,'w+');

code1=1;
for ii=1:460
    % Determine code drops
    if rem(ii,stim.nmax)==0,
        code2=stim.nmax; 
    else
        code2=rem(ii,stim.nmax);
    end
    if (rem(ii,stim.nmax)==1 && ii>1), code1=code1+1; end    
    
    % First let's get the mappings
    stim_data{length(stim.head)+2}.col(ii)=code1+stim.noffset;
    stim_data{length(stim.head)+3}.col(ii)=code2+stim.noffset;
    
    % Setup the arrays
    if ii==1
        fprintf(fid_img1,'array {\n');
        fprintf(fid_img2,'array {\n');
        fprintf(fid_txt,'array {\n');
        fprintf(fid_let1,'array {\n');
        fprintf(fid_let2,'array {\n');
    end

    fprintf(fid_img1,['\tbitmap { filename = "' name1dir stim_data{file_idx}.col{ii} '";};\n']);
    fprintf(fid_img2,['\tbitmap { filename = "' name2dir stim_data{xmpl_idx}.col{ii} '";};\n']);
    fprintf(fid_txt,['\ttext { caption = "' stim_data{name_idx}.col{ii} '"; font_size=$L3;};\n']);
    fprintf(fid_let1,['\ttext { caption = "' stim_data{lett_idx}.col{ii} '"; font_size=$L4;};\n']);
    fprintf(fid_let2,['\ttext { caption = "' stim_data{lett2_idx}.col{ii} '"; font_size=$L4;};\n']);

    if ii==430;
        fprintf(fid_img1,'} pic1_stimuli;\n');
        fprintf(fid_img2,'} pic2_stimuli;\n');
        fprintf(fid_txt,'} txt_stimuli;\n');
        fprintf(fid_let1,'} let1_stimuli;\n');
        fprintf(fid_let2,'} let2_stimuli;\n');
    end
end
fclose('all');

write_struct(stim_data,fullfile(stim.out_dir,'stimuli_list_codes.csv'));
        
        
    
    

