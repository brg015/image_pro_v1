function transparent_text_overlay(design,rando)

if ~exist(design.out_dir,'dir'), mkdir(design.out_dir); end

Alphabet={'a' 'b' 'c' 'd' 'e' 'f' 'g' 'h' 'i' 'j' 'k' 'l' 'm' 'n' 'o' ...
    'p' 'q' 'r' 's' 't' 'u' 'v' 'w' 'x' 'y' 'z'};
write_images=0;
%=========================================================================%
% Load in stimuli data
%=========================================================================%
data=excel_reader(design.file,design.col,design.head);
% Collect index of "good" stimuli
item_idx=strcmp(data{strcmp(design.head,'set')}.col,'ok');
category=data{strcmp(design.head,'category')}.col(item_idx);
category_u=unique(category);
for ii=1:length(category_u)
    category_c(ii)=sum(strcmp(category,category_u{ii}));
end
keyboard
% Trial splits
enco=[22 18 29 23 20 25 36 22 23 29 22 31]; % ~75
retr=[8  7  9  8  8  8  10 8  8  9  8  9];  % ~25
ctch=[2  1  3  2  1  2  7  1  2  3  2  4];  % N=2
N_enco=sum(enco);
N_ctch=sum(ctch);
N_retr=sum(retr);
N_trials=N_enco+N_ctch+N_retr;
Assignment=nan(1,length(category));

code1_idx=strcmp(design.head,'code1');
code2_idx=strcmp(design.head,'code2');
ID_idx=strcmp(design.head,'ID');
file_idx=strcmp(design.head,'filename');
let_idx=strcmp(design.head,'letter');
for ii=1:N_trials % Trim extra spaces from letters
    if isspace(data{let_idx}.col{ii}(end))
        data{let_idx}.col{ii}(end)='';
    end
end
%=========================================================================%
% Determine splits of data 
%=========================================================================%
for ii=1:rando.nfiles

% Where to save files
out_dir=fullfile(design.out_dir,['set' num2str(ii)],filesep);
if ~exist(out_dir,'dir'), mkdir_tree(out_dir); end

stim_data{1}.header='File Idx';
stim_data{2}.header='Position';
stim_data{3}.header='Phase';      stim_data{3}.col(1:N_trials)=1;
stim_data{4}.header='Code1';
stim_data{5}.header='Code2';
stim_data{6}.header='Assingment'; 
stim_data{7}.header='jpg';
stim_data{8}.header='File Idx';
stim_data{9}.header='Position';
stim_data{10}.header='Phase';      stim_data{10}.col(1:400)=2;
stim_data{11}.header='Code1';
stim_data{12}.header='Code2';
stim_data{13}.header='Assignment';
stim_data{14}.header='jpg';
stim_data{15}.header='Letter';

% Sort out ER based upon category assignment
for jj=1:length(category_u)
    I=find(strcmp(category,category_u{jj}));
    A=randperm(sum(strcmp(category,category_u{jj})));
    enco_set=A(1:enco(jj));
    retr_set=A(enco(jj)+1:enco(jj)+retr(jj));
    ctch_set=A(enco(jj)+retr(jj)+1:end);
    Assignment(I(enco_set))=1;
    Assignment(I(retr_set))=3;
    Assignment(I(ctch_set))=2;
end

%=========================================================================%
% Determine presentation order
%=========================================================================%
Eorder1=randperm(N_enco+N_ctch);
Rorder2=random_lag(Eorder1,rando.lag,100);

c_ovr=N_enco+N_ctch+1; c_enc=1;
for jj=1:N_trials
    % Encoding
    if Assignment(jj)==1 || Assignment(jj)==2
        Eorder2(jj)=Eorder1(c_enc); c_enc=c_enc+1;
    else
        Eorder2(jj)=c_ovr; c_ovr=c_ovr+1;
    end
end

% Modify R order so that the catch trials are last and thus never presented
CR=1;
for jj=1:N_trials
    
    stim_data{1}.col(jj)=str2double(data{ID_idx}.col{(Eorder2==jj)});
    stim_data{2}.col(jj)=find(Eorder2==jj);
    stim_data{4}.col(jj)=str2double(data{code1_idx}.col{(Eorder2==jj)});
    stim_data{5}.col(jj)=str2double(data{code2_idx}.col{(Eorder2==jj)});   
    stim_data{6}.col(jj)=Assignment(Eorder2==jj);
    
    if Assignment(Rorder2==jj)~=2
        stim_data{8}.col(CR)=str2double(data{ID_idx}.col{(Rorder2==jj)});
        stim_data{9}.col(CR)=find(Rorder2==jj);
        stim_data{11}.col(CR)=str2double(data{code1_idx}.col{(Rorder2==jj)});
        stim_data{12}.col(CR)=str2double(data{code2_idx}.col{(Rorder2==jj)});
        stim_data{13}.col(CR)=Assignment(Rorder2==jj);
    end
    
    found=1;
    File_Read=fullfile(design.raw_dir,data{file_idx}.col{(Eorder2==jj)});
    save_name=fullfile(out_dir,data{file_idx}.col{(Eorder2==jj)});
    if ~exist(File_Read,'file')
        [pre,name,suf]=fileparts(File_Read);
        idx_under=strfind(name,'_');
        if ~isempty(idx_under), name=name(1:idx_under-1); end
        File_Read=fullfile(pre,[name,suf]);
    end
    if ~exist(File_Read,'file'),
        display(['DNE: ' File_Read]); found=1;
        File_Read=fullfile(design.raw_dir,'X.jpg');
    end

    Letter=data{let_idx}.col{(Eorder2==jj)};
    if Assignment(Eorder2==jj)==2 % Catch trials, compare single letters
        Lop=randperm(26);
        if ~strcmp(Alphabet{Lop(1)},Letter(1))
            Letter=Alphabet{Lop(1)};
        else
            Letter=Alphabet{Lop(2)};
        end
    end
    stim_data{15}.col{jj}=Letter;
    % Jpg files
    stim_data{7}.col{jj}=data{file_idx}.col{(Eorder2==jj)};
    
    if Assignment(Rorder2==jj)~=2
        stim_data{14}.col{CR}=data{file_idx}.col{(Rorder2==jj)};
        CR=CR+1;
    end

    % Load in the image file
    if write_images==1
        if (found==1 && ~exist(save_name,'file'))
            try
                C=imread(File_Read); 
                % Find center 
                [x,y,~]=size(C);
                x_center=round(x/2); y_center=round(y/2);
                [~,~,~,~]=AddTextToImage(C,Letter,[x_center,y_center],...
                    [255,255,255],'Serif',54,save_name,1);
                close all;
            catch err
                smart_err(err);
                display(['Failed: ' File_Read]);
                imwrite(C,save_name);
            end
        end
    end
end

write_struct_txt(stim_data,fullfile(out_dir,'stim_map.txt'));
write_struct(stim_data,fullfile(out_dir,'stim_map.csv'));
% Write out MTurk files
mturk_file=fopen(fullfile(out_dir,'mturk_var.txt'),'w+');
var1_str=['var order = new Array('];
var2_str=['var trialType = new Array('];
var3_str=['var blockStart = new Array('];
var4_str=['var catchTrial = new Array('];
O=[stim_data{2}.col(1:330)'; stim_data{9}.col(:)]; % Order
C=[stim_data{6}.col(1:330), zeros(1,400)];         % Catch
P=[ones(1,330), ones(1,400)*2];                    % Phase
S=zeros(1,760);                                     % Stop
S(1)=1;
S(165)=1;
S(331)=1;
S(431)=1;
S(531)=1;
S(631)=1;
for kk=1:length(O)
    var1_str=[var1_str num2str(O(kk)) ','];
    var2_str=[var2_str num2str(P(kk)) ','];
    var3_str=[var3_str num2str(S(kk)) ','];
    var4_str=[var4_str num2str(C(kk)) ','];
end
var1_str=[var1_str(1:end-1) ');\n'];
var2_str=[var2_str(1:end-1) ');\n'];
var3_str=[var3_str(1:end-1) ');\n'];
var4_str=[var4_str(1:end-1) ');\n'];
fprintf(mturk_file,var1_str);
fprintf(mturk_file,var2_str);
fprintf(mturk_file,var3_str);
fprintf(mturk_file,var4_str);
fclose(mturk_file);

end % Random file loop

end % File end

function R2=random_lag(R,lag,addto)
    lag_fit=0; c=0; MAX_Iterations=1000; L=length(R)+addto;
    while lag_fit==0,
        c=c+1;
        if c>MAX_Iterations,
            % Prevent endless repitition as problem can be impossible...
            error('Specified lag is too long!\n');
        end
        R2=randperm(L); 
        for ii=1:L,
            idx=find(R2(ii)==R);        % Find matching element in R1
            if ~isempty(idx)
                distance_diff(ii)=ii+L-idx; % Find distance to element
                % If distance is less then required lag, try again
                if distance_diff(ii)<=lag, break; end
            end
        end
        if min(distance_diff)>=lag, lag_fit=1; end % Distance always exceeds lag, R2 is found
    end
    fprintf('Randomization Facts\n')
    fprintf(['\tIterations: ' num2str(c) '\n']);
    fprintf(['\tMin Dstnce: ' num2str(min(distance_diff)) '\n']);
    fprintf(['\tAvg Dstnce: ' num2str(mean(distance_diff)) '\n']);
    fprintf(['\tMax Dstnce: ' num2str(max(distance_diff)) '\n']);
    
end
    