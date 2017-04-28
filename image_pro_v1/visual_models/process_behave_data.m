function process_behave_data(behav,design,ctch)
% Future Updates:
% Look at creating compiled measures w/ output
%=========================================================================%
%% Declare Variables
%=========================================================================%
run_filters={'base'};

switch behav.type
    case 'behav'
        col{1}={'B' 'F' 'G' 'J' 'L' 'M'};
        headers{1}={'block' 'code1' 'code2' 'RT' 'r_N' 'button'};
        col{2}={'B' 'F' 'G' 'K'  'M' 'N'};
        headers{2}={'block' 'code1' 'code2' 'RT'  'r_N' 'button'};
        col{3}={'B' 'F' 'G' 'K'  'M' 'N' 'D'};
        headers{3}={'block' 'code1' 'code2' 'RT' 'r_N' 'button' 'ans_afc' };    
    case 'eeg'
        % Output jitters as well, thus columns change
        col{1}={'B' 'F' 'G' 'L' 'N' 'O'};
        headers{1}={'block' 'code1' 'code2' 'RT' 'r_N' 'button'};
        col{2}={'B' 'F' 'G' 'M' 'O' 'P'};
        headers{2}={'block' 'code1' 'code2' 'RT'  'r_N' 'button'};
        col{3}={'B' 'F' 'G' 'M' 'O' 'P' 'D'};
        headers{3}={'block' 'code1' 'code2' 'RT' 'r_N' 'button' 'ans_afc' };
end

key_cor=0;
key_inc=1;
key_NR=-1;
key_ER=-2;
key_catch=[-2 5 6 9 10];
key_hi=[8 4 7 11];
key_lo=[9 5 6 10];

Ef_idx=1; 
R1f_idx=2; 
R2f_idx=3;
    
bugger=0;
%=========================================================================%
%% Declare Header Information
%=========================================================================%
N=430;
data=excel_reader(fullfile(design.file),design.col,design.head);
c1_idx=strcmp(design.head,'code1');
c2_idx=strcmp(design.head,'code2');
old_idx=strcmp(design.head,'Old');
catch_idx=strcmp(design.head,'Catch');
left_idx=strcmp(design.head,'L');
cat_idx=strcmp(design.head,'category');

for ii=1:N
    code_mat(ii,1)=str2double(data{c1_idx}.col{ii});
    code_mat(ii,2)=str2double(data{c2_idx}.col{ii});
    ans_ctc(ii)=str2double(data{catch_idx}.col{ii});
    ans_old(ii)=str2double(data{old_idx}.col{ii});
    ans_afc(ii)=str2double(data{left_idx}.col{ii});
end

category_u=unique(data{cat_idx}.col);
for ii=1:length(category_u),
    F{ii}.vect=strcmp(category_u{ii},data{cat_idx}.col);   
    F{ii}.name=category_u{ii}; F{ii}.d=1;
end
% Load in the following variables - I think these are fully sufficient for
% any analysis that I'd wish to conduct.
% code 1
% code 2
% RT
% r_N
% button

save_file{1}='enc_save.csv'; % Encoding file
save_file{2}='ret_save.csv'; % Retrieval file
save_file{3}='afc_save.csv'; % AFC file
sess{1}='1'; sess{2}='2'; sess{3}='2';

M_headers={'E1_catch_R' 'E1_catch_hit' 'E1_catch_RT' ...
        'R1_hc' 'R1_hit' 'R1_miss' 'R1_cr' 'R1_fa' 'R1_RT' ...
        'R2_hc' 'R2_hit' 'R2_miss' 'R2_cr' 'R2_fa' 'R2_RT' 'R2_s' ...
        'R1_e' 'R2_e' 'R1R2_ON' 'R1R2_NO' 'E_block'};
Mc={0,0,0,...
    0,[5 6],[5 6],[7 8],[7 8],0,...
    0,[11 12],[11 12],[13 14],[13 14],0,0,0,0,[19,20,22,23],[19,20,22,23],0,[19,20,22,23],[19,20,22,23]};
ML=length(M_headers);


%=========================================================================%
%% Subject specific processing
%=========================================================================%
subj_str='subject';
if ctch==1
    b1_idx=5; 
    for s=1:length(behav.subjects) 
        f=1; % Encoding file only
        load_file=fullfile(behav.dir,[subj_str behav.subjects{s} sess{f}],save_file{f});
        if exist(load_file,'file')
            temp=excel_reader(load_file,col{f},headers{f});
            Nt=length(temp{1}.col);
            c=1;
            for jj=1:length(temp)
                if ~strcmp(temp{jj}.header,'x')
                    for kk=1:Nt, data_mat{f}(kk,c)=str2double(temp{jj}.col{kk}); end
                    c=c+1;
                end
            end
        else
            data_mat{f}=[];
        end
        % We now have data_mat
        for jj=1:N
            [~,E_idx]=ismember(code_mat(jj,:),data_mat{Ef_idx}(:,2:3),'rows');
            if E_idx~=0
                M(jj,strcmp('E1_catch_R',M_headers))=~...
                    isempty(find(key_catch==data_mat{Ef_idx}(E_idx,b1_idx),1));
                % Check if actual check trial
                if (M(jj,strcmp('E1_catch_R',M_headers))==1 && ans_ctc(jj)==1), 
                    M(jj,strcmp('E1_catch_hit',M_headers))=1; 
                else
                    M(jj,strcmp('E1_catch_hit',M_headers))=0; 
                end
            end
        end
        C{s}=M;
    end
    
    catch_data{1}=data{4};
    catch_data{2}=data{8};
    catch_data{3}=data{9};
    catch_data{4}=data{13};
    for ii=1:length(C)
        catch_data{4+ii}.header=behav.subjects{ii};
        catch_data{4+ii}.col=C{ii}(:,1);
        if ii==1
            V=C{ii}(:,1);
        else
            V=C{ii}(:,1)+V;
        end
    end
    catch_data{4+length(C)+1}.header='Group';
    catch_data{4+length(C)+1}.col=V;
   % data{4}=name
   % data{8}=letter
   % data{9}=catch
   write_struct(catch_data,fullfile(behav.save_dir,'Catch_Report.csv'));
   return;
end
    
for s=1:length(behav.subjects) 
    switch behav.hc
        case 0
            save_name=fullfile(behav.save_dir,'pre',[subj_str behav.subjects{s} '.csv']);
            save_name_RT=fullfile(behav.save_dir,'pre',[subj_str behav.subjects{s} '_RT.csv']);
        case 1 % Only HC count as hits
            save_name=fullfile(behav.save_dir,'pre_hc',[subj_str behav.subjects{s} '.csv']);
            save_name_RT=fullfile(behav.save_dir,'pre_hc',[subj_str behav.subjects{s} '_RT.csv']);
        case 2 % Only HC count as hits
            save_name=fullfile(behav.save_dir,'pre_hc_only',[subj_str behav.subjects{s} '.csv']);
            save_name_RT=fullfile(behav.save_dir,'pre_hc_only',[subj_str behav.subjects{s} '_RT.csv']);
        case 3
            save_name=fullfile(behav.save_dir,'pre_hc_toss',[subj_str behav.subjects{s} '.csv']);
            save_name_RT=fullfile(behav.save_dir,'pre_hc_toss',[subj_str behav.subjects{s} '_RT.csv']);
        case 4
            save_name=fullfile(behav.save_dir,'pre_lc_only',[subj_str behav.subjects{s} '.csv']);
            save_name_RT=fullfile(behav.save_dir,'pre_lc_only',[subj_str behav.subjects{s} '_RT.csv']);
    end
    if ~exist(fullfile(behav.save_dir,'pre_lc_only'),'dir'),
        mkdir(fullfile(behav.save_dir,'pre_lc_only'));
    end
        
%     if exist(save_name,'file'), continue; end
    display(['  Preprocessing: ' behav.subjects{s}]);
    % Creating data_mats [Nt X (code1,code2,RT,button)] from the subjects
    % individual save files
    
    block_idx=1;
    RT_idx=4; 
    b1_idx=5; 
    b2_idx=6; 
    ans_afc_idx=7;

    for f=1:3
        load_file=fullfile(behav.dir,[subj_str behav.subjects{s} sess{f}],save_file{f});
        if exist(load_file,'file')
            temp=excel_reader(load_file,col{f},headers{f});
            Nt=length(temp{1}.col);
            c=1;
            for jj=1:length(temp)
                if ~strcmp(temp{jj}.header,'x')
                    for kk=1:Nt, data_mat{f}(kk,c)=str2double(temp{jj}.col{kk}); end
                    c=c+1;
                end
            end
        else
            data_mat{f}=[];
        end
    end
    %---------------------------------------------------------------------%
    % HC only analysis 
    %---------------------------------------------------------------------%
    if behav.hc==1
       % in data_mat{2} & data_mat{3}
       %   col 5 : 9 -> 6
       %   col 6 : 0 -> 1 (of same index)
       II=find(data_mat{2}(:,5)==9);
       data_mat{2}(II,5)=6;
       data_mat{2}(II,6)=1;
       II=find(data_mat{3}(:,5)==9);
       data_mat{3}(II,5)=6;
       data_mat{3}(II,6)=1;
    elseif behav.hc==2 % HC items only
        % in data_mat{2} & data_mat{3}
        %   col 5 : 9 -> 6
        %   col 6 : 0 -> 1 (of same index)
        % LC responses i.e. 5 6 9 10 are all NR now
        for uu=2:3
            II=find(data_mat{uu}(:,5)==5);
            data_mat{uu}(II,5)=-1;
            data_mat{uu}(II,6)=-1;
            II=find(data_mat{uu}(:,5)==6);
            data_mat{uu}(II,5)=-1;
            data_mat{uu}(II,6)=-1;
            II=find(data_mat{uu}(:,5)==9);
            data_mat{uu}(II,5)=-1;
            data_mat{uu}(II,6)=-1;
            II=find(data_mat{uu}(:,5)==10);
            data_mat{uu}(II,5)=-1;
            data_mat{uu}(II,6)=-1;
        end
    elseif behav.hc==3
       % in data_mat{2} & data_mat{3}
       %   col 5 : 9 -> 6
       %   col 6 : 0 -> 1 (of same index)
       % Here we're just throwing out instead of tossing
       II=find(data_mat{2}(:,5)==9);
       data_mat{2}(II,5)=-1;
       data_mat{2}(II,6)=-11;
       II=find(data_mat{3}(:,5)==9);
       data_mat{3}(II,5)=-1;
       data_mat{3}(II,6)=-1; 
    elseif behav.hc==4 % HC items only
        % in data_mat{2} & data_mat{3}
        %   col 5 : 9 -> 6
        %   col 6 : 0 -> 1 (of same index)
        % HC respones ie. 7 8 11 12 are all NR now
        for uu=2:3
            II=find(data_mat{uu}(:,5)==7);
            data_mat{uu}(II,5)=-1;
            data_mat{uu}(II,6)=-1;
            II=find(data_mat{uu}(:,5)==8);
            data_mat{uu}(II,5)=-1;
            data_mat{uu}(II,6)=-1;
            II=find(data_mat{uu}(:,5)==11);
            data_mat{uu}(II,5)=-1;
            data_mat{uu}(II,6)=-1;
            II=find(data_mat{uu}(:,5)==12);
            data_mat{uu}(II,5)=-1;
            data_mat{uu}(II,6)=-1;
        end
    end
    % setup a save file
    L=length(data); save_data=data; save_data_RT=data;
    for ii=1:length(M_headers), 
        save_data{L+ii}.header=M_headers{ii}; 
        save_data_RT{L+ii}.header=M_headers{ii}; 
    end
    for jj=1:ML, 
        save_data{L+jj}.col=nan(N,1); 
        save_data_RT{L+jj}.col=nan(N,1); 
    end; 
    M=nan(N,ML); M_RT=nan(size(M));
    
    for jj=1:N
        [~,E_idx]=ismember(code_mat(jj,:),data_mat{Ef_idx}(:,2:3),'rows');
        if ~isempty(data_mat{R1f_idx})
            [~,R1_idx]=ismember(code_mat(jj,:),data_mat{R1f_idx}(:,2:3),'rows');
            [~,R2_idx]=ismember(code_mat(jj,:),data_mat{R2f_idx}(:,2:3),'rows');
        else
            R1_idx=0; R2_idx=0; % Just says Ret unrun as of now.
        end
        
        if E_idx~=0
            % Check if subject responded
            M(jj,strcmp('E1_catch_R',M_headers))=~isempty(find(key_catch==data_mat{Ef_idx}(E_idx,b1_idx),1));
            % Check if actual check trial
            if (M(jj,strcmp('E1_catch_R',M_headers))==1 && ans_ctc(jj)==1), 
                M(jj,strcmp('E1_catch_hit',M_headers))=1; 
            else
                M(jj,strcmp('E1_catch_hit',M_headers))=0; 
            end
            % Save RT if available
            if M(jj,strcmp('E1_catch_R',M_headers))==1, 
                M(jj,strcmp('E1_catch_RT',M_headers))=data_mat{Ef_idx}(E_idx,RT_idx); 
            end
            M(jj,strcmp('E_block',M_headers))=data_mat{Ef_idx}(E_idx,block_idx);
        end

        if R1_idx~=0
            R=data_mat{R1f_idx}(R1_idx,b2_idx);
            if (R==0 || R==1) % Response, check for hi 
                M(jj,strcmp('R1_hc',M_headers))=~isempty(find(key_hi==data_mat{R1f_idx}(R1_idx,b1_idx),1));
                M(jj,strcmp('R1_RT',M_headers))=data_mat{R1f_idx}(R1_idx,RT_idx);
            end
            tf=0;
            if (ans_old(jj)==1 && R==key_cor), M(jj,strcmp('R1_hit',M_headers))=1; tf=1; else M(jj,strcmp('R1_hit',M_headers))=0; end
            if (ans_old(jj)==1 && R==key_inc), M(jj,strcmp('R1_miss',M_headers))=1; tf=1; else M(jj,strcmp('R1_miss',M_headers))=0; end
            if (ans_old(jj)==0 && R==key_cor), M(jj,strcmp('R1_cr',M_headers))=1; tf=1; else M(jj,strcmp('R1_cr',M_headers))=0; end
            if (ans_old(jj)==0 && R==key_inc), M(jj,strcmp('R1_fa',M_headers))=1; tf=1; else M(jj,strcmp('R1_fa',M_headers))=0; end
            if tf==0 
                M(jj,strcmp('R1_e',M_headers))=1; 
                M(jj,strcmp('R1_hit',M_headers))=NaN;
                M(jj,strcmp('R1_miss',M_headers))=NaN;
                M(jj,strcmp('R1_cr',M_headers))=NaN;
                M(jj,strcmp('R1_fa',M_headers))=NaN;
            end
        end
        if R2_idx~=0
            R=data_mat{R2f_idx}(R2_idx,b2_idx);
            if (R==0 || R==1) % Response, check for hi 
                M(jj,strcmp('R2_hc',M_headers))=~isempty(find(key_hi==data_mat{R2f_idx}(R2_idx,b1_idx),1));
                M(jj,strcmp('R2_RT',M_headers))=data_mat{R2f_idx}(R2_idx,RT_idx);
            end
            M(jj,strcmp('R2_hit',M_headers))=0; 
            M(jj,strcmp('R2_miss',M_headers))=0; 
            M(jj,strcmp('R2_cr',M_headers))=0; 
            M(jj,strcmp('R2_fa',M_headers))=0; 
            % We're reading from the data file as compared to the master as
            % this has changed between subjects.
            afc_resp=data_mat{R2f_idx}(R2_idx,ans_afc_idx); %ans_afc(jj)
            M(jj,strcmp('R2_s',M_headers))=afc_resp;
            tf=0;
            if ((ans_old(jj)==1 && afc_resp==1) && R==key_cor), M(jj,strcmp('R2_hit',M_headers))=1; tf=1; end % Hit
            if ((ans_old(jj)==1 && afc_resp==0) && R==key_cor), M(jj,strcmp('R2_cr',M_headers))=1; tf=1; end % CR
            if ((ans_old(jj)==1 && afc_resp==1) && R==key_inc), M(jj,strcmp('R2_miss',M_headers))=1; tf=1; end % Miss
            if ((ans_old(jj)==1 && afc_resp==0) && R==key_inc), M(jj,strcmp('R2_fa',M_headers))=1; tf=1; end % FA
            
            if ((ans_old(jj)==0 && afc_resp==1) && R==key_cor), M(jj,strcmp('R2_cr',M_headers))=1; tf=1; end % CR
            if ((ans_old(jj)==0 && afc_resp==0) && R==key_cor), M(jj,strcmp('R2_cr',M_headers))=1; tf=1; end % CR
            if ((ans_old(jj)==0 && afc_resp==1) && R==key_inc), M(jj,strcmp('R2_fa',M_headers))=1; tf=1; end % FA
            if ((ans_old(jj)==0 && afc_resp==0) && R==key_inc), M(jj,strcmp('R2_fa',M_headers))=1; tf=1; end % FA
            if tf==0 
                M(jj,strcmp('R2_e',M_headers))=1; 
                M(jj,strcmp('R2_hit',M_headers))=NaN;
                M(jj,strcmp('R2_miss',M_headers))=NaN;
                M(jj,strcmp('R2_cr',M_headers))=NaN;
                M(jj,strcmp('R2_fa',M_headers))=NaN;
            end
        end
        %=================================================================%
        % Incongruency Check
        %=================================================================%
        try
            R1_old=M(jj,strcmp('R1_hit',M_headers)) || M(jj,strcmp('R1_fa',M_headers));
            R1_new=M(jj,strcmp('R1_miss',M_headers)) || M(jj,strcmp('R1_cr',M_headers));
            R2_old=M(jj,strcmp('R2_hit',M_headers)) || M(jj,strcmp('R2_fa',M_headers));
            R2_new=M(jj,strcmp('R2_miss',M_headers)) || M(jj,strcmp('R2_cr',M_headers));
            M(jj,strcmp('R1R2_ON',M_headers))=(R1_old && R2_new);
            M(jj,strcmp('R1R2_NO',M_headers))=(R1_new && R2_old);
        catch err
            M(jj,strcmp('R1R2_ON',M_headers))=NaN;
            M(jj,strcmp('R1R2_NO',M_headers))=NaN;
        end
    end
    %=====================================================================%
    % Create an RT compatible template as well... it has the exact same
    % format as M, but contains RTs instead
    %=====================================================================%
    RT_idx=~cellfun('isempty',(strfind(M_headers,'RT')));
    E1_idx=setdiff(find(~cellfun('isempty',(strfind(M_headers,'E1')))),find(RT_idx));
    R1_idx=setdiff(find(~cellfun('isempty',(strfind(M_headers,'R1')))),find(RT_idx));
    R2_idx=setdiff(find(~cellfun('isempty',(strfind(M_headers,'R2')))),find(RT_idx));
    
    E1_RTs=M(:,strcmp('E1_catch_RT',M_headers));
    R1_RTs=M(:,strcmp('R1_RT',M_headers));
    R2_RTs=M(:,strcmp('R2_RT',M_headers));
    
    for ii=1:length(E1_idx), M_RT(:,E1_idx(ii))=M(:,E1_idx(ii)).*E1_RTs; end
    for ii=1:length(R1_idx), M_RT(:,R1_idx(ii))=M(:,R1_idx(ii)).*R1_RTs; end
    for ii=1:length(R2_idx), M_RT(:,R2_idx(ii))=M(:,R2_idx(ii)).*R2_RTs; end
    M_RT(M_RT==0)=nan;

    for jj=1:ML, 
        save_data{L+jj}.col=M(:,jj); 
        save_data_RT{L+jj}.col=M_RT(:,jj); 
    end
    %=====================================================================%
    % Short mod to remove catch trials and fix ON of template files
    %=====================================================================%
    % Need to fix .mat and .csv files, lets fix the .mat and use it to
    % update the .csv after
    rm_catch=(M(:,1)+M(:,2))>0;
    rm_columns=[4:8 10:14];
    % 1) gotta kill catch trials
    for ii=rm_columns
        M(rm_catch,ii)=NaN;
    end
    % 2) Then update measures
    R1_hit=M(:,5)==1;
    R1_miss=M(:,6)==1;
    R2_hit=M(:,11)==1;
    R2_miss=M(:,12)==1;
    R2_s=M(:,16)==1;
    M(:,19)=and(and(R1_hit,R2_miss),R2_s);
    M(:,20)=and(and(R1_miss,R2_hit),R2_s);
    M(:,22)=and(and(R1_hit,R2_hit),R2_s);
    M(:,23)=and(and(R1_miss,R2_miss),R2_s);
    % 4) Update save_data
    save_data{19}.col=M(:,5);
    save_data{20}.col=M(:,6);
    save_data{21}.col=M(:,7);
    save_data{22}.col=M(:,8);
    save_data{25}.col=M(:,11);
    save_data{26}.col=M(:,12);
    save_data{27}.col=M(:,13);
    save_data{28}.col=M(:,14);
    save_data{33}.col=M(:,19);
    save_data{34}.col=M(:,20);
    save_data{36}.header='R1R2_OO'; save_data{36}.col=M(:,22);
    save_data{37}.header='R1R2_NN'; save_data{37}.col=M(:,23);
    clear R1_hit R1_miss R2_hit R2_miss R2_s;
    write_struct(save_data,save_name);
    write_struct(save_data_RT,save_name_RT);
    switch behav.hc
        case 0
            save(fullfile(behav.save_dir,'pre',[subj_str behav.subjects{s} '.mat']),'M');
            save(fullfile(behav.save_dir,'pre',[subj_str behav.subjects{s} '_RT.mat']),'M_RT');
        case 1
            save(fullfile(behav.save_dir,'pre_hc',[subj_str behav.subjects{s} '.mat']),'M');
            save(fullfile(behav.save_dir,'pre_hc',[subj_str behav.subjects{s} '_RT.mat']),'M_RT');
        case 2
            save(fullfile(behav.save_dir,'pre_hc_only',[subj_str behav.subjects{s} '.mat']),'M');
            save(fullfile(behav.save_dir,'pre_hc_only',[subj_str behav.subjects{s} '_RT.mat']),'M_RT');    
        case 3
            save(fullfile(behav.save_dir,'pre_hc_toss',[subj_str behav.subjects{s} '.mat']),'M');
            save(fullfile(behav.save_dir,'pre_hc_toss',[subj_str behav.subjects{s} '_RT.mat']),'M_RT');  
    end
    clear save_data save_data_RT  M_RT M data_mat
end
%=========================================================================%
%% Save group data
%=========================================================================%
% 3) add new headers
M_headers{22}='R1R2_OO';
M_headers{23}='R1R2_NN';
ML=ML+2;

L=length(data); save_data=data; 
for ii=1:length(M_headers), 
    save_data{L+ii}.header=M_headers{ii}; 
end
for jj=1:ML, 
    save_data{L+jj}.col=nan(N,1);  
end; 
Mgroup=nan(N,ML,length(behav.subjects)); 

for s=1:length(behav.subjects)
    switch behav.hc
        case 0
            load(fullfile(behav.save_dir,'pre',[subj_str behav.subjects{s} '.mat']));
        case 1
            load(fullfile(behav.save_dir,'pre_hc',[subj_str behav.subjects{s} '.mat']));
        case 2
            load(fullfile(behav.save_dir,'pre_hc_only',[subj_str behav.subjects{s} '.mat']));
        case 3
            load(fullfile(behav.save_dir,'pre_hc_toss',[subj_str behav.subjects{s} '.mat']));
    end
    try
    Mgroup(:,:,s)=M;
    catch err
        keyboard
    end
end

M=nansum(Mgroup,3);

 for jj=1:ML, 
    save_data{L+jj}.col=M(:,jj); 
 end
switch behav.hc
    case 0
        write_struct(save_data,fullfile(behav.save_dir,'pre',[subj_str '00.csv']));
        save(fullfile(behav.save_dir,'pre',[subj_str '00.mat']),'M');
    case 1
        write_struct(save_data,fullfile(behav.save_dir,'pre_hc',[subj_str '00.csv']));
        save(fullfile(behav.save_dir,'pre_hc',[subj_str '00.mat']),'M');
    case 2
        write_struct(save_data,fullfile(behav.save_dir,'pre_hc_only',[subj_str '00.csv']));
        save(fullfile(behav.save_dir,'pre_hc_only',[subj_str '00.mat']),'M');
    case 3
        write_struct(save_data,fullfile(behav.save_dir,'pre_hc_toss',[subj_str '00.csv']));
        save(fullfile(behav.save_dir,'pre_hc_toss',[subj_str '00.mat']),'M');
    case 4
        write_struct(save_data,fullfile(behav.save_dir,'pre_lc_only',[subj_str '00.csv']));
        save(fullfile(behav.save_dir,'pre_lc_only',[subj_str '00.mat']),'M');
end
% behav.subjects{length(behav.subjects)+1}='00';

%=========================================================================%
%% Subject filters and save templates
%=========================================================================%
FL=length(F);
% run_filters={'Fo'};
% 'FcFsFR1FR2onFR2s' 'FcFsFR1R2noFR2s' 

% run_filters=[];



% for ii=1:length(M_headers), template_RT{ii}.header=M_headers{ii}; end
% template_RT{ML+1}.header='ID';

P=1;
plot_set{P}.v=[5:8;11:14];
plot_set{P}.name='dprime';
plot_set{P}.legend={'Wdprime(O)' 'Pdprime(O)'};
plot_set{P}.type='bar';
plot_set{P}.i=2;
plot_set{P}.dprime=1;
plot_set{P}.set='values';
P=P+1;

plot_set{P}.v=[5:8,11:14];
plot_set{P}.name='Percent';
plot_set{P}.legend={'Whit' 'Wmiss' 'Wcr' 'Wfa' 'Phit' 'Pmiss' 'Pcr' 'Pfa'};
plot_set{P}.type='bar';
plot_set{P}.i=2;
plot_set{P}.dprime=0;
plot_set{P}.set='values';
P=P+1;

plot_set{P}.v=[5:8,11:14];
plot_set{P}.name='N';
plot_set{P}.legend={'Whit' 'Wmiss' 'Wcr' 'Wfa' 'Phit' 'Pmiss' 'Pcr' 'Pfa'};
plot_set{P}.type='bar';
plot_set{P}.i=4;
plot_set{P}.dprime=0;
plot_set{P}.set='values';
P=P+1;

plot_set{P}.v=[5:8,11:14];
plot_set{P}.name='RTs';
plot_set{P}.legend={'Whit' 'Wmiss' 'Wcr' 'Wfa' 'Phit' 'Pmiss' 'Pcr' 'Pfa'};
plot_set{P}.type='bar';
plot_set{P}.i=2;
plot_set{P}.dprime=0;
plot_set{P}.set='RT';
P=P+1;

% plot_set{P}.v=[1:2];
% plot_set{P}.name='Encoding_N';
% plot_set{P}.legend={'R' 'RR'};
% plot_set{P}.type='bar';
% plot_set{P}.i=4;
% plot_set{P}.dprime=0;
% plot_set{P}.set='values';
% P=P+1;
% 
% plot_set{P}.v=[17:18];
% plot_set{P}.name='Errors';
% plot_set{P}.legend={'R1_e' 'R2_e'};
% plot_set{P}.type='bar';
% plot_set{P}.i=4;
% plot_set{P}.dprime=0;
% plot_set{P}.set='values';
% P=P+1;

category_u{13}='';
for ii=1:13 
    model{ii,1}={category_u{ii},'Fc','Fs'};
    model{ii,2}={category_u{ii},'Fc','Fs','R1R2ee','FR2s'};   
    model{ii,3}={category_u{ii},'Fc','Fs','FR1R2on','FR2s'};
    model{ii,4}={category_u{ii},'Fc','Fs','FR1R2no','FR2s'};
    model{ii,5}={category_u{ii},'Fc','Fs','FR2s'};
    model{ii,6}={category_u{ii}};
    model{ii,7}={category_u{ii},'b1'};
    model{ii,8}={category_u{ii},'b2'};
    model{ii,9}={category_u{ii},'Fc','Fs','Fo'};
    model{ii,10}={category_u{ii},'Fo'};
end
category_u{13}='All';




File1_h={'Subject' 'HC_O1' 'HC_FA1' 'HC_CR1' 'HC_M1' 'dp_HC_O1' 'dp_HC_CR1' ...
    'HC_O2' 'HC_FA2' 'HC_CR2' 'HC_M2' 'dp_HC_O2' 'dp_HC_CR2' ...
    'LC_O1' 'LC_FA1' 'LC_CR1' 'LC_M1' 'dp_LC_01' 'dp_LC_CR1' ...
    'LC_O2' 'LC_FA2' 'LC_CR2' 'LC_M2' 'dp_LC_02' 'dp_LC_CR2'};
%=========================================================================%
%% Filter and save data
%=========================================================================%
for s=1:length(behav.subjects)
    sdisp(['Analyzing: ' behav.subjects{s}],1); 
    template_RT{ML+1}.col{1}=behav.subjects{s};
    %================%
    % Loads in M
    %================%
    switch behav.hc
        case 0
            load(fullfile(behav.save_dir,'pre',[subj_str behav.subjects{s} '.mat']));
            load(fullfile(behav.save_dir,'pre',[subj_str behav.subjects{s} '_RT.mat']));
        case 1
            load(fullfile(behav.save_dir,'pre_hc',[subj_str behav.subjects{s} '.mat']));
            load(fullfile(behav.save_dir,'pre_hc',[subj_str behav.subjects{s} '_RT.mat']));
        case 2
            load(fullfile(behav.save_dir,'pre_hc_only',[subj_str behav.subjects{s} '.mat']));
            load(fullfile(behav.save_dir,'pre_hc_only',[subj_str behav.subjects{s} '_RT.mat']));
        case 3
            load(fullfile(behav.save_dir,'pre_hc_toss',[subj_str behav.subjects{s} '.mat']));
            load(fullfile(behav.save_dir,'pre_hc_toss',[subj_str behav.subjects{s} '_RT.mat']));
          case 4
            load(fullfile(behav.save_dir,'pre_lc_only',[subj_str behav.subjects{s} '.mat']));
            load(fullfile(behav.save_dir,'pre_lc_only',[subj_str behav.subjects{s} '_RT.mat']));
    end
    %==========================%
    % Subject specific filters
    %==========================%
    % M_headers -> 23 IDs
    % M         -> All IDs
    % M_RT      -> Doesn't have final two IDs
    %================%
    % Memory Measures A'
    %================%
    % HCO / (HCO+HCFA)    => A
    % HCFA / (HCO+HCFA)   => B
    % HCCR / (HCN+HCCR)    => C
    % HCM / (HCN+HCCR)   => D
   
    M(isnan(M))=0;
    R1_HCO =and(logical(M(:,strcmp('R1_hc',M_headers))),...
        logical(M(:,strcmp('R1_hit',M_headers))));
    R1_HCFA=and(logical(M(:,strcmp('R1_hc',M_headers))),...
        logical(M(:,strcmp('R1_fa',M_headers))));
    R1_HCM =and(logical(M(:,strcmp('R1_hc',M_headers))),...
        logical(M(:,strcmp('R1_miss',M_headers))));
    R1_HCCR=and(logical(M(:,strcmp('R1_hc',M_headers))),...
        logical(M(:,strcmp('R1_cr',M_headers))));

    R2_HCO =and(logical(M(:,strcmp('R2_hc',M_headers))),...
        logical(M(:,strcmp('R2_hit',M_headers))));
    R2_HCFA=and(logical(M(:,strcmp('R2_hc',M_headers))),...
        logical(M(:,strcmp('R2_fa',M_headers))));
    R2_HCM =and(logical(M(:,strcmp('R2_hc',M_headers))),...
        logical(M(:,strcmp('R2_miss',M_headers))));
    R2_HCCR=and(logical(M(:,strcmp('R2_hc',M_headers))),...
        logical(M(:,strcmp('R2_cr',M_headers))));
    
    A_R1=sum(R1_HCO)/(sum(R1_HCO)+sum(R1_HCFA));
    B_R1=sum(R1_HCFA)/(sum(R1_HCO)+sum(R1_HCFA));
    C_R1=sum(R1_HCCR)/(sum(R1_HCCR)+sum(R1_HCM));
    D_R1=sum(R1_HCM)/(sum(R1_HCCR)+sum(R1_HCM));
    % AB_R1_dprime=dprime(A_R1,B_R1);
    % CD_R1_dprime=dprime(C_R1,D_R1);
    AB_R1_dprime=dprime_calc(R1_HCO,R1_HCFA);
    CD_R1_dprime=dprime_calc(R1_HCCR,R1_HCM);
    
    A_R2=sum(R2_HCO)/(sum(R2_HCO)+sum(R2_HCFA));
    B_R2=sum(R2_HCFA)/(sum(R2_HCO)+sum(R2_HCFA));
    C_R2=sum(R2_HCCR)/(sum(R2_HCCR)+sum(R2_HCM));
    D_R2=sum(R2_HCM)/(sum(R2_HCCR)+sum(R2_HCM));
    % AB_R2_dprime=dprime(A_R2,B_R2);
    % CD_R2_dprime=dprime(C_R2,D_R2);
    AB_R2_dprime=dprime_calc(R2_HCO,R2_HCFA);
    CD_R2_dprime=dprime_calc(R2_HCCR,R2_HCM);
    
     %================%
    % Memory Measures B'
    %================%
    % LCO / (LCO+LCFA)    => E
    % LCFA / (LCO+LCFA)   => F
    % LCCR / (LCN+LCCR)   => G
    % LCM / (LCN+LCCR)    => H
    R1_LCO =and(~logical(M(:,strcmp('R1_hc',M_headers))),...
        logical(M(:,strcmp('R1_hit',M_headers))));
    R1_LCFA=and(~logical(M(:,strcmp('R1_hc',M_headers))),...
        logical(M(:,strcmp('R1_fa',M_headers))));
    R1_LCM =and(~logical(M(:,strcmp('R1_hc',M_headers))),...
        logical(M(:,strcmp('R1_miss',M_headers))));
    R1_LCCR=and(~logical(M(:,strcmp('R1_hc',M_headers))),...
        logical(M(:,strcmp('R1_cr',M_headers))));
    
    R2_LCO =and(~logical(M(:,strcmp('R2_hc',M_headers))),...
        logical(M(:,strcmp('R2_hit',M_headers))));
    R2_LCFA=and(~logical(M(:,strcmp('R2_hc',M_headers))),...
        logical(M(:,strcmp('R2_fa',M_headers))));
    R2_LCM =and(~logical(M(:,strcmp('R2_hc',M_headers))),...
        logical(M(:,strcmp('R2_miss',M_headers))));
    R2_LCCR=and(logical(M(:,strcmp('R2_hc',M_headers))),...
        logical(M(:,strcmp('R2_cr',M_headers))));
    
    E_R1=sum(R1_LCO)/(sum(R1_LCO)+sum(R1_LCFA));
    F_R1=sum(R1_LCFA)/(sum(R1_LCO)+sum(R1_LCFA));
    G_R1=sum(R1_LCCR)/(sum(R1_LCCR)+sum(R1_LCM));
    H_R1=sum(R1_LCM)/(sum(R1_LCCR)+sum(R1_LCM));
    % EF_R1_dprime=dprime(E_R1,F_R1);
    % GH_R1_dprime=dprime(G_R1,H_R1);
    EF_R1_dprime=dprime_calc(R1_LCO,R1_LCFA);
    GH_R1_dprime=dprime_calc(R1_LCCR,R1_LCM);
    
    E_R2=sum(R2_LCO)/(sum(R2_LCO)+sum(R2_LCFA));
    F_R2=sum(R2_LCFA)/(sum(R2_LCO)+sum(R2_LCFA));
    G_R2=sum(R2_LCCR)/(sum(R2_LCCR)+sum(R2_LCM));
    H_R2=sum(R2_LCM)/(sum(R2_LCCR)+sum(R2_LCM));
    % EF_R2_dprime=dprime(E_R2,F_R2);
    % GH_R2_dprime=dprime(G_R2,H_R2);
    EF_R2_dprime=dprime_calc(R2_LCO,R2_LCFA);
    GH_R2_dprime=dprime_calc(R2_LCCR,R2_LCM);

    File1_v(s,:)=[A_R1 B_R1 C_R1 D_R1 AB_R1_dprime CD_R1_dprime ...
        A_R2 B_R2 C_R2 D_R2 AB_R2_dprime CD_R2_dprime ...
        E_R1 F_R1 G_R1 H_R1 EF_R1_dprime GH_R1_dprime ...
        E_R2 F_R2 G_R2 H_R2 EF_R2_dprime GH_R2_dprime];
    
    fi=FL+1; 

    F{fi}.vect=~ans_ctc;   F{fi}.name='Fc'; F{fi}.d=2;   fi=fi+1;

    F{fi}.vect=M(:,strcmp('E1_catch_R',M_headers));    F{fi}.name='Fs'; F{fi}.d=2;   
    F{fi}.vect(isnan(F{fi}.vect))=0; F{fi}.vect=~F{fi}.vect; fi=fi+1;
    
    F{fi}.vect=ans_old;   F{fi}.name='Fo'; F{fi}.d=2;   fi=fi+1;
    F{fi}.vect=~ans_old;  F{fi}.name='Fn'; F{fi}.d=2;   fi=fi+1;
    
    F{fi}.vect=M(:,strcmp('R1_hc',M_headers));    
    F{fi}.name='FR1c'; F{fi}.d=2; 
    F{fi}.vect(isnan(F{fi}.vect))=0; fi=fi+1;
    
    F{fi}.vect=(M(:,strcmp('R1_hit',M_headers))+M(:,strcmp('R1_fa',M_headers)));    
    F{fi}.vect(isnan(F{fi}.vect))=0; F{fi}.name='FR1o'; F{fi}.d=2; fi=fi+1;
    
    F{fi}.vect=(M(:,strcmp('R1_cr',M_headers))+M(:,strcmp('R1_miss',M_headers)));    
    F{fi}.vect(isnan(F{fi}.vect))=0; F{fi}.name='FR1n'; F{fi}.d=2; fi=fi+1;
    
    F{fi}.vect=M(:,strcmp('R2_s',M_headers));   
    F{fi}.name='FR2s'; F{fi}.d=2; 
    F{fi}.vect(isnan(F{fi}.vect))=0; fi=fi+1;
    
    F{fi}.vect=~F{fi-1}.vect;  
    F{fi}.name='FR2e'; F{fi}.d=2; fi=fi+1;
    
    F{fi}.vect=M(:,strcmp('R1R2_ON',M_headers));
    F{fi}.vect(isnan(F{fi}.vect))=0;
    F{fi}.name='FR1R2on'; F{fi}.d=2; fi=fi+1;
    
    F{fi}.vect=M(:,strcmp('R1R2_NO',M_headers));
    F{fi}.vect(isnan(F{fi}.vect))=0;
    F{fi}.name='FR1R2no'; F{fi}.d=2; fi=fi+1;

    Fsum=M(:,strcmp('R1R2_NO',M_headers))+M(:,strcmp('R1R2_ON',M_headers));
    Fsum(isnan(Fsum))=2;
    Fsum(Fsum==2)=1;

    F{fi}.vect=~Fsum;
    F{fi}.name='R1R2ee'; F{fi}.d=2; fi=fi+1;

    F{fi}.vect=(M(:,strcmp('E_block',M_headers))==1);
    F{fi}.name='b1'; F{fi}.d=2; fi=fi+1;
    
    F{fi}.vect=(M(:,strcmp('E_block',M_headers))==2);
    F{fi}.name='b2'; F{fi}.d=2; fi=fi+1;
    
    for jj=1:length(F), 
        f_names{jj}=F{jj}.name; 
        F{jj}.vect(isnan(F{jj}.vect))=0;
    end
    %==========================%
    % Calculate values for each filter, based the counts upon the Mc
    % variable
    %==========================%
    for jj=1:size(model,2)
        for ii=1:size(model,1)
            Mtemp=M;
            MtempRT=M_RT;
            Mmask=zeros(size(M));
            t=0;
            for kk=1:length(model{ii,jj})
                if ~isempty(model{ii,jj}{kk})
                    I=strcmp(model{ii,jj}{kk},f_names);
                    Mmask(logical(F{I}.vect),:)=Mmask(logical(F{I}.vect),:)+1;   
                    t=t+1;
                end
            end
            Mmask(Mmask<t)=0; Mmask(Mmask==t)=1;
            save_rows{ii}=logical(mean(Mmask,2));
            Mtemp=Mtemp(save_rows{ii},:);
            MtempRT=MtempRT(save_rows{ii},:);

            % Create value var
            for kk=1:ML
%                 valueRT{ii,jj}(1,kk)=nansum(MtempRT(:,kk),1); value{ii,jj}(1,kk)=0;  
                if Mc{kk}==0
                    value{ii,jj}(1,kk)=nansum(Mtemp(:,kk),1);
                else      
                    for zz=1:length(Mc{kk})
                        if zz==1,
                            value{ii,jj}(1,kk)=nansum(Mtemp(:,Mc{kk}(zz)),1);
                        else
                            value{ii,jj}(1,kk)=nansum(Mtemp(:,Mc{kk}(zz)),1)+value{ii,jj}(1,kk);
                        end
                    end
                end 
            end
           
            value{ii,jj}(2,:)=nansum(Mtemp,1)./value{ii,jj}(1,:);
%             valueRT{ii,jj}(2,:)=nanmean(MtempRT,1);
            
            value{ii,jj}(3,:)=nanstd(Mtemp,1); % Really only makes sense for RT calculations
%             valueRT{ii,jj}(3,:)=nanstd(MtempRT,1);
            
            value{ii,jj}(4,:)=nansum(Mtemp,1);
%             valueRT{ii,jj}(4,:)=nansum(MtempRT,1);
        end
    end
    %==========================%
    % Save the data and save the resulting plots - may be worth integrating
    % these at some point.
    %==========================%
    % model [category X filter]
    for jj=1:size(model,2)
        
        model_name=[model{1,jj}{2:end}];
        sdisp(model_name,2);
        if isempty(model_name), model_name='base'; end
        switch behav.hc
            case 0, model_dir=fullfile(behav.save_dir,'pro',model_name); 
            case 1, model_dir=fullfile(behav.save_dir,'pro_hc',model_name); 
            case 2, model_dir=fullfile(behav.save_dir,'pro_hc_only',model_name); 
            case 3, model_dir=fullfile(behav.save_dir,'pro_hc_toss',model_name); 
                 case 4, model_dir=fullfile(behav.save_dir,'pro_lc_only',model_name); 
        end
        if ~exist(model_dir,'dir'), mkdir(model_dir); end
        
        if (~isempty(run_filters) && sum(strcmp(model_name,run_filters))==0), continue; end
        display(['     Running']);
            

        for ii=1:size(model,1)
            %=============================================================%
            % Save Data 
            %=============================================================%
            
            for oo=1:length(M_headers), template{ii}{oo}.header=M_headers{oo}; end
            template{ii}{ML+1}.header='ID';

            template{ii}{ML+1}.col{1}=behav.subjects{s};
            for zz=1:4
                for kk=1:length(template{ii})-1
                    template{ii}{kk}.col(1)=value{ii,jj}(zz,kk);
%                     template_RT{kk}.col(1)=valueRT{ii,jj}(zz,kk);
                end
                switch zz
                    case 1, sa=[model{ii,jj}{1} 'Ntype'];
                    case 2, sa=[model{ii,jj}{1} 'mean'];
                    case 3, sa=[model{ii,jj}{1} 'std'];
                    case 4, sa=[model{ii,jj}{1} 'Nevents'];
                end
                 
                % This only saves the group model
%                 if isempty([model{ii,jj}{1}]),
                    if s==1,
                        write_struct(template{ii},fullfile(model_dir,[sa '_summary.csv']));
%                         write_struct(template_RT,fullfile(model_dir,[sa '_RT_summary.csv']));
                    else
                        write_struct(template{ii},fullfile(model_dir,[sa '_summary.csv']),'a');
%                         write_struct(template_RT,fullfile(model_dir,[sa '_RT_summary.csv']),'a');
                    end
%                 end
            end
        end
        %=============================================================%
        % Make interesting plots
        %=============================================================%
        my_cmap=jet(13);
        for kk=1:length(plot_set)
           X=[]; Xe=[];
 
           switch plot_set{kk}.set
               case 'values', plot_values=value;
               case 'RT', plot_values=value;
           end
           
           switch plot_set{kk}.type
               case 'line' % Not used
                   for pp=1:size(model,1)
                       plot(plot_values{pp,jj}(2,plot_set{kk}.v),'color',my_cmap(pp,:),...
                           'linewidth',3); 
                       hold on;     
                   end
                   legend(category_u); 
               case 'bar'
                   for pp=1:size(model,1) 
                       if plot_set{kk}.dprime==1
                           % Should be [hit miss cr fa]
                           for gg=1:size(plot_set{kk}.v,1)
                               A=plot_values{pp,jj}(:,plot_set{kk}.v(gg,:));
                               HR=(A(1,1)*A(2,1)+A(1,3)*A(2,3))/(A(1,1)+A(1,3));
                               FA=(A(1,2)*A(2,2)+A(1,4)*A(2,4))/(A(1,2)+A(1,4));
                                if HR==1, 
                                    HR=2^-(1/(A(1,1)+A(1,3))); 
                                elseif HR==0
                                    HR=1-2^-(1/(A(1,1)+A(1,3)));
                                end

                                if FA==1, 
                                    FA=2^-(1/(A(1,2)+A(1,4))); 
                                elseif FA==0,
                                    FA=1-2^-(1/(A(1,2)+A(1,4)));
                                end
                                X(pp,(gg-1)*3+1)=norminv(HR)-norminv(FA);
                                HR=A(2,1); if HR==1, HR=2^-(1/A(1,1)); elseif HR==0, HR=2^-(1/A(1,1)); end
                                FA=A(2,4); if FA==1, FA=2^-(1/A(1,4)); elseif FA==0, FA=2^-(1/A(1,4)); end
                                X(pp,(gg-1)*3+2)=norminv(HR)-norminv(FA);
                                HR=A(2,3); if HR==1, HR=2^-(1/A(1,3)); elseif HR==0, HR=2^-(1/A(1,3)); end
                                FA=A(2,2); if FA==1, FA=2^-(1/A(1,2)); elseif FA==0, FA=2^-(1/A(1,2)); end
                                X(pp,(gg-1)*3+3)=norminv(HR)-norminv(FA);
                           end
                           Xe(pp,:)=X(pp,[2,5]); % Only Hit-FA calc
                           if pp==size(model,1)
                               dprime_dir{jj}=model_dir;
                               dprime_data{jj}{1}.header='Word';
                               dprime_data{jj}{2}.header='Picture';
                               dprime_data{jj}{1}.col(s)=Xe(end,1);
                               dprime_data{jj}{2}.col(s)=Xe(end,2);
                           end
                       else
                            Xe(pp,:)=plot_values{pp,jj}(plot_set{kk}.i,plot_set{kk}.v);
                       end
                   end

                   if ~exist(fullfile(model_dir,[plot_set{kk}.name '_' behav.subjects{s} '.png']),'file')
                      h=figure(kk);
                       bar(Xe'); 
                       title(plot_set{kk}.name);
                       legend(category_u);
                       ylabel('%');
                       set(gca,'XTickLabel',plot_set{kk}.legend);
                       set(gcf,'Position',[50,50,1200,800]);
                       set(gcf,'color','w');
                       export_fig(fullfile(model_dir,[plot_set{kk}.name '_' behav.subjects{s}]));
                       close(h);
                   end
           end
           clear Xe X HR FA
        end
    end
    clear M M_RT
end

clear save_data;
save_data{1}.header='ID';
for ii=1:length(behav.subjects)
    save_data{1}.col{ii}=behav.subjects{ii};
    for jj=1:length(File1_v)
        save_data{jj+1}.col(ii)=File1_v(ii,jj);
        if ii==1, save_data{jj+1}.header=File1_h{jj+1}; end
    end
end
write_struct(save_data,fullfile(behav.save_dir,'Confidence_Measures.csv'));
sdisp('Save dprime',2);
for ii=1:length(dprime_dir),
    if ~isempty(dprime_dir{ii})
        write_struct(dprime_data{ii},fullfile(dprime_dir{ii},'dprime.csv'));
    end
end





