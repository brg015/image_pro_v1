% BRG2014(WINTER)
function create_randomization_files(rando,save_dir)

image_mapping=fullfile(rando.stim_dir,'image_mapping.csv');
if ~exist(image_mapping,'file'),
    display([image_mapping ': Doesn''t exist - quitting']); return;
end

% Reading mapping file
din=excel_reader(image_mapping,{'A' 'B' 'C'},{'FileNames' 'c1' 'c2'});
L=length(din{1}.col); % Because we have exemplars

% Define the fourth column to have living information...
LYN=rand(L,1); LYN(LYN<0.5)=0; LYN(LYN>=0.5)=1;
for ii=1:L, din{4}.col{ii}=LYN(ii); end

for n=1:length(rando.nfiles)    
    % Thus create an offset matrix of random 0s and 1s
    % R_offset => determines which exemplar to use 0/1
    R_offset=rand(1,L/2);
    R_offset(R_offset<0.5)=0;
    R_offset(R_offset>=0.5)=1;

    % Create lagged randomized list for unique items
    R1=randperm(L/2);
    R2=random_lag(R1,rando.lag);
    
    % Randomization procedure highly similar to the previous ones, except
    % R3 has even more constraints.
    lag_fit=0; c=0;
    while lag_fit==0;
        c=c+1;
        R3=random_lag(R2,int32(rando.lag*1.5)); % Artifical lag increase to help fit
        
        % Sort R3 so that first items match encoded items - this is really
        % only change that happens here...
        for ii=1:rando.enc_L, 
            [~,AFC_order(ii)]=find(R3==R1(ii)); 
        end
        [~,I]=sort(AFC_order);
        R3_reorder=[AFC_order(I) setdiff(1:L/2,AFC_order(I))];
        R3=R3(R3_reorder);
    
        L2=length(R2);
        for ii=1:length(R2),
            idx=find(R3(ii)==R2);        % Find matching element in R1
            distance_diff(ii)=ii+L2-idx; % Find distance to element
            % If distance is less then required lag, try again
            if distance_diff(ii)<=rando.lag, break; end
        end
        if min(distance_diff)>=rando.lag, lag_fit=1; end % Distance always exceeds lag, R2 is found
    end
    fprintf('R3 Randomization Facts\n')
    fprintf(['\tIterations: ' num2str(c) '\n']);
    fprintf(['\tMin Dstnce: ' num2str(min(distance_diff)) '\n']);
    fprintf(['\tAvg Dstnce: ' num2str(mean(distance_diff)) '\n']);
    fprintf(['\tMax Dstnce: ' num2str(max(distance_diff)) '\n']);
   
    % Before corrupting indices, let's check for old/new assesments
    for ii=1:length(R2)
        % If location of retrieval item ii is greater than the number of
        % items shown at encoding then it is marked as new (==1)
        % 
        % Also determine the offset used in R1 and apply to R2 as well
        Match_Idx=find(R1==R2(ii));
        OldNew(ii)=Match_Idx>rando.enc_L;
        R2_offset(ii)=R_offset(Match_Idx);
    end
    
    % Now add the potential offset. We multiply the list of potential items
    % by two to match the length of the actual list. We then subtract by
    % one such that the first exemplar is always referenced. We then add an
    % array of 0/1s which speficies that exemplars should be randomly
    % presented.
    %
    % In general, these represent presentation order
    R1=(R1*2-1)+R_offset;
    R2=(R2*2-1)+R2_offset;  % Forces R2 codes to equal R1 codes 
    R4=(R3*2-1)+~R_offset;             % Other exemplar
    R3=(R3*2-1)+R_offset;              % Presented item
    
    % Specify headers of output file
    data{1}.header='file_idx';    % Encoding Set
    data{2}.header='phase';       data{2}.col(1:L/2)=1;
    data{3}.header='code1';
    data{4}.header='code2';
    data{5}.header='NonLiving';
    data{6}.header='file_idx';    % Word Retrieval Set
    data{7}.header='phase';       data{7}.col(1:L/2)=2;
    data{8}.header='code1';
    data{9}.header='code2';
    data{10}.header='New';        data{10}.col(1:L/2)=OldNew; % 0==Old, 1==New
    data{11}.header='file_idx1';  % 2AFC Left Side Present
    data{12}.header='phase';      data{12}.col(1:L/2)=3;
    data{13}.header='code1';
    data{14}.header='code2';
    data{15}.header='file_idx2';  % 2AFC Right Side Present
    data{16}.header='phase';      data{16}.col(1:L/2)=3;
    data{17}.header='code1';
    data{18}.header='code2';
    data{19}.header='Right';      % 0==L, 1==R is original
    
    for ii=1:(L/2)
        data{1}.col{ii}=R1(ii);
        data{3}.col(ii)=str2num(din{2}.col{R1(ii)});
        data{4}.col(ii)=str2num(din{3}.col{R1(ii)});
        data{5}.col(ii)=din{4}.col{R1(ii)};
        
        data{6}.col{ii}=R2(ii);
        data{8}.col{ii}=str2num(din{2}.col{R2(ii)});
        data{9}.col{ii}=str2num(din{3}.col{R2(ii)});
        
        data{11}.col{ii}=R3(ii);
        data{13}.col(ii)=str2num(din{2}.col{R3(ii)});
        data{14}.col(ii)=str2num(din{3}.col{R3(ii)});
        
        data{15}.col{ii}=R4(ii);
        data{17}.col(ii)=str2num(din{2}.col{R4(ii)});
        data{18}.col(ii)=str2num(din{3}.col{R4(ii)});      
    end
    
    % Compare original encoding data to left presented data, if a match is
    % found then we know correct item is on the left.
    Original_Present=[data{3}.col(:),data{4}.col(:)];
    L_Present=[data{13}.col(:),data{14}.col(:)];
    for ii=1:(L/2)
        data{19}.col(ii)=~ismember(L_Present(ii,:),Original_Present,'rows');
    end
  
    write_struct_txt(data,fullfile(save_dir,[rando.save_name num2str(n) '.txt']));
    write_struct(data,fullfile(save_dir,[rando.save_name num2str(n) '.csv']));
end
    
end
    

