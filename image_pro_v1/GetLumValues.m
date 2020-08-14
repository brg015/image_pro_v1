% BRG 1/27/20
%% Set image path
image_path='\\ccn-cabezaserv1.win.duke.edu\Data\SAM\Exemplar1_Mar2016_square\';
%% Read in relevant files
%
wrk_dir='C:\Users\brg13\Google Drive\Papers (manuscripts)\SEA\Cerebral Cortex Submission\';
addpath(genpath('\\ccn-cabezaserv1.win.duke.edu\Data\Geib\Scripts\Public\function_files\'));
data1=excel_reader(fullfile(wrk_dir,'PresentationStimList.csv'));
data2=excel_reader(fullfile(wrk_dir,'SAM_enc_pro_3500.csv'));
%% Compute luminance
% see: 
% https://stackoverflow.com/questions/596216/formula-to-determine-brightness-of-rgb-color
for iStim=1:954
    X=double(imread(fullfile(image_path,data1{3}.col{iStim})));
    % Step 1: Normalize
    R=reshape(X(:,:,1)/255,1,[]);
    G=reshape(X(:,:,2)/255,1,[]);
    B=reshape(X(:,:,3)/255,1,[]);
    % Step 2: Convert
    Rnew=NaN(1,length(R));
    Rnew(R<=0.04045)=R(R<=0.04054)/12.92;
    Rnew(R>0.04045)=((R(R>0.04045)+0.055)/1.055).^2.4;
    
    Gnew=NaN(1,length(G));
    Gnew(G<=0.04045)=G(G<=0.04054)/12.92;
    Gnew(G>0.04045)=((G(G>0.04045)+0.055)/1.055).^2.4;
    
    Bnew=NaN(1,length(B));
    Bnew(B<=0.04045)=B(B<=0.04054)/12.92;
    Bnew(B>0.04045)=((B(B>0.04045)+0.055)/1.055).^2.4;
    
    Y1=Rnew*0.299 + Gnew*0.7152 + Bnew*0.0722;
    % Perception
    Y2=NaN(1,length(Y1));
    Y2(Y1<=0.008856)=Y1(Y1<=0.008856)*903.3;
    Y2(Y1>0.008856)=(Y1(Y1>0.008856).^(1/3))*116-16;
    
    ImInfo(iStim).name=data1{3}.col{iStim};
    ImInfo(iStim).lum=mean2(Y1);
    ImInfo(iStim).perlum=mean2(Y2);
    clear Y1 X R G B;
    clear Rnew Gnew Bnew Y2;
end
%% Compare to subject data
data3{1}.header='StimulusName';
data3{2}.header='LuminanceValue';
data3{3}.header='PerLuminanceValue';
for iTrial=1:640
    Ind=cell2num(data2{4}.col(iTrial));
    if Ind ~= 0
        data3{1}.col{iTrial}=ImInfo(Ind).name;
        data3{2}.col{iTrial}=ImInfo(Ind).lum;
        data3{3}.col{iTrial}=ImInfo(Ind).perlum;
    else
        data3{1}.col{iTrial}='Catch';
        data3{2}.col{iTrial}='Nan';
        data3{3}.col{iTrial}='Nan';
    end
end

write_struct(data3,fullfile(wrk_dir,'LumValues.csv'));