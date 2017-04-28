function randomize_stimuli(design,rando)
keyboard;
data=excel_reader(design.file,design.col,design.head);
L=length(data);
data{L+1}.header='E_order';
data{L+2}.header='E_pos';
data{L+3}.header='R1_order';
data{L+4}.header='R1_pos';
data{L+5}.header='R2_order';
data{L+6}.header='R2_pos';

Eorder1=randperm(460); 
Rorder1=random_lag(Eorder1,rando.lag,0);
Rorder2=random_lag(Rorder1,rando.lag,0);

Old_idx=strcmp(design.head,'Old');
catch_idx=strcmp(design.head,'Catch');
for ii=1:460, c(ii)=str2num(data{catch_idx}.col{ii}); end
for ii=1:460, o(ii)=str2num(data{Old_idx}.col{ii}); end
Eorder1(~logical(o))=461;
Rorder1(logical(c))=461;
Rorder2(logical(c))=461;
Rorder2(~logical(o))=461;

for jj=1:460    
    data{L+3}.col(jj)=Rorder1(jj);
    data{L+1}.col(jj)=Eorder1(jj);
    data{L+5}.col(jj)=Rorder2(jj);
end

[~,Epos1]=sort(data{L+1}.col);
[~,Rpos1]=sort(data{L+3}.col);
[~,Rpos2]=sort(data{L+5}.col);
for jj=1:460
    data{L+2}.col(jj)=Epos1(jj);        
    data{L+4}.col(jj)=Rpos1(jj);
    data{L+6}.col(jj)=Rpos2(jj);
end
write_struct_txt(data,[rando.save_file '.txt']);
write_struct(data,[rando.save_file '.csv']);