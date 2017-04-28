function hmax_c3(stim)

switch stim.study
    case 'ER'
        % Generates random sampling of 16 exemplar images
        sets=1;
        for jj=1:sets
            train_set(jj,:)=jj:3:288;
        end
end

for kk=1:sets
    
clear p g m 
lib = struct;
% 1) Choose model parameters.
% p = hmax_cvpr06_params_base;
% p = hmax_cvpr06_params_full;k
% p = hmax_cvpr06_params_full_ac;
[p, g] = hmax_pnas07_params;
%-----------------------------------------------------------------------------------------------------------------------

gs = [g.s2b, g.s2, g.s3]; % Each of these stages needs a dictionary.
% ns = [2048 , 2048, 1024]; % Size of each dictionary.
ns = [84 84 42];
display(' Training...');

for i = 1 : numel(gs)
    m = hmax.Model(p, lib);
    d = hmax_s.EmptyDict(m, gs(i), ns(i));
    cns('init', m, 'cpu'); tic
    for image=train_set(kk,:);
        fn = stim.epl{image};
        im = fancy_read(stim,fn);
        hmax.LoadImage(m, im);
        cns('run');
        d = hmax_s.SampleFeatures(m, gs(i), d, ns(i));
        fprintf('sampled %u "%s" features from "%s"\n', ns(i), m.groups{gs(i)}.name, fn);
    end
    cns('done');
    d = hmax_s.SortFeatures(d);
    if cns_istype(m, -gs(i), 'ss')
        d = hmax_ss.SparsifyDict(d);
    end
    lib.groups{gs(i)} = d;
    toc
end
clear i d j fn im;
m = hmax.Model(p, lib);

cns('init', m, 'cpu');
display(' Building...');
for image=1:length(stim.exp)
    tic
    fn = stim.exp{image};
    im = fancy_read(stim,fn);
    
    hmax.LoadImage(m, im);
    cns('run');
    c2b = cns('get', -m.c2b, 'val'); c2b = cat(1, c2b{:});
    c3  = cns('get', -m.c3 , 'val'); c3  = cat(1, c3 {:});
    c2b(c2b == cns_fltmin) = 0;
    c3 (c3  == cns_fltmin) = 0;
    cmat{1}(image,:)=c2b;
    cmat{2}(image,:)=c3;
    
    clear j fn im time;
    toc
end
cns('done');

mkdir_tree(fullfile(stim.save,'hmaxc3'));
mkdir_tree(fullfile(stim.save,'hmaxc3_nan'));
mkdir_tree(fullfile(stim.save,'hmaxc2b'));
mkdir_tree(fullfile(stim.save,'hmaxc2b_nan'));
stim_ID_num=stim.stim_ID_num;

R=corr(cmat{2}');
x=eye(size(R));
save(fullfile(stim.save,'hmaxc3','model.mat'),'R','stim_ID_num');
R(x==1)=nan;
save(fullfile(stim.save,'hmaxc3_nan','model.mat'),'R','stim_ID_num');
R=corr(cmat{1}');
save(fullfile(stim.save,'hmaxc2b','model.mat'),'R','stim_ID_num');
R(x==1)=nan;
save(fullfile(stim.save,'hmaxc2b_nan','model.mat'),'R','stim_ID_num');





end