function vision_shell(stim)
% stim.model
%   GaborJet
%   RadonLum
%   hmax_c3
%   hmax_c2
%=========================================================================%
% HMAX setup (retained as reference)
%=========================================================================%
% Need functional mex compilers in order to do this.
% Initial Setup - this must be run a windows 64bit machine, akin to the
% serv1, in order for the code to work on serv1. Compilers unfortunately
% don't work on Windows Server 2012 w/ Matlab
% ERROR: Files\Microsoft was unexpected at this time -> is expected, just
% says that gpu processing failed
install=0;
if install==1
    cns_install; cns_build demopkg; demopkg_run; cns_build hmax;
end

switch stim.model
    case 'GaborJet'
        GaborJet(stim);  % Complete
    case 'RadonLum'
        radon_lum(stim); % Complete
    case 'hmaxc3'
        hmax_c3(stim);
    case 'hmaxc2'
        hmax_c2(stim); % Complete
end