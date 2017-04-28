% Master_Image (Summer 2014)
% BR Geib
%
% Clear all and setup code directries
clc; clear all; close all;
serv='serv1';
switch serv
    case 'local', root_dir='C:\';
    case 'serv1', root_dir='D:\Data\';
end
code_dir{1}=fullfile(root_dir,'Geib\Scripts\Matlab_Scripts\visual_models');
code_dir{2}=fullfile(root_dir,'Geib\Scripts\Matlab_Scripts\function_files\');

for ii=1:length(code_dir), addpath(genpath(code_dir{ii})); end
%=========================================================================%
% Input Parameters
%=========================================================================%
% Notes:
% 1) Images should be RGB uint8 - other formats may work, but this is what
%    I've been working with atm.
% 2) The flow of information follows the order of the variables listed.
%    That is images are 1) cropped, 2) squared by adding white padding to
%    the shorter dimension, 3) Resized to specified dimensions, 4) thrown
%    onto a white background.
% 3) If you don't want to write to a directory, e.g. wrte==0, then you
%    don't need to define a write directory.
% 4) The screen will alert as to the processing of images and if specified
%    images could not be found

% Image array setup (example): all I'm doing here is looking for all the
% jpgs in a folder and assigning the image names to a cell array.
X=dir('F:\SAM\PresentationScripts\Exemplar1_Mar2016\*.jpg');
image_array=cell(1,length(X));
for ii=1:length(X),
    image_array{ii}=X(ii).name;
end

params.raw_dir='F:\SAM\PresentationScripts\Exemplar1_Mar2016\';   % Where your images are located

params.crop.dir='';
params.crop.wrte=0;

params.sqre.dir='';
params.sqre.pad=100;     % Optional padding around image
params.sqre.wrte=0;

params.resz.dir='F:\SAM\PresentationScripts\Exemplar1_Mar2016_altered\';
params.resz.size=[300 300];
params.resz.wrte=1;

image_modifier(image_array,params)