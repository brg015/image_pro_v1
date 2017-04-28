% B.R. Geib Winter 2012 - br.geib33@gmail.com
%
% Reads excel files that have been stripped of commas. Discards the top row
% as suspected headers.
%
% Inputs
%   file => .csv file
%       e.g. file.csv
%   varargin{1} => an array of excel columns to use 
%       e.g. {'A' 'CX'} etc.
%   varargin{2} => corresponding names for the loaded columns
%       e.g. {'column A' 'column B'} etc.
% Outputs
%   data => a structure with header and column components. 
%
% Additional Notes: This script works well with the write_struct.m script.
% The write_struct.m script writes out .csv files based upon the
% information in the data structure.
%
function [data]=excel_reader(file,varargin)

data={};
if exist(file,'file')
	fid=fopen(file); show=0;
else
	data={};
	win_print(['ERROR: ' file ' DNE\n']);
	fprintf('  Continue [Y/N]: ');
	reply=input('','s'); 
	if strcmpi(reply,'y'),
		return;
	else
		error('Error: Correct pathways\n');
	end
end

% Read in the data
count=1;
while 1
    line=fgetl(fid);
    if ~ischar(line), break, fclose(fid); end
    split{count}=strread(line,'%s','delimiter',',');
    count=count+1;
end

file_data={}; header=split{1}; c=1;

% Allow for dication of columns
if length(varargin)==2
    file_data=let2num(varargin{1});
else
    file_data=1:size(header,1);
end

% Organize and save as a structure
for n = file_data
    % if headers are defined
    if length(varargin)==2
        data{c}.header=varargin{2}{c}; 
        % First is header, so start with 2nd structure
        for i=2:length(split)
            if ~isempty(split{i})
                try
                    data{c}.col(i-1) = split{i}(n);
                catch err
                    data{c}.col(i-1) = {''};
                end
            end
        end
        c=c+1;
    % if headers aren't defined
    else
        data{n}.header = header(n);
       % First is header, so start with 2nd structure
        for i=2:length(split)
            try
                data{n}.col(i-1) = split{i}(n);
            catch err
                 data{c}.col(i-1) = {''};
            end
                
        end
        c=c+1;
    end
 
end

% Print out the first 5 examples to make sure looks okay
if show==1
    for i=1:length(data), fprintf([char(data{i}.header) '\t\t' ]); end
    fprintf('\n');
    for i=1:5
        for j=1:length(data), fprintf([char(data{j}.col(i)) '\t\t' ]); end
        fprintf('\n');
    end
end

fclose(fid);

end

%=========================================================================%
% Internal Function: let2num
%=========================================================================%
%
% Function converts letter columns to numbers, e.g. AA = 27
function num_array=let2num(array)

for i=1:length(array)
    sum=0;
    N=length(array{i});
    for j=1:length(array{i})
        switch array{i}(j)
            case 'A', L=1;
            case 'B', L=2;
            case 'C', L=3;
            case 'D', L=4;
            case 'E', L=5;
            case 'F', L=6;
            case 'G', L=7;
            case 'H', L=8;
            case 'I', L=9;
            case 'J', L=10;
            case 'K', L=11;
            case 'L', L=12;
            case 'M', L=13;
            case 'N', L=14;
            case 'O', L=15;
            case 'P', L=16;
            case 'Q', L=17;
            case 'R', L=18;
            case 'S', L=19;
            case 'T', L=20;
            case 'U', L=21;
            case 'V', L=22;
            case 'W', L=23;
            case 'X', L=24;
            case 'Y', L=25;
            case 'Z', L=26;
        end
        % Add the base 26
        if j==N, sum=sum+L; else sum=sum+L*26*(N-j); end
    end
    num_array(i)=sum;
    clear sum N;
end

end
    
