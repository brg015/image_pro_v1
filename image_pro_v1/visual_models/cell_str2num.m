% BRG 2014 Winter
%
% Description: Convert string cell arrays to numbers, anything that can't
% be changed is left as NaN
%
% Notes: Add more robust checks at some date
function num_array=cell_str2num(cell_array)

L=length(cell_array);
num_array=nan(L,1);
for ii=1:L
    try num_array(ii)=str2num(cell_array{ii}); end
end