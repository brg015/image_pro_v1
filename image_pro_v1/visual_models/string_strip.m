function num_r=string_strip(string)
    string(strfind(string,'i'))='o'; % Avoid imaginary numbers
    num=''; num_r=[];
    for ii=1:length(string)
        try
            num=[num,num2str(str2num(string(ii)))];
        end
    end
    try num_r=str2num(num); end
end