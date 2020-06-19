clr
cd('E:\SiteResp_GMSelect\Python_SourceCode')

tmp = dir('*.AT2');

% k(1).name

len = length(tmp);

for ii = 1:len

% K = importdata(tmp(ii).name);
% name = K{index};

fid = fopen(tmp(ii).name);
i = 1;
tline = fgetl(fid);
A{i} = tline;

while ischar(tline)
    i = i+1;
    tline = fgetl(fid);
    A{i} = tline;
end
fclose(fid);


req = A{4};
% g = textscan(fid,'%s','delimiter','\n');
% lenreq = length(g{:});

temp = req(6:13);
temp1 = req(21:26);
if strcmp(temp(length(temp)),',') == 1
    numGMPTS = str2double(temp(1:length(temp)-1));
    temp1 = req(21:26);
    dT = str2double(temp1);
else
    numGMPTS = str2double(temp(1:length(temp)));
    dT = str2double(temp1);
end

mod = strcat(num2str(numGMPTS),{'  '},num2str(dT),{'  '},'NPTS, DT');

% A{4} = num2cell(mod);

fid = fopen(strcat('E:\SiteResp_GMSelect\Python_SourceCode\',tmp(ii).name), 'w');

for jj = 1:numel(A)
    if jj == 4
        fprintf(fid,'%s\n', char(mod));
    else
    fprintf(fid,'%s\n', A{jj});
    end
end

fclose all;
end