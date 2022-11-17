function [head, tag, tcspc, IRF, timname, AP, DC] = load_data(name)

head    = [];
tcspc   = [];
IRF     = [];
tag     = [];
AP      = [];
DC      = [];
timname = '';

tmp = strfind(name,'\');
pathname = name(1:tmp(end));
filename = name(tmp(end)+1:end);
name1 = [filename(1:end-4) '_DATA.mat'];

names1 = dir([pathname name1]);

if (numel(names1)==1)
    load([pathname name1], 'head', 'tag', 'tcspc', 'IRF', 'AP', 'DC');
    tag(tag<0)=0;
    timname = [filename(1:end-4) '_FLIM.mat'];
    
else
    fprintf('\n Data-files were not found in current directory.');
    fprintf('\n Please call "Process_file" to create necessary files.');
end
     