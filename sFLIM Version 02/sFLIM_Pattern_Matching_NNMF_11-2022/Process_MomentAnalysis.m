
[filenames, pathname]  =  uigetfile( ...
    {'*.ptu','MAT-files (*.ptu)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Pick a file', ...
    'MultiSelect', 'on');
fprintf('\nTotal Number of file to be processed: %d \n',numel(filenames));
tic;

parfor k = 1:numel(filenames)
    
    filename  = [pathname cell2mat(filenames(k))];
    MomentAnalysis(filename);
    fprintf('\n Finished Processing file: %d',k);
end

toc