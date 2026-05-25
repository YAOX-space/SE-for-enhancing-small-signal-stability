function data = import_reproduction_data(dataDir)
%IMPORT_REPRODUCTION_DATA Load all machine-readable paper tables.

data = struct();
files = dir(fullfile(dataDir, '*.csv'));
for k = 1:numel(files)
    [~, name] = fileparts(files(k).name);
    path = fullfile(files(k).folder, files(k).name);
    opts = detectImportOptions(path, 'Delimiter', ',');
    data.(matlab.lang.makeValidName(name)) = readtable(path, opts);
end
end
