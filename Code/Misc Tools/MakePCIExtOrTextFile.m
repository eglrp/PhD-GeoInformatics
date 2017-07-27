function MakePCIExtOrTextFile(inFile, idx, outFile, imDir)
fid = fopen(inFile, 'r');
for i = 1:5
    fgets(fid);
end
% data = textscan(fid, '%d %n %n %n %n %n %n', -1);
data = fscanf(fid, '%g', [7 inf])';
fclose(fid);

%3321C_2010_318_01_0001_RGB
if isempty(idx)
    idx = 1:size(data,1);
end
fid = fopen(outFile, 'w');
for i = 1:length(idx) %size(data,1)
%     imageName = sprintf('%s%04.0f_RGB', imageNamePrefix, data(i, 1));
    imageNum = sprintf('%04.0f', data(idx(i), 1));
    dirList = dir([imDir '\\*_*_*_*_' imageNum '_*.tif']);
    if (length(dirList) > 1)
        error('File index %d not unique in %s', i, imDir);
    end
    fprintf(fid, '%s %f %f %f %f %f %f\n', dirList(1).name(1:end-4), data(idx(i), 2:end));
end
fclose(fid);
end