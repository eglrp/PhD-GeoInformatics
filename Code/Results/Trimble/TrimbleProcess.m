%%
% Rename and copy files
% NB: Matjiesvlei.* raw (now renamed on device) is actually Matjiesvlei5.*
% Groenefontein raw was renamed as Rooiberg

% There are Marjiesvlei and Marjiesvlei1 on the HD - these map the same
% area(s).  Actually there are 2 polygons in each of these files (the same
% polygons that is)

baseDir = 'D:\Data\Development\Projects\MSc GeoInformatics\Docs\Trimble\';

dirList = dir([baseDir '*']);
dirList = dirList(3:end);
dirList = dirList([dirList.isdir]);

gdalCmdLine = '"C:\Program Files\GDAL\gdal_translate.exe" -mo "BitsPerSample=12" -a_nodata 0' ;

for i = 1:length(dirList)
    dirName = [baseDir dirList(i).name '\'];
    subDirList = dir([dirName '*_gen.*']);
    for j = 1:length(subDirList)
        sourceFileName = [dirName subDirList(j).name];
        fn = dirList(i).name;
%         fn = [fn(1) (fn(2:end))];
        destFileName = [dirName fn subDirList(j).name(end-3:end)];
        fprintf('Moving %s to %s\n', sourceFileName, destFileName);
        success = movefile(sourceFileName, destFileName);
        if (success)
            destFileName2 = [dirName];
            fprintf('Copying %s to %s\n', destFileName, baseDir);
            copyfile(destFileName, baseDir);
        else
            error(lasterr);
        end
    end
end
%%
clear all
%Rename corrected files
% srcNames = {'GROENEFONTEIN', 'GROENFONTEIN', 'GROOTKOP', 'MATJIESVLEI1'};
% destNames = {'RooiBerg', 'GroenFontein', 'GrootKop', 'MatjiesVlei'};
baseDir = 'D:\Data\Development\Projects\MSc GeoInformatics\Docs\Trimble\Corrected\';
dirList = dir([baseDir '*']);
dirList = dirList(3:end);
dirList = dirList([dirList.isdir]);
% 
% for i = 1:length(dirList)
%     matchIndex = strcmp(dirList(i).name(1:end-1), srcNames);
%     matchIndex = find(matchIndex);
%     
%     if ~isempty(matchIndex)
%         matchIndex = matchIndex(1);
%         srcDirName = [baseDir dirList(i).name];
%         newDirName = [baseDir, destNames{matchIndex}, dirList(i).name(end)];
%         fprintf('Renaming %s\n to %s\n', srcDirName, newDirName);
%         ok = input('OK? ', 's');
%         if strcmpi(ok, 'y')
%             movefile(srcDirName, newDirName)
%         end
%     end
% end
%


for i = 1:length(dirList)
    dirName = [baseDir dirList(i).name '\'];
    subDirList = dir([dirName '*_gen.*']);
    for j = 1:length(subDirList)
        sourceFileName = [dirName subDirList(j).name];
        fn = dirList(i).name;
%         fn = [fn(1) lower(fn(2:end))];
        destFileName = [dirName fn subDirList(j).name(end-3:end)];
        fprintf('Moving %s to %s\n', sourceFileName, destFileName);
        success = movefile(sourceFileName, destFileName);
        if (success)
%             destFileName2 = [dirName];
            fprintf('Copying %s to %s\n', destFileName, baseDir);
            copyfile(destFileName, baseDir);
        else
            error(lasterr);
        end
    end
end


%%

%%
%Combine corrected shape files into 1
clear all;
baseDir = 'D:\Data\Development\Projects\MSc GeoInformatics\Docs\Trimble\Corrected\';
dirList = dir([baseDir '*.shp']);
metaFile = 'D:\Data\Development\Projects\MSc GeoInformatics\Docs\My Docs\Field Trip 2012\2012FieldTripNotes.xls';
[n metaText] = xlsread(metaFile);
metaFields = metaText(1,:);
metaText = metaText(2:end,:);

res = [];
for i = 1:length(dirList)
    name = dirList(i).name(1:end-4);
    shape = shaperead([baseDir dirList(i).name]);
    xlen = arrayfun(@(x)(length(x.X)), shape);
    shape = shape(xlen>50); %remove mistake polygons
    for j = 1:length(shape)
        if (length(shape) > 1)
            shape(j).Name = [name char(97 + j - 1)];
        else
            shape(j).Name = name;
        end
        shape(j).Geometry = 'Polygon';
        shape(j).FileName = [baseDir dirList(i).name];
        
        metaIdx = strcmpi(shape(j).Name, metaText(:,1));
        for k = 2:length(metaFields)
             val = str2double(metaText{metaIdx,k}); %convert numeric if possible
             if isnan(val)
                val = metaText{metaIdx,k};
             end
            shape(j).(metaFields{k}) = val;
        end
    end
    res = [res;shape];
end
shapewrite(res, [baseDir '..\GroundTruthCombinedCorrected2.shp']);
shapewrite(res(11), [baseDir '..\MatjiesVlei2Corrected.shp']);

%%
D:\Data\Development\Projects\MSc GeoInformatics\Docs\Trimble