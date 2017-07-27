%%
% Set nodata and 12 bits per pixel for import to arcmap

baseDir = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\3321D_2010_319\';
dirList = dir([baseDir '*RGBIR.tif']);
gdalCmdLine = '"C:\Program Files\GDAL\gdal_translate.exe" -mo "BitsPerSample=12" -a_nodata 0' ;

for i = 1:length(dirList)
    infileName = [baseDir dirList(i).name];
    outfileName = [infileName '.new'];
    fprintf('Processing %d of %d, %s\n', i, length(dirList), dirList(i).name);
    cmdLine = [gdalCmdLine ' "' infileName '" "' outfileName '"'];
    [status, result] = dos(cmdLine, '-echo');
    if (status == 0)
        r = recycle;
        recycle('off');
        delete(infileName);
%         delete([outfileName '.aux.xml']);
        delete([infileName '.aux.xml']);
        delete([infileName '.ovr']);
%         delete([infileName(1:end-4) '.pox']);
        recycle(r);
        status = movefile(outfileName, infileName, 'f');
        if (status == 0)
            fprintf('ERROR: could not rename %s to %s, %s\n', outfileName, infileName, lasterr);
        end
    end
end