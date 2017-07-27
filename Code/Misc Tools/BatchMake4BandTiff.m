function BatchMake4BandTiff(inBaseDirName, outDirName)
% gdal_translate -ot UInt16 -of GTiff -co "TILED=YES" -co "COMPRESS=DEFLATE" 3321A_2010_316_01_0003_RGB.tif 3321A_2010_316_01_0003_RGB_16b.tif
%     path1 = getenv('PATH')
%     path1 = [path1 ':/usr/local/bin']
%     setenv('PATH', path1)
%     !echo $PATH 
    doSkipExisting = true;
    gdalCmdLine = '"C:\Program Files\GDAL\gdal_translate.exe" -ot UInt16 -of GTiff -co "TILED=YES" ' ;
    rgbDirList = dir([inBaseDirName '\RGB\*RGB.tif']);
    for i = 1:length(rgbDirList)
        fprintf('>Processing image %d of %d\n', i, length(rgbDirList));
        rgbFileName = [inBaseDirName '\RGB\' rgbDirList(i).name];
        cirFileName = [inBaseDirName '\CIR\' rgbDirList(i).name(1:end-7) 'CIR.tif'];
%         rgbOutFileName = [outDirName '\' rgbDirList(i).name(1:end-4) '_16.tif'];
%         cirOutFileName = [outDirName '\' rgbDirList(i).name(1:end-7) 'CIR_16.tif'];
        rgbOutFileName = [outDirName '\RGB_16.tif'];
        cirOutFileName = [outDirName '\CIR_16.tif'];
        rgbIrOutFileName = [outDirName '\' rgbDirList(i).name(1:end-4) 'IR.tif'];
        
        if (doSkipExisting && exist(rgbIrOutFileName, 'file'))
            continue;
        end
        
        if (~exist(cirFileName, 'file'))
            fprintf('ERROR: %s does not exist\n', cirFileName);
            continue;
        end        

        %run these 2 commands in parallel
        fprintf('Writing RGB %s ...\n', rgbOutFileName);
        rgbCmdLine = [gdalCmdLine ' "' rgbFileName '" "' rgbOutFileName]; % '" &']; %exit closes the window
        [status, result] = dos(rgbCmdLine, '-echo');
        fprintf('Writing CIR %s ...\n', cirOutFileName);
        cirCmdLine = [gdalCmdLine ' "' cirFileName '" "' cirOutFileName]; % '" &'];
        [status, result] = dos(cirCmdLine, '-echo');
        
        fprintf('Reading RGB %s ...\n', rgbOutFileName);
        rgbTiff = Tiff(rgbOutFileName, 'r+');
        rgbIm = rgbTiff.read();

        fprintf('Reading CIR %s ...\n', cirOutFileName);
        cirTiff = Tiff(cirOutFileName, 'r+');
        cirIm = cirTiff.read();
        
        fprintf('Writing RGBIR %s ...\n', rgbIrOutFileName);
        rgbIrTiff = Tiff(rgbIrOutFileName, 'w');
        try
            rgbIrIm = cat(3, cirIm(:,:,[2 3]), rgbIm(:,:,3), cirIm(:,:,1));

            %copy tags to output file
            tagId = Tiff.TagID;
            fn = fieldnames(tagId);
    %         tagStruct = tagId;
            for j = 1:length(fn)
                try
                    val = rgbTiff.getTag(fn{j});
                    rgbIrTiff.setTag(fn{j}, val);
                catch ME
                    x = 1;
                end
            end
            rgbIrTiff.setTag('BitsPerSample', 16);
            rgbIrTiff.setTag('SamplesPerPixel', 4);
            % rgbIrTiff.setTag('MaxSampleValue', 4095);
            % tout.setTag('MaxSampleValue', 4095);
            % tout.setTag('NoDataValue', 4095);
            rgbIrTiff.setTag('Compression', Tiff.Compression.None); %DEFLATE IS THE PNG ALGORITHM
            
            rgbIrTiff.write(rgbIrIm);
            
        catch ME
            rgbTiff.close();
            cirTiff.close();
            rgbIrTiff.close();
            
            state = recycle;
            recycle('off');
%             delete(rgbOutFileName);
%             delete(cirOutFileName);
            recycle(state);
            
            rethrow(ME);
        end
        rgbTiff.close();
        cirTiff.close();
        rgbIrTiff.close();

%         state = recycle;
%         recycle('off');
%         delete(rgbOutFileName);
%         delete(cirOutFileName);
%         recycle(state);
    end

end

