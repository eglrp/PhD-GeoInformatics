function res = ParseZippsLog(logFileName)

fid = fopen(logFileName);

tline = fgets(fid);
redcordCount = 1;
while ischar(tline)
    if strncmp(tline, '------------------------ Image', 30)
        redcordCount = redcordCount + 1;
        res(redcordCount).Name = sscanf(tline, '------------------------ Image %s ');
        tmp = textscan(res(redcordCount).Name, '%*s%d', 'delimiter', '~');
        res(redcordCount).Id = tmp{1};
    elseif strncmp(tline, 'Exposure Settings:', 18)
        tline = fgets(fid);
        tline = fgets(fid); %fstop
        tline = fgets(fid); %exposure
        res(redcordCount).Exp = sscanf(tline, 'Exposure [ms] %f %f %f %f %f');
    elseif strncmp(tline, 'MS Platform calib.:', 19)
        tmp = deblank(sscanf(tline, 'MS Platform calib.: %s'));
        if (strcmpi(tmp, 'successful'))
            res(redcordCount).MsPlatformCalibOk = true;
        else
            res(redcordCount).MsPlatformCalibOk = false;
        end
    elseif strncmp(tline, 'Radiometric platform calibration:', 33)
        tmp = deblank(sscanf(tline, 'Radiometric platform calibration: %s'));
        if (strcmpi(tmp, 'successful'))
            res(redcordCount).RadiomPlatformCalibOk = true;
        else
            res(redcordCount).RadiomPlatformCalibOk = false;
        end
    end
    tline = fgets(fid);
end
fclose(fid);

end

