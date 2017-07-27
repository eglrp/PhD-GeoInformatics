function res = ParseZippsTxt(logFileName)
%parses text file produced by "gpplog2txt *.log"

fid = fopen(logFileName);

tline = fgets(fid);
fields = textscan(tline, '%s');
fields = fields{1};
recordCount = 1;
tline = fgets(fid);
while ischar(tline) & ~feof(fid)
    fieldData = textscan(tline, '%s');
    fieldData = fieldData{1};
    %hack for date time with space
    res(recordCount).(fields{1}) = fieldData{1};
    res(recordCount).(fields{2}) = [fieldData{2} fieldData{3}];
    for i=4:length(fieldData);
        res(recordCount).(fields{i-1}) = fieldData{i};
    end
    tline = fgets(fid);
    recordCount = recordCount + 1;
end
fclose(fid);

end

