function [invLut lut] = ReadIntergraphLutFile(lutFileName)

fid = fopen(lutFileName);
for i = 1:4
    line = fgetl(fid);
end
prevPos = ftell(fid);
while (~feof(fid))
    line = fgetl(fid);
    if (length(line) > 3 && line(1) ~= '#')
        fseek(fid, prevPos, -1);
        lut = fscanf(fid, '%u');
        break;
    end
    prevPos = ftell(fid);
end
lut = uint16(reshape(lut, 4, []))';

%invert lut
invLut = lut;
for i = 2:4
%     uniquer
    %we should do something better with repeated vals like take the mean
    invLut(:, i) = interp1q(double(lut(:,i)), double(lut(:,1)), double(lut(:,1)));
end
%hack the last line
tmp = invLut(end-1, 2:4);
tmp = tmp + (4095.-tmp)/2;
invLut(end, 2:4) = tmp; 

%DH HACK!!! ASSUME BGR ORDERING!!!!
% lut(:,2:4) = lut(:,4:-1:2);
% invLut(:,2:4) = invLut(:,4:-1:2);




