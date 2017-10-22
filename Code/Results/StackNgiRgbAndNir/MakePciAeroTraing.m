%%
%make PCI aero-tri file
%Note: we should rather make 4 band ims (and invert lut?) before rectifying

%TO DO
%- Compress rectified tiffs same as source ones 
%- Make 4 band images
%- Set 12-16bit ot scaling or whatever is needed to get arcmap to display
%them properly

close all;clear all;
% idx = [2 3 68 69 72 73 433 434 436 437 440 441 ...
%     138 139 142 143 208 209 212 213 278 289 282 283 348 349 352 353 444 445 448 449 ...
%     452 453 455 456 459 460 462 463 465 466 469 470];

ngiOrFileName = {...
	'E:\Triangulation\3323D_2015_1001\extori3323D_2015_1001_lo23wgs84n_e_rect.or',...
    'E:\Triangulation\3324C_2015_1004\extori3324C_2015_1004_lo25wgs84n_e_rect.or',...
    };

outFile = {...
	'C:\Data\Development\Projects\PhD GeoInformatics\Code\Results\StackNgiRgbAndNir\extori3323D_2015_1001_lo23wgs84n_e_rect.txt',...
    'C:\Data\Development\Projects\PhD GeoInformatics\Code\Results\StackNgiRgbAndNir\extori3324C_2015_1004_lo25wgs84n_e_rect.txt',...
    };

imDir = {...
    'E:\Unrectified_Aerials\3323D_2015_1001\',...
    'E:\Unrectified_Aerials\3324C_2015_1004\',...
    };

% 3324C_2015_1004_01_0002_RGB

for i = 1:length(ngiOrFileName)
    MakePCIExtOrTextFile(ngiOrFileName{i}, [], outFile{i}, imDir{i});
end


%%