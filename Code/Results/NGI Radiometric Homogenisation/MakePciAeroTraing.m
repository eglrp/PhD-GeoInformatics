%%
%make PCI aero-tri file
%Note: we should rather make 4 band ims (and invert lut?) before rectifying

%TO DO
%- Compress rectified tiffs same as source ones 
%- Make 4 band images
%- Set 12-16bit ot scaling or whatever is needed to get arcmap to display
%them properly

close all; clear all;
% idx = [2 3 68 69 72 73 433 434 436 437 440 441 ...
%     138 139 142 143 208 209 212 213 278 289 282 283 348 349 352 353 444 445 448 449 ...
%     452 453 455 456 459 460 462 463 465 466 469 470];

ngiOrFileName = {...
	'E:\Triangulation\3318B_2016_1142\extori3318B_2016_1142_lo19wgs84n_e_rect.or',...
    'E:\Triangulation\3318D_2016_1143\extori3318D_2016_1143_lo19wgs84n_e_rect.or',...
    };

outFile = {...
	'C:\Data\Development\Projects\PhD GeoInformatics\Docs\PCI\NGI Orthorectification\extori3318B_2016_1142_lo19wgs84n_e_rect.txt',...
    'C:\Data\Development\Projects\PhD GeoInformatics\Docs\PCI\NGI Orthorectification\extori3318D_2016_1143_lo19wgs84n_e_rect.txt',...
    };

imDir = {...
    'E:\Raw\3318B_2016_1142',...
    'E:\Raw\3318D_2016_1143',...
    };

% 3324C_2015_1004_01_0002_RGB

for i = 1:length(ngiOrFileName)
    MakePCIExtOrTextFile(ngiOrFileName{i}, [], outFile{i}, imDir{i});
end


%%