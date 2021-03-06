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
    'F:\MSc GeoInformatics\Data\NGI\Aerial Triangulation\3321A_2010_316\extori3321A_2010_316_lo21wgs84n_e_rect.or',... 
    'F:\MSc GeoInformatics\Data\NGI\Aerial Triangulation\3321B_2010_317\extori3321B_2010_317_lo21wgs84n_e_rect.or',...
    'F:\MSc GeoInformatics\Data\NGI\Aerial Triangulation\3321C_2010_318\extori3321C_2010_318_lo21wgs84n_e_rect.or',...
    'F:\MSc GeoInformatics\Data\NGI\Aerial Triangulation\3321D_2010_319\extori3321D_2010_319_lo21wgs84n_e_rect.or',...
    'F:\MSc GeoInformatics\Data\NGI\Aerial Triangulation\3322A_2010_320\extori3322A_2010_320_lo23wgs84n_e_rect.or',...
    'F:\MSc GeoInformatics\Data\NGI\Aerial Triangulation\3322B_2010_321\extori3322B_2010_321_lo23wgs84n_e_rect.or',...
    'F:\MSc GeoInformatics\Data\NGI\Aerial Triangulation\3322C_2010_322\extori3322C_2010_322_lo23wgs84n_e_rect.or',...
    'F:\MSc GeoInformatics\Data\NGI\Aerial Triangulation\3322D_2010_323\extori3322D_2010_323_lo23wgs84n_e_rect.or',...
    };

outFile = {...
    'D:\Data\Development\Projects\MSc GeoInformatics\Docs\PCI\Batch\extori3321A_2010_316_lo21wgs84n_e_rect.txt';...
    'D:\Data\Development\Projects\MSc GeoInformatics\Docs\PCI\Batch\extori3321B_2010_317_lo21wgs84n_e_rect.txt';...
    'D:\Data\Development\Projects\MSc GeoInformatics\Docs\PCI\Batch\extori3321C_2010_318_lo21wgs84n_e_rect.txt';...
    'D:\Data\Development\Projects\MSc GeoInformatics\Docs\PCI\Batch\extori3321D_2010_319_lo21wgs84n_e_rect.txt';...
    'D:\Data\Development\Projects\MSc GeoInformatics\Docs\PCI\Batch\extori3322A_2010_320_lo23wgs84n_e_rect.txt';...
    'D:\Data\Development\Projects\MSc GeoInformatics\Docs\PCI\Batch\extori3322B_2010_321_lo23wgs84n_e_rect.txt';...
    'D:\Data\Development\Projects\MSc GeoInformatics\Docs\PCI\Batch\extori3322C_2010_322_lo23wgs84n_e_rect.txt';...
    'D:\Data\Development\Projects\MSc GeoInformatics\Docs\PCI\Batch\extori3322D_2010_323_lo23wgs84n_e_rect.txt';...
    };

imDir = {...
    'F:\MSc GeoInformatics\Data\NGI\Calibrated\3321A_2010_316\RGB\';...
    'F:\MSc GeoInformatics\Data\NGI\Calibrated\3321B_2010_317\RGB\';...
%     'F:\MSc GeoInformatics\Data\NGI\Calibrated\3321C_2010_318\RGB\';...
    'G:\MSc GeoInformatics\Data\NGI\My Calibrated\3321C_2010_318\';...
%     'F:\MSc GeoInformatics\Data\NGI\Calibrated\3321D_2010_319\RGB\';...
    'G:\MSc GeoInformatics\Data\NGI\My Calibrated\3321D_2010_319\';...
    'F:\MSc GeoInformatics\Data\NGI\Calibrated\3322A_2010_320\RGB\';...
    'F:\MSc GeoInformatics\Data\NGI\Calibrated\3322B_2010_321\RGB\';...
    'F:\MSc GeoInformatics\Data\NGI\Calibrated\3322C_2010_322\RGB\';...
    'F:\MSc GeoInformatics\Data\NGI\Calibrated\3322D_2010_323\RGB\';...
    };


for i = 1:length(ngiOrFileName)
    MakePCIExtOrTextFile(ngiOrFileName{i}, [], outFile{i}, imDir{i});
end


%%