close all;clear all;
%%
%use gdalwarp to warp modis to ngi format

cd 'D:\Data\Development\Projects\MSc GeoInformatics\Code\Results\MappingToolboxExp'

srcFileName = 'D:\Data\Development\Projects\MSc GeoInformatics\modisSa.tif';
destFileName__ = 'F:\MSc GeoInformatics\Data\NGI\Rectified\RGB\3320AA\3320AA_01_2010_312_RGB_RECT.tif';
outFileName = 'D:\Data\Development\Projects\MSc GeoInformatics\Code\Results\MappingToolboxExp\modisSaOut.tif';

dos(sprintf('gdalinfo "%s"', srcFileName)) %Y["EPSG","4326"]] 
dos(sprintf('gdalinfo "%s"', destFileName__))

dos(sprintf('gdalwarp "%s" "%s"', ...    
    srcFileName, outFileName));

dos(sprintf('gdalwarp -multi -r bilinear -tr 460 460 -t_srs "%s" "%s" "%s"', ...
    'D:\Data\Development\Projects\MSc GeoInformatics\Code\Results\MappingToolboxExp\dest.prf', ...
    srcFileName, outFileName));

si = geotiffinfo(srcFileName)
geotiffinfo(outFileName)
dos(sprintf('gdalinfo "%s"', outFileName))

newExtents = [-456000 -3960600 (-456000+3700*460) (-3960600+3300*460)];

dos(sprintf('gdalwarp -multi -r bilinear -tr 460 460 -te %d %d %d %d -t_srs "%s" "%s" "%s"', newExtents, ...
    'D:\Data\Development\Projects\MSc GeoInformatics\Code\Results\MappingToolboxExp\dest.prf', ...
    srcFileName, outFileName));

dos(sprintf('gdalinfo "%s"', outFileName))

%note that there are some cool options to gdalwarp, eg you can warp into a
%mosaic and use multiple threads

%%
%make a test mosaic

mosaicFolder = 'H:\MSc GeoInformatics\Data\NGI\Temp\';
outFile = [ 'GdalMosaicTest.tif'];

dos('gdal_merge.py')

fileList = dir([mosaicFolder '3320*.tif']);
fileStr = '';
for i = 1:length(fileList)
    fileStr = sprintf('%s "%s%s"', fileStr, '', fileList(i).name);
end

%delete output file first
cd(mosaicFolder);
delete(outFile)
dos(sprintf('gdal_merge.py -o %s -ps 10 10 %s', outFile, fileStr));

im = geotiffread(outFile);

figure;
% \imshow(im(:,:,[4 1 2])*16)
imshow(im)

%%
%Test gdal2tiles
%see also http://www.maptiler.org/
dos('gdalbuildvrt -tr 1 1 "H:\MSc GeoInformatics\Data\NGI\Temp\out.vrt" "H:\MSc GeoInformatics\Data\NGI\Temp\3320*.tif"');
%for google earth
dos('gdal2tiles.py -k -p raster -r near "H:\MSc GeoInformatics\Data\NGI\Temp\out.vrt" "H:\MSc GeoInformatics\Data\NGI\Temp\Tiles"');
%for google maps
% dos('gdal2tiles.py -k -r near "H:\MSc GeoInformatics\Data\NGI\Temp\out.vrt" "H:\MSc GeoInformatics\Data\NGI\Temp\Tiles"');
%Note - it works with 8bit RGB imnages but not 16bit RGBN images

%%
% convert HDF to ???

%%
%use gdalwarp to warp modis HDF to ngi projection

cd 'D:\Data\Development\Projects\MSc GeoInformatics\Code\Results\MappingToolboxExp'
MCD43A4.A2010033.h19v12.005.2010052115603.hdf
srcFileName = 'G:\MSc GeoInformatics\Data\MODIS\MCD43A4.A2010025.h19v12.005.2010043064233.hdf';
destFileName__ = 'G:\MSc GeoInformatics\Data\NGI\Rectified\RGB\3320AA\3320AA_01_2010_312_RGB_RECT.tif';
outFileName = 'D:\Data\Development\Projects\MSc GeoInformatics\Code\Results\MappingToolboxExp\modisSaOut.tif';

dos(sprintf('gdalinfo "%s"', srcFileName))
dos(sprintf('gdalinfo "%s"', destFileName__))

dos(sprintf('gdalwarp "%s" "%s"', ...    
    srcFileName, outFileName));

dos(sprintf('gdalwarp -multi -r bilinear -tr 460 460 -t_srs "%s" "%s" "%s"', ...
    'D:\Data\Development\Projects\MSc GeoInformatics\Code\Results\MappingToolboxExp\dest.prf', ...
    srcFileName, outFileName));

si = geotiffinfo(srcFileName)
geotiffinfo(outFileName)
dos(sprintf('gdalinfo "%s"', outFileName))

newExtents = [-456000 -3960600 (-456000+3700*460) (-3960600+3300*460)];

dos(sprintf('gdalwarp -multi -r bilinear -tr 460 460 -te %d %d %d %d -t_srs "%s" "%s" "%s"', newExtents, ...
    'D:\Data\Development\Projects\MSc GeoInformatics\Code\Results\MappingToolboxExp\dest.prf', ...
    srcFileName, outFileName));

dos(sprintf('gdalinfo "%s"', outFileName))

%%
%Test OSGeo4W with HDF support etc
% cd 'D:\Data\Development\Projects\MSc GeoInformatics\Code\Results\MappingToolboxExp'
cd 'C:\OSGeo4W\bin'
% MCD43A4.A2010033.h19v12.005.2010052115603.hdf
srcFileName = 'G:\MCD43A4.A2010025.h19v12.005.2010043064233.hdf';
destFileName__ = 'G:\MSc GeoInformatics\Data\NGI\Rectified\RGB\3320AA\3320AA_01_2010_312_RGB_RECT.tif';
outFileName = 'D:\Data\Development\Projects\MSc GeoInformatics\Code\Results\MappingToolboxExp\OsGeo4Wtst.tif';

dos(sprintf('gdalinfo -proj4 "%s"', srcFileName))
dos(sprintf('gdalinfo -proj4 "%s"', destFileName__))

% The below works in an OSGeo4W shell
% C:\>gdaltest HDF4_EOS:EOS_GRID:"G:\MCD43A4.A2010025.h19v12.005.2010043064233.hdf
% ":MOD_Grid_BRDF:Nadir_Reflectance_Band2 g:\osgeotst.tif

dos(sprintf('gdalwarp -t_srs "%s" HDF4_EOS:EOS_GRID:"%s":MOD_Grid_BRDF:Nadir_Reflectance_Band1 "%s"', ...    
    '+proj=tmerc +lat_0=0 +lon_0=21 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs',srcFileName, outFileName));

dos(sprintf('gdalinfo -proj4 "%s"', outFileName))

t = Tiff(outFileName);
im = t.read();
figure;
imagesc(im);
t.close;

%%--------------------------------------------------------------------------
% Look at the difference in PCS between NGI provided WGS84 and PCI generated
% Hartebeeshoek94 files

ngiFileName = 'G:\MSc GeoInformatics\Data\NGI\Rectified\RGB\3320AA\3320AA_01_2010_312_RGB_RECT.tif';
pciFileName = 'H:\MSc GeoInformatics\Data\NGI\My Rectified\3321C_2010_318\o3321C_2010_318_01_0004_RGBIR.tif';
pciFileName2 = 'H:\MSc GeoInformatics\Data\NGI\My Rectified\3321D_2010_319\o3321D_2010_319_01_0004_RGBIR.tif';

dos(sprintf('gdalinfo -proj4 "%s"', ngiFileName)) %'+proj=tmerc +lat_0=0 +lon_0=21 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ' 
dos(sprintf('gdalinfo -proj4 "%s"', pciFileName)) %+proj=tmerc +lat_0=0 +lon_0=21 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs '
dos(sprintf('gdalinfo -proj4 "%s"', pciFileName2)) %+proj=tmerc +lat_0=0 +lon_0=21 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs '

%- Note that only difference in the gen proj.4 is that PCI files have
%ellps=wgs84 whereas the NGI files have datum=wgs84
%- The PROJCS[...] strings are very different though - then it assumes
%defaults and generates the proj.4 strings

%Generate a file using PCI string as opposed to NGI string
outFileNamePci = [outFileName(1:end-4) '2.tif'];
dos(sprintf('gdalwarp -t_srs "%s" HDF4_EOS:EOS_GRID:"%s":MOD_Grid_BRDF:Nadir_Reflectance_Band1 "%s"', ...        
    '+proj=tmerc +lat_0=0 +lon_0=21 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs',srcFileName, outFileNamePci));

%now compare NGI and PCI string gen files
if false
    t = Tiff(outFileName);
    im = t.read();
    t.close;

    tPci = Tiff(outFileNamePci);
    imPci = tPci.read();
    tPci.close;

    figure;
    imagesc([im imPci]);

    figure;
    imagesc(log10(double(abs(im-imPci))));
else
    gi = geotiffinfo(outFileName);
    giPci = geotiffinfo(outFileNamePci);
    [im, r] = geotiffread(outFileName);
    [imPci, rPci] = geotiffread(outFileNamePci);

    figure;
    subplot(1,2,1)
    mapshow(uint16(im)*4, r)
    subplot(1,2,2)
    mapshow(uint16(imPci)*4, rPci)
end
%Not identical but v similar
%looking at the theory, I dont expect them to be identical though, the
%datum includes the origin and alignment which is not necessarily the same
%as the defaults with the the ellipse
%%

%--------------------------------------------------------------------------
% Try gdalwarp with "target aligned pixels" -tap option.
close all;clear all;
cd C:\OSGeo4W\bin
modisFileName = 'G:\MCD43A4.A2010025.h19v12.005.2010043064233.hdf';
ngiFileNames = {'G:\MSc GeoInformatics\Data\NGI\Rectified\RGB\3320AA\3320AA_02_2010_312_RGB_RECT.tif';...
    'G:\MSc GeoInformatics\Data\NGI\Rectified\RGB\3320AA\3320AA_03_2010_312_RGB_RECT.tif';...
    'G:\MSc GeoInformatics\Data\NGI\Rectified\RGB\3320AA\3320AA_04_2010_312_RGB_RECT.tif';...
    'G:\MSc GeoInformatics\Data\NGI\Rectified\RGB\3320AA\3320AA_05_2010_312_RGB_RECT.tif';...
    'G:\MSc GeoInformatics\Data\NGI\Rectified\RGB\3320AA\3320AA_05_2010_312_RGB_RECT.tif';...
    'G:\MSc GeoInformatics\Data\NGI\Rectified\RGB\3320AA\3320AA_06_2010_312_RGB_RECT.tif';...
    'G:\MSc GeoInformatics\Data\NGI\Rectified\RGB\3320AA\3320AA_07_2010_312_RGB_RECT.tif';...
    'G:\MSc GeoInformatics\Data\NGI\Rectified\RGB\3320AA\3320AA_08_2010_312_RGB_RECT.tif';...
    'G:\MSc GeoInformatics\Data\NGI\Rectified\RGB\3320AA\3320AA_09_2010_312_RGB_RECT.tif'};

%reproject modis to ngi co-ords
dos(sprintf('gdalinfo -proj4 "%s"', modisFileName)) %'+proj=tmerc +lat_0=0 +lon_0=21 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ' 
modisProjFileName = 'G:\modisProj.tif';
dos(sprintf('gdalwarp -r cubic -tr 250 250 -tap -t_srs "%s" HDF4_EOS:EOS_GRID:"%s":MOD_Grid_BRDF:Nadir_Reflectance_Band1 "%s"', ...        
    '+proj=tmerc +lat_0=0 +lon_0=21 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs', modisFileName, modisProjFileName));

[modisProjIm, modisProjR] = geotiffread(modisProjFileName);
figure
mapshow(uint16(abs(modisProjIm)*8), modisProjR)
hold on

% [x y] = modisProjR.intrinsicToWorld(1,1)

%resample ngi image to modis reprojected co-ords
for i = 7:length(ngiFileNames)
    dos(sprintf('gdalinfo -proj4 "%s"', ngiFileNames{i})) %'+proj=tmerc +lat_0=0 +lon_0=21 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ' 
    ngiProjFileName = sprintf('G:\\ngiProj%d.tif', i);
    dos(sprintf('gdalwarp -r bilinear -tr 250 250 -tap -t_srs "%s" "%s" "%s"', ...        
        '+proj=tmerc +lat_0=0 +lon_0=21 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs', ngiFileNames{i}, ngiProjFileName));
    [ngiProjIm, ngiProjR] = geotiffread(ngiProjFileName);
    mapshow(ngiProjIm, ngiProjR)
    rectangle('Position', [ngiProjR.XLimWorld(1) ngiProjR.YLimWorld(1) ...
        ngiProjR.RasterWidthInWorld ngiProjR.RasterHeightInWorld], 'EdgeColor', 'r')
end

%upsample back to orig res
% lanczos
upSampledFileName = [ngiProjFileName(1:end-4) '_.tif'];
dos(sprintf('gdalwarp -r cubicspline -tr .5 .5 -tap -t_srs "%s" "%s" "%s"', ...        
    '+proj=tmerc +lat_0=0 +lon_0=21 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs', ngiProjFileName, upSampledFileName));

[im, r] = geotiffread(upSampledFileName);

imshow(im)
%Note - tap does work, but it can put a black single pixel border around the
%NGI image.  Adjacent images can then end up without in between border
%pixels

%%---------------------------------------------------------------------------

%% BUDGET
%
% satrix 1000000
% absa mm 600000
% absa share 900000
% aims 470000
% investecbond 1500000
% investec ra 100000
% liberty 100000
% dbpf 550000

%retirement 1.2
%other 3.5

%%
% budget
% med aid 1500
% medication 2000
% fritz 1500
% rent 1000
% food 3000
% petrol+maint 1000
% rates + elec 500
% insurance 1000
% internet 500
% cell phone 200
% entertainment 1000
% grace+goodman 500

% 16000 %checked on bank statements (excl rent and RA)
% which can be covered by 10% growth on 2e6 of my investments
% then there is 2.7 left for retirement

% 16*2.7 @ 60

rs = 2.7*((1.1)^20)*1000000
ri = (16000)*((1.05)^30)

rs/(12*ri)


t=1:60;
x = t+37;
rst_ = 2.7e6*(1.1.^t);
% 
rit_ = 20000*(1.05.^t);

ra = 20;

rit = 20000;
rst = 2.7e6;
for i = 2:length(t)
    rit(i) = rit(i-1)*1.06; %inflation
    rst(i) = rst(i-1)*1.09; %retirement savings investment income
    if (t(i) > ra) %after retirement
        rst(i) = rst(i)-(12*rit(i)); %income spent after retirement
    else
        rst(i) = rst(i)+(12*rit(i)*0); %savings before retirement
    end
end

figure;
semilogy(x,rst)
hold all;
semilogy(x,rit)
legend({'savings','income'})
grid on

figure;
plot(x,rst)
hold all;
plot(x,rit)
legend({'savings','income'})
grid on


%%
clear all;

t = 1:15;
x = t+37;

rit = 15000;
rst = 2e6;
ra = 0;
for i = 2:length(t)
    rit(i) = rit(i-1)*1.06; %inflation
    rst(i) = rst(i-1)*1.1; %retirement savings investment income
    if (t(i) > ra) %after retirement
        rst(i) = rst(i)-(12*rit(i)); %income spent after retirement
    else
        rst(i) = rst(i)+(12*rit(i)*0); %savings before retirement
    end
end

figure;
semilogy(x,rst)
hold all;
semilogy(x,rit)
legend({'savings','income'})
grid on

figure;
plot(x,rst)
hold all;
plot(x,rit)
legend({'savings','income'})
grid on


%%
%Notes
%
% - The msg is more about how much return I get on what I have now than how
% much I save per month.
