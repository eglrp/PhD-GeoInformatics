%%
fileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\3321C_2010_318\o3321C_2010_318_01_0003_RGBIR.tif';
fileName = 'F:\MSc GeoInformatics\Data\NGI\Rectified\RGB\3320AA\3320AA_01_2010_312_RGB_RECT.tif';
fileName = 'F:\MSc GeoInformatics\Data\NGI\Rectified\RGB\3320AA\3320AA_01_2010_312_RGB_RECT.tfw';
fileName = 'F:\MSc GeoInformatics\Data\CGA\DEMs\SUDEM\x3321c_1_L2.tif';

gi = geotiffinfo(fileName)

%%
h = worldmap('Africa');
getm(h, 'MapProjection')
geoshow('landareas.shp', 'FaceColor', [0.15 0.5 0.15])
geoshow('worldlakes.shp', 'FaceColor', 'cyan')
geoshow('worldrivers.shp', 'Color', 'blue')
geoshow('worldcities.shp', 'Marker', '.',...
                           'Color', 'red')
%%
cd(fullfile(matlabroot,'toolbox','map','mapdata'))
mapview
boston_roads = shaperead('boston_roads.shp');
surveyFeetPerMeter = unitsratio('survey feet','meter');
for k = 1:numel(boston_roads)
    boston_roads(k).X = surveyFeetPerMeter * boston_roads(k).X;
    boston_roads(k).Y = surveyFeetPerMeter * boston_roads(k).Y;
end 
roadcolors = makesymbolspec('Line', ...
{'CLASS',1,'Color',[1 1 1]}, {'CLASS',2,'Color',[1 1 0]}, ...
{'CLASS',3,'Color',[0 1 0]}, {'CLASS',4,'Color',[0 1 1]}, ...
{'CLASS',5,'Color',[1 0 1]}, {'CLASS',6,'Color',[0 0 1]})

%%
figure
mapshow boston.tif
axis image off

% The orthophoto is in survey feet and the roads are in meters.
% Convert the road units to feet before overlaying them.
S = shaperead('boston_roads.shp');
surveyFeetPerMeter = unitsratio('sf','meter');
x = surveyFeetPerMeter * [S.X];
y = surveyFeetPerMeter * [S.Y]; 
mapshow(x,y)

[Z, R] = sdtsdemread('9129CATD.DDF');

% View the Mount Washington terrain data as a mesh.
figure
mapshow(Z, R, 'DisplayType', 'mesh');
demcmap(Z)

%%
load korea
S = shaperead('landareas', 'UseGeoCoords', true);

figure;
worldmap(map, refvec)
geoshow(map, refvec, 'DisplayType', 'texturemap');
demcmap(map)
axis off

geoshow([S.Lat], [S.Lon], 'Color', 'black');

%%
clear all

[ortho, cmap] = imread('concord_ortho_w.tif');
R = worldfileread('concord_ortho_w.tfw', 'planar', size(ortho));
figure
mapshow(ortho, cmap, R)

xLimits = get(gca,'XLim');
yLimits = get(gca,'YLim');
pond = shaperead('concord_hydro_area.shp', 'RecordNumbers', 14);
hold on
mapshow(pond, 'FaceColor', [0.3 0.5 1], 'EdgeColor', 'black')
mapshow('concord_roads.shp', 'Color', 'red', 'LineWidth', 1);
xlabel('easting in meters')
ylabel('northing in meters')

set(gca,'XLim',xLimits,'YLim',yLimits)

%%
clear all
[X, cmap] = imread('concord_ortho_w.tif');
I_orig = ind2gray(X, cmap);

R_orig = worldfileread('concord_ortho_w.tfw','planar',size(X));

currentFormat = get(0,'format');
format short g
R_orig

I_half = imresize(I_orig, size(I_orig)/2, 'bicubic');

R_half = R_orig;
R_half.RasterSize = size(I_half)

figure
h1 = mapshow(I_orig,R_orig);
ax1 = get(h1,'Parent');
set(ax1, 'XLim', [208000 208250], 'YLim', [911800 911950],'TickDir','out')

figure
h2 = mapshow(I_half,R_half);
ax2 = get(h2,'Parent');
set(ax2, 'XLim', [208000 208250], 'YLim', [911800 911950],'TickDir','out')

x = 208202.21;
y = 911862.70;
line(x, y, 'Parent', ax1, 'Marker', '+', 'MarkerEdgeColor', 'r');
line(x, y, 'Parent', ax2, 'Marker', '+', 'MarkerEdgeColor', 'r');

[xIntrinsic1, yIntrinsic1] = R_orig.worldToIntrinsic(x, y)
[xIntrinsic2, yIntrinsic2] = R_half.worldToIntrinsic(x, y)

format(currentFormat);

%%
clear all

[baseImage1,cmap1] = imread('concord_ortho_w.tif');
[baseImage2,cmap2] = imread('concord_ortho_e.tif');

currentFormat = get(0,'format');
format short g
R1 = worldfileread('concord_ortho_w.tfw','planar',size(baseImage1))
R2 = worldfileread('concord_ortho_e.tfw','planar',size(baseImage2))

close all
mapshow(baseImage1,cmap1,R1)
ax1 = gca;
mapshow(ax1,baseImage2,cmap2,R2)
title('Map View, Massachusetts State Plane Coordinates');
xlabel('Easting (meters)');
ylabel('Northing (meters)');

inputImage = imread('concord_aerial_sw.jpg');
figure, imshow(inputImage)
title('Unregistered Aerial Photograph')

%%
clear all; close all;

imageFileName = 'D:\Data\Development\Projects\MSc GeoInformatics\modisSa.tif';

if true
    shapeFileName = 'D:\Data\Development\Projects\MSc GeoInformatics\Docs\ArcGIS\karoo_veg_bound_DD_WGS84.shp';
    shp = shaperead(shapeFileName, 'UseGeoCoords', true);
    symbolSpec = makesymbolspec('Polygon', {'Default', 'FaceColor', '', 'FaceAlpha', 0,...
            'EdgeColor', 'r', 'EdgeAlpha', 1});
else
    shapeFileName = 'D:\Data\Development\Projects\MSc GeoInformatics\Docs\ArcGIS\karoo_veg_poly20_DD_WGS84.shp';
    si = shapeinfo(shapeFileName);
    shp = shaperead(shapeFileName, 'UseGeoCoords', true, 'Attributes', {'BIOME'});
    classNames = unique({shp.BIOME});
    c = hsv(10);
    rule = {};
    for i = 1:length(classNames)
%         rule{i} = {'BIOME', classNames{i}, 'FaceColor', c(i,:), 'FaceAlpha', 0.5,...
%             'EdgeColor', c(i,:), 'EdgeAlpha', 0.5};
%         rule{i} = {'BIOME', classNames{i}, 'EdgeColor', 'k', 'EdgeAlpha', 1, ...
%             'FaceColor', c(i,:), 'FaceAlpha', 0.6};
        rule{i} = {'BIOME', classNames{i}, 'FaceColor', c(i,:), 'FaceAlpha', 0};
    end
    rule{end + 1} = {'Default', 'Visible', 'off'};
    symbolSpec = makesymbolspec('Polygon', rule{:});
end

shpProj = shp;
shpProj = rmfield(shpProj, {'Lat', 'Lon'});
proj = geotiffinfo(imageFileName);
for i = 1:length(shpProj)
    [shpProj(i).X shpProj(i).Y] = projfwd(proj, ...
        shp(i).Lat, shp(i).Lon);
    shpProj(i).BoundingBox = [min(shpProj(i).X) min(shpProj(i).Y); max(shpProj(i).X) max(shpProj(i).Y)];
end
%     figure;
%     mapshow(shpProj, 'SymbolSpec', symbolSpec);

% [x, y] = projfwd(proj, lat, lon);
% mapshow(im, imRef);

gi = geotiffinfo(imageFileName)
[im, imRef] = geotiffread(imageFileName);
imRef.RasterInterpretation = 'cell';

imSmall = imresize(im, [size(im,1) size(im,2)]/4, 'nearest');

imRefSmall = imRef;
imRefSmall.RasterSize = size(imSmall)

figure
mapshow(imSmall, imRefSmall); %, 'DisplayType', 'surface');
% mapshow(im, imRef); %, 'DisplayType', 'surface');
% mapshow(shpProj, 'FaceColor', [0.3 0.5 1], 'EdgeColor', 'black')
mapshow(shpProj, 'SymbolSpec', symbolSpec)
hold on
axis image off

%this shows the alternative of using geoshow to do the projection with 
%correctly configures map axes

figure
ax = axesm('mapprojection', 'eqaconicstd', 'geoid', almanac('earth', 'wgs84', proj.UOMLength), 'MapParallels', ...
    [proj.ProjParm(1:2)]);

mapshow(imSmall, imRefSmall); 
setm(ax, 'MapLatLimit', [-35 -20], 'MapLonLimit', [15 35])
axis off; gridm; mlabel on; plabel on; framem on
geoshow(shp, 'SymbolSpec', symbolSpec)
setm(ax, 'MapLatLimit', [-35 -20], 'MapLonLimit', [15 35])
setm(ax, 'MLabelLocation', 5)
setm(ax, 'PLabelLocation', 5)


%%

imageFileName = 'D:\Data\Development\Projects\MSc GeoInformatics\Code\Results\BatchProcess\3321A_2010_316_01_0003_RGB_16b2.tif';
% imagePyrFileName = 'D:\Data\Development\Projects\MSc GeoInformatics\Code\Results\BatchProcess\3321A_2010_316_01_0003_RGB_Pyr.tif';

imagePyrFileName = rsetwrite(imageFileName);

imtool(imagePyrFileName)

%%
clear all;close all;
cd 'D:\Data\Development\Projects\MSc GeoInformatics\Code\Results\LutTest\'
imageFileName = '3321D_2010_319_01_0033_CIR1.tif';

worldFileName = getworldfilename(imageFileName);

im = imread(imageFileName);
im = uint8(im/16);

r = worldfileread(worldFileName, 'planar', [size(im,1) size(im,2)])

figure
mapshow(im, r); %, 'DisplayType', 'surface');

%%
clear all; close all;
imageFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\3321D_2010_319\o3321D_2010_319_01_0003_RGBIR.tif';

proj = geotiffinfo(imageFileName);

% worldFileName = getworldfilename(imageFileName);

[im, R] = geotiffread(imageFileName);


ellipsoidName = 'clarke66';

mapUnits = proj.HorizontalUnits;

clarke66 = almanac('earth', ellipsoidName, 'metres')
type 9129.txt

%NOTES
% - usually possible to recover geographic coordinates if the projection parameters
%   and datum are known. Using this information, you can perform an inverse projection, 
%   running the projection backward to solve for latitude and longitude. The toolbox can 
%   perform accurate inverse projections for any of its projection functions as long as the 
%   original projection parameters and reference ellipsoid (or spherical radius) are provided to it.

