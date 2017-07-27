%%
%Testing how to plot exported ping with meter scale
imFileName = 'D:\Data\Development\Projects\MSc GeoInformatics\Code\Results\Report1\MixedCloseupRgb.png';
refImFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\3321D_2010_319\o3321D_2010_319_04_0141_RGBIR.tif';
%'F:\MSc GeoInformatics\Data\NGI\Rectified\CIR\3322DA\3322DA_01_2010_323_CIR_RECT.tif';
im = imread(imFileName);
R = worldfileread(getworldfilename(rgbFilenames{end}));

ms = geotiff2mstruct(geotiffinfo(refImFileName));

figure
axesm(ms)
mapshow(im, R);
axis on
%%
%Sanity tests 2013

close all; clear all;
cd 'D:\Data\Development\Projects\MSc GeoInformatics\Code\Results\Report3';

rgbFilenames = {'MixedCloseupRgb.png';...
    'SpekboomThicketCloseupRgb3.png';...
    'SpekboomThicketCloseupRgb4.png';...
    'SpekboomThicketCloseupRgb5.png'};

cirFilenames = {'MixedCloseupCir.png';...
    'SpekboomThicketCloseupCir3.png';...
    'SpekboomThicketCloseupCir4.png';...
    'SpekboomThicketCloseupCir5.png'};

for i = 1:length(rgbFilenames)
    [combinedImRgbIr cirIm R] = CombineRgbCir(rgbFilenames{i}, cirFilenames{i});
    [data imData] = ExtractFeatures(combinedImRgbIr);

    % load(dataFilename);

    feats = [6 9];
    load('w1');
    out = im2feat(imData(:,:,feats))*w*classc;
    VisualiseClassification(out, combinedImRgbIr(:, :, 1:3), cirIm, R);
end
