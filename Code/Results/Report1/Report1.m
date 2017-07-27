%%
close all; clear all;
cd 'D:\Data\Development\Projects\MSc GeoInformatics\Code\Results\Report1'
rgbFilename = 'F:\MSc GeoInformatics\Data\NGI\Temp\SpekboomThicketCloseupRgb1.png';
cirFilename = 'F:\MSc GeoInformatics\Data\NGI\Temp\SpekboomThicketCloseupCir1.png';
dataFilename = 'F:\MSc GeoInformatics\Data\NGI\Temp\SpekboomThicketCloseupDataset1.mat';

[combinedImRgbIr cirIm] = CombineRgbCir(rgbFilename, cirFilename);

% [data imData] = ExtractFeatures(combinedImRgbIr);

load(dataFilename);

data = setprior(data, 0);
[tr ts] = gendat(data, 0.5);
feats = [6 9];
w = tr(:, feats)*qdc;

confmat(ts(:, feats)*w, 'disagreement')
c = confmat(ts(:, feats)*w);
cn = c./repmat(sum(c, 2), 1, size(c, 2));
[cnr] = ReduceConfMat(cn, {[1 3] [2]}, true);
err = 1-mean(diag(cnr))

load('w1');

out = im2feat(imData(:,:,feats))*w*classc;

VisualiseClassification(out, combinedImRgbIr(:, :, 1:3), cirIm);


figure;
imshow(combinedImRgbIr(:, :, 1:3))
title('Training scene')

%%
close all; clear all;
cd 'D:\Data\Development\Projects\MSc GeoInformatics\Code\Results\Report1';
rgbFilename2 = 'F:\MSc GeoInformatics\Data\NGI\Temp\SpekboomThicketCloseupRgb2.png';
cirFilename2 = 'F:\MSc GeoInformatics\Data\NGI\Temp\SpekboomThicketCloseupCir2.png';
dataFilename2 = 'F:\MSc GeoInformatics\Data\NGI\Temp\SpekboomThicketCloseupDataset2.mat';

[combinedImRgbIr2 cirIm2] = CombineRgbCir(rgbFilename2, cirFilename2);

[data2 imData2] = ExtractFeatures(combinedImRgbIr2);

% load(dataFilename);

feats = [6 9];
load('w1');

out2 = im2feat(imData2(:,:,feats))*w*classc;

VisualiseClassification(out2, combinedImRgbIr2(:, :, 1:3), cirIm2);

%%
close all; clear all;
cd 'D:\Data\Development\Projects\MSc GeoInformatics\Code\Results\Report1';
rgbFilename = 'F:\MSc GeoInformatics\Data\NGI\Temp\NonSpekboomThicketCloseupRgb3.png';
cirFilename = 'F:\MSc GeoInformatics\Data\NGI\Temp\NonSpekboomThicketCloseupCir3.png';
dataFilename = 'F:\MSc GeoInformatics\Data\NGI\Temp\NonSpekboomThicketCloseupDataset3.mat';

[combinedImRgbIr cirIm] = CombineRgbCir(rgbFilename, cirFilename);

[data imData] = ExtractFeatures(combinedImRgbIr);

save(dataFilename, 'data', 'imData');
% load(dataFilename);


feats = [6 9];
load('w1');

out = im2feat(imData(:,:,feats))*w*classc;

VisualiseClassification(out, combinedImRgbIr(:, :, 1:3), cirIm);


%%
%Sanity tests 2013

close all; clear all;
cd 'D:\Data\Development\Projects\MSc GeoInformatics\Code\Results\Report1';

rgbFilenames = {'MixedCloseupRgb.png';...
    'SpekboomThicketCloseupRgb3.png';...
    'SpekboomThicketCloseupRgb4.png';...
    'SpekboomThicketCloseupRgb5.png'};

cirFilenames = {'MixedCloseupCir.png';...
    'SpekboomThicketCloseupCir3.png';...
    'SpekboomThicketCloseupCir4.png';...
    'SpekboomThicketCloseupCir5.png'};

for i = 1:length(rgbFilenames)
    [combinedImRgbIr cirIm] = CombineRgbCir(rgbFilenames{i}, cirFilenames{i});
    [data imData] = ExtractFeatures(combinedImRgbIr);

    % load(dataFilename);

    feats = [6 9];
    load('w1');
    out = im2feat(imData(:,:,feats))*w*classc;
    VisualiseClassification(out, combinedImRgbIr(:, :, 1:3), cirIm);
end
