%%
close all; clear all;
cd 'F:\MSc GeoInformatics\Data\NGI\Temp'
rgbFilename = 'F:\MSc GeoInformatics\Data\NGI\Temp\SpekboomThicketCloseupRgb1.png';
cirFilename = 'F:\MSc GeoInformatics\Data\NGI\Temp\SpekboomThicketCloseupCir1.png';

[combinedImRgbIr  cir ndvi] = CombineRgbCir(rgbFilename, cirFilename);

figure;
imshow(combinedImRgbIr(:, :, [1:3]))

figure;
imagesc(cir)


thr = graythresh(ndvi);
figure;
imagesc(ndvi>0.2)
colormap gray

figure;
imagesc(ndvi)
colormap gray

irNorm = cir(:, :, 1)./sum(cir, 3);
figure;
imagesc(irNorm)
colormap gray

t = entropyfilt(mean(combinedImRgbIr, 3));
figure;
imagesc(t)
colormap gray

t = stdfilt(mean(combinedImRgbIr, 3), true(9,9));
figure;
imagesc(t)
colormap gray

cirMask = DrawMask(cir);

%%
close all; clear all;
cd 'F:\MSc GeoInformatics\Data\NGI\Temp'
rgbFilename = 'F:\MSc GeoInformatics\Data\NGI\Temp\SpekboomThicketCloseupRgb1.png';
cirFilename = 'F:\MSc GeoInformatics\Data\NGI\Temp\SpekboomThicketCloseupCir1.png';
maskFilename = 'F:\MSc GeoInformatics\Data\NGI\Temp\SpekboomThicketCloseupMask1.mat';
bgMaskFilename = 'F:\MSc GeoInformatics\Data\NGI\Temp\SpekboomThicketCloseupBgMask1.mat';

dataFilename = 'F:\MSc GeoInformatics\Data\NGI\Temp\SpekboomThicketCloseupDataset1.mat';

[combinedImRgbIr cirIm] = CombineRgbCir(rgbFilename, cirFilename);

mask1 = DrawMask(combinedImRgbIr(:,:,1:3), maskFilename);
bgMask1 = DrawMask(combinedImRgbIr(:,:,1:3), bgMaskFilename);
bushMask1 = DrawMask(combinedImRgbIr(:,:,1:3), bgMaskFilename);
bgMask = bgMask1;
save(bgMaskFilename, 'bgMask');

load(maskFilename)
load(bgMaskFilename)

[data imData] = ExtractFeatures(combinedImRgbIr, mask, bgMask);

save(dataFilename, 'data', 'mask', '*');

load(dataFilename)

scatterdui(data, ['ro'; 'bx'], 'legend')

data = setprior(data, 0);
[tr ts] = gendat(data, 0.5)

feats = [6 9];
w = qdc(tr(:, feats))

save('w1', 'w');
confmat(ts(:, feats)*w, 'disagreement')

out = im2feat(imData(:,:,feats))*w*classc;

outIm = reshape(+out(:,2), size(mask));

figure;
subplot(1,2,1)
imshow(outIm);
subplot(1,2,2)
imshow(cirIm);
%%

rgbFilename2 = 'F:\MSc GeoInformatics\Data\NGI\Temp\SpekboomThicketCloseupRgb2.png';
cirFilename2 = 'F:\MSc GeoInformatics\Data\NGI\Temp\SpekboomThicketCloseupCir2.png';

[combinedImRgbIr2 cirIm2] = CombineRgbCir(rgbFilename2, cirFilename2);

[data2 imData2] = ExtractFeatures(combinedImRgbIr2);

load 'w1'

feats = [6 9];

out = im2feat(imData2(:,:,feats))*w*classc;

outIm = reshape(+out(:,2), [size(combinedImRgbIr2, 1) size(combinedImRgbIr2, 2)]);

figure;
subplot(1,2,1)
imshow(outIm);
subplot(1,2,2)
imshow(cirIm2);

%%

close all; clear all;
cd 'F:\MSc GeoInformatics\Data\NGI\Temp'
rgbFilename = 'F:\MSc GeoInformatics\Data\NGI\Temp\SpekboomThicketCloseupRgb1.png';
cirFilename = 'F:\MSc GeoInformatics\Data\NGI\Temp\SpekboomThicketCloseupCir1.png';
maskFilename = 'F:\MSc GeoInformatics\Data\NGI\Temp\SpekboomThicketCloseupMask1.mat';
bgMaskFilename = 'F:\MSc GeoInformatics\Data\NGI\Temp\SpekboomThicketCloseupBgMask1.mat';
bushMaskFilename = 'F:\MSc GeoInformatics\Data\NGI\Temp\SpekboomThicketCloseupBushMask1.mat';

dataFilename = 'F:\MSc GeoInformatics\Data\NGI\Temp\SpekboomThicketCloseupDataset1.mat';

[combinedImRgbIr cirIm] = CombineRgbCir(rgbFilename, cirFilename);

load(maskFilename)
load(bgMaskFilename)
load(bushMaskFilename)

% bushMask = DrawMask(combinedImRgbIr(:,:,1:3), bushMaskFilename, mask | bgMask);

% bushMask = mask;
% save(bushMaskFilename, 'bushMask');

combinedMask = zeros(size(mask));
combinedMask(mask & (~bushMask)) = 1;
combinedMask(bgMask & (~bushMask)) = 2;
combinedMask(bushMask) = 3;

classLabels = {'Spekboom', 'Background', 'Tree'};
[data imData] = ExtractFeatures(combinedImRgbIr, combinedMask, classLabels);

save(dataFilename, 'data', 'imData', 'mask', 'combinedMask');
% load(dataFilename)

scatterdui(data(1:5:end,:), ['ro'; 'bx'; 'g^'], 'legend')

data = setprior(data, [0.4 0.2 0.4]);
idx = ~any(isinf(+data), 2);
data = data(idx, :);

% [wfs,rfs] = featselo(tr(1:5:end,:), scalem([], 'variance')*naivebc, 3, ts(1:5:end,:)); %NIR G / bN NIR / entropyIrRat

[tr ts] = gendat(data, 0.5)

feats = [6 9];
% feats = [4 7];
w = qdc(tr(:, feats))
figure;
scatterd(tr(:, feats), 'legend')
plotc(w)

save('w1', 'w');
confmat(ts(:, feats)*w, 'disagreement');
c = confmat(ts(:, feats)*w);
cn = c./repmat(sum(c, 2), 1, size(c, 2))

[CRn] = ReduceConfMat(cn, {[1 3], [2]}, true)
% CRn = CR./repmat(sum(CR, 2), 1, size(CR, 2))
% CRn =
%    0.960600308777462   0.039399691222538
%    0.050822427050975   0.949177572949025

out = im2feat(imData(:,:,feats))*w*classc;

VisualiseClassification(out, combinedImRgbIr(:, :, 1:3), cirIm);

figure;
h1 = subplot(2,2,1);
imshow(combinedImRgbIr(:,:,1))
h2 = subplot(2,2,2);
imshow(combinedImRgbIr(:,:,2))
h3 = subplot(2,2,3);
imshow(combinedImRgbIr(:,:,3))
h4 = subplot(2,2,4);
imshow(combinedImRgbIr(:,:,4))
linkaxes([h1 h2 h3 h4])


figure;
h1 = subplot(2,2,1);
imagesc(imData(:,:,4))
colormap('gray')
h2 = subplot(2,2,2);
imagesc(imData(:,:,9))
colormap('gray')
h3 = subplot(2,2,3);
imagesc(imData(:,:,10))
colormap('gray')
h4 = subplot(2,2,4);
imagesc(imData(:,:,13))
colormap('gray')
linkaxes([h1 h2 h3 h4])

PlotFeatureImage(imData, cellstr(getfeatlab(data)))

%%
close all; clear all;
cd 'F:\MSc GeoInformatics\Data\NGI\Temp'
rgbFilename = 'F:\MSc GeoInformatics\Data\NGI\Temp\SpekboomThicketCloseupRgb2.png';
cirFilename = 'F:\MSc GeoInformatics\Data\NGI\Temp\SpekboomThicketCloseupCir2.png';
dataFilename = 'F:\MSc GeoInformatics\Data\NGI\Temp\SpekboomThicketCloseupDataset2.mat';

[combinedImRgbIr cirIm] = CombineRgbCir(rgbFilename, cirFilename);

[data imData] = ExtractFeatures(combinedImRgbIr);

save(dataFilename, 'data', 'imData');

feats = [6 9];
load('w1');

out = im2feat(imData(:,:,feats))*w*classc;

VisualiseClassification(out, combinedImRgbIr(:, :, 1:3), cirIm);
