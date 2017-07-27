function [data imData] = ExtractFeatures(imRgbIr, mask, bgMask)

if (nargin == 1)
    mask = true(size(imRgbIr));
    bgMask = false(size(imRgbIr));
end

%when you normalise by intensity, you really want to be normalising by the
%real spectrum, whereas here, we are normalising by a weighted version (esp IR)

rgG = imRgbIr(:,:,1:3)./repmat(sum(imRgbIr(:,:,1:3), 3), [1 1 3]);

rgGIR = imRgbIr./repmat(sum(imRgbIr, 3), [1 1 4]);

%scaling here? suspect...
scale = 1;
ndvi = (imRgbIr(:,:,4)*scale - imRgbIr(:,:,1))./(imRgbIr(:,:,1) + ...
    imRgbIr(:,:,4)*scale);

irRat = imRgbIr(:,:,4)./imRgbIr(:,:,1);

entropyI = entropyfilt(mean(imRgbIr, 3));
stdI = stdfilt(mean(imRgbIr, 3), true(9,9));

entropyIrRat = entropyfilt(irRat);
stdIrRat = stdfilt(irRat, true(9,9));

features = [];
bgFeatures = [];
labels = {};
labels = [labels {'R','G','B','NIR'}];
for i = 1:size(imRgbIr, 3)
    band = imRgbIr(:,:,i);
    features(:, end+1) = band(mask); 
    bgFeatures(:, end+1) = band(bgMask); 
end

labels = [labels {'rN','gN','bN','nirN'}];
for i = 1:size(rgGIR, 3)
    band = rgGIR(:,:,i);
    features(:, end+1) = band(mask); 
    bgFeatures(:, end+1) = band(bgMask); 
end

labels = [labels {'NDVI'}];
features(:, end+1) = ndvi(mask); 
bgFeatures(:, end+1) = ndvi(bgMask); 

labels = [labels {'irRat'}];
features(:, end+1) = irRat(mask); 
bgFeatures(:, end+1) = irRat(bgMask); 

labels = [labels {'entropyI'}];
features(:, end+1) = entropyI(mask); 
bgFeatures(:, end+1) = entropyI(bgMask); 

labels = [labels {'stdI'}];
features(:, end+1) = stdI(mask); 
bgFeatures(:, end+1) = stdI(bgMask); 

labels = [labels {'entropyIrRat'}];
features(:, end+1) = entropyIrRat(mask); 
bgFeatures(:, end+1) = entropyIrRat(bgMask); 

labels = [labels {'stdIrRat'}];
features(:, end+1) = stdIrRat(mask); 
bgFeatures(:, end+1) = stdIrRat(bgMask); 

classLabels = [repmat({'Spekboom'}, size(features, 1), 1); repmat({'Background'}, size(bgFeatures, 1), 1)];
data = dataset([features; bgFeatures], classLabels);
data = setfeatlab(data, labels);
