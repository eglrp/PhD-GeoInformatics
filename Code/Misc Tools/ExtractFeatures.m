function [data dataIm] = ExtractFeatures(imRgbIr, combinedMask, classLabels)

% if (nargin == 1)
%     mask = true(size(imRgbIr));
%     bgMask = false(size(imRgbIr));
% end

%when you normalise by intensity, you really want to be normalising by the
%real spectrum, whereas here, we are normalising by a weighted version (esp IR)
% mask = mask & (~bushMask);


% rgG = imRgbIr(:,:,1:3)./repmat(sum(imRgbIr(:,:,1:3), 3), [1 1 3]);

rgGIR = uint32(imRgbIr)./uint32(repmat(sum(imRgbIr, 3), [1 1 4]));

%scaling here? suspect...
scale = 1;
ndvi = (imRgbIr(:,:,4)*scale - imRgbIr(:,:,1))./(imRgbIr(:,:,1) + ...
    imRgbIr(:,:,4)*scale);

irRat = imRgbIr(:,:,4)./(imRgbIr(:,:,1) + eps);

entropyI = entropyfilt(mean(imRgbIr, 3));
stdI = stdfilt(mean(imRgbIr, 3), true(9,9));

entropyIrRat = entropyfilt(irRat);
stdIrRat = stdfilt(irRat, true(9,9));

features = [];
bgFeatures = [];
labels = {};
labels = [labels {'R','G','B','NIR'}];
dataIm = zeros(size(imRgbIr, 1), size(imRgbIr, 2), 14); %NB feat len must be right here
dataIm(:, :, 1:4) = imRgbIr;

labels = [labels {'rN','gN','bN','nirN'}];
dataIm(:, :, 5:8) = rgGIR;

labels = [labels {'NDVI'}];
dataIm(:, :, 9) = ndvi;

labels = [labels {'irRat'}];
dataIm(:, :, 10) = irRat;

labels = [labels {'entropyI'}];
dataIm(:, :, 11) = entropyI;

labels = [labels {'stdI'}];
dataIm(:, :, 12) = stdI;

labels = [labels {'entropyIrRat'}];
dataIm(:, :, 13) = entropyIrRat;

labels = [labels {'stdIrRat'}];
dataIm(:, :, 14) = stdIrRat;

if (nargin > 1 && ~isempty(combinedMask))
    data = [];
    classLabelsV = {};
    featuresV = [];
    for i = 1:max(combinedMask(:))
        mask  = combinedMask == i;
        mask = repmat(mask, [1 1 size(dataIm, 3)]);
        features = reshape(dataIm(mask), [], size(dataIm, 3));
        classLabelsV = [classLabelsV; repmat(classLabels(i), size(features, 1), 1)];
        featuresV = [featuresV; features];
    end
    data = dataset(featuresV, classLabelsV);
    data = setfeatlab(data, labels);
else
    data = [];
end
