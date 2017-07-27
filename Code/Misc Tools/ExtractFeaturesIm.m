function [dataIm, labels] = ExtractFeaturesIm(subIm)
    if isstruct(subIm)
        subIm = subIm.data;
    end
    rgGIR = subIm./repmat(sum(subIm, 3), [1 1 4]);

    %scaling here? suspect...
    scale = 1;
    ndvi = (subIm(:,:,4)*scale - subIm(:,:,1))./(subIm(:,:,1) + ...
        subIm(:,:,4)*scale);

    irRat = subIm(:,:,4)./(subIm(:,:,1) + eps);

    entropyI = entropyfilt(mean(subIm, 3)/(2^12), true(7,7));
    stdI = stdfilt(mean(subIm, 3), true(7,7));

    entropyIrRat = entropyfilt(irRat/(10), true(7,7));
    stdIrRat = stdfilt(irRat, true(7,7));

    labels = {};
    labels = [labels {'R','G','B','NIR'}];

    dataIm = zeros(size(subIm, 1), size(subIm, 2), 14); %NB feat len must be right here
    dataIm(:, :, 1:4) = subIm;

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

end

