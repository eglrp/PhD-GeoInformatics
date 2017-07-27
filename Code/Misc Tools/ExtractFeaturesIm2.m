function [dataIm, labels] = ExtractFeaturesIm2(subIm, varargin)
    %% 
    wPcaSpekboom = 'F:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\PcaSpekBoom.mat';
    wPcaAll = 'F:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\PcaAll.mat';
    wPcaRgG = 'F:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\PcaRgG.mat';
    lbpFiltSize = 4;
    lbpFiltRadius = 1;
    
    winSize = [5 5];
    feats = 1:50;
    ModifyDefaultArgs(varargin);

    if (ischar(wPcaSpekboom))
        prload(wPcaSpekboom);
    end

    if (ischar(wPcaAll))
        prload(wPcaAll);
    end

    if (ischar(wPcaRgG))
        prload(wPcaRgG);
    end
    
    %convert and scale
    %note that some of the algorithms require the image to be in the range
    %0-1
    if isstruct(subIm)
        subIm = double(subIm.data)/(2^14);
    elseif (max(subIm(:))>1)
        subIm = subIm/(2^14); %TO DO: bug - should be 2^12!!! Is this why our image classifier doesn't work?  Or is it a bug?  How the hell does the above work then?????
    end
    
    noDataMask = all(subIm==0, 3);

    featIm = reshape(subIm, [], size(subIm, 3));
    subImPcaSpekboom = featIm*wPcaSpekboom; %something like tasselled cap
    subImPcaAll = featIm*wPcaAll; %for texture extraction

    %scaling here? suspect...
    rgGIR = subIm./repmat(sum(subIm, 3), [1 1 4]);

    rgGIRIm = reshape(rgGIR, [], 4);
    subImPcaRgG = rgGIRIm*wPcaRgG;

    %% per-pixel features
    scale = 1;
    ndvi = (subIm(:,:,4)*scale - subIm(:,:,1))./(subIm(:,:,1) + ...
        subIm(:,:,4)*scale);

    irRat = subIm(:,:,4)./(10*subIm(:,:,1) + eps); %try scale 0-1

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

    labels = [labels {'tc1','tc2','tc3','tc4'}];
    dataIm(:, :, 11:14) = reshape(subImPcaSpekboom, size(subIm));

    labels = [labels {'pc1','pc2','pc3','pc4'}];
    dataIm(:, :, 15:18) = reshape(subImPcaAll, size(subIm));

    labels = [labels {'rc1','rc2','rc3','rc4'}];
    dataIm(:, :, 19:22) = reshape(subImPcaRgG, size(subIm));
%     dataIm = int16(dataIm(:,:,[6 9 11])*(2^12));
%     dataIm = double(dataIm(:,:,[6 9 11 12]));
%     return
    %% sliding win features
    subImPcaAll = reshape(subImPcaAll(:,1), size(subIm, 1), size(subIm, 2));
    
    [dataIm_ labels_] = SlidingWinFeatures(subImPcaAll, varargin{:}, 'feats', feats-23+1);
    labels = [labels CatLabels(labels_, 'Pc1')];
%     dataIm(:, :, 15:21) = dataIm_;

if true
    dataIm(:, :, 23:29) = dataIm_;

    if (max(feats) < 30)
        dataIm(repmat(noDataMask, [1 1 size(dataIm, 3)])) = 0;        
        return;
    end
    
    [dataIm_ labels_] = SlidingWinFeatures(irRat, varargin{:}, 'feats', feats-30+1);
    labels = [labels CatLabels(labels_, 'IrRat')];
    dataIm(:, :, 30:36) = dataIm_;

    [dataIm_ labels_] = SlidingWinFeatures(ndvi, varargin{:}, 'feats', feats-37+1);
    labels = [labels CatLabels(labels_, 'Ndvi')];
    dataIm(:, :, 37:43) = dataIm_;

    [dataIm_ labels_] = SlidingWinFeatures(rgGIR(:,:,2), varargin{:}, 'feats', feats-44+1);
    labels = [labels CatLabels(labels_, 'Gn')];
    dataIm(:, :, 44:50) = dataIm_;

else
    dataIm(:, :, 23:25) = dataIm_;
    dataIm = double(dataIm(:, :, [5 6 7 9 11 12 23:25])); %FEATS
%     dataIm = double(dataIm(:, :, [16 23:25])); %FEATS2
end
    dataIm(repmat(noDataMask, [1 1 size(dataIm, 3)])) = 0;
    
end

%%
function [dataIm labels] = SlidingWinFeatures(subIm, varargin)
% tmp = internal.images.isFigureAvailable;
% internal.images.isFigureAvailable = false;
%TO DO: we need to scale the features better to get them 0-1 for
%entropyfilt, also we need to offset them!!!

    wPcaSpekboom = 'F:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\PcaSpekBoom.mat';
    wPcaAll = 'F:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\PcaAll.mat';
    wPcaRgG = 'F:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\PcaRgG.mat';
    lbpFiltSize = 4;
    lbpFiltRadius = 1;
    
    winSize = [5 5];

    feats = 1:7;
    ModifyDefaultArgs(varargin);
    
%     filtR = generateRadialFilterLBP(lbpFiltSize, lbpFiltRadius); % DH 8/2016
    featFnParams = {@entropyfilt, {true(winSize)};...
        @stdfilt, {true(winSize)};...
        @nlfilter, {winSize, @(x) mean(x(:))} ;...
        @nlfilter, {winSize, @(x) median(x(:))} ;...
        @nlfilter, {winSize, @(x) kurtosis(x(:), 0)} ;...
        @nlfilter, {winSize, @(x) skewness(x(:), 0)} ;...
        % @efficientLBP, {'filtR', filtR, 'isRotInv', true,
        % 'isChanWiseRot', false} ;...  % DH 8/2016
        };
    labels = {'Entropy','Std','Mean','Median','Kurtosis','Skewness','Lbp'};

    dataIm = zeros(size(subIm, 1), size(subIm, 2), 7);
    for i = 1:size(featFnParams, 1)
        if any(feats == i)
            dataIm(:,:,i) = featFnParams{i,1}(subIm, featFnParams{i,2}{:});
        end
    end

%     labels = {};
%     if (any(feats==1))
%         dataIm(:,:,1) = entropyfilt(subIm, true(winSize));
%     else
%         dataIm(:,:,1) = zeros(size(subIm, 1), size(subIm, 2));
%     end
%     labels = [labels {'Entropy'}];
%     if (any(feats==2))
%         dataIm(:,:,2) = stdfilt(subIm, true(winSize));
%     else
%         dataIm(:,:,2) = zeros(size(subIm, 1), size(subIm, 2));
%     end
%     labels = [labels {'Std'}];
% if true    
%     dataIm(:,:,3) = nlfilter(subIm , winSize, @(x) mean(x(:)));
%     labels = [labels {'Mean'}];
%     dataIm(:,:,4) = nlfilter(subIm , winSize, @(x) median(x(:)));
%     labels = [labels {'Median'}];
%     dataIm(:,:,5) = nlfilter(subIm , winSize, @(x) kurtosis(x(:),0));
%     labels = [labels {'Kurtosis'}];
%     dataIm(:,:,6) = nlfilter(subIm , winSize, @(x) skewness(x(:),0));
%     labels = [labels {'Skewness'}];
%     
%     filtR = generateRadialFilterLBP(lbpFiltSize, lbpFiltRadius);
%     dataIm(:,:,7) = efficientLBP(subIm, 'filtR', filtR, 'isRotInv', true, 'isChanWiseRot', false);
%     labels = [labels {'Lbp'}];
% %     uniqueRotInvLBP = findUniqValsRILBP(lbpFiltSize);  
% %     tightValsRILBP = 1:length(uniqueRotInvLBP);
% %     effTightRILBP = tightHistImg(dataIm(:,:,7), 'inMap', uniqueRotInvLBP, 'outMap', tightValsRILBP);
% %     hlbp = hist(single(effTightRILBP(:)), tightValsRILBP);
% 
% else
%     filtR = generateRadialFilterLBP(lbpFiltSize, lbpFiltRadius);
% %     dataIm(:,:,7) = efficientLBP(subIm, 'filtR', filtR, 'isRotInv', true, 'isChanWiseRot', false);
%     dataIm(:,:,3) = efficientLBP(subIm, 'filtR', filtR, 'isRotInv', true, 'isChanWiseRot', false);
%     labels = [labels {'Lbp'}];
% end
% internal.images.isFigureAvailable = tmp;
end

%%
function labels = CatLabels(labels, postFix)
% labels = char(labels);
    for i = 1:length(labels)
        labels{i} = [labels{i} postFix];
    end
end