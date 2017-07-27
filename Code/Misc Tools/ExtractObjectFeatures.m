function [data dataIm] = ExtractObjectFeatures(imRgbIr, shapeArray, spatialRef, varargin)
    wPcaSpekboom = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\PcaSpekBoom2.mat';
    wPcaAll = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\PcaAll2.mat';
    wPcaRgG = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\PcaRgG2.mat';
    winSize = [5 5];
    border = ceil(winSize(1)/2) - 1;
    lbpFiltSize = 4;
    lbpFiltRadius = 1;

    ModifyDefaultArgs(varargin);

    
    % if (nargin == 1)
%     mask = true(size(imRgbIr));
%     bgMask = false(size(imRgbIr));
% end

    %when you normalise by intensity, you really want to be normalising by the
    %real spectrum, whereas here, we are normalising by a weighted version (esp IR)
    % mask = mask & (~bushMask);


    % rgG = imRgbIr(:,:,1:3)./repmat(sum(imRgbIr(:,:,1:3), 3), [1 1 3]);
    data = [];
    classLabelsV = {};
    subClassLabelsV = {};
    featuresV = [];
    areaLabelsV = {};
    classAreaLabelsV = {};
    habitatLabelsV = {};

    figure;
    for i = 1:length(shapeArray)
        fprintf('.');        
        [shapeArray(i).Y, shapeArray(i).X] = map2pix(spatialRef, [shapeArray(i).X], [shapeArray(i).Y]);
        idx = find(isnan(shapeArray(i).X));
        idx = 1:idx(1)-1;
        shapeArray(i).X = shapeArray(i).X(idx);
        shapeArray(i).Y = shapeArray(i).Y(idx);

%         xIdx = floor((min(shapeArray(i).X)-border)):ceil((max(shapeArray(i).X)+border));
%         yIdx = floor((min(shapeArray(i).Y)-border)):ceil((max(shapeArray(i).Y)+border));
        xIdx = floor((min(shapeArray(i).X))):ceil((max(shapeArray(i).X)));
        yIdx = floor((min(shapeArray(i).Y))):ceil((max(shapeArray(i).Y)));
        subIm = double(imRgbIr(yIdx, xIdx, :));
        maskOrig =  poly2mask(shapeArray(i).X(idx), shapeArray(i).Y(idx), size(imRgbIr,1), size(imRgbIr,2));
        maskOrig = maskOrig(yIdx, xIdx, :);
        mask = bwmorph(maskOrig, 'erode', border);
        if (sum(mask(:)) < 9)
            mask = maskOrig;
        end

        subplot(1,3,1);
        imshow(maskOrig);
        subplot(1,3,2);
        imshow(mask);
        subplot(1,3,3);
        imshow(subIm(:,:,[4 1 2])/(2^12));
        
        [dataIm, labels] = ExtractFeaturesIm2(subIm, 'wPcaSpekboom', wPcaSpekboom, ...
            'wPcaAll', wPcaAll, 'wPcaRgG', wPcaRgG, 'winSize', winSize, ...
            'lbpFiltSize', lbpFiltSize, 'lbpFiltRadius', lbpFiltRadius);
        
        if (isempty(shapeArray(i).Class) || isempty(shapeArray(i).Habitat))
            fprintf('Skipping unlabelled blob\n');
            continue;
        end
        if (isempty(shapeArray(i).SubClass))
            shapeArray(i).SubClass = shapeArray(i).Class;
        end
        mask = repmat(mask, [1 1 size(dataIm, 3)]);
        features = reshape(dataIm(mask), [], size(dataIm, 3));

        %Find LBP histograms 
        lbpFeatIdx = [29 43 50];
        uniqueRotInvLBP = findUniqValsRILBP(lbpFiltSize);
        tightValsRILBP = 1:length(uniqueRotInvLBP);
        hlbp = [];
        for j = 1:length(lbpFeatIdx)
            effTightRILBP = tightHistImg(features(:, lbpFeatIdx(j)), 'inMap', uniqueRotInvLBP, 'outMap', tightValsRILBP); %CHECK INDICES
            h = hist(single(effTightRILBP(:)), tightValsRILBP);
            hlbp = [hlbp h./sum(h)];            
            for k = 1:length(h)
                labels{end+1} = sprintf('%sHist%d', labels{lbpFeatIdx(j)}, k);
            end
        end
        features = [mean(features, 1) hlbp];

        %Find GLCM features
        texImIdx = [6 9 15];
        for j = 1:length(lbpFeatIdx)
            tmpIm = dataIm(:,:,texImIdx(j));
            tmpIm(~maskOrig) = nan;

            glcms = graycomatrix(tmpIm, 'NumLevels', 20, 'GrayLimits', [0 0.5], 'Symmetric', true);
            stats = graycoprops(glcms, 'all'); %'Contrast', 'Correlation', 'Energy', 'Homogeneity'
            
            if (isnan(stats.Correlation))
                stats.Correlation = 0;
            end
            
            features(end+1) = stats.Contrast;
            labels{end+1} = sprintf('%s%s', labels{texImIdx(j)}, 'Contrast');
            features(end+1) = stats.Correlation;
            labels{end+1} = sprintf('%s%s', labels{texImIdx(j)}, 'Correlation');
            features(end+1) = stats.Energy;
            labels{end+1} = sprintf('%s%s', labels{texImIdx(j)}, 'Energy');
            features(end+1) = stats.Homogeneity;
            labels{end+1} = sprintf('%s%s', labels{texImIdx(j)}, 'Homogeneity');
        end

        if (any(isnan(features)))
            fprintf('Skipping NAN feature\n');
            labels(isnan(features))
            continue;
        end
        
        if (isempty(features))
            fprintf('Skipping empty blob\n');
            continue;
        end
        
        classLabelsV = [classLabelsV; cellstr(repmat(shapeArray(i).Class, size(features, 1), 1))];
        subClassLabelsV = [subClassLabelsV; cellstr(repmat(shapeArray(i).SubClass, size(features, 1), 1))];
        featuresV = [featuresV; features];
        areaLabelsV = [areaLabelsV; cellstr(repmat(shapeArray(i).Area, size(features, 1), 1))];
        classAreaLabel = sprintf('%s - %s', shapeArray(i).Area, shapeArray(i).Class);
        classAreaLabelsV = [classAreaLabelsV; cellstr(repmat(classAreaLabel, size(features, 1), 1))];
        havitatLabel = sprintf('%s', shapeArray(i).Habitat);
        habitatLabelsV = [habitatLabelsV; cellstr(repmat(shapeArray(i).Habitat, size(features, 1), 1))];
        %end
    end
    fprintf('\n');

    data = prdataset(featuresV);
    data = setfeatlab(data, labels);
    data = addlabels(data, char(classLabelsV), 'Default');
    data = addlabels(data, char(areaLabelsV), 'Area');
    data = addlabels(data, char(subClassLabelsV), 'SubClass');
    data = addlabels(data, char(classAreaLabelsV), 'ClassArea');
    data = addlabels(data, char(habitatLabelsV), 'Habitat');

