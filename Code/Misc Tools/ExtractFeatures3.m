function [data dataIm] = ExtractFeatures3(imRgbIr, shapeArray, spatialRef, varargin)
    wPcaSpekboom = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\PcaSpekBoom.mat';
    wPcaAll = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\PcaAll.mat';
    wPcaRgG = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\PcaRgG.mat';
    winSize = [5 5];
    border = ceil(winSize(1)/2) - 1;

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
    slopeLabelsV = {};

    figure;
    for i = 1:length(shapeArray)
        fprintf('.');        
        [shapeArray(i).Y, shapeArray(i).X] = map2pix(spatialRef, [shapeArray(i).X], [shapeArray(i).Y]);
        idx = find(isnan(shapeArray(i).X));
        idx = 1:idx(1)-1;
        shapeArray(i).X = shapeArray(i).X(idx);
        shapeArray(i).Y = shapeArray(i).Y(idx);
        
        offset = ceil(winSize(1)/2);  %add extra to allow window to extend to boundary

        xIdx = floor((min(shapeArray(i).X)))-offset:ceil((max(shapeArray(i).X)))+offset;
        yIdx = floor((min(shapeArray(i).Y)))-offset:ceil((max(shapeArray(i).Y)))+offset;
        subIm = double(imRgbIr(yIdx, xIdx, :));
        maskOrig =  poly2mask(shapeArray(i).X(idx), shapeArray(i).Y(idx), size(imRgbIr,1), size(imRgbIr,2));
        maskOrig = maskOrig(yIdx, xIdx, :);

        if (border > 0)
            mask = bwmorph(maskOrig, 'erode', border); %avoid edges of objects for texture features
        else
            mask = maskOrig;
        end

        subplot(1,2,1);
        imshow(mask);
        subplot(1,2,2);
        imshow(subIm(:,:,[4 1 2])/(2^12));

        [dataIm, labels] = ExtractFeaturesIm2(subIm, 'wPcaSpekboom', wPcaSpekboom, 'wPcaAll', wPcaAll, 'wPcaRgG', wPcaRgG, 'winSize', winSize);
        %for i = 1:max(combinedMask(:))
    %         mask  = combinedMask == i;
            if (isempty(shapeArray(i).Class) || isempty(shapeArray(i).Habitat) || (shapeArray(i).Slope ~= 'N' && shapeArray(i).Slope ~= 'S'))
                fprintf('Skipping unlabelled blob\n');
                continue;
            end
            if (isempty(shapeArray(i).SubClass))
                shapeArray(i).SubClass = shapeArray(i).Class;
            end
            mask = repmat(mask, [1 1 size(dataIm, 3)]);
            features = reshape(dataIm(mask), [], size(dataIm, 3));
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
            slopeLabel = sprintf('%s', shapeArray(i).Slope);
            slopeLabelsV = [slopeLabelsV; cellstr(repmat(shapeArray(i).Slope, size(features, 1), 1))];
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
    data = addlabels(data, char(slopeLabelsV), 'Slope');

