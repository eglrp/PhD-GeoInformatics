function [data dataIm] = ExtractFeatures2(imRgbIr, shapeArray, spatialRef)

% if (nargin == 1)
%     mask = true(size(imRgbIr));
%     bgMask = false(size(imRgbIr));
% end

%when you normalise by intensity, you really want to be normalising by the
%real spectrum, whereas here, we are normalising by a weighted version (esp IR)
% mask = mask & (~bushMask);


% rgG = imRgbIr(:,:,1:3)./repmat(sum(imRgbIr(:,:,1:3), 3), [1 1 3]);
border = 5;
data = [];
classLabelsV = {};
subClassLabelsV = {};
featuresV = [];
areaLabelsV = {};
classAreaLabelsV = {};
habitatLabelsV = {};

figure;
for i = 1:length(shapeArray)
    [shapeArray(i).Y, shapeArray(i).X] = map2pix(spatialRef, [shapeArray(i).X], [shapeArray(i).Y]);
    idx = find(isnan(shapeArray(i).X));
    idx = 1:idx(1)-1;
    shapeArray(i).X = shapeArray(i).X(idx);
    shapeArray(i).Y = shapeArray(i).Y(idx);

    xIdx = floor((min(shapeArray(i).X)-border)):ceil((max(shapeArray(i).X)+border));
    yIdx = floor((min(shapeArray(i).Y)-border)):ceil((max(shapeArray(i).Y)+border));
    subIm = double(imRgbIr(yIdx, xIdx, :));
    mask =  poly2mask(shapeArray(i).X(idx), shapeArray(i).Y(idx), size(imRgbIr,1), size(imRgbIr,2));
    mask = mask(yIdx, xIdx, :);

    subplot(1,2,1);
    imshow(mask);
    subplot(1,2,2);
    imshow(subIm(:,:,[4 1 2])/(2^12));

    [dataIm, labels] = ExtractFeaturesIm(subIm);
    %for i = 1:max(combinedMask(:))
%         mask  = combinedMask == i;
        if (isempty(shapeArray(i).Class) || isempty(shapeArray(i).Habitat))
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
    %end
    fprintf('.');
end
fprintf('\n');

data = dataset(featuresV);
data = setfeatlab(data, labels);
data = addlabels(data, char(classLabelsV), 'Default');
data = addlabels(data, char(areaLabelsV), 'Area');
data = addlabels(data, char(subClassLabelsV), 'SubClass');
data = addlabels(data, char(classAreaLabelsV), 'ClassArea');
data = addlabels(data, char(habitatLabelsV), 'Habitat');

