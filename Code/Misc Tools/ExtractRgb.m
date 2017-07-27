function [data dataIm] = ExtractRgb(imRgbIr, shapeArray, spatialRef)


% rgG = imRgbIr(:,:,1:3)./repmat(sum(imRgbIr(:,:,1:3), 3), [1 1 3]);
% border = 5;
border = 0; 
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
    subIm = double(imRgbIr(yIdx, xIdx, :))/(2^14);
    mask =  poly2mask(shapeArray(i).X(idx), shapeArray(i).Y(idx), size(imRgbIr,1), size(imRgbIr,2));
    mask = mask(yIdx, xIdx, :);

%     if isstruct(subIm)
%         subIm = subIm.data;
%     end
    dataIm = subIm;
    labels = {'R','G','B','NIR'};

%     dataIm = zeros(size(subIm, 1), size(subIm, 2), 14); %NB feat len must be right here
%     dataIm(:, :, 1:4) = subIm;
    
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

data = prdataset(featuresV);
data = setfeatlab(data, labels);
data = addlabels(data, char(classLabelsV), 'Default');
data = addlabels(data, char(areaLabelsV), 'Area');
data = addlabels(data, char(subClassLabelsV), 'SubClass');
data = addlabels(data, char(classAreaLabelsV), 'ClassArea');
data = addlabels(data, char(habitatLabelsV), 'Habitat');

