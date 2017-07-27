function res = ValidateClfrAgainstJvGroundTruth(imageFileName, gtShapeFileName, clfrFileName, varargin)
%% Classify image only in areas of ground truth to save time
winSize = [5 5];
border = 0;
ModifyDefaultArgs(varargin);
global w feats
if (~isempty(clfrFileName))
    prload(clfrFileName); %prload for new mappings, load for old (eg dd tools)
end
%     w = prmapping(w);

%     out = imread(outImageFileNames{g});
%     [dum outC] = max(+out, [], 3);

%     fprintf('Processing %s\n', imageFileName);
shapeArray = shaperead(gtShapeFileName);
gi = geotiffinfo(imageFileName);
imRgbIr = imread(imageFileName);
count = 1;
figure;
res = []; %NB, res is saved in clfrFileName
for i = 1:length(shapeArray)
    [shapeArray(i).Y, shapeArray(i).X] = map2pix(gi.SpatialRef, [shapeArray(i).X], [shapeArray(i).Y]);
    if (min(shapeArray(i).Y) > 0 && max(shapeArray(i).Y) <= size(imRgbIr,1) && min(shapeArray(i).X) > 0 && max(shapeArray(i).X) <= size(imRgbIr,2))
        
        %% Extract features for this region (Taken from ExtractFeatures3...)
        %             fprintf('.');
        %             [shapeArray(i).Y, shapeArray(i).X] = map2pix(gi.SpatialRef, [shapeArray(i).X], [shapeArray(i).Y]);
        %             idx = find(isnan(shapeArray(i).X));
        %             idx = 1:idx(1)-1;
        %             shapeArray(i).X = shapeArray(i).X(idx);
        %             shapeArray(i).Y = shapeArray(i).Y(idx);
        
        offset = ceil(winSize(1)/2);  %add extra to allow window to extend to boundary
        
        xIdx = floor((min(shapeArray(i).X)))-offset:ceil((max(shapeArray(i).X)))+offset;
        yIdx = floor((min(shapeArray(i).Y)))-offset:ceil((max(shapeArray(i).Y)))+offset;
        subIm = double(imRgbIr(yIdx, xIdx, :));
        
        maskOrig = false(size(imRgbIr,1), size(imRgbIr,2));
        startIdx = 1;
        idx = find(isnan(shapeArray(i).X));
        for j = 1:length(idx)
            polyIdx = startIdx:idx(j)-1;
            maskOrig = maskOrig | poly2mask(shapeArray(i).X(polyIdx), shapeArray(i).Y(polyIdx), size(imRgbIr,1), size(imRgbIr,2));
            startIdx = idx(j)+1;
        end
        %             maskOrig =  poly2mask(shapeArray(i).X(idx), shapeArray(i).Y(idx), size(imRgbIr,1), size(imRgbIr,2));
        mask = maskOrig(yIdx, xIdx, :);
        
        
        %             subplot(1,2,1);
        %             imshow(mask);
        %             subplot(1,2,2);
        %             imshow(subIm(:,:,[4 1 2])/(2^12));
        
        if (true)
            [featIm, labels] = ExtractFeaturesIm2(subIm, 'feats', feats);
            
            %% Classify region
            
            out = im2feat(featIm(:,:,feats))*w;
            %     outIm = reshape(+out(:,2), [size(imstruct.data, 1) size(imstruct.data, 2)]);
            outIm = im2uint8(reshape(+(out*classc), [size(subIm, 1) size(subIm, 2) size(w, 2)]));
            if (size(outIm, 3) == 2)
                outIm(:, :, 3) = 0;
            end
            outC = out*nlabeld;
        else
            outIm = blockproc(subIm, [128 128], @(x)ClassifyIm(x, w, feats));
            %                     , 'Destination', ...
            %                         fn(i).outImageFileName, 'UseParallel', true);
            [dum outC] = max(+outIm, [], 3);
        end
        outC = outC==2; %NB check index corresponds to SB
        outC = reshape(outC, size(outIm, 1), size(outIm, 2));
        se = strel('disk', 1);
        outC = imopen(outC,se);
        outC = imclose(outC,se);
        %             Mat strEl = getStructuringElement(MORPH_ELLIPSE, Size(3, 3), Point(-1, -1));
        % 			cv::morphologyEx(outIm, outIm, CV_MOP_OPEN, strEl);
        % 			cv::morphologyEx(outIm, outIm, CV_MOP_CLOSE, strEl);
        
        
        %%
        %             rowIdx = floor(min(shapeArray(i).Y)):ceil(max(shapeArray(i).Y));
        %             colIdx = floor(min(shapeArray(i).X)):ceil(max(shapeArray(i).X));
        
        %         idx = 1:idx(1)-1;
        fprintf('%s %d-%d, %f\n', shapeArray(i).Name, shapeArray(i).CoverMin, shapeArray(i).CoverMax, 100*sum(outC(mask))/sum(mask(:)));
        
        res(count).name = shapeArray(i).Name;
        res(count).coverMin = shapeArray(i).CoverMin;
        res(count).coverMax = shapeArray(i).CoverMax;
        res(count).gtCover = mean([shapeArray(i).CoverMin shapeArray(i).CoverMax]);
        res(count).clfCover = 100*sum(outC(mask))/sum(mask(:));
        count = count + 1;
        
        if true
            %                 figure;
            
            h1 = subplot(2,2,1);
            imshow(mask);
            h2 = subplot(2,2,2);
            imshow(outIm);
            h3 = subplot(2,2,3);
            imshow(subIm(:,:, [1 2 3])/2^13);
            h4 = subplot(2,2,4);
            imshow(subIm(:,:, [4 1 2])/2^13);
            linkaxes([h1 h2 h3 h4], 'xy')
        end
    end
end

end

