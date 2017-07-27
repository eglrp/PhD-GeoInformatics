function res = MakeGroundTruthImage(imageFileName, gtShapeFileName, varargin)
%% Classify image only in areas of ground truth to save time
    winSize = [5 5];
    ModifyDefaultArgs(varargin);

    shapeArray = shaperead(gtShapeFileName);
    gi = geotiffinfo(imageFileName);
    imRgbIr = imread(imageFileName);
    count = 1;
    figure;
    res = []; %NB, res is saved in clfrFileName

    stretchValsCir = stretchlim(imRgbIr(:, :, [4 1 2]), [0.02 0.98])*(2^16-1);
    stretchValsRgb = stretchlim(imRgbIr(:, :, [1 2 3]), [0 0.99])*(2^16-1);

    for i = 1:length(shapeArray)
        [shapeArray(i).Y, shapeArray(i).X] = map2pix(gi.SpatialRef, [shapeArray(i).X], [shapeArray(i).Y]);
        if (min(shapeArray(i).Y) > 0 && max(shapeArray(i).Y) <= size(imRgbIr,1) && min(shapeArray(i).X) > 0 && max(shapeArray(i).X) <= size(imRgbIr,2))
            
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
            maskBoundary = logical(imdilate(mask, strel('disk', 2)) - mask);

            rgbIm = subIm(:,:,1:3);
%             rgbIm(rgbIm > 1) = 1;
            cirIm = subIm(:,:,[4 1 2]);

            cirIm = imadjust(cirIm/5000, stretchValsCir/5000, [0 1]);
            rgbIm = imadjust(rgbIm/5000, stretchValsRgb/5000, [0 1]);

            color = [1 1 0];
            for k = 1:3
                band = rgbIm(:,:,k);
                band(maskBoundary) = color(k);
                rgbIm(:,:,k) = band;

                band = cirIm(:,:,k);
                band(maskBoundary) = color(k);
                cirIm(:,:,k) = band;
            end

%             subShapeX = round(shapeArray(i).X - xIdx(1) + 1);
%             subShapeY = round(shapeArray(i).Y  - yIdx(1) + 1);
%             subShapeX(isnan(subShapeX)) = [];
%             subShapeY(isnan(subShapeX)) = [];
%             for k = 1:length(subShapeX)
%                 rgbIm(subShapeY(k), subShapeX(k),:) = [1 1 0];
%                 cirIm(subShapeY(k), subShapeX(k),:) = [1 1 0];
%             end

            res(count).rgbIm = rgbIm;
            res(count).cirIm = cirIm;
            res(count).name = shapeArray(i).Name;
            res(count).coverMin = shapeArray(i).CoverMin;
            res(count).coverMax = shapeArray(i).CoverMax;
            count = count + 1;
            
            
%%
            if true
%                 figure;
                h1 = subplot(1,2,1);
                imshow(rgbIm);
                h2 = subplot(1,2,2);
                imshow(cirIm);
                linkaxes([h1 h2], 'xy')
                drawnow;
            end
        end
    end

end

