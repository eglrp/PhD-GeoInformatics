function [combinedImRgbIr cirIm R] = CombineRgbCir(rgbFilename, cirFilename)
    rgbIm = im2double(imread(rgbFilename));
    cirIm = im2double(imread(cirFilename)); %IR,R,G

    combinedImRgbIr = cat(3, rgbIm, cirIm(:, :, 1));
    scale = 1.5; %mean(mean(rgbIm(:,:,1)))/mean(mean(cirIm(:,:,2)));

    if (nargout == 3)
        R = worldfileread(getworldfilename(rgbFilename));
    end
    %note that the images have been scaled differently so just
    %substituting bands doesn't produce a realistic image. bands seem to
    %have been scaled individually
%     combinedImR2g2bIr = scale*cat(3, cirIm(:, :, 2:3), combinedImRgbIr(:,:,3)/scale);

%     ndvi = (combinedImRgbIr(:,:,4)*scale - combinedImRgbIr(:,:,1))./(combinedImRgbIr(:,:,1) + ...
%         combinedImRgbIr(:,:,4)*scale);
% 
%     rgGIR = combinedImRgbIr./repmat(sum(combinedImRgbIr, 3), [1 1 4]);
%     rgG = combinedImRgbIr(:, :, 1:3)./repmat(sum(combinedImRgbIr(:, :, 1:3), 3), [1 1 3]);
end

