function VisualiseClassification(out, visIm, cirIm, R)

outIm = reshape(+out(:,2) - (+out(:,3) + +out(:,1)), size(visIm(:,:,1)));

outMask = reshape((out*nlabeld) == 2, size(visIm(:,:,1)));

% outMask = bwareaopen(outMask, 10);
CC = bwconncomp(outMask);
numPixels = cellfun(@numel,CC.PixelIdxList);
for i = 1:length(numPixels)
    if (numPixels(i) < 20)
        outMask(CC.PixelIdxList{i}) = false;
    end
end

B = bwboundaries(outMask);
figure;
% h1 = subplot(1,2,1);
% imshow(outIm);
% h2 = subplot(1,2,2);
imshow(visIm);
hold on;
for i = 1:length(B)
    boundary = B{i};
    plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2);
end
% linkaxes([h1 h2]);
if (nargin ==4)
    figure;
    h1 = subplot(1,2,1);
    mapshow(outIm, R);
%     xlabel('Lon')
%     ylabel('Lat')
    title('P. Afra Probability Map')
    h2 = subplot(1,2,2);
    mapshow(cirIm, R);
    title('Scene')
%     xlabel('Lon')
%     ylabel('Lat')
    hold on
    for i = 1:length(B)
        boundary = B{i};
        [lat, lon] = pix2latlon(R, boundary(:,1), boundary(:,2));
        mapshow(lon, lat, 'DisplayType' ,'line', 'LineWidth', 2, 'Color', [0 0 1])
%         plot(boundary(:,2)*R(2,1)+R(3,1), boundary(:,1)*R(2,2)+R(3,2), 'b', 'LineWidth', 2);
    end
    
else
    figure;
    h1 = subplot(1,2,1);
    imshow(outIm);
    h2 = subplot(1,2,2);
    imshow(cirIm);
    hold on
    for i = 1:length(B)
        boundary = B{i};
        plot(boundary(:,2), boundary(:,1), 'b', 'LineWidth', 2);
    end
end
% linkaxes([h1 h2]);
fontsize(12)
set(h1, 'FontSize', 10);
set(h2, 'FontSize', 10);

if (nargin ==4)
    figure
    h3 = subplot(1,2,1);
    mapshow(visIm, R);
    title('Scene')
%     xlabel('Lon')
%     ylabel('Lat')
    h4 = subplot(1,2,2);
    mapshow(outIm, R);
%     xlabel('Lon')
%     ylabel('Lat')
    title('P. Afra Probability Map')
else
    figure
    h3 = subplot(1,2,1);
    imshow(visIm);
    title('Scene')
    h4 = subplot(1,2,2);
    imshow(outIm);
    title('P. Afra Probability Map')
end
fontsize(12)
set(h3, 'FontSize', 10);
set(h4, 'FontSize', 10);
linkaxes([h1 h2 h3 h4]);
