function [Spectra IjLocs] = PlotMultibandImage(Image)

% figure('WindowButtonMotionFcn',@MouseCallback)
disp('"+", "-" cycles bands, "1".. "N" sets band, "m" shows mean, "d" deletes last pt, "c" clears all pts, "q" exits');

%remove borders
% m = median(Image(:));
% Image([1:5 end-4:end], :, :) = m;
% Image(:, [1:5 end-4:end], :) = m;
% Image = Image([6:end-5], :, :);
% Image = Image(:, [6:end-5], :);

dImage = double(Image - min(Image(:)));
dImage = dImage./max(dImage(:));
% loHi = prctile(dImage(1:100:end)', [1 99]);
loHi = stretchlim(dImage(1:100:end)', [0.05 0.95]);

for i = 1:size(Image, 3)
    dImage(:, :, i) = imadjust(dImage(:, :, i), loHi);
end

f1 = figure; 
% iptaddcallback(f1, 'WindowButtonMotionFcn', @MouseMotionCallback);
% iptaddcallback(f1, 'WindowButtonDownFcn', @MouseDownCallback);
% set(f1, 'pointer', 'fullcrosshair');
% return
if (size(Image, 1) > size(Image, 2))
    subplot(1, 2, 1);
else
    subplot(2, 1, 1);
end
imH = imagesc(mean(dImage, 3));
set(imH, 'EraseMode', 'none');
colormap('gray');
% colorbar;
axis image;
hold on;
title('Mean(Bands) Image')
imAxH = get(f1, 'CurrentAxes');

if (size(Image, 1) > size(Image, 2))
    subplot(1, 2, 2);
else
    subplot(2, 1, 2);
end

plot(squeeze(Image(1, 1, :)));
% loHi = prctile(Image(:)', [0.1 99.9]);
% axis([1 size(Image, 3) loHi(1) loHi(2)]);
title('Spectra')
xlabel('Band')
% get(f1, 'CurrentAxes');
hold on;
grid on;
specAxH = gca;
cla;
set(specAxH, 'PlotBoxAspectRatio', get(imAxH, 'PlotBoxAspectRatio'));

but = 0;
markerSymbols = {'o', 'x', '+', '*', 's', 'd', '^', 'v', '<', '>', 'p', 'h','o', 'x', '+', '*', 's', 'd', '^', 'v', '<', '>', 'p', 'h'};
markerColours = {'r', 'b', 'g', 'y', 'm', 'c', 'k'};
markers = {};
markerCount = 0;
markerH = [];
specH = [];
band = 0;
Spectra = [];
IjLocs = [];

for i = 1:length(markerSymbols)
    for j = 1:length(markerColours)
        markers = [markers, [markerSymbols{i} markerColours{j}]];
    end
end

while (but ~= 3)
    [x, y, but] = ginput(1);
    if (but == 1)
        markerCount = markerCount + 1;
        IjLocs(markerCount, :) = [round(y), round(x)];
        markerH(markerCount) = plot(imAxH, x, y, markers{markerCount});
        Spectra(markerCount, :) = Get4ConnSpectra(Image, [round(y), round(x)]);
%         Spectra(markerCount, :) = squeeze(Image(round(y), round(x), :))';
        specH(markerCount) = plot(specAxH, Spectra(markerCount, :),  [markers{markerCount} '-']);
    elseif (but >= uint8('1') && but <= uint8(num2str(min([size(Image, 3) 9]))))
        band = str2num(char(but));
        set(imH, 'CData', dImage(:, :, band));
        title(imAxH, sprintf('Band %d Image', band));
    elseif (but == uint8('+'))
        band = band + 1;
        if (band > size(Image, 3))
            band = size(Image, 3);
        end
        set(imH, 'CData', dImage(:, :, band));
        title(imAxH, sprintf('Band %d Image', band));
    elseif (but == uint8('-'))
        band = band - 1;
        if (band < 1)
            band = 1;
        end
        set(imH, 'CData', dImage(:, :, band));
        title(imAxH, sprintf('Band %d Image', band));
    elseif (but == uint8('m'))
        set(imH, 'CData', mean(dImage, 3));
        title(imAxH, 'Mean(Bands) Image');
        band = 0;
    elseif (but == uint8('c'))
        for i = 1:markerCount
            delete(markerH(i));
        end
        markerCount = 0;
        markerH = [];
        specH = [];
        Spectra = [];
        IjLocs = [];
        cla(specAxH);
    elseif (but == uint8('d'))
        delete(markerH(markerCount));
        delete(specH(markerCount));
        markerCount = markerCount - 1;
        Spectra = Spectra(1:markerCount, :);
        IjLocs = IjLocs(1:markerCount, :);
    elseif ((but == 3) || (but == uint8('q')))
        break;
    end
end
% 
% function MouseMotionCallback(Src, EventData, Image)
%     disp(Src);
%     disp(EventData);
% return
% function MouseDownCallback(Src, EventData, Image)
%     disp(Src);
%     disp(EventData);
% return
function Spectrum = Get4ConnSpectra(Image, CoOrd)
if 0
    CoOrd = repmat(round(CoOrd), 5, 1) + [0 1;0 -1;1 0;-1 0;0 0];
%     Spectrum = squeeze(Image(CoOrd(:, 1), CoOrd(:, 2), :))';
    for i=1:size(CoOrd, 1)
        Spectrum(i, :) = squeeze(Image(CoOrd(i, 1), CoOrd(i, 2), :))';
    end
    Spectrum = mean(Spectrum, 1);
else
    CoOrd = round(CoOrd);
    Spectrum = squeeze(Image(CoOrd(1), CoOrd(2), :))';
end
return