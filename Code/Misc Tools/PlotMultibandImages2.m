function [Spectra IjLocs] = PlotMultibandImages(Images, IjLocs)

% figure('WindowButtonMotionFcn',@MouseCallback)
disp('"+", "-" cycles bands, "1".. "N" sets band, "m" shows mean, "d" deletes last pt, "c" clears all pts, "q" exits');

%remove borders
% m = median(Image(:));
% Image([1:5 end-4:end], :, :) = m;
% Image(:, [1:5 end-4:end], :) = m;
%setup limits based on first image i.e. scale all ims the same
for i = 1:length(Images)
%     Images{i} = Images{i}([6:end-5], :, :);
%     Images{i} = Images{i}(:, [6:end-5], :);
    dImages{i} = double(Images{i} - min(Images{i}(:)));
    dImages{i} = dImages{i}./max(dImages{i}(:));
%     loHi = prctile(dImages{i}(1:100:end)', [1 99]);

    for j = 1:size(Images{1}, 3)
        band = dImages{i}(:, :, j);
        loHi(j,:) = stretchlim(band(1:100:end)', [0.03 0.97]);
        dImages{i}(:, :, j) = imadjust(dImages{i}(:, :, j), loHi(j, :));
    end
end

f1 = figure; 
% iptaddcallback(f1, 'WindowButtonMotionFcn', @MouseMotionCallback);
% iptaddcallback(f1, 'WindowButtonDownFcn', @MouseDownCallback);
% set(f1, 'pointer', 'fullcrosshair');
% return
% if (size(Image, 1) > size(Image, 2))
%     subplot(1, 1+length, 1);
% else
%     subplot(2, 1, 1);
% end
for i = 1:length(Images)
    ax(i) = subplot(2, length(Images), i);
    imH{i} = imagesc(mean(dImages{i}, 3));
    set(imH{i}, 'EraseMode', 'none');
    colormap('gray');
    % colorbar;
    axis image;
    hold on;
    title('Mean(Bands) Image')    
    imAxH{i} = get(f1, 'CurrentAxes');

    subplot(2, length(Images), length(Images)+i);
    
    plot(squeeze(Images{i}(1, 1, :)));
    % loHi = prctile(Image(:)', [0.1 99.9]);
    % axis([1 size(Image, 3) loHi(1) loHi(2)]);
    title('Spectra')
    xlabel('Band')
    % get(f1, 'CurrentAxes');
    hold on;
    grid on;
    specAxH{i} = gca;
    cla;
    set(specAxH{i}, 'PlotBoxAspectRatio', get(imAxH{i}, 'PlotBoxAspectRatio'));
end
linkaxes(ax, 'xy')
but = 0;
markerSymbols = {'o', 'x', '+', '*', 's', 'd', '^', 'v', '<', '>', 'p', 'h','o', 'x', '+', '*', 's', 'd', '^', 'v', '<', '>', 'p', 'h'};
markerColours = {'r', 'b', 'g', 'y', 'm', 'c', 'k'};
markers = {};
markerCount = 0;
markerH = {[]};
specH = {[]};
band = 0;
Spectra = {[]};
%IjLocs = [];

for i = 1:length(markerSymbols)
    for j = 1:length(markerColours)
        markers = [markers, [markerSymbols{i} markerColours{j}]];
    end
end

if nargin > 1
    for i = 1:size(IjLocs, 1)
        markerCount = markerCount + 1;
        %IjLocs(markerCount, :) = [round(y), round(x)];
        x = IjLocs(markerCount, 2);
        y = IjLocs(markerCount, 1);
        for i = 1:length(Images)
            markerH{i}(markerCount) = plot(imAxH{i}, x, y, markers{markerCount});
            Spectra{i}(markerCount, :) = Get4ConnSpectra(Images{i}, [round(y), round(x)]);
            specH{i}(markerCount) = plot(specAxH{i}, Spectra{i}(markerCount, :),  [markers{markerCount} '-']);
        end
        % hack for msc stuff
        if length(Images)==2
            d = abs(Spectra{1}-Spectra{2});
            error = mean(d(:));
            title(specAxH{i}, sprintf('Spectra - MAD: %2.3f', error))
        end
    end
end

zoomPanOn = false;
while (but ~= 3)
    if ~zoomPanOn
        [x, y, but] = ginput(1);
    else
        k = waitforbuttonpress;
        if k ==1
            but = 'z';
        end
    end
    if (but == 1)
        markerCount = markerCount + 1;
        IjLocs(markerCount, :) = [round(y), round(x)];
        for i = 1:length(Images)
            markerH{i}(markerCount) = plot(imAxH{i}, x, y, markers{markerCount});
            Spectra{i}(markerCount, :) = Get4ConnSpectra(Images{i}, [round(y), round(x)]);
            specH{i}(markerCount) = plot(specAxH{i}, Spectra{i}(markerCount, :),  [markers{markerCount} '-']);
        end
        % hack for msc stuff
        if length(Images)==2
            d = abs(Spectra{1}-Spectra{2});
            error = mean(d(:));
            title(specAxH{i}, sprintf('Spectra - MAD: %2.3f', error))
        end
%         Spectra(markerCount, :) = squeeze(Image(round(y), round(x), :))';
    elseif (but >= uint8('1') && but <= uint8(num2str(min([size(Images{1}, 3) 9]))))
        band = str2num(char(but));
        for i = 1:length(Images)
            set(imH{i}, 'CData', dImages{i}(:, :, band));
            title(imAxH{i}, sprintf('Band %d Image', band));
        end
    elseif (but == uint8('+'))
        band = band + 1;
        if (band > size(Images{1}, 3))
            band = size(Images{1}, 3);
        end
        for i = 1:length(Images)
            set(imH{i}, 'CData', dImages{i}(:, :, band));
            title(imAxH{i}, sprintf('Band %d Image', band));
        end
    elseif (but == uint8('-'))
        band = band - 1;
        if (band < 1)
            band = 1;
        end
        for i = 1:length(Images)
            set(imH{i}, 'CData', dImages{i}(:, :, band));
            title(imAxH{i}, sprintf('Band %d Image', band));
        end
    elseif (but == uint8('m'))
        for i = 1:length(Images)
            set(imH{i}, 'CData', mean(dImages{i}, 3));
            title(imAxH{i}, 'Mean(Bands) Image');
        end
        band = 0;
    elseif (but == uint8('c'))
        for i = 1:markerCount
            for j = 1:length(markerH)
                delete(markerH{j}(i));
            end
        end
        markerCount = 0;
        markerH = {[]};
        specH = {[]};
        Spectra = {[]};
        IjLocs = [];
        for i = 1:length(specAxH)
            cla(specAxH{i});
        end
    elseif (but == uint8('d'))
        for j = 1:length(markerH)
            delete(markerH{j}(markerCount));
            delete(specH{j}(markerCount));
        end
        markerCount = markerCount - 1;
        for j = 1:length(Spectra)
            Spectra{j} = Spectra{j}(1:markerCount, :);
        end
        IjLocs = IjLocs(1:markerCount, :);
    elseif (but == uint8('x'))
        for i = 1:length(Images)
            set(imH{i}, 'CData', dImages{i}(:,:,1:3));
            title(imAxH{i}, 'Bands 1-3');
        end
        band = 0;
    elseif (but == uint8('z'))
        zoomPanOn = ~zoomPanOn;
        but = ' ';
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
    CoOrd = repmat(round(CoOrd), 5, 1) + [0 1;0 -1;1 0;-1 0;0 0];
%     Spectrum = squeeze(Image(CoOrd(:, 1), CoOrd(:, 2), :))';
    for i=1:size(CoOrd, 1)
        Spectrum(i, :) = squeeze(Image(CoOrd(i, 1), CoOrd(i, 2), :))';
    end
    Spectrum = mean(Spectrum, 1);
return