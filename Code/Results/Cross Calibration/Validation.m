%%
clear all;
close all; 
spotImFilenames = {'D:\Data\Development\Projects\PhD GeoInformatics\Code\Results\Cross Calibration\AtcorSpotCIR1.png';...
    'D:\Data\Development\Projects\PhD GeoInformatics\Code\Results\Cross Calibration\AtcorSpotCIR2.png';...
    'D:\Data\Development\Projects\PhD GeoInformatics\Code\Results\Cross Calibration\AtcorSpotCIR3.png';...    
    };

xcalibImFilenames = {'D:\Data\Development\Projects\PhD GeoInformatics\Code\Results\Cross Calibration\XCalib2CIR1.png';...
    'D:\Data\Development\Projects\PhD GeoInformatics\Code\Results\Cross Calibration\XCalib2CIR2.png';...
    'D:\Data\Development\Projects\PhD GeoInformatics\Code\Results\Cross Calibration\XCalib2CIR3.png';...    
    };

modisImFilenames = {'D:\Data\Development\Projects\PhD GeoInformatics\Code\Results\Cross Calibration\ModisCIR1.png';...
    'D:\Data\Development\Projects\PhD GeoInformatics\Code\Results\Cross Calibration\ModisCIR2.png';...
    'D:\Data\Development\Projects\PhD GeoInformatics\Code\Results\Cross Calibration\ModisCIR3.png';...    
    };

whiteBal = [];
scaleFactor = 50/255;
cd 'D:\Data\Development\Projects\PhD GeoInformatics\Code\Results\Cross Calibration\';

for i = 1:length(spotImFilenames)
    
    spotIm = (imread(spotImFilenames{i}));
    xcalibIm = (imread(xcalibImFilenames{i}));
    modisIm = (imread(modisImFilenames{i}));
    R = worldfileread(getworldfilename(spotImFilenames{i}));
    
    [p f] = fileparts(xcalibImFilenames{i});
    if false
        imFigFilename = [p '\' f 'ValidColBalImFig.png'];
        densityFigFilename = [p '\' f 'ValidColBalDensityFig.png'];
        %only find white bal for 1st im, then apply that to others
        [spotIm whiteBal] = ColourBalImage(xcalibIm, spotIm, 'whiteBal', whiteBal);
    else
        imFigFilename = [p '\' f 'ValidImFig.png'];
        densityFigFilename = [p '\' f 'ValidDensityFig.png'];
    end

    figure;
    h1 = subplot(2,2,1);
%     imshow(spotIm);
    mapshow(spotIm, R);
%     scaleruler on
    axis tight
    title('SPOT')
    h2 = subplot(2,2,2);
%     imshow(xcalibIm);
    mapshow(xcalibIm, R);
    axis tight
    title('Cross Calibrated NGI')
    h3 = subplot(2,2,3);
%     imshow(modisIm);
    mapshow(modisIm, R);
    axis tight
    title('MODIS')
    h4 = subplot(2,2,4);
    diffIm = abs(double(xcalibIm) - double(spotIm))*scaleFactor; %/double(max([spotIm(:); xcalibIm(:)]));
%    imshow(log10(diffIm)/log10(50));
%    mapshow(log10(diffIm)/log10(50), R);
    mapshow(uint8(abs(double(xcalibIm) - double(spotIm))), R);
    axis tight
    linkaxes([h1 h2 h3 h4], 'xy')
    lbl = sprintf('Error=|SPOT-XCalib| (mean(Error): %.2f%%, Std(Error): %.2f%%)', nanmean(diffIm(:)), nanstd(diffIm(:)));
%     lbl = sprintf('Error');
%     fprintf(lbl)
    title(lbl)
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 30 20])
    drawnow;

    print('-dpng', imFigFilename, '-r150');
    
    icons = {'m','r','g', 'k'};
    titles = {'IR Density', 'Red Density', 'Green Density', 'Error Density'};
    figure;
    xmesh = linspace(0,255,50);
    for j = 1:size(spotIm, 3)
        spotBand = spotIm(:,:,j);
%         spotDensity(j, :) = hist(spotBand(:), xmesh);
        [bandwidth, spotDensity, spotXmesh, cdf] = kde(double(spotBand(:)), 100, 0, 255);
        xcalibBand = xcalibIm(:,:,j);
        [bandwidth, xcalibDensity, xcalibXmesh, cdf] = kde(double(xcalibBand(:)), 100, 0, 255);
        modisBand = modisIm(:,:,j);
        [bandwidth, modisDensity, modisXmesh, cdf] = kde(double(modisBand(:)), 100, 0, 255);
        diffBand = diffIm(:,:,j);
        [bandwidth, diffDensity, diffXmesh, cdf] = kde(double(diffBand(:)), 100, 0, 255);
        subplot(1,3,j)
        h(j) = plot(spotXmesh*scaleFactor, spotDensity);
        hold all;
        plot(xcalibXmesh*scaleFactor, xcalibDensity);
        hold all;
        plot(modisXmesh*scaleFactor, modisDensity);
        hold all;
%         plot(diffXmesh*scaleFactor, diffDensity, [icons{j} '-+']);
%         hold on;
        title(titles{j});
        legend({'SPOT', 'Cross Calib', 'MODIS'});
        grid on;
        axis tight;
        axis square;
        xlabel('Reflectance (%)')
        ylabel('Density')
    end
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 30 20])
    drawnow;

    print('-dpng', densityFigFilename, '-r150');
    %saveas(gcf, densityFigFilename, 'png');

end

delete('ValidationColBal.zip');
delete('Validation.zip');
zip('ValidationColBal.zip', {'*ValidColBalImFig.png', '*ValidColBalDensityFig.png'});
zip('Validation.zip', {'*ValidImFig.png', '*ValidDensityFig.png'});



%%
% Same as above but histograms instead of kde
clear all;
%% 
% 
close all; 

spotImFilenames = {'D:\Data\Development\Projects\PhD GeoInformatics\Code\Results\Cross Calibration\AtcorSpotCIR1.png';...
    'D:\Data\Development\Projects\PhD GeoInformatics\Code\Results\Cross Calibration\AtcorSpotCIR2.png';...
    'D:\Data\Development\Projects\PhD GeoInformatics\Code\Results\Cross Calibration\AtcorSpotCIR3.png';...    
    };

xcalibImFilenames = {'D:\Data\Development\Projects\PhD GeoInformatics\Code\Results\Cross Calibration\XCalib2CIR1.png';...
    'D:\Data\Development\Projects\PhD GeoInformatics\Code\Results\Cross Calibration\XCalib2CIR2.png';...
    'D:\Data\Development\Projects\PhD GeoInformatics\Code\Results\Cross Calibration\XCalib2CIR3.png';...    
    };

modisImFilenames = {'D:\Data\Development\Projects\PhD GeoInformatics\Code\Results\Cross Calibration\ModisCIR1.png';...
    'D:\Data\Development\Projects\PhD GeoInformatics\Code\Results\Cross Calibration\ModisCIR2.png';...
    'D:\Data\Development\Projects\PhD GeoInformatics\Code\Results\Cross Calibration\ModisCIR3.png';...    
    };
whiteBal = [];
scaleFactor = 50/255;
cd 'D:\Data\Development\Projects\PhD GeoInformatics\Code\Results';


for i = 1:length(spotImFilenames)
    
    spotIm = (imread(spotImFilenames{i}));
    xcalibIm = (imread(xcalibImFilenames{i}));
    modisIm = (imread(modisImFilenames{i}));
    
    [p f] = fileparts(xcalibImFilenames{i});
    if true
        imFigFilename = [p '\' f 'ValidColBalImFig.png'];
        densityFigFilename = [p '\' f 'ValidColBalDensityFig.png'];
        %only find white bal for 1st im, then apply that to others
        [spotIm whiteBal] = ColourBalImage(xcalibIm, spotIm, 'whiteBal', whiteBal);
    else
        imFigFilename = [p '\' f 'ValidImFig.png'];
        densityFigFilename = [p '\' f 'ValidDensityFig.png'];
    end

    figure;
    subplot(2,2,1)
    imshow(spotIm);
    title('SPOT')
    subplot(2,2,2)
    imshow(xcalibIm);
    title('Cross Calib')
    subplot(2,2,3)
    imshow(modisIm);
    title('MODIS')
    subplot(2,2,4)
    diffIm = abs(double(xcalibIm) - double(spotIm))*scaleFactor; %/double(max([spotIm(:); xcalibIm(:)]));
    imshow(log10(diffIm)/log10(50));
    lbl = sprintf('Error=|SPOT-XCalib| (mean(Error): %.2f%%, Std(Error): %.2f%%)', nanmean(diffIm(:)), nanstd(diffIm(:)));
%     fprintf(lbl)
    title(lbl)
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    saveas(gcf, imFigFilename, 'png');
    
    icons = {'m','r','g', 'k'};
    titles = {'IR Density', 'Red Density', 'Green Density', 'Error Density'};
    figure;
    xmesh = linspace(0,255,50);
    for j = 1:size(spotIm, 3)
        spotBand = spotIm(:,:,j);
%         spotDensity(j, :) = hist(spotBand(:), xmesh);
        binCentres = linspace(0.5, 254.5, 50);
        [spotDensity, spotXmesh] = hist(double(spotBand(:)), binCentres);
        spotDensity = spotDensity./(numel(spotBand).*(binCentres(2)-binCentres(1)));
        
        xcalibBand = xcalibIm(:,:,j);
        [xcalibDensity, xcalibXmesh] = hist(double(xcalibBand(:)), binCentres);
        xcalibDensity = xcalibDensity./(numel(xcalibBand).*(binCentres(2)-binCentres(1)));
        
        binCentres = linspace(0.5, 254.5, 20);
        modisBand = modisIm(:,:,j);
        [modisDensity, modisXmesh] = hist(double(modisBand(:)), binCentres);
        modisDensity = modisDensity./(numel(modisBand).*(binCentres(2)-binCentres(1)));
        
        subplot(1,3,j)
        h(j) = plot(spotXmesh*scaleFactor, spotDensity, [icons{j} 'o-' ]);
        hold on;
        plot(xcalibXmesh*scaleFactor, xcalibDensity, [icons{j} 'x-']);
        hold on;
        plot(modisXmesh*scaleFactor, modisDensity, [icons{j} '^-']);
        hold on;
%         plot(diffXmesh*scaleFactor, diffDensity, [icons{j} '-+']);
%         hold on;
        title(titles{j});
        legend({'SPOT', 'Cross Calib', 'MODIS'});
        grid on;
        axis tight;
        axis square;
        xlabel('Reflectance (%)')
        ylabel('Density')
    end
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    saveas(gcf, densityFigFilename, 'png');

end


%% Round 2: Full georef images
close all; clear all;
spotFileName = 'D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\My Rectified\AATCORCorrected_oS131022114824832b_278400201_Uint16.tif';
ngiFileName = 'D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\My Rectified\StudyAreaXCalibMosaicSpotExtent.tif';

spotIm = int16(imread(spotFileName));
ngiIm = int16(imread(ngiFileName));
spotIm = spotIm(:,:,1:3);

%%
errorIm = (spotIm - ngiIm);
% errorIm(spotIm==0)=nan;
% imwrite(errorIm, 'G:\PhD GeoInformatics\Data\NGI\My Rectified\errorIm.tif');

fprintf('\tRMS\n');
fprintf('\tBand\tMean\tStd\tMedian\tMad\n');
fprintf('-------------------------------\n');

for i = 1:3
    errorBand = errorIm(:,:,i);
    errorBand(spotIm(:,:,i)==0) = [];
    errorBand = (single(errorBand)/50).^2;    
    fprintf('\t%d\t%2.4f\t%2.4f\t%2.4f\t%2.4f\n', i, sqrt(mean(errorBand(:))), sqrt(std(errorBand(:))), sqrt(median(errorBand(:))), sqrt(mad(errorBand(:), 1)));
end

errorIm_ = (single(errorIm)./50).^2;
errorIm_(spotIm==0) = [];

fprintf('\tAll\t%2.4f\t%2.4f\t%2.4f\t%2.4f\n', sqrt(mean(errorIm_(:))), sqrt(std(errorIm_(:))), sqrt(median(errorIm_(:))), sqrt(mad(errorIm_(:), 1)));
fprintf('-------------------------------\n');

fprintf('\n\n\tABS\n');
fprintf('\tBand\tMean\tStd\tMedian\tMad\n');
fprintf('-------------------------------\n');

for i = 1:3
    errorBand = errorIm(:,:,i);
    errorBand(spotIm(:,:,i)==0) = [];
    errorBand = abs(single(errorBand))/50;    
    fprintf('\t%d\t%2.4f\t%2.4f\t%2.4f\t%2.4f\n', i, mean(errorBand(:)), std(errorBand(:)), median(errorBand(:)), mad(errorBand(:), 1));
end

errorIm_ = abs(single(errorIm))./50;
errorIm_(spotIm==0) = [];

fprintf('\tAll\t%2.4f\t%2.4f\t%2.4f\t%2.4f\n', mean(errorIm_(:)), std(errorIm_(:)), median(errorIm_(:)), mad(errorIm_(:), 1));
fprintf('-------------------------------\n');

% figure;
% subplot(1, 3, 1)
% imshow(spotIm/5000)
% subplot(1, 3, 2)
% imshow(ngiIm/5000)
% subplot(1, 3, 3)
% imshow(errorIm/5000)

%% Look at relationship
close all; clear all;
spotFileName = 'G:\PhD GeoInformatics\Data\NGI\My Rectified\AATCORCorrected_oS131022114824832b_278400201_Uint16.tif';
ngiFileName = 'G:\PhD GeoInformatics\Data\NGI\My Rectified\StudyAreaXCalibMosaicSpotExtent.tif';

spotIm = int16(imread(spotFileName));
ngiIm = int16(imread(ngiFileName));
spotIm = spotIm(:,:,1:3);

figure
for i = 1:3
    spotBand = spotIm(:,:,i);
    ngiBand = ngiIm(:,:,i);
    subplot(1,3,i)
    plot(spotBand(1:1000:end), ngiBand(1:1000:end), 'kx');
end

%% 2017
% Check min vals of raw ngi im to justify C=0
ngiRawFileName = 'F:\PhD GeoInformatics\Data\NGI\Cross Calibration\Mosaics\StudyAreaUncalibratedMosaicSpotMask.tif'
im = imread(ngiRawFileName);
im(im<0)=0;
imd = double(im)/4096;
imshow(imd)
mask = imd(:,:,1)>0;
for i = 1:3
    b = im(:,:,i);
    %max(b(:))
    prctile(b(mask),99.99)
    prctile(b(mask),0.001)
    min(b(mask))
end

%% compare SPOT and DMC spectra in invariant locations
% spotFileName = 'F:\PhD GeoInformatics\Data\NGI\My Rectified\spotInt16_2.tif';
% ngiRawFileName = 'F:\PhD GeoInformatics\Data\NGI\Cross Calibration\Mosaics\StudyAreaUncalibratedMosaicSpotMask.tif';
% ngiCalibFileName = 'F:\PhD GeoInformatics\Data\NGI\My Rectified\StudyAreaXCalibMosaicSpotMask.tif';
spotFileName = 'D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\My Rectified\spotInt16_2.tif';
ngiRawFileName = 'D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\Cross Calibration\Mosaics\StudyAreaUncalibratedMosaicSpotMask.tif';
ngiCalibFileName = 'D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\My Rectified\StudyAreaXCalibMosaicSpotMask.tif';

spotIm = imread(spotFileName);
ngiCalibIm = imread(ngiCalibFileName);
ngiRawIm = imread(ngiRawFileName);

mask = ~(ngiCalibIm < 0);
ngiCalibIm(~mask)=0;
mask = mask(:,:,1);
%imshow(mask)

% do global adjustment of images
spotImN = zeros(size(spotIm), 'int16');
if true
    for i = 1:3
        bs = single(spotIm(:,:,i));
        bn = single(ngiCalibIm(:,:,i));
        bsn = bs*mean(bn(mask))/mean(bs(mask));
        spotImN(:,:,i) = int16(bsn);
    end
end

clear error*
errorImN = abs(spotImN(:,:,1:3) - ngiCalibIm);
errorIm = abs(spotIm(:,:,1:3) - ngiCalibIm);

for i = 1:3
    ben = single(errorImN(:, :, i));
    be = single(errorIm(:, :, i));
    errorImNF(:,i) = ben(mask);
    errorImF(:,i) = be(mask);
end

reflFactor = 5000;
mean(errorImNF, 1)/reflFactor
mean(errorImF, 1)/reflFactor
%%
[Spectra, IjLocs] = PlotMultibandImages2({single(spotIm(:, :, 1:3))/reflFactor, single(ngiCalibIm)/reflFactor}, IjLocs);
labels = {'Water', 'Shadow', 'Vegetation 1', 'Bright Sand 1', 'Bright Sand 2', 'Bright Sand 3', 'Bare Ground 1', 'Bare Ground 2', 'Vegetation 3', 'Vegetation 2'}
save 'C:\Data\Development\Projects\PhD GeoInformatics\Code\Results\Cross Calibration\SpotDmcSpectra2.mat' Spectra IjLocs labels

save 'C:\Data\Development\Projects\PhD GeoInformatics\Code\Results\Cross Calibration\SpotDmcSpectra.mat' Spectra SpectraN IjLocs IjLocsN
% SpectraN = Spectra
% IjLocsN = IjLocs

%%
figure
for i = 1:size(Spectra{1}, 1)
    subplot(4,3, i)
    plot(Spectra{1}(i,:),'k:')
    hold on
    plot(Spectra{2}(i,:),'k')
    grid on
    title(labels{i})
end
legend({'Spot', 'DMC'})

%% plot for paper
load 'C:\Data\Development\Projects\PhD GeoInformatics\Code\Results\Cross Calibration\SpotDmcSpectra2.mat'
labels = {'Water', 'Shadow', 'Vegetation 1', 'Bright Sand', 'Bright Sand 2', 'Bright Sand 3', 'Bare Ground', 'Bare Ground', 'Vegetation 3', 'Vegetation 2'};

% get raw dn spectra for ij locs
ngiRawIm = imread(ngiRawFileName);
mask = ~(ngiRawIm < 0);
ngiRawIm(~mask)=0;
SpectraDn = [];
[SpectraDn, IjLocsDn] = PlotMultibandImages2({single(spotIm(:, :, 1:3))/reflFactor, ngiRawIm}, IjLocs);
save 'C:\Data\Development\Projects\PhD GeoInformatics\Code\Results\Cross Calibration\SpotDmcSpectra2.mat' Spectra IjLocs labels SpectraDn
% %4 conn spectrum
% for i = 1:size(IjLocs, 1)
%     CoOrd = repmat(round(IjLocs(i, :)), 5, 1) + [0 1;0 -1;1 0;-1 0;0 0];
% %     Spectrum = squeeze(Image(CoOrd(:, 1), CoOrd(:, 2), :))';
%     SpectraDn(i, :) = squeeze(ngiRawIm(CoOrd(1), CoOrd(2), :))';
%     SpectraDn(i, :) = mean(SpectraDn(i, :), 1);
% end
% close all;
% figure
% for i = 1:length(specIdx)
%     subplot(2, 3, i)
%     plot(SpectraDn{2}(specIdx(i), bandIdx), 'k-x');
% end
%%
load 'C:\Data\Development\Projects\PhD GeoInformatics\Code\Results\Cross Calibration\SpotDmcSpectra2.mat'
labels = {'Water', 'Shadow', 'Vegetation 1', 'Bright sand', 'Bright sand 2', 'Bright sand 3', 'Bare ground', 'Bare ground', 'Vegetation 3', 'Vegetation 2'};
bandIdx = [3 2 1];
bandLabels = {'Green', 'Red', 'NIR'};
specIdx = [1 4 8 3 10];
fontSize = 12;
titles = {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)'};
% close all;
fig = figure;
left_color = zeros(1,3);
right_color = zeros(1,3);
set(fig,'defaultAxesColorOrder',[left_color; right_color]);

for i = 1:length(specIdx)
    subplot(2, 3, i)

    yyaxis left
    plot(Spectra{1}(specIdx(i), bandIdx), 'k-x', 'LineWidth', 1, 'MarkerSize', 7);
    hold on
    plot(Spectra{2}(specIdx(i), bandIdx), 'k-.x', 'LineWidth', 1, 'MarkerSize', 7);
    ytickformat('%.1f')
    yyaxis right
    plot(SpectraDn{2}(specIdx(i), bandIdx), 'k--x', 'LineWidth', 1, 'MarkerSize', 7);
    ytickformat('%g')
    ylim([0 2048]);
    ylabel('DN')
    %set(gca, 'Font', 'arial')
%     set(gca, 'FontSize', 8)
    yyaxis left
    %grid on
    ylim([0 1]);
    ylabel('Surface reflectance')
    ax = gca;
    ax.YTick = [0:0.2:1];
    ax.XTick = [1 2 3];
    ax.XTickLabels = bandLabels;
%     title(labels(specIdx(i)), 'FontWeight', 'Normal', 'FontSize', fontSize)
    title(titles{i}, 'FontWeight', 'Normal', 'FontSize', fontSize)
    set(gca, 'FontName', 'Arial')
    set(gca, 'FontSize', fontSize-1)
    ax = gca;
    ax.TickDir = 'in';
    ax.TickLength = ax.TickLength*2;
    set(gca, 'box', 'off')
end
legend({['SPOT 5 surface' newline 'reflectance'], ['DMC surface' newline 'reflectance'], 'DMC DN'}, 'Location', 'NorthWest', 'FontSize', fontSize)
legend boxoff
subplot(2,3,6)
axis off
mad = mean(mean(abs(Spectra{2}(specIdx,:)-Spectra{1}(specIdx,:))))
rms = sqrt(mean(mean((Spectra{2}(specIdx,:)-Spectra{1}(specIdx,:)).^2)))
s = sprintf('MAD = %1.2f%%\nRMS = %1.2f%%', mad*100, rms*100)
% text(0, 0, s, 'FontName', 'Arial', 'FontSize', 14)
t = annotation('textbox');
t.String = s;
s = t.FontSize;
t.FontSize = fontSize;
t.FontName = 'Arial';
t.BackgroundColor = 'white';
t.Position(1:2) = [0.7 0.1];
t.LineStyle = 'none';
% t.LineWidth = 0;

legend boxoff

print('C:\Data\Development\Projects\PhD GeoInformatics\Docs\My Docs\Thesis\Retrieval of Surface Reflectance from Aerial Imagery\Figure 16 - Comparison of DMC and SPOT Spectra.eps',... 
    '-depsc', '-r600')
%%

mean(mean(abs(SpectraN{2}-SpectraN{1})))
mean(mean(abs(Spectra{2}-Spectra{1})))

