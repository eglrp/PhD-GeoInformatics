close all; clear all;

imFileNames = { 'F:\MSc GeoInformatics\Data\PCI\Radiometric Check\o3321C_2010_318_01_0002_RGB1.tif', ...
    'F:\MSc GeoInformatics\Data\PCI\Radiometric Check\o3321D_2010_319_01_0033_RGB1.tif'};
% imFileNames = { 'F:\MSc GeoInformatics\Data\NGI\Calibrated\3321C_2010_318\RGB\3321C_2010_318_01_0001_RGB.tif'};

% lutFileNames = {'F:\MSc GeoInformatics\Data\NGI\LUT\3321D\rgb_12.lut',...
%     'F:\MSc GeoInformatics\Data\NGI\LUT\3321C\rgb_12.lut'};
lutFileNames = {'F:\MSc GeoInformatics\Data\NGI\Aerial Triangulation\3321A_2010_316\ZI_Final_Files\LUT\rgb_12.lut',...
    'F:\MSc GeoInformatics\Data\NGI\Aerial Triangulation\3321B_2010_317\ZI_Final_Files\LUT\rgb_12.lut'};
icons = {'-x', '-o'};
for i = 1:length(imFileNames)
    t{i} = Tiff(imFileNames{i}, 'r');
    im{i} = t{i}.read();
    
    [invLut{i} lut{i}] = ReadIntergraphLutFile(lutFileNames{i});
    imLut{i} = ApplyIntergraphLut(im{i}, invLut{i});
    imLutLut{i} = ApplyIntergraphLut(imLut{i}, lut{i});

    figure(1)
    plot(lut{i}, icons{i});
    hold on;
    grid on;

    figure(2)
    subplot(1,2,i)
    imagesc(bitshift(im{i}, 4))
    axis image

    figure(3)
    subplot(1,2,i)
%     imagesc(bitshift(imLut{i}, 6))
    imagesc(imLut{i}*50)
%     imagesc(sqrt(float(imLut{i})*16)*256)
    axis image
end
%%
%make PCI aero-tri file
close all;clear all;
idx = [2 3 68 69 72 73 433 434 436 437 440 441 ...
    138 139 142 143 208 209 212 213 278 289 282 283 348 349 352 353 444 445 448 449 ...
    452 453 455 456 459 460 462 463 465 466 469 470];

ngiOrFileName = 'F:\MSc GeoInformatics\Data\NGI\Aerial Triangulation\3321C_2010_318\extori3321C_2010_318_lo21wgs84n_e_rect.or';
outFile = 'D:\Data\Development\Projects\MSc GeoInformatics\Docs\PCI\Acquisition Time Effect\extori3321D_2010_318_lo21wgs84n_e_rect.txt';
imDir = 'F:\MSc GeoInformatics\Data\NGI\Calibrated\3321C_2010_318\RGB\';

MakePCIExtOrTextFile(ngiOrFileName, idx, outFile, imDir);
%%
%make PCI aero-tri file
close all;clear all;
idx = [33 34 39 40 105 106 111 112 178 179 184 185 286 287 292 293 358 359 ...
    364 365 401 402 ...
    437 438 441 442 444 445 447 448 450 451 453 454 456 457 459 460 464 465 ...
    468 469 471 472 474 475];

% idx = 1:518;

ngiOrFileName = 'F:\MSc GeoInformatics\Data\NGI\Aerial Triangulation\3321D_2010_319\extori3321D_2010_319_lo21wgs84n_e_rect.or';
outFile = 'D:\Data\Development\Projects\MSc GeoInformatics\Docs\PCI\Acquisition Time Effect\extori3321D_2010_319_lo21wgs84n_e_rect.txt';
imDir = 'F:\MSc GeoInformatics\Data\NGI\Calibrated\3321D_2010_319\RGB\';

MakePCIExtOrTextFile(ngiOrFileName, [], outFile, imDir);
%%
%Make eg images for bernhard

close all; clear all;

imFileNames = { 'D:\Data\Development\Projects\MSc GeoInformatics\Code\Results\LutTest\3321D_2010_319_01_0033_CIR1.tif'};
% imFileNames = { 'F:\MSc GeoInformatics\Data\NGI\Calibrated\3321C_2010_318\RGB\3321C_2010_318_01_0001_RGB.tif'};

lutFileNames = {'F:\MSc GeoInformatics\Data\NGI\Aerial Triangulation\3321D_2010_319\ZI_Final_Files\LUT\cir_12.lut'};
icons = {'-x', '-o'};
for i = 1:length(imFileNames)
    t{i} = Tiff(imFileNames{i}, 'r');
    im{i} = t{i}.read();

    [invLut{i} lut{i}] = ReadIntergraphLutFile(lutFileNames{i});
    imLut{i} = ApplyIntergraphLut(im{i}, invLut{i});
    imLutLut{i} = ApplyIntergraphLut(imLut{i}, lut{i});

    figure(1)
    plot(lut{i}, icons{i});
    hold on;
    grid on;

    figure(2)
%     subplot(1,2,i)
    subplot(1,2,1)
    imagesc(bitshift(im{i}, 4))
    axis image
    axis off

% figure(3)
%     subplot(1,2,i)
%     imagesc(bitshift(imLut{i}, 6))
%     imLut{i}(:,:,3) = imLut{i}(:,:,3)/1.5;
    subplot(1,2,2)
    imagesc(imLut{i}*32)
%     imagesc(sqrt(float(imLut{i})*16)*256)
    axis image
    axis off
end

lutd = sum(lut{1}(:, 2:4), 2);
lutn = double(lut{1}(:,2:4))./repmat(lutd, 1, 3);
figure;
plot(lutn) 
%NOTE this shows that normalised colour is damaged by non-linear LUT

%%
%Compare vis and ir luts
close all; clear all;
lutFileNames = {'F:\MSc GeoInformatics\Data\NGI\Aerial Triangulation\3321D_2010_319\ZI_Final_Files\LUT\rgb_12.lut';...
    'F:\MSc GeoInformatics\Data\NGI\Aerial Triangulation\3321D_2010_319\ZI_Final_Files\LUT\cir_12.lut'};
icons = {'-', '--'};

for i = 1:length(lutFileNames)

    [invLut{i} lut{i}] = ReadIntergraphLutFile(lutFileNames{i});
    figure(1)
    plot(lut{i}, icons{i}); 
    hold on;
    grid on;
end

% function F = myfun(x,xdata)
% F = x(1)*exp(x(2)*xdata);
%see how well a power law fits the LUT data
x0 = [1 2];
xdata = double(invLut{1}(:, 1));
ydata = double(invLut{1}(:, 3));
[x,resnorm] = lsqcurvefit(@(x,xdata)x(1)*xdata.^x(2), x0, xdata, ydata);
yhat = x(1)*xdata.^x(2);

figure;
plot(xdata, ydata,'r-');
hold on;
plot(xdata, yhat,'b--');

%
%NOTE this shows LUT is not a power law and that same bands from RGB and
%CIR LUTs are weirdly unrelated


%%
%Compare inverted LUT R & G bands from RGB and CIR images for LUT accuracy

close all; clear all;

imFileNames = { 'D:\Data\Development\Projects\MSc GeoInformatics\Code\Results\LutTest\3321D_2010_319_01_0033_RGB1.tif'; ...
    'D:\Data\Development\Projects\MSc GeoInformatics\Code\Results\LutTest\3321D_2010_319_01_0033_CIR1.tif'};

lutFileNames = {'F:\MSc GeoInformatics\Data\NGI\Aerial Triangulation\3321D_2010_319\ZI_Final_Files\LUT\rgb_12.lut';...
    'F:\MSc GeoInformatics\Data\NGI\Aerial Triangulation\3321D_2010_319\ZI_Final_Files\LUT\cir_12.lut'};
    
icons = {'-x', '-o'};
for i = 1:length(imFileNames)
    t{i} = Tiff(imFileNames{i}, 'r');
    tmp = t{i}.read();
    im{i} = tmp(1:10:end,1:10:end,:);
    
    [invLut{i} lut{i}] = ReadIntergraphLutFile(lutFileNames{i});
    imLut{i} = ApplyIntergraphLut(im{i}, invLut{i});
    imLutLut{i} = ApplyIntergraphLut(imLut{i}, lut{i});
    
    im{i} = single(im{i})/4096;
    imLut{i} = single(imLut{i})/4096;
    imLutLut{i} = single(imLutLut{i})/4096;

    figure(1)
    plot(lut{i}, icons{i});
    hold on;
    grid on;

    figure(2)
    subplot(1,2,i)
    imagesc(im{i})
    axis image
    axis off

    figure(3)
    subplot(1,2,i)
    imagesc(imLut{i})
%     imLut{i}(:,:,3) = imLut{i}(:,:,3)/1.5;
%     subplot(1,2,2)
%     imagesc(imLut{i}*32)
%     imagesc(sqrt(float(imLut{i})*16)*256)
    axis image
    axis off
end

figure;
subplot(1,2,1)
imagesc(abs(im{1}(:,:,1) - im{2}(:,:,2)));
subplot(1,2,2)
dr = imLut{2}(:,:,2) - imLut{1}(:,:,1);
imagesc(abs(dr));

figure;
subplot(1,2,1)
imagesc(abs(im{2}(:,:,3) - im{1}(:,:,2)));
subplot(1,2,2)
dg = imLut{2}(:,:,3) - imLut{1}(:,:,2);
imagesc(abs(dg));
%NOTE - we see there are significant differences between the same bands
%inverted from different sources with their appropriate LUT

im4band = cat(3, imLut{2}(:,:,2), imLut{2}(:,:,3), imLut{1}(:,:,3), imLut{2}(:,:,1));

figure;
subplot(1,2,1)
imshow(2*im4band(:,:,1:3));
subplot(1,2,2)
imshow(2*im4band(:,:,[4 1 2]));

figure;
imshow(2*imLut{1});

%%
%Find a white balance to make nice 4 band images - no LUT
close all;


wb = m{1}([1 2])./m{2}([2 3])
% wb = [1 1]
im4band2 = cat(3, wb(1)*im{2}(:,:,2), wb(2)*im{2}(:,:,3), im{1}(:,:,3), im{2}(:,:,1));

figure;
subplot(1,2,1)
imshow(im4band2(:,:,1:3))
subplot(1,2,2)
imshow(im{1})

%%
% Check the histograms of inverted lut r,g bands from rgb and cir images
% to see if we can make sense of what is missing from the LUT inversion to
% get these equal


close all;
figure;
idx = [1 2 3; 2 3 1];
for i = 1:2
    for j = 1:size(idx, 2)
        data = im{i}(:, :, idx(i, j));
%         [bandwidth, density{i,j}, xmesh, cdf] = kde(data(:), 512, 0, 1);
        subplot(size(idx,2),1,j)
        xmesh = linspace(0,1,512);
        density{i,j} = hist(data(:), xmesh);
        plot(xmesh,density{i,j});
        hold all;
        grid on
        min(data(:))
        sum(density{i,j})
    end
end
%NOTE the histograms of the r,g bands are tantalisingly close

dataR1 = imLut{1}(:,:,2);
dataR1 = dataR1(1:end)';
dataR2 = imLut{2}(:,:,3);
dataR2 = dataR2(1:end);

figure;
plot(dataR1,dataR2,'x')

dataR2 = [dataR2(:) ones(numel(dataR1),1)];

A = dataR2\dataR1
mean(dataR1)-mean(dataR2)


dataR1 = im{1}(:,:,2);
dataR1 = dataR1(1:end)';
dataR2 = im{2}(:,:,3);
dataR2 = dataR2(1:end);

figure;
plot(dataR1,dataR2,'x')

dataR2 = [dataR2(:) ones(numel(dataR1),1)];

A = dataR2\dataR1
mean(dataR1)-mean(dataR2)

%%
%check rN, gN from rgG between LUT and inv LUT images
rN1 = im{1}(:,:,2)./sum(im{1}(:,:,[1 2]), 3);
rN2 = im{2}(:,:,3)./sum(im{2}(:,:,[2 3]), 3);

mean(mean(abs(rN1-rN2)))/mean(mean(rN1))

figure;
imagesc(log(abs(rN1-rN2)))

rN1 = imLut{1}(:,:,2)./sum(imLut{1}(:,:,[1 2]), 3);
rN2 = imLut{2}(:,:,3)./sum(imLut{2}(:,:,[2 3]), 3);

mean(mean(abs(rN1-rN2)))/mean(mean(rN1))

figure;
imagesc(log(abs(rN1-rN2)))
%NOTE: these are closish

%%
%Check histograms of "calibrated" images across images inside jobs and
%between jobs
clear t
fileNames = {...
    'G:\MSc GeoInformatics\Data\NGI\My Calibrated\3321C_2010_318\3321C_2010_318_01_0003_RGBIR.tif',...
    'G:\MSc GeoInformatics\Data\NGI\My Calibrated\3321C_2010_318\3321C_2010_318_07_0221_RGBIR.tif',...
    'G:\MSc GeoInformatics\Data\NGI\My Calibrated\3321C_2010_318\3321C_2010_318_13_0443_RGBIR.tif',...
    'G:\MSc GeoInformatics\Data\NGI\My Calibrated\3321D_2010_319\3321D_2010_319_01_0003_RGBIR.tif',...
    'G:\MSc GeoInformatics\Data\NGI\My Calibrated\3321D_2010_319\3321D_2010_319_09_0308_RGBIR.tif',...
    'G:\MSc GeoInformatics\Data\NGI\My Calibrated\3321D_2010_319\3321D_2010_319_14_0499_RGBIR.tif',...
    'G:\MSc GeoInformatics\Data\NGI\My Calibrated\3321D_2010_319\3321D_2010_319_01_0033_RGBIR.tif',...
    };
% lutFileNames = {...
%     'F:\MSc GeoInformatics\Data\NGI\Aerial Triangulation\3321D_2010_319\ZI_Final_Files\LUT\rgb_12.lut';...
%     'F:\MSc GeoInformatics\Data\NGI\Aerial Triangulation\3321D_2010_319\ZI_Final_Files\LUT\cir_12.lut'
%     };

for i = 1:length(fileNames)
    t(i) = Tiff(fileNames{i});
    
end

%check different tiles in 1 image
figure;
icons = {'r','g','b','m'};
for i = 1:10:t(1).numberOfTiles
    tile = t(1).readEncodedTile(i);
    for b = 1:4
        band = tile(:,:,b);
        [bandwidth, density{i,b}, xmesh, cdf] = kde(band(:), 255, 0, 4096);
        subplot(4,1,b)
        plot(xmesh, density{i,b});
        hold on;
    end
end

%check over many images
figure;
icons = {'r','g','b','m'};
for i = 1:length(t)
      tmp = t(i).read();
      im{i} = tmp(1:10:end,1:10:end,:);
    
    for b = 1:4
          band = im{i}(:,:,b);
          [bandwidth, density{i,b}, xmesh, cdf] = kde(band(:), 255, 0, 4096);
        subplot(length(t),1,i)
        plot(xmesh, density{i,b});
        hold all;
    end
end
legend(fileNames);

%NOTE: I thought dodging may be aiming at achieving a specific histogram
%but the above plots would indicate not.  Some images have similar
%histograms between bands but not all and there are definitely not similar
%histograms between images.

%NOTES
%------------------------------------------------------------------------
%- Colour (rgG) remains the same under any scaling of intensity i.e. all 
%channels scaled the same.
%- Colour (rgG) does not remain the same under arbitrary gamma / LUT
%therefore we have to be able to undo the LUT which seems to be a problem
%as we get blue bias as opposed to green bias.
%- LUT followed by dodging (dodgning assumed=intensity scale): if we invert
%the LUT we then get a scaled version of the original value.  Only if the
%LUT is identical for all channels will we end up with the same colour rgG.
%NB The above only applies for a power law LUT which is apparently not what
%we have.
%- So we are a bit screwed: we don't seem to be inverting the LUT correctly
%and even if we could we would have a problem as the channels each have
%different LUTs
%- The last graph tells us that the LUT kills colour esp for vals below
%~1500
%- Given that the CIR inv lut seems closer to the real thing, perhaps we
%should make our 4 band image from CIR + B from RGB
%- A wb gain between R,G from CIR to R,G from RGB does very little and does
%not fix the blueness

%MORE NOTES ON HISTOGRAMS
%--------------------------------------------------------------------------
%- The "calibrated" images have equal histograms for RGB bands and
%different one for IR
%- The above suggests that either the LUT is there to achieve this and the
%dodging has no effect or the dodging achieves this. It seems more likely
%that it is the dodging.
%- Inverting the LUT gives different but similarly shaped histograms
%between R,G,B and between RG in RGB and RG in CIR.  There is a roughly
%linear but noisy (incl offset) relationship between R in RGB and R in CIR (ditto for G).
%- WHat is the reason for the noise??? Dodging?
%- The inv LUT CIR channels have a significant offset which is a concern.
%We dont see that the inv LUT channels are scaled versions of each other as
%would be the case if the LUT was an id power law for all channels and
%dodging was simply an intensity scale.
%- The histogram shapes of the inv LUT channels look suspiciously similar
%which cant be natural.  This implies either the dodging is doing some kind
%of histogram shaping or this is an artifact of inv LUT on dodged vals.
%- How does it happen that RG from RGB and RG from CIR have = histograms
%but not equal in a per pixel sense????
%- Basically the conclusion is that dodging does not seem to be merely a
%simple scaling which makes LUT inversion a pretty meaningless task.  It
%also means there will no longer be a linear relationship between DN and
%surface refl i.e. we can't do any kind of theoretical calibration to
%reflectance.  Out best bet would probably be to adjust the mosaic using
%histagram transformations found from the overlapping areas.
%- I dont think we can actually expect the inv LUT bands to be equal unless
%its a v special case of LUT and dodge
%- We should see how close rN (=R/(R+G)?) from RGB and CIR are
%- Do all images in a job have the same histograms?  Then perhaps we just
%need to adjust histograms between jobs to be = for a sensible calibration.

