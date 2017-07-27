%%
%Look at logs to understand light section in 3172, dark patch in 322, calib
%failures etc
clear all;close all;
cd 'D:\Data\Development\Projects\MSc GeoInformatics\Code\Results\ZIPPS'

fn = {'GPP_3321b_3172_20Apr2013_1033.log', ...
    'GPP_3321d_319_11Apr2013_0941.log', ...
    'GPP_3322a_320_15Apr2013_1556.log', ...
    'GPP_3322c_322_10Apr2013_2234.log'};

for i = 1:length(fn)
    res = ParseZippsLog(fn{i});

    idx = ~cellfun(@isempty, {res.MsPlatformCalibOk}) & ~cellfun(@isempty, {res.RadiomPlatformCalibOk});
    res = res(idx);

    figure;
    h1 = subplot(3,1,1);
    plot([res.Id], [res.Exp], '--')
    legend({'Pan', 'NIR', 'Red', 'Green', 'Blue'});
    title([fn{i} ' -  Exposure'])
    h2 = subplot(3,1,2);
    plot([res.Id], [res.MsPlatformCalibOk])
    title([fn{i} ' -  MsPlatformCalibOk'])
    h3 = subplot(3,1,3);
    plot([res.Id], [res.RadiomPlatformCalibOk])
    title([fn{i} ' -  RadiomPlatformCalibOk'])

    linkaxes([h1 h2 h3], 'x');
end



%%
%Following the above, analyse the more comprehensive txt files produced by
%gplog2txt from above log files
close all;clear all;
fn = {'GPP_3321b_3172_20Apr2013_1033.txt', ...
    'GPP_3321d_319_11Apr2013_0941.txt', ...
    'GPP_3322a_320_15Apr2013_1556.txt', ...
    'GPP_3322c_322_10Apr2013_2234.txt'};

for i = 1:length(fn)
    res = ParseZippsTxt(fn{i});

%     idx = ~cellfun(@isempty, {res.MsPlatformCalibOk}) & ~cellfun(@isempty, {res.RadiomPlatformCalibOk});
%     res = res(idx);
    ims = char({res.Image});
%     delimIdx = strfind(ims, '~');
    ids = str2double(cellstr(ims(:,15:end)));
    
    [ids, sortIdx] = sort(ids);
    res = res(sortIdx);
    MSplatformcalib = {res.MSplatformcalib};
    MSplatformcalib = strcmpi(MSplatformcalib, 'successful');
    
    PANplatformcalib = {res.PANplatformcalib};
    PANplatformcalib = strcmpi(PANplatformcalib, 'successful');
    
    GraduateNIR = {res.GraduateNIR};
    GraduateNIR = strcmpi(GraduateNIR, 'good');
    
    GraduatePan = {res.GraduatePan};
    GraduatePan = strcmpi(GraduatePan, 'good');

    GraduateRGB = {res.GraduateRGB};
    GraduateRGB = strcmpi(GraduateRGB, 'good');
    
    figure;
    h1 = subplot(3,1,1);
    plot(ids, MSplatformcalib)
    hold all;
    plot(ids, PANplatformcalib)
    plot(ids, GraduateNIR)
    plot(ids, GraduatePan)
    plot(ids, GraduateRGB)
    legend({'MSplatformcalib', 'PANplatformcalib', 'GraduateNIR', 'GraduatePan', 'GraduateRGB'})
    title(fn{i})
    axis tight
    
    h2 = subplot(3,1,2);
    plot(ids, str2double({res.GNRMS}))
    hold all
    plot(ids, str2double({res.GBRMS}))
    plot(ids, str2double({res.GVRMS}))
    plot(ids, str2double({res.GVDist})/10)
    legend({'GNRMS', 'GBRMS', 'GVRMS', 'GVDist'})
    axis tight

    h3 = subplot(3,1,3);
    plot(ids, str2double({res.FStopRed}), 'r')
    hold all
    plot(ids, str2double({res.FStopNIR}), 'g')
    plot(ids, str2double({res.FStopPan}), 'b')
    plot(ids, str2double({res.ExposureRed}), 'r--')
    plot(ids, str2double({res.ExposureNIR}), 'g--')
    plot(ids, str2double({res.ExposurePan}), 'b--')
    evRed = str2double({res.ExposureRed})./(str2double({res.FStopRed}).*str2double({res.FStopRed}));
    evNir = str2double({res.ExposureNIR})./(str2double({res.FStopNIR}).*str2double({res.FStopNIR}));
    evPan = str2double({res.ExposurePan})./(str2double({res.FStopPan}).*str2double({res.FStopPan}));
    plot(ids, 100*evRed, 'k');
    plot(ids, 100*evNir, 'k--');
    plot(ids, 1000*evPan, 'k:');
    
    legend({'FStopRed', 'FStopNIR', 'FStopPan', 'ExposureRed', 'ExposureNIR', 'ExposurePan', 'EVRed', 'EVNIR', 'EVPan'})
    axis tight

    linkaxes([h1 h2 h3], 'x');
end
tilefigs

figure;
plot(evPan./evRed)
hold all
plot(evNir./evRed)
hold all

%%
%NOTES
%-------------------------------------------------------------------------
%- In 3172 the images show a change in brightness around 211.  The
%exposures also show a change around here. 
%- Each channel has different exposures - is this compensated for ????
%- In 322 aorund 264, there is a dark patch that looks like it could be a
%cloud (an irregular dark section with overexposed surrounds).  This
%corresponds to a sudden drop in exposure time which doesn't quite make
%sense - would expect exposure to have gone up. All *calibOk=true round here
%though.
%- In general there is a lot of jumping around of exposures.  This is not
%correlated betw the different bands.
%- Checking out the images from before/after an IR exposure jump.  There
%does seem to a a relative IR difference in the calibrated CIR ims so I dont think
%ZIPPS is calibrating for exposure.
%- Checking out the images from before/after an RGB exposure jump.  There
%is no obvious RGB difference as expected (the R,G & B exposures all change
%together).  There is however PERHAPS a CIR difference as the relative IR/RGB
%exposure changes!!!
%- NB Summary of above: RGB exposure change together so RGB exposure
%compensation not necessary but IR-RGB exposure change independently and
%dont seem to be compensated for.  Also PAN-RGB change independently and
%PAN is used for PAN sharpening, does "fixed histogram option" cater for
%this?  See more below
% - MsPlatformCalibOK: this is how accurately the bands are registered to
% each other.  It fails occasionally.  It seems to always fail with NIR.
% When it fails, it seems to mostly fail with still small errors so is
% probably not serious.  After failure, ZIPPS tries to reprocess (robust
% registration?).  This mostly seems to fail again.
% - Fstop: assume f number is the denom std eg f/16 (because PAN has higher
% f nums - we would xpect the PAN aperture to be narrower because it is
% spectrally wider)
% - EV: amount of light is prop to exp./(fstop^2).  Looking at EV plots
% sheds some more light on what is happening with the auto fstop and exp
% vals.  The RGB, Pan and NIR EV vals vary roughly proportionally to each
% other i.e. the auto fstop/exp is compensating for changes in the ttl
% amount of light more than changes in colour. BUT it is not exactly
% proportional so it is also compensating for changes in colour.  This is
% NB for calibration. 
%- I assume the auto exp/fstop is achieving a mean val in each colour band.
% The question is how much of this value is due to varying ambient
% conditions (BRDF, sun angle, haze) and how much is due to the actual
% scene.  There are still some unknowns here as it does not seem to be
% making all images "white" on avg which is my assumption for the auto alg
% above.  
%- Looking at more before/after EV change images does seem to show a slight IR/RGB
%relative change.  I dont think this can be purely compensated for based on
%EV scaling of bands or something similar as the EV adjustments are
%obviously occurring due to real colour/brightness changes in the scene -
%so you could actually make things worse by asjusting for EV.
%- Looking at the MS registration rms errors/distances: there doesn't seem
%to be anything untoward happening here so I wouldn't worry about
%MSPlatformCalibOk=false too much - they seem to be only just outside spec
%(although I can't see what threshold spec this is from the graphs)
%RadiomPlatformCalibOk=false is still a mystery to me - the .sum files say
%all the images are OK ?!
%- Checking between NGI RGB and my rectified RGB: there is little
%difference in sharpness of edges.  Band registration and focus accuracy
%looks similar.  There are more shadows (& hence texture) apparent in my imagery.  
%IR is in general not as well focused as RGB (expected due to wavelen?).
%There are registration differences between the 2 of at least 1 meter (sometimes a lot more).  Who
%is better/worse and what is the source of the issue?
%- The Brovey transformation uses the fact that the PAN camera nearly has 
%the same sensitivity as the sum of all multi-spectral camera heads. 
%Therefore, the NIR correction and blue correction for the RGB and CIR 
%image are not necessary. (The settings in the Pan Sensitivity tab are 
%automatically ignored)

%ZIPPS PROC CHECKING TO DO:
%- Check is ZIPPS compensates for exposures
%- (1) Can you see spekboom? yes - checked through all the ground truth
%polygons and it pretty much makes sense.  The rooiberg ground truth looked
%a bit out but I think this is human factor
%- (2) are bands registered? yes - errors in logs mostly ok and ims look
%more or less ok
%- (3) is rectification OK? yes - difficult to choose between NGI and mine,
%usually mine actually looks more convincing with the ground truth polygons
%- Can PAN be recovered from sum(RGB,IR)? Yes
%- Relook at doc for pan sharpening options. - done - we have selected ok
%- Have we chosen a good format - what about RGB&PAN at orig res? PAN can
%be recovered so not necessary.  our format is convenient for both
%visualisation and classification