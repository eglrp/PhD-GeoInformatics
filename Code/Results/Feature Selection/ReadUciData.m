%% Indian Pines (Hyperspec data http://www.ehu.eus/ccwintco/index.php?title=Hyperspectral_Remote_Sensing_Scenes)
% 200 bands, 16 class, vegetation class\
% class separation is poor and there looks to be some weird quantisation/noise
% issue with it, so rather leave this one out
close all; clear all;
load 'C:\Data\Development\Projects\MSc GeoInformatics\Data\Feature Selection\Hyperspectral\Indian_pines_corrected.mat'
load 'C:\Data\Development\Projects\MSc GeoInformatics\Data\Feature Selection\Hyperspectral\Indian_pines_gt.mat'
classNames = {'Alfalfa', ... 
    'Corn-notill', ...
    'Corn-mintill', ...
    'Corn', ...
    'Grass-pasture', ...
    'Grass-trees', ...
    'Grass-pasture-mowed', ...
    'Hay-windrowed', ...
    'Oats', ...
    'Soybean-notill', ...
    'Soybean-mintill', ...
    'Soybean-clean', ...
    'Wheat', ...
    'Woods', ...
    'Buildings-Grass-Trees-Drives', ...
    'Stone-Steel-Towers', ...
    };
unique(indian_pines_gt)

s = size(indian_pines_corrected);
data = reshape(indian_pines_corrected, [prod(s(1:2)), s(3)]);
lbl = indian_pines_gt(:);

data(lbl==0, :) = []; % remove background
data = prdataset(data, lbl(lbl>0));
data = remclass(data);
data = setlablist(data, classNames);
data = setprior(data, 0);
data = setfeatlab(data, num2str([1:size(data, 2)]'));

datan = data./repmat(sum(+data, 2), 1, size(data, 2));
PrPlotSpec(data(1:10:end, :))

[tr, ts] = gendat(data, 0.5);
feats = feast('jmi', 20, +tr, getnlab(tr))

w = tr(:, feats)*(pcam([],0.95)*opencvknnc([], 5))
w = tr(:, feats)*(scalem([], 'variance')*libsvc([], proxm([], 'r', 0.3)))
out = ts(:, feats)*w
out*testc
c = confmat(out)

%% KSC (Hyperspec data http://www.ehu.eus/ccwintco/index.php?title=Hyperspectral_Remote_Sensing_Scenes)
% 176 bands, 13 class (vegetation)

clear all;
load 'C:\Data\Development\Projects\MSc GeoInformatics\Data\Feature Selection\Hyperspectral\KSC.mat'
load 'C:\Data\Development\Projects\MSc GeoInformatics\Data\Feature Selection\Hyperspectral\KSC_gt.mat'
cl = {'Scrub', 'Willow swamp', 'Cabbage palm hammock', 'Cabbage palm/oak hammock', 'Slash pine', ...
    'Oak/broadleaf hammock', 'Hardwood swamp', 'Graminoid marsh', 'Spartina marsh', 'Cattail marsh',...
    'Salt marsh', 'Mud flats', 'Water'};

unique(KSC_gt)
s = size(KSC);
data = reshape(KSC, [prod(s(1:2)), s(3)]);
lbl = KSC_gt(:);

%data(lbl==0, :) = []; % remove background
data = prdataset(data(lbl>0, :), lbl(lbl>0)); % remove background
% remove outliers
data(any(+data > 1e4, 2), :) = [];
data = setprior(data, 0);
data = setprior(data, 0);
data = setfeatlab(data, num2str([1:size(data, 2)]'));
data = setlablist(data, cl);

% datan = data./repmat(sum(+data, 2), 1, size(data, 2));
PrPlotSpec(data(1:10:end, :))
%PlotMultibandImage(KSC)
c = corr(+data);
figure; imagesc(c);
clear KSC*

%%
% quick clfr test
[tr, ts] = gendat(data, 0.6);
tr_ = tr(1:end, :);
feats = feast('jmi', 10, +tr_, getnlab(tr_))

[err, cerr, nlabOut] = prcrossval(data, scalem([], 'variance')*opencvsvc([], [], {'SVMType', 'C_SVC', 'KernelType', 'RBF', 'Gamma', 1, 'C', 1, 'ClassWeights', ones(1, getsize(data, 3))}), 10)

w = tr(:, feats)*(scalem([], 'variance')*opencvsvc([], [], {'SVMType', 'C_SVC', 'KernelType', 'RBF', 'Gamma', 1, 'C', 1, 'ClassWeights', ones(1, getsize(data, 3))}))
w = tr(:, feats)*(scalem([], 'variance')*opencvknnc([], 3))
w = tr(:, feats)*(scalem([], 'variance')*libsvc([], proxm([], 'r', .3)))
w = tr(:, feats)*opencvdtreec([], 12, {'Priors', ones(1, getsize(data, 3))/getsize(data, 3), 'MaxDepth', 12, 'Use1seRule', false, 'UseSurrogates', false, 'CVFolds', 5, 'TruncatePrunedTree', true, 'MinSampleCount', min(classsizes(data))/100})
w = tr(:, feats)*libsvc([], proxm([], 'r', 30))
out = ts(:, feats)*w
out*testc
c = confmat(out)

%% num features
[tr, ts] = gendat(data, 0.6);
feats = feast('jmi', 25, +tr, getnlab(tr));

for f = 3:25
    w = tr(:, feats(1:f))*(scalem([], 'variance')*opencvknnc([], 3));
    resF(f) = ts(:, feats(1:f))*w*testc
end
figure
plot(resF, 'k-x')


%% save file
myPreferredFeatures = [];
myNumFeatures = 7;
myClusterThresh = 0.07; %0.05;
save('C:\Data\Development\Projects\MSc GeoInformatics\Data\Feature Selection\Hyperspectral\KscPr.mat', ...
    'data', 'myPreferredFeatures', 'myNumFeatures', 'myClusterThresh');


%% figure out best params & test different fs algorithms
m = FeatSelClusterRankM([], 'mi', 0, [], 'preferredFeatures', ...
        [], 'clusterThresh', 0.05, 'showFigures', false, 'jmiFormulation', false);

resFcr = BootstrapFsEval(data, m, 'numBootStraps', 5, ...
    'numFeatures', 7);

% m = FeatSelClusterRankM([], 'mi', 0, [], 'preferredFeatures', ...
%         [], 'clusterThresh', 0.02, 'showFigures', false, 'jmiFormulation', false);
% 
% resFcr2 = BootstrapFsEval(data, m, 'numBootStraps', 5, ...
%     'numFeatures', 10);

m = FeatSelFeastM([], 'jmi', 0);

resJmi = BootstrapFsEval(data, m, 'numBootStraps', 5, ...
    'numFeatures', 7);


%% Botswana (Hyperspec data http://www.ehu.eus/ccwintco/index.php?title=Hyperspectral_Remote_Sensing_Scenes)
% 145 bands, 14 class (vegetation)
clear all;
load 'C:\Data\Development\Projects\MSc GeoInformatics\Data\Feature Selection\Hyperspectral\Botswana.mat'
load 'C:\Data\Development\Projects\MSc GeoInformatics\Data\Feature Selection\Hyperspectral\Botswana_gt.mat'
cl = {'Water', 'Hippo grass', 'Floodplain grasses1', 'Floodplain grasses2', ...
    'Reeds1', 'Riparian', 'Firescar2', 'Island interior', 'Acacia woodlands', 'Acacia shrublands', ...
    'Acacia grasslands', 'Short mopane', 'Mixed mopane', 'Exposed Soils'};
unique(Botswana_gt);
s = size(Botswana);
data = reshape(Botswana, [prod(s(1:2)), s(3)]);
lbl = Botswana_gt(:);

%data(lbl==0, :) = []; % remove background
data = prdataset(data(lbl>0, :), lbl(lbl>0));
data = setprior(data, 0);
data = setfeatlab(data, num2str([1:size(data, 2)]'));
data = setlablist(data, cl);

datan = data./repmat(sum(+data, 2), 1, size(data, 2));
PrPlotSpec(data(1:10:end, :))
c = corr(+data);
figure; imagesc(c);
clear Botswana*
% PlotMultibandImage(Botswana)

% quick clfr test
[tr, ts] = gendat(data);
feats = feast('jmi', 25, +tr, getnlab(tr));

w = tr(:, feats)*(scalem([], 'variance')*opencvsvc([], [], {'SVMType', 'C_SVC', 'KernelType', 'RBF', 'Gamma', 2, 'C', 1, 'ClassWeights', ones(1, getsize(data, 3))}))
w = tr(:, feats)*(scalem([], 'variance')*opencvdtreec([], 12, {'Priors', ones(1, getsize(data, 3))/getsize(data, 3), 'MaxDepth', 12, 'Use1seRule', false, 'UseSurrogates', false, 'CVFolds', 5, 'TruncatePrunedTree', true, 'MinSampleCount', min(classsizes(data))/100}));
w = tr(:, feats)*(scalem([], 'variance')*opencvrtreec([], [], {'Priors', ones(1, getsize(data, 3))/getsize(data, 3), 'MaxNumOfTreesInTheForest', 50, 'NActiveVars', 4, 'CalcVarImportance', true, 'MaxDepth', 10, 'ForestAccuracy', 0.01}))

w = tr(:, feats)*(scalem([], 'variance')*opencvknnc([], 3))
w = tr(:, feats)*(scalem([], 'variance')*libsvc([], proxm([], 'r', 0.5)))
out = ts(:, feats)*w
out*testc
c = confmat(out)

%% best num features

for f = 3:20
    w = tr(:, feats(1:f))*(scalem([], 'variance')*opencvknnc([], 3));
    resF(f) = ts(:, feats(1:f))*w*testc
end
figure
plot(resF, 'k-x')

%% save file
myPreferredFeatures = [];
myNumFeatures = 7;
myClusterThresh = 0.03;
save('C:\Data\Development\Projects\MSc GeoInformatics\Data\Feature Selection\Hyperspectral\BotswanaPr.mat', ...
    'data', 'myPreferredFeatures', 'myNumFeatures', 'myClusterThresh');


%% figure out best params & test different fs algorithms
m = FeatSelClusterRankM([], 'mi', 0, [], 'preferredFeatures', ...
        [], 'clusterThresh', 0.03, 'showFigures', false, 'jmiFormulation', false);...

resFcr = BootstrapFsEval(data, m, 'numBootStraps', 5, ...
    'numFeatures', 7);

% m = FeatSelClusterRankM([], 'mi', 0, [], 'preferredFeatures', ...
%         [], 'clusterThresh', 0.017, 'showFigures', false, 'jmiFormulation', false);...
% 
% resFcr2 = BootstrapFsEval(data, m, 'numBootStraps', 5, ...
%     'numFeatures', 7);

m = FeatSelFeastM([], 'jmi', 0);

resJmi = BootstrapFsEval(data, m, 'numBootStraps', 5, ...
    'numFeatures', 7);


%% Cover Type (cartographic incl categorical vars)
% note that the last 44 columns are binary feats that require special
% interpretation.  running fcr on this shows that there is not much
% correlation between features
raw = csvread('C:\Data\Development\Projects\MSc GeoInformatics\Data\Feature Selection\UCI\Cover Type\covtype.data');
% 
data = prdataset(raw(:, 1:end-1), raw(:, end));
fl = {'Elevation', 'Aspect', 'Slope', 'Horizontal_Distance_To_Hydrology', 'Vertical_Distance_To_Hydrology', ...
'Horizontal_Distance_To_Roadways', 'Hillshade_9am', 'Hillshade_Noon', 'Hillshade_3pm', 'Horizontal_Distance_To_Fire_Points', ...
'Wilderness_Area (4 binary columns)', 'Soil_Type (40 binary columns)'};
cl = {'1 -- Spruce/Fir', '2 -- Lodgepole Pine', '3 -- Ponderosa Pine', '4 -- Cottonwood/Willow', '5 -- Aspen',...
    '6 -- Douglas-fir', '7 -- Krummholz'};

data = setfeatlab(data, [fl(1:10), repmat(fl(11),1,4), repmat(fl(12),1,40)]);
data = setlablist(data, cl);
scatterdui(data, 'legend')
save('C:\Data\Development\Projects\MSc GeoInformatics\Data\Feature Selection\UCI\Cover Type\Cover Type.mat', 'data');

%% Forest Types
% Not all the attributes may count as features - check the paper
%% Initialize variables.
filename = 'C:\Data\Development\Projects\MSc GeoInformatics\Data\Feature Selection\UCI\Forest Types\training.csv';
delimiter = ',';
startRow = 2;
formatSpec = '%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';

fileID = fopen(filename,'r');
dataArray1 = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);

filename = 'C:\Data\Development\Projects\MSc GeoInformatics\Data\Feature Selection\UCI\Forest Types\testing.csv';
fileID = fopen(filename,'r');
dataArray2 = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);

for i = 1:length(dataArray1)
    dataArray{i} = [dataArray1{i}; dataArray2{i}];
end

cl = dataArray{1};
data = [dataArray{2:end-1}];
fl = {'b1','b2','b3','b4','b5','b6','b7','b8','b9','pred\_minus\_obs\_H\_b1','pred\_minus\_obs\_H\_b2','pred\_minus\_obs\_H\_b3','pred\_minus\_obs\_H\_b4','pred\_minus\_obs\_H\_b5','pred\_minus\_obs\_H\_b6','pred\_minus\_obs\_H\_b7','pred\_minus\_obs\_H\_b8','pred\_minus\_obs\_H\_b9','pred\_minus\_obs\_S\_b1','pred\_minus\_obs\_S\_b2','pred\_minus\_obs\_S\_b3','pred\_minus\_obs\_S\_b4','pred\_minus\_obs\_S\_b5','pred\_minus\_obs\_S\_b6','pred\_minus\_obs\_S\_b7','pred\_minus\_obs\_S\_b8','pred\_minus\_obs\_S\_b9'};
cll =  {'s - Sugi forest', 'h - Hinoki forest', 'd - Mixed deciduous forest', 'o - Other non-forest land'};

data = prdataset(data, cl);
data = setfeatlab(data, fl);
data = setlablist(data, cll);

save('C:\Data\Development\Projects\MSc GeoInformatics\Data\Feature Selection\UCI\Forest Types\Forest Types.mat', 'data');

%% LRS data set
% 10 main classes with 10 sub-classes each - usage suspect

%'C:\Data\Development\Projects\MSc GeoInformatics\Data\Feature Selection\UCI\Low Resolution Spectrometer\lrs.data'

%% Statlog Landsat
% 	The database is a (tiny) sub-area of a scene, consisting of 82 x 100
% 	pixels. Each line of data corresponds to a 3x3 square neighbourhood
% 	of pixels completely contained within the 82x100 sub-area. Each line
% 	contains the pixel values in the four spectral bands 
% 	(converted to ASCII) of each of the 9 pixels in the 3x3 neighbourhood
% 	and a number indicating the classification label of the central pixel. 
% 	The number is a code for the following classes:
fn = 'C:\Data\Development\Projects\MSc GeoInformatics\Data\Feature Selection\UCI\Statlog (Landsat Satellite)\sat.trn';
data = dlmread(fn, ' ');

cl = {'1 - red soil', '2 - cotton crop', '3 - grey soil', '4 - damp grey soil', '5 - soil with vegetation stubble', ...
    	'6 - mixture class (all types present)', '7 - very damp grey soil'};
    
tr = prdataset(data(:, 1:end-1), data(:, end));
tr = setlablist(tr, cl);
tr = remclass(tr);
tr = setfeatlab(tr, num2str([1:size(tr, 2)]'));
tr = setprior(tr, 0);

fn = 'C:\Data\Development\Projects\MSc GeoInformatics\Data\Feature Selection\UCI\Statlog (Landsat Satellite)\sat.tst';
data = dlmread(fn, ' ');

ts = prdataset(data(:, 1:end-1), data(:, end));
ts = setlablist(ts, cl);
ts = remclass(ts);
ts = setfeatlab(ts, num2str([1:size(ts, 2)]'));
ts = setprior(ts, 0);

data = [tr; ts];

%myPreferredFeatures = 5*4+1:5*4+4;
myPreferredFeatures = [];
myNumFeatures = 4;
myClusterThresh = 0.2;
save('C:\Data\Development\Projects\MSc GeoInformatics\Data\Feature Selection\UCI\Statlog (Landsat Satellite)\Statlog Landsat.mat', ...
    'tr', 'ts', 'data', 'myPreferredFeatures', 'myNumFeatures', 'myClusterThresh');
%% Figure out the best params (see above)

method = FeatSelClusterRankM([], naivebc, 0, [], 'preferredFeatures', ...
        myPreferredFeatures, 'clusterThresh', myClusterThresh, 'showFigures', false, 'jmiFormulation', false);
w = tr*method;

res = BootstrapFsEval(data, method, 'numBootStraps', 10, 'numFeatures', myNumFeatures);

%% num features
[tr, ts] = gendat(data, 0.6);
feats = feast('jmi', 25, +tr, getnlab(tr));

for f = 3:25
    w = tr(:, feats(1:f))*(scalem([], 'variance')*opencvknnc([], 3));
    resF(f) = ts(:, feats(1:f))*w*testc
end
figure
plot(resF, 'k-x')



%% Urban land cover
% Contains training and testing data for classifying a high resolution
% aerial image into 9 types of urban land cover. Multi-scale spectral,
% size, shape, and texture information are used for classification. There
% are a low number of training samples for each class (14-30) and a high
% number of classification variables (148), so it may be an interesting
% data set for testing feature selection methods. The testing data set is
% from a random sampling of the image.
% 
% Class is the target classification variable. The land cover classes are:
% trees, grass, soil, concrete, asphalt, buildings, cars, pools, shadows.
clear all;
filenames = {'C:\Data\Development\Projects\MSc GeoInformatics\Data\Feature Selection\UCI\Urban Land Cover\training.csv',...
    'C:\Data\Development\Projects\MSc GeoInformatics\Data\Feature Selection\UCI\Urban Land Cover\testing.csv'};
delimiter = ',';
startRow = 2;
formatSpec = '%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
for i = 1:length(filenames)
    fileID = fopen(filenames{i},'r');
    dataArray{i} = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
    fclose(fileID);
    fileID = fopen(filenames{i},'r');
    header = textscan(fileID, '%s', length(dataArray{i})-1, 'Delimiter', delimiter);
    fclose(fileID);    
end
fl = header{1}(2:end);
tr = prdataset([dataArray{1}{2:end-1}], dataArray{1}{1});
tr = setfeatlab(tr, header{1}(2:end));
tr = setprior(tr, 0);
ts = prdataset([dataArray{2}{2:end-1}], dataArray{2}{1});
ts = setfeatlab(ts, header{1}(2:end));
ts = setprior(ts, 0);
data = [tr; ts];
%     19
%     29
%     75
%      8
%     30
myPreferredFeatures = [19, 29, 75, 8]; %from below
myNumFeatures = 4;
myClusterThresh = 0.2; %0.3; %25;

save('C:\Data\Development\Projects\MSc GeoInformatics\Data\Feature Selection\UCI\Urban Land Cover\Urban Land Cover.mat', ...
    'tr', 'ts', 'data', 'myPreferredFeatures', 'myNumFeatures', 'myClusterThresh');
%% Figure out the best params (see above)

method = FeatSelClusterRankM([], 'mi', 0, [], 'preferredFeatures', ...
        [], 'clusterThresh', myClusterThresh, 'showFigures', false, 'jmiFormulation', false);
w = tr*method;

for i = 1:10
    res{i} = BootstrapFsEval(data, method, 'numBootStraps', 1, 'numFeatures', i);
end
clfAcc = [];

% see what best num features is
for i = 1:10
    clfAcc(i,:) = res{i}.ClfMeanAcc;
end
figure
plot(clfAcc)

% see what preferredFeatures is
res = BootstrapFsEval(data, FeatSelClusterRankM([], naivebc, 0, [], 'preferredFeatures', ...
        [], 'clusterThresh', 0.3, 'showFigures', true, 'jmiFormulation', false), 'numBootStraps', 1, 'numFeatures', 5);

selFeats = unique(res.FeatIdx(:))
NselFeats = hist(res.FeatIdx(:), selFeats)
[dum sortedFeats] = sort(-NselFeats)
myPreferredFeatures = selFeats(sortedFeats(1:5))

res = BootstrapFsEval(data, FeatSelClusterRankM([], naivebc, 0, [], 'preferredFeatures', ...
        myPreferredFeatures, 'clusterThresh', myClusterThresh, 'showFigures', false, 'jmiFormulation', false), ...
        'numBootStraps', 10, 'numFeatures', myNumFeatures);

%% num features
[tr, ts] = gendat(data, 0.5);
feats = feast('jmi', 25, +tr, getnlab(tr));
feats = +featself(tr, naivebc, 15, 5)

for f = 1:length(feats)
    w = tr(:, feats(1:f))*(scalem([], 'variance')*opencvknnc([], 3));
    resF(f) = ts(:, feats(1:f))*w*testc
end
figure
plot(resF, 'k-x')

    
    
%% make synthetic data
%  gendat      - Random sampling of datasets for training and testing
%  gensubsets  - Generation of a consistent series of subsets of a dataset
%  gendatgauss - Generation of multivariate Gaussian distributed data
%  gendatb     - Generation of banana shaped classes
%  gendatc     - Generation of circular classes
%  gendatd     - Generation of two difficult classes
%  gendath     - Generation of Highleyman classes
%  gendati     - Generation of random windows from images
%  gendatk     - Nearest neighbour data generation
%  gendatl     - Generation of Lithuanian classes
%  gendatm     - Generation of 8 2d classes
%  gendatmm    - Generation of 4 multi-modal 2d classes
%  gendatp     - Parzen density data generation
%  gendatr     - Generate regression dataset from data and target values
%  gendats     - Generation of two Gaussian distributed classes
%  gendatw     - Sample dataset by given weigths
%  gendatv     - Generation of a very large dataset
%  gentrunk    - Generation of Trunk's example
%  genmdat     - Generation of a multi-dimensional dataset
close all; clear my*;
data = gendatv(10000, 2, 5);
% scatterdui(data)
datan = +data + randn(size(data)) * 0.25;
idx = [1:5];
datan2 = +data(:, idx) + randn(size(data(:, idx))) * 0.5;
datas = randn(size(data, 1), 2);
data = [data, datan, datan2, datas];
% scatterdui(data)
data = setfeatlab(data, num2str([1:size(data, 2)]'));

myNumFeatures = 5;  % the last one is spurious
myPreferredFeatures = 1:5;
myClusterThresh = 0.15;

save('C:\Data\Development\Projects\MSc GeoInformatics\Data\Feature Selection\Synthetic.mat', ...
    'data', 'myPreferredFeatures', 'myNumFeatures', 'myClusterThresh');

%% figure out params

res = BootstrapFsEval(data, FeatSelClusterRankM([], naivebc, 0, [], 'preferredFeatures', ...
        myPreferredFeatures, 'clusterThresh', myClusterThresh, 'showFigures', false, 'jmiFormulation', false), ...
        'numBootStraps', 10, 'numFeatures', myNumFeatures);

%% num features
[tr, ts] = gendat(data, 0.5);
feats = feast('jmi', 17, +tr, getnlab(tr));
feats = +featself(tr, naivebc, 17, 5)

for f = 1:length(feats)
    w = tr(:, feats(1:f))*(scalem([], 'variance')*opencvknnc([], 3));
    resF(f) = ts(:, feats(1:f))*w*testc
end
figure
plot(resF, 'k-x')

    
    
%% Scratch patch
clear all;
load 'C:\Data\Development\Projects\MSc GeoInformatics\Data\Feature Selection\UCI\Statlog (Landsat Satellite)\Statlog Landsat.mat'
load 'C:\Data\Development\Projects\MSc GeoInformatics\Data\Feature Selection\UCI\Urban Land Cover\Urban Land Cover.mat';
load 'C:\Data\Development\Projects\MSc GeoInformatics\Data\Feature Selection\UCI\Forest Types\Forest Types.mat';
load 'C:\Data\Development\Projects\MSc GeoInformatics\Data\Feature Selection\UCI\Cover Type\Cover Type.mat';

data = remclass(data);
data = setprior(data, 0);
cs = classsizes(data);
subData = gendat(data, min(cs)*ones(1, length(cs)));

resj = BootstrapFsEval(subData(:, 1:10), FeatSelFeastM([], 'jmi', 0), 'numBootStraps', 10, 'numFeatures', 5);
resfcr = BootstrapFsEval(subData(:, 1:10), FeatSelClusterRankM([], naivebc, 5, [], 'clusterThresh', 0.25, ...
    'showFigures', true), 'numBootStraps', 10, 'numFeatures', 5);

resfs = BootstrapFsEval(data, featself([], naivebc, 5), 'numBootStraps', 10, 'numFeatures', 5);

% resf = FeatureClusterRank([tr;ts(1:400,:)], 'clusterThresh', 0.3, 'showFigures', true, 'criterion', naivebc);

%%
load 'C:\Data\Development\Projects\MSc GeoInformatics\Data\Feature Selection\UCI\Statlog (Landsat Satellite)\Statlog Landsat.mat'

data = remclass(data);
data = setprior(data, 0);
data = setfeatlab(data, num2str([1:size(data, 2)]'));

cs = classsizes(data);
subData = gendat(data, min(cs)*ones(1, length(cs)));
fl = cellstr(getfeatlab(data));

tr_ = data(1:size(tr, 1), :);
resj = BootstrapFsEval(tr_, FeatSelFeastM([], 'jmi', 0), 'numBootStraps', 10, 'numFeatures', 5);
resfcr = BootstrapFsEval(tr_, FeatSelClusterRankM([], 'nmi', 0, [], 'clusterThresh', 0.11, ...
    'showFigures', false, 'jmiFormulation', true), 'numBootStraps', 10, 'numFeatures', 5);

%% does/can new jmi formulation select my original features?
% no - it does not...
clear all; close all;
%load('F:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\DataAllWin5NoBorder2.mat');
load('C:\Data\Development\Projects\MSc GeoInformatics\Data\Feature Selection\DataAllWin5NoBorder2.mat')
dataAll = changelablist(dataAll, 'Default');
fl = cellstr(getfeatlab(dataAll));
%remove LBP
idx = strmatch('Lbp', fl);
dataAll(:,idx)=[];

close all hidden
cs = classsizes(dataAll);
%cs(1) = cs(2);
cs(1:2) = cs(3);
randreset;
subData = gendat(gendat(dataAll), cs); %if N spec'd no sample with repl ??
subData = setprior(subData, 0);
fl = cellstr(getfeatlab(dataAll));

preferredFeatures = [9 5:8 10 1:4 19:22 15:18];

resfcr = BootstrapFsEval(subData, FeatSelClusterRankM([], naivebc, 0, [], 'preferredFeatures', ...
    [], 'clusterThresh', 0.175, 'showFigures', false, 'jmiFormulation', false), ...
    'numBootStraps', 10, 'numFeatures', 6);
resfcr2 = BootstrapFsEval(subData, FeatSelClusterRankM([], 'nmi', 0, [], 'preferredFeatures', ...
    [], 'clusterThresh', 0.175, 'showFigures', false, 'jmiFormulation', false), ...
    'numBootStraps', 10, 'numFeatures', 6);
resfcr3 = BootstrapFsEval(subData, FeatSelClusterRankM([], 'nmi', 0, [], 'preferredFeatures', ...
    [], 'clusterThresh', 0.175, 'showFigures', false, 'jmiFormulation', true), ...
    'numBootStraps', 10, 'numFeatures', 6);
resj = BootstrapFsEval(subData, FeatSelFeastM([], 'jmi', 0), 'numBootStraps', 10, 'numFeatures', 6);

