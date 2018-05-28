%% compare FCR and std approaches on different data sets
% first load / create the data sets to test on

%% Spekboom data
close all hidden; clear all;
load('D:\Data\Development\Projects\PhD GeoInformatics\Data\Feature Selection\DataAllWin5NoBorder2.mat')
dataAll = changelablist(dataAll, 'Default');
dataAll = setprior(dataAll, 0);
fl = cellstr(getfeatlab(dataAll));
idx = strmatch('Lbp', fl);
dataAll(:, idx)=[];
fl = cellstr(getfeatlab(dataAll));

feats = [9 15 20 23 7 6]; %ranked cluster
fl(feats)

cs = classsizes(dataAll);
if true
    cs(1) = cs(2);
else
    cs = min(cs)*ones(1,3);
end
% randreset(2);
randreset
subData = gendat(dataAll, cs);
subData = changelablist(subData, 'Default');
subData = setprior(subData, 0);

cdata{1} = subData;
numFeatures = 6;
cdataNames = {'Spekboom'};
preferredFeatures = {[9 15 20 23 7 6]};
clusterThresh = 0.2;

%% Synthetic data
clear my* data* tr ts subData
load('D:\Data\Development\Projects\PhD GeoInformatics\Data\Feature Selection\Synthetic.mat');

cdata{end+1} = data;
numFeatures(end+1) = myNumFeatures;
cdataNames{end+1} = 'Synthetic';
preferredFeatures{end+1} = myPreferredFeatures;
clusterThresh(end+1) = myClusterThresh;

%% UCI data
close all; clear my*;
uciFiles = { 'D:\Data\Development\Projects\PhD GeoInformatics\Data\Feature Selection\UCI\Statlog (Landsat Satellite)\Statlog Landsat.mat'; ...
     'D:\Data\Development\Projects\PhD GeoInformatics\Data\Feature Selection\UCI\Urban Land Cover\Urban Land Cover.mat';...
%     'D:\Data\Development\Projects\PhD GeoInformatics\Data\Feature Selection\UCI\Forest Types\Forest Types.mat';...
%     'D:\Data\Development\Projects\PhD GeoInformatics\Data\Feature Selection\UCI\Cover Type\Cover Type.mat';...
    };
uciNames = {'Statlog Landsat', 'Urban Land Cover', 'Forest Types', 'Cover Type'};
for i = 1:length(uciFiles)
    load(uciFiles{i});
    data = remclass(data);
    data = setprior(data, 0);
    cs = classsizes(data);
    subData = gendat(data, min(cs)*ones(1, length(cs)));
    cdata{end+1} = subData;
    cdataNames{end+1} = uciNames{i};
    numFeatures(end+1) = myNumFeatures;
    preferredFeatures{end+1} = myPreferredFeatures;
    clusterThresh(end+1) = myClusterThresh;
end

%% Hyperspectral data
close all; clear my*;
hsFiles = { 'D:\Data\Development\Projects\PhD GeoInformatics\Data\Feature Selection\Hyperspectral\BotswanaPr.mat'; ...
     'D:\Data\Development\Projects\PhD GeoInformatics\Data\Feature Selection\Hyperspectral\KscPr.mat';...
    };
hsNames = {'Botswana', 'KSC'};
for i = 1:length(hsFiles)
    load(hsFiles{i});
    data = remclass(data);
    data = setprior(data, 0);
    cs = classsizes(data);
    subData = gendat(data, min(cs)*ones(1, length(cs)));
    cdata{end+1} = subData;
    cdataNames{end+1} = hsNames{i};
    numFeatures(end+1) = myNumFeatures;
    preferredFeatures{end+1} = myPreferredFeatures;
    clusterThresh(end+1) = myClusterThresh;
end


%%  Compare the methods
load 'D:\Data\Development\Projects\PhD GeoInformatics\Data\Feature Selection\CompareFsMethodsHs4.mat'
% cdata = cdata(1)
close all
clear my*
clear i m tr ts subData uci* hs* data feats cs idx fl

% clusterThresh?
% res = {};
% methods = {...
%         FeatSelClusterRankM([], naivebc, 0, [], 'preferredFeatures', ...
%         preferredFeatures{i}, 'clusterThresh', 0.2, 'showFigures', false, 'jmiFormulation', false);...
%         featseli([], naivebc, numFeatures(i));...
%         featseli([], 'mi', numFeatures(i));...
%         featself([], naivebc, numFeatures(i));...
%         featself([], 'mi', numFeatures(i));...
%         FeatSelFeastM([], 'jmi', 0);...
%         FeatSelClusterRankM([], 'mi', 0, [], 'preferredFeatures', ...
%             preferredFeatures{i}, 'clusterThresh', 0.2, 'showFigures', false, 'jmiFormulation', true);...
%         FeatSelClusterRankM([], 'mi', 0, [], 'preferredFeatures', ...
%             preferredFeatures{i}, 'clusterThresh', 0.2, 'showFigures', false, 'jmiFormulation', false);...
%         FeatSelClusterRankM([], 'mi', 0, [], 'preferredFeatures', ...
%             preferredFeatures{i}, 'clusterThresh', 0.2, 'showFigures', false, 'jmiFormulation', false, 'useCorrelation', false);...
%         FeatSelClusterRankM([], 'mi', 0, [], 'preferredFeatures', ...
%             preferredFeatures{i}, 'clusterThresh', 0.2, 'showFigures', false, 'jmiFormulation', true, 'useCorrelation', false);...
%     };
numMethods = 7;
numFeatures(5)=6;
%cdata = cdata([1,5])
res = cell(numMethods, length(cdata));
delete(gcp('nocreate'))
parpool(4)
% ms = [1,8,9];
% is = [4,6];
for m = 1:numMethods
%for mi = 1:3
%    m = ms(mi);
    parfor i = 1:length(cdata)
%    parfor ii = 1:2
%        i = is(ii);
%         fprintf('Method %d, Data %d -----------------------------------------------\n', m, i);
        fprintf('Method %d, Data %d -----------------------------------------------\n', m, i);
        innerClusterThresh = clusterThresh(i);
        innerNumFeatures = numFeatures(i);
        innerMethods = {...
             FeatSelClusterRankM([], naivebc, 0, [], 'clusterThresh', innerClusterThresh, 'showFigures', false, ...
                    'jmiFormulation', false, 'clusterMethod', 'heirarchical', 'preferredFeatures', preferredFeatures{i});...
             FeatSelClusterRankM([], 'mi', 0, [], 'clusterThresh', innerClusterThresh, 'showFigures', false, ...
                    'jmiFormulation', false, 'clusterMethod', 'heirarchical', 'preferredFeatures', preferredFeatures{i});...
             FeatSelClusterRankM([], naivebc([], 25), 0, [], 'clusterThresh', innerClusterThresh, 'showFigures', false, ...
                    'jmiFormulation', false, 'clusterMethod', 'heirarchical', 'preferredFeatures', preferredFeatures{i});...
             FeatSelClusterRankM([], naivebc([], 25), 0, [], 'clusterThresh', innerClusterThresh, 'showFigures', false, ...
                    'jmiFormulation', false, 'clusterMethod', 'ap', 'preferredFeatures', preferredFeatures{i});...
             FeatSelClusterRankM([], 'mi', 0, [], 'clusterThresh', innerClusterThresh, 'showFigures', false, ...
                    'jmiFormulation', false, 'clusterMethod', 'ap', 'preferredFeatures', preferredFeatures{i});...
             FeatSelClusterRankM([], naivebc([], 25), 0, [], 'clusterThresh', innerClusterThresh, 'showFigures', false, ...
                    'jmiFormulation', false, 'clusterMethod', 'ap');...
             FeatSelClusterRankM([], 'mi', 0, [], 'clusterThresh', innerClusterThresh, 'showFigures', false, ...
                    'jmiFormulation', true, 'clusterMethod', 'ap');...
%                 featseli([], naivebc, innerNumFeatures, 5);...
%                 featself([], naivebc, innerNumFeatures, 5);...
%                 featseli([], 'mi', innerNumFeatures);...
%                 featself([], 'mi', innerNumFeatures);...
%                 featself([], naivebc([], 30), innerNumFeatures);...
            };

%         innerMethods = {...
%             FeatSelClusterRankM([], naivebc, 0, [], 'preferredFeatures', ...
%                 preferredFeatures{i}, 'clusterThresh', innerClusterThresh, 'showFigures', true, 'jmiFormulation', false);...
%                 featseli([], naivebc, innerNumFeatures, 5);...
%                 featself([], naivebc, innerNumFeatures, 5);...
% %             featseli([], 'mi', innerNumFeatures);...
% %             featself([], 'mi', innerNumFeatures);...
% %             featself([], 'nmi', innerNumFeatures);...
% %             FeatSelFeastM([], 'jmi', 0);...
%             FeatSelClusterRankM([], 'mi', 0, [], 'preferredFeatures', ...
%                 preferredFeatures{i}, 'clusterThresh', innerClusterThresh, 'showFigures', false, 'jmiFormulation', false);...
% %             FeatSelClusterRankM([], 'mi', 0, [], 'preferredFeatures', ...
% %                 preferredFeatures{i}, 'clusterThresh', innerClusterThresh, 'showFigures', false, 'jmiFormulation', true);...
% % %             FeatSelClusterRankM([], 'mi', 0, [], 'preferredFeatures', ...
% % %                 preferredFeatures{i}, 'clusterThresh', innerClusterThresh, 'showFigures', true, 'jmiFormulation', false, 'useCorrelation', false);...
% % %             FeatSelClusterRankM([], 'mi', 0, [], 'preferredFeatures', ...
% % %                 preferredFeatures{i}, 'clusterThresh', innerClusterThresh, 'showFigures', true, 'jmiFormulation', true, 'useCorrelation', false);...
% %             featselb([], naivebc, 0);...
% %             featselb([], 'nmi', 0);...
%         };
% 
        fprintf('Data: %s, Method: %s, clusterThresh: %d, numFeatures: %d\n', ...
            cdataNames{i}, struct(innerMethods{m}).name, innerClusterThresh, innerNumFeatures);
        res{m, i} = BootstrapFsEval(cdata{i}, innerMethods{m}, 'numBootStraps', 10, ...
            'numFeatures', innerNumFeatures);
    end
end
delete(gcp('nocreate'))

%%
if false
    res_ = res; 
%     load 'D:\Data\Development\Projects\PhD GeoInformatics\Data\Feature Selection\CompareFsMethodsHs3.mat'
    load 'D:\Data\Development\Projects\PhD GeoInformatics\Data\Feature Selection\CompareFsMethodsHs4.mat'
    
    ms = [1,8,9];
    is = [4,6];

    for mi = 1:3
        m = ms(mi);
    %    parfor i = 1:length(cdata)
        for ii = 1:2
            i = is(ii);
            res{m,i} = res_{mi,ii}
        end
    end
end

methodMedianFsDuration = zeros(numMethods, length(cdata));
methodMedianStability = zeros(numMethods, length(cdata));
methodMedianAcc = zeros(numMethods, length(cdata));

methodNames = {'FCR-H-naivebc', 'FCR-H-mi', 'FCR-H-naivebc25', 'FCR-ap-naivebc25', 'FCR-ap-mi', 'FCR-ap-naivebc25-nopref', 'FCR-ap-mi-jmi-nopref'};
methodNames = {'featseli-naivebc', 'featself-naivebc', 'featseli-mi', 'featself-mi', 'featself-naivebc30'};


% methodNames = {'featself-naivebc',...
%     'FCR-mi-jmi off-corr on'};

% methodNames = {'FCR-naivebc-jmi off-corr on', 'featseli-naivebc', 'featself-naivebc',...
%     'featseli-mi', 'featself-mi', 'featself-nmi', 'JMI', 'FCR-mi-jmi off-corr on', 'FCR-mi-jmi on-corr on', ...
%     'FCR-mi-jmi off-corr off', 'FCR-mi-jmi on-corr off', 'featselb-naivebc', 'featselb-nmi'};
% methodNames = {'FCR-naivebc-jmi off-corr on', 'featseli-naivebc', 'featself-naivebc',...
%     'featseli-mi', 'featself-mi', 'featself-nmi', 'JMI', 'FCR-mi-jmi off-corr on', 'FCR-mi-jmi on-corr on', ...
%     'featselb-naivebc', 'featselb-nmi'};
% for i = 1:size(res, 2)
%     table = {'Method', 'Tanimoto', 'Consistency', 'SpearmanCC', 'FsDuration', 'ClfMeanAcc'};
%     for m = 1:size(res, 1) % use cluster index rather than feature index for FCR
%         table(m+1,:) = {methodNames{m}, res{m, i}.TanimotoStability, res{m, i}.Consitency, ...
%             res{m, i}.SpearmanRankCorrCoeffStab, mean(res{m, i}.FsDuration), res{m, i}.ClfMeanAcc(1)};
%         methodMedianStability(m, i) = res{m, i}.Consitency; %(res{m, i}.TanimotoStability + res{m, i}.Consitency)/2;
%         methodMedianAcc(m, i) = mean(res{m, i}.ClfMeanAcc(end));  %3nn accuracy only
%         methodMedianFsDuration(m, i) = mean(res{m, i}.FsDuration);
%     end
%     disp(cdataNames{i})
%     disp(table)
% end
% 


for i = 1:size(res, 2)
    table = {'Method', 'Tanimoto', 'Consistency', 'SpearmanCC', 'FsDuration', 'ClfMeanAcc'}; %, 'ClfMeanAcc', 'ClfMeanAcc', 'ClfMeanAcc'};
    for m = 1:size(res, 1) % use cluster index rather than feature index for FCR
        table(m+1,:) = {methodNames{m}, res{m, i}.TanimotoStability, res{m, i}.Consitency, ...
            res{m, i}.SpearmanRankCorrCoeffStab, mean(res{m, i}.FsDuration), res{m, i}.ClfMeanAcc(end)}; %, res{m, i}.ClfMeanAcc(2), res{m, i}.ClfMeanAcc(3), res{m, i}.ClfMeanAcc(4)};
        methodMedianStability(m, i) = res{m, i}.Consitency; %(res{m, i}.TanimotoStability + res{m, i}.Consitency)/2;
        methodMedianAcc(m, i) = mean(res{m, i}.ClfMeanAcc(end));  %3nn accuracy only
        methodMedianFsDuration(m, i) = mean(res{m, i}.FsDuration);
    end
    disp(cdataNames{i})
    disp(table)
end

% for i = 1:size(res, 2)
%     table = {'Method', 'Tanimoto', 'Consistency', 'SpearmanCC', 'FsDuration', 'ClfMeanAcc', 'ClfMeanAcc', 'ClfMeanAcc', 'ClfMeanAcc'};
%     for m = 1:size(res, 1) 
%         table(m+1,:) = {methodNames{m}, res{m, i}.TanimotoStability, res{m, i}.Consitency, ...
%             res{m, i}.SpearmanRankCorrCoeffStab, mean(res{m, i}.FsDuration), res{m, i}.ClfMeanAcc(1), res{m, i}.ClfMeanAcc(2), res{m, i}.ClfMeanAcc(3), res{m, i}.ClfMeanAcc(4)};
%         methodMedianStability(m, i) = (res{m, i}.TanimotoStability + res{m, i}.Consitency)/2;
%         methodMedianAcc(m, i) = mean(res{m, i}.ClfMeanAcc);
%         methodMedianFsDuration(m, i) = mean(res{m, i}.FsDuration);
%     end
%     disp(cdataNames{i})
%     disp(table)
% end

table = {'Method', 'Stability', 'Accuracy', 'FsDuration', 'Overall'};
for m = 1:size(res, 1) % use cluster index rather than feature index for FCR
    table(m+1,:) = {methodNames{m}, mean(methodMedianStability(m,:)),...
        mean(methodMedianAcc(m, :)), mean(methodMedianFsDuration(m,:)), mean(methodMedianAcc(m, :).*methodMedianStability(m, :))};
end
disp(table)

%%NNB create innerMethods above
if false
    % do it again but use cluster index instead of feature index for all FCRs
    resFci = res;
    for i = 1:size(res, 2)
        table = {'Method', 'Tanimoto', 'Consistency', 'SpearmanCC', 'Duration', 'ClfMeanAcc'}; % , 'ClfMeanAcc', 'ClfMeanAcc', 'ClfMeanAcc'};
        for m = 1:size(res, 1) % use cluster index rather than feature index for FCR
            % re-eval stability using cluster indices
            if (strcmpi(innerMethods{m}.name, 'Feature Clustering and Ranking'))
                %DH 2018 - NNB - RenumClustAcrossBootstraps must only be called
                %once!!! and it is already called in BootstrapFsEval
                %resFci{m, i} = RenumClustAcrossBootstraps(resFci{m, i});
                featIdx = resFci{m, i}.FeatIdx;
                try
                    for b = 1:size(resFci{m, i}.FeatIdx, 2)
                        featClustIdx = resFci{m, i}.FeatClustIdx(:, b);
                        resFci{m, i}.FeatIdx(:, b) = featClustIdx(resFci{m, i}.FeatIdx(:, b));
                    end
                    resFci{m, i} = FsStabilityEval(resFci{m, i});
                catch ex
                    disp(ex.message);
                end
    % DH 2018 - the below needs attention - we are re-assigning the orig
    % FeatIdx after FsStabilityEval
                if true
                    resFci{m, i}.FeatIdx_ = resFci{m, i}.FeatIdx;
                    resFci{m, i}.FeatIdx = featIdx;
                else
                    resFci{m, i}.FeatIdxOrig = featIdx;
                end
            end
            table(m+1,:) = {methodNames{m}, resFci{m, i}.TanimotoStability, resFci{m, i}.Consitency, ...
                resFci{m, i}.SpearmanRankCorrCoeffStab, mean(res{m, i}.FsDuration), resFci{m, i}.ClfMeanAcc(1)}; %, resFci{m, i}.ClfMeanAcc(2), resFci{m, i}.ClfMeanAcc(3), resFci{m, i}.ClfMeanAcc(4)};
            methodMedianStability(m, i) = resFci{m, i}.Consitency; %(resFci{m, i}.TanimotoStability + resFci{m, i}.Consitency)/2;
            methodMedianAcc(m, i) = mean(resFci{m, i}.ClfMeanAcc(end));  % 3nn only
            methodMedianFsDuration(m, i) = mean(res{m, i}.FsDuration);
        end
        disp(cdataNames{i})
        disp(table)
    end
else       % Old Oct 2016 code used to gen paper results in CompareFsMethodsHs4.mat
    resFci = res;
    for i = 1:size(res, 2)
        table = {'Method', 'Tanimoto', 'Consistency', 'SpearmanCC', 'Duration', 'ClfMeanAcc'}; %, 'ClfMeanAcc', 'ClfMeanAcc', 'ClfMeanAcc'};
        for m = 1:size(res, 1) % use cluster index rather than feature index for FCR
            % re-eval stability using cluster indices
            if (strcmpi(innerMethods{m}.name, 'Feature Clustering and Ranking'))
                featIdx = resFci{m, i}.FeatIdx;
                try
                    resFci{m, i}.FeatIdx = resFci{m, i}.FeatClustIdx(resFci{m, i}.FeatIdx);
                    resFci{m, i} = FsStabilityEval(resFci{m, i});
                catch ex
                    disp(ex.message);
                    rethrow(ex)
                end
                resFci{m, i}.FeatIdx = featIdx;
            end
            table(m+1,:) = {methodNames{m}, resFci{m, i}.TanimotoStability, resFci{m, i}.Consitency, ...
                resFci{m, i}.SpearmanRankCorrCoeffStab, mean(res{m, i}.FsDuration), resFci{m, i}.ClfMeanAcc(end)};            
%                 resFci{m, i}.SpearmanRankCorrCoeffStab, mean(res{m, i}.FsDuration), resFci{m, i}.ClfMeanAcc(1), resFci{m, i}.ClfMeanAcc(2), resFci{m, i}.ClfMeanAcc(3), resFci{m, i}.ClfMeanAcc(4)};            methodMedianStability(m, i) = (resFci{m, i}.TanimotoStability + resFci{m, i}.Consitency)/2;
            methodMedianStability(m, i) = resFci{m, i}.Consitency; %(resFci{m, i}.TanimotoStability + resFci{m, i}.Consitency)/2;
            methodMedianAcc(m, i) = mean(resFci{m, i}.ClfMeanAcc);
            methodMedianFsDuration(m, i) = mean(res{m, i}.FsDuration);
        end
        disp(cdataNames{i})
        disp(table)
    end
end

table = {'Method', 'Stability', 'Accuracy', 'FsDuration', 'Overall'};
for m = 1:size(res, 1) % use cluster index rather than feature index for FCR
    table(m+1,:) = {methodNames{m}, mean(methodMedianStability(m,:)),...
        mean(methodMedianAcc(m, :)), mean(methodMedianFsDuration(m,:)), mean(methodMedianAcc(m, :).*methodMedianStability(m, :))};
end
disp(table)

table = {'Data', 'Most Stable', 'Most Accurate', 'Fastest', 'Best Overall'};
for i = 1:size(res, 2) % use cluster index rather than feature index for FCR
    [ms msi] = max(methodMedianStability(:, i));
    [ma mai] = max(methodMedianAcc(:, i));
    [ma mfi] = min(methodMedianFsDuration(:, i));
    [mo moi] = max(methodMedianStability(:, i).*methodMedianAcc(:, i));
    table(i+1,:) = {cdataNames{i}, methodNames{msi}, methodNames{mai}, methodNames{mfi}, methodNames{moi}};
end
disp(table)

[o oi] = sort(-mean(methodMedianStability.*methodMedianAcc, 2))
% methods ranked overall
disp('Methods ranked overall:'); 
disp(methodNames(oi)')

[a ai] = sort(-mean(methodMedianAcc, 2));
disp('Methods ranked by accuracy:'); 
disp(methodNames(ai)')

[s si] = sort(-mean(methodMedianStability, 2));
disp('Methods ranked by stability:'); 
disp(methodNames(si)')

[s si] = sort(-mean(methodMedianFsDuration, 2));
disp('Methods ranked by speed:'); 
disp(methodNames(si)')

%%
% resFeatIdx = res; % FCR res using preferred feats and featIdx for stability calcs 
% save 'D:\Data\Development\Projects\PhD GeoInformatics\Data\Feature Selection\resFeatIdx.mat' resFeatIdx res
save 'D:\Data\Development\Projects\PhD GeoInformatics\Data\Feature Selection\resHs4.mat' res resFci
save 'D:\Data\Development\Projects\PhD GeoInformatics\Data\Feature Selection\CompareFsMethodsHs4.mat'

save 'D:\Data\Development\Projects\PhD GeoInformatics\Data\Feature Selection\CompareFsMethodsApNoPrefFeat.mat'
save 'D:\Data\Development\Projects\PhD GeoInformatics\Data\Feature Selection\CompareFsMethodsFCR_MVFS.mat'

save 'D:\Data\Development\Projects\PhD GeoInformatics\Data\Feature Selection\CompareFsMethodsFCR_AP.mat'

%% Make eg dendgrogram (for paper)
close all;clear all;

load('D:\Data\Development\Projects\PhD GeoInformatics\Data\Feature Selection\Synthetic.mat');

res = BootstrapFsEval(data, FeatSelClusterRankM([], naivebc, 0, [], 'preferredFeatures', ...
        myPreferredFeatures, 'clusterThresh', myClusterThresh, 'showFigures', true, 'jmiFormulation', false), ...
        'numBootStraps', 1, 'numFeatures', myNumFeatures);

figure(2)
print('C:\Data\Development\Projects\PhD GeoInformatics\Docs\My Docs\Thesis\Feature Clustering and Ranking\Figure 1 - Example dendrogram.eps', ...
    '-deps', '-r1200')

%% Visualise results (for paper)
% NB remember to replace feat idx with cluster idx for FCR
clear all; close all;

methodNames = {'FCR-naivebc-jmi off-corr on', 'featseli-naivebc', 'featself-naivebc',...
    'featseli-mi', 'featself-mi', 'featself-nmi', 'JMI', 'FCR-mi-jmi off-corr on', 'FCR-mi-jmi on-corr on', ...
    'featselb-naivebc', 'featselb-nmi'};

methodNamesAbbr = methodNames

methodNamesAbbr = {'FCR-NaiveBC', 'Rank-NaiveBC', 'FS-NaiveBC',...
    'Rank-MI', 'FS-MI', 'FS-MI', 'JMI', 'FCR-MI', 'FCR-MI-JMI', ...
    'BE-NaiveBC', 'BE-MI'};
cdataNamesAbbr = {'Spekboom', 'Synthetic', 'Landsat', 'Urban', 'Botswana', 'KSC'};

    %'FCR-mi-jmi off-corr off', 'FCR-mi-jmi on-corr off', 'featselb-naivebc', 'featselb-nmi'};
load 'D:\Data\Development\Projects\PhD GeoInformatics\Data\Feature Selection\CompareFsMethodsHs4.mat'
%methodIdx = [7, 1, 2, 3, 12, 8, 4, 5, 13];
methodIdx = 1:length(methodNames);
methodIdx = [7, 1, 2, 3, 10, 8, 4, 5, 11];

methodNamesAbbr(methodIdx)'
res = resFci;  %%NNB

knnAcc = cellfun(@(x) x.ClfMeanAcc(end), res, 'UniformOutput', true);
clfMeanAcc = cellfun(@(x) mean(x.ClfMeanAcc), res, 'UniformOutput', true);
consistency = cellfun(@(x) x.Consitency, res, 'UniformOutput', true);
stability = cellfun(@(x) x.TanimotoStability, res, 'UniformOutput', true);
duration = cellfun(@(x) mean(x.FsDuration), res, 'UniformOutput', true);

overall = knnAcc.*consistency;
knnAccMean = mean(knnAcc, 2);
consistencyMean = mean(consistency, 2);
overallMean = mean(overall, 2);
clfMeanAccMean = mean(clfMeanAcc, 2);
[dum knnAccIdx] = sort(knnAccMean(methodIdx));
[dum clfrMeanAccIdx] = sort(clfMeanAccMean(methodIdx));
[dum consistencyIdx] = sort(consistencyMean(methodIdx));
[dum overallIdx] = sort(overallMean(methodIdx));

% consistencyIdx(methodIdx)
%find data set difficulties\
for i = 1:length(cdata)
    s = size(cdata{i});
    s(3) = getsize(cdata{i}, 3);
    df(i) = s(1)/(s(2)*s(3));
end
disp('Difficulty ratio')
disp([cdataNamesAbbr' num2cell(df(:))])

durationTtl = sum(duration, 2);
[durationTtlSort, durationIndex] = sort(durationTtl(methodIdx));
disp('Cumularive execution time')
disp([methodNamesAbbr(methodIdx(durationIndex))' num2cell(durationTtl(methodIdx(durationIndex)))])

figure;
bar(consistency(methodIdx(consistencyIdx), :))
h = legend(cdataNamesAbbr, 'Location', 'NorthEastOutside');
v = get(h, 'title');
set(v, 'string','Data Sets');
grid on;
ax = gca;
ax.XTick = 1:length(methodIdx);
ax.XTickLabel = methodNamesAbbr(methodIdx(consistencyIdx));
ax.XTickLabelRotation = 45;
ax.XLim(2) = ax.XLim(2);  %make space for legend
xlabel('Method', 'FontWeight', 'Bold')
ylabel('Stability', 'FontWeight', 'Bold')
%title('Consistency by Method')
fontsize(9)
print('C:\Data\Development\Projects\PhD GeoInformatics\Docs\My Docs\Thesis\Feature Clustering and Ranking\Figure 2 - Method stability per data set.eps', ...
    '-depsc', '-r1200')


figure;
bar(consistency(methodIdx(consistencyIdx), :)')
legend(methodNamesAbbr(methodIdx(consistencyIdx)), 'Location', 'NorthEastOutside')
grid on;
ax = gca;
ax.XTick = 1:length(cdataNamesAbbr);
ax.XTickLabel = cdataNamesAbbr;
ax.XTickLabelRotation = 45;
xlabel('Data')
ylabel('Stability')

figure;
bar(100*knnAcc(methodIdx(knnAccIdx), :))
h = legend(cdataNamesAbbr, 'Location', 'NorthEastOutside');
v = get(h, 'title');
set(v, 'string','Data Sets');
grid on
ax = gca;
ax.XTick = 1:length(methodIdx);
ax.XTickLabel = methodNamesAbbr(methodIdx(knnAccIdx));
ax.XTickLabelRotation = 45;
xlabel('Method', 'FontWeight', 'Bold')
ylabel('Accuracy (%)', 'FontWeight', 'Bold')
%title('Consistency by Method')
fontsize(9)
print('C:\Data\Development\Projects\PhD GeoInformatics\Docs\My Docs\Thesis\Feature Clustering and Ranking\Figure 3 - Method accuracy per data set.eps', ...
    '-depsc', '-r1200')

figure;
bar(knnAcc(methodIdx(knnAccIdx), :)')
legend(methodNamesAbbr(methodIdx(knnAccIdx)), 'Location', 'NorthEastOutside')
grid on
ax = gca;
ax.XTick = 1:length(cdataNamesAbbr);
ax.XTickLabel = cdataNamesAbbr;
ax.XTickLabelRotation = 45;
xlabel('Data')
ylabel('Accuracy')

figure;
bar(clfMeanAcc(methodIdx(clfrMeanAccIdx), :))
legend(cdataNamesAbbr, 'Location', 'NorthEastOutside')
grid on
ax = gca;
ax.XTick = 1:length(methodIdx);
ax.XTickLabel = methodNamesAbbr(methodIdx(clfrMeanAccIdx));
ax.XTickLabelRotation = 45;
xlabel('Method')
ylabel('Accuracy')
title('Mean Clfr Acc by Method')

% figure;
% plot(clfMeanAcc', 'x-')
% legend(methodNamesAbbr, 'Location', 'NorthEastOutside');
% ax = gca;
% ax.XTick = 1:size(consistency, 2);
% ax.XTickLabel = cdataNamesAbbr;
% ax.XTickLabelRotation = 45;
% title('Clf Acc by Data');

figure;
bar(overall(methodIdx(overallIdx), :))
legend(cdataNamesAbbr, 'Location', 'NorthEastOutside');
grid on
ax = gca;
ax.XTick = 1:length(methodIdx);
ax.XTickLabel = methodNamesAbbr(methodIdx(overallIdx));
ax.XTickLabelRotation = 45;
xlabel('Method')
ylabel('Accuracy x Stability')

figure;
bar(overall(methodIdx(overallIdx), :)')
legend(methodNamesAbbr(methodIdx(overallIdx)), 'Location', 'NorthEastOutside');
grid on
ax = gca;
ax.XTick = 1:length(cdataNamesAbbr);
ax.XTickLabel = cdataNamesAbbr;
ax.XTickLabelRotation = 45;
xlabel('Data','FontWeight', 'bold')
ylabel('Accuracy x Stability','FontWeight', 'bold')

% figure;
% plot(overall(methodIdx, :)', 'x-')
% legend(methodNamesAbbr);
% ax = gca;
% ax.XTick = 1:size(consistency, 2);
% ax.XTickLabel = cdataNamesAbbr;
% ax.XTickLabelRotation = 45;
% title('Overall by Data');


%%
% get data for all bootstraps
close all;
knnAccBs = [];
for m = 1:size(res, 1)
    rr = [];
    ss = [];
    for d = 1:size(res, 2)
        rr = [rr, res{m, d}.ClfAcc(:, end)'];
    end
    knnAccBs(m, :) = rr;
end
overallBs = (consistency.*knnAcc);

%find non dominant rank of each method (includes both acc and stability) as
%in Brown et al
for i = 1:size(consistency, 2)
    R(:, i) = NonDominantRank(-[knnAcc(methodIdx, i) consistency(methodIdx, i)]);
end
Rm = mean(R, 2);
[Rs Ridx] = sort(Rm);
fprintf('Combined accuracy and stability ranking\n')
disp([methodNamesAbbr(methodIdx(Ridx))' num2cell(Rm(Ridx))])
disp([methodNamesAbbr(methodIdx(Ridx))' num2cell(R(Ridx,:))])

for i = 1:size(consistency, 2)
    R(:, i) = NonDominantRank(-[knnAcc(methodIdx, i)]);
end
Rm = mean(R, 2);
[Rs Ridx] = sort(Rm);
fprintf('Accuracy ranking\n')
disp([methodNamesAbbr(methodIdx(Ridx))' num2cell(Rm(Ridx))])

for i = 1:size(consistency, 2)
    R(:, i) = NonDominantRank(-[consistency(methodIdx, i)]);
end
Rm = mean(R, 2);
[Rs Ridx] = sort(Rm);
fprintf('Stability ranking\n')
disp([methodNamesAbbr(methodIdx(Ridx))' num2cell(Rm(Ridx))])


%plot consistency vs accuracy
figure
icons = {'x','o','+','s','*','^','<','>','.'};
for i = 1:size(consistency(methodIdx, :), 2)
    subplot(2, 3, i)
    for j = 1:size(consistency(methodIdx, :), 1)
        text(consistency(methodIdx, i)', knnAcc(methodIdx, i)', methodNamesAbbr(methodIdx), 'FontSize', 8);
    end
    axis([min(consistency(methodIdx, i)) max(consistency(methodIdx, i)) ...
        min(knnAcc(methodIdx, i)) max(knnAcc(methodIdx, i))]) 
    grid on
    title(cdataNamesAbbr(i))
    xlabel('Stability')
    ylabel('Accuracy')
end


% test for significance of non-dominant ranks
for i = 1:size(consistency, 2)
    R(:, i) = NonDominantRank(-[knnAcc(methodIdx, i) consistency(methodIdx, i)]);
end
Rm = mean(R, 2);
[Rs Ridx] = sort(Rm);
fprintf('Combined accuracy and stability ranking\n')
disp([methodNamesAbbr(methodIdx(Ridx))' num2cell(Rm(Ridx))])
disp([methodNamesAbbr(methodIdx(Ridx))' num2cell(R(Ridx,:))])

[p,tbl,stats] = friedman(R(Ridx, :)', 1);
figure
multcompare(stats);
ax = gca;
ax.YTick = 1:length(methodNamesAbbr((methodIdx(Ridx))));
ax.YTickLabel = flipud(methodNamesAbbr((methodIdx(Ridx)))');
title('Non Dominant Ranks');
%nemenyi test as in brown et al
clear nemenyi
RR = [R R]*0; %hack to see effect of repeats and more data sets
RR(:, 1:2:11) = R;
RR(:, 2:2:12) = R+randn(size(R));
RR = R;

figure
[p, n, meanrank, CDa, rankmean] = nemenyi(RR(Ridx, :)', 1, 'labels', ...
    methodNamesAbbr((methodIdx(Ridx))), 'alpha', 0.05, 'ploton', 'odetail');
figure
[p, n, meanrank, CDa, rankmean] = nemenyi(R(Ridx, :)', 1, 'labels', ...
    methodNamesAbbr((methodIdx(Ridx))), 'alpha', 0.05, 'ploton', 'voline');


[p,tbl,stats] = friedman(1-knnAcc');
figure
multcompare(stats);
ax = gca;
ax.YTick = 1:length(methodNamesAbbr);
ax.YTickLabel = flipud(methodNamesAbbr');
title('3NN Acc');
figure
[p, n, meanrank, CDa, rankmean] = nemenyi(1-knnAcc', 1, 'labels', ...
    methodNamesAbbr, 'alpha', 0.05, 'ploton', 'odetail');


[p,tbl,stats] = friedman(1-consistency');
figure
multcompare(stats);
ax = gca;
ax.YTick = 1:length(methodNamesAbbr);
ax.YTickLabel = flipud(methodNamesAbbr');
title('Consistency');
figure
[p, n, meanrank, CDa, rankmean] = nemenyi(1-consistency', 1, 'labels', ...
    methodNamesAbbr, 'alpha', 0.05, 'ploton', 'odetail');

[p,tbl,stats] = friedman(1-overallBs');
multcompare(stats);
ax = gca;
ax.YTick = 1:length(methodNamesAbbr);
ax.YTickLabel = flipud(methodNamesAbbr');
title('Overall');


figure;
boxplot(knnAcc(methodIdx(knnAccIdx), :)')
ax = gca;
ax.XTick = 1:size(knnAcc, 1);
ax.XTickLabel = methodNamesAbbr(methodIdx(knnAccIdx));
ax.XTickLabelRotation = 45;
title('Acc over data')
grid on

figure;
boxplot(knnAccBs(methodIdx(knnAccIdx), :)')
ax = gca;
ax.XTick = 1:size(knnAccBs, 1);
ax.XTickLabel = methodNamesAbbr(methodIdx(knnAccIdx));
ax.XTickLabelRotation = 45;
title('Acc over data and bootstraps')
grid on

figure;
boxplot(consistency')
ax = gca;
ax.XTick = 1:size(consistency, 1);
ax.XTickLabel = methodNamesAbbr;
ax.XTickLabelRotation = 45;
title('Consistency over data and bootstraps')
grid on

figure;
boxplot(overallBs(methodIdx, :)')
ax = gca;
ax.XTick = 1:size(consistency(methodIdx, :), 1);
ax.XTickLabel = methodNamesAbbr(methodIdx);
ax.XTickLabelRotation = 45;
title('Overall')
grid on

%% see what classification results we get using different FS feats on jv ground truth
% note that the new FCR algorithm does not select the same features as the
% old one.  we want to see if these new feats are actual improvements on
% the old ones...

%% load and init my spekboom data
close all hidden; clear all;
load('D:\Data\Development\Projects\PhD GeoInformatics\Data\Feature Selection\DataAllWin5NoBorder2.mat')
dataAll = changelablist(dataAll, 'Default');
dataAll = setprior(dataAll, 0);
fl = cellstr(getfeatlab(dataAll));
idx = strmatch('Lbp', fl);
dataAll(:, idx) = [];
fl = cellstr(getfeatlab(dataAll));

feats = [9 15 20 23 7 6]; %ranked cluster
fl(feats)

cs = classsizes(dataAll);
if true
    cs(1) = cs(2);
else
    cs = min(cs)*ones(1,3);
end
randreset;
subData = gendat(dataAll, cs);
subData = changelablist(subData, 'Default');
subData = setprior(subData, 0);

%% do the feature selection
preferredFeatures = [9 15 20 23 7 6]; % 5 8 10 1:4 19 21:22 16:18];

res = FeatureClusterRank(subData, 'criterion', naivebc, 'preferredFeatures', ...
    preferredFeatures, 'clusterThresh', 0.175, 'showFigures', true, 'jmiFormulation', true);

fl(res.Feats(1:6))

resfcr = BootstrapFsEval(subData, FeatSelClusterRankM([], naivebc, 0, [], 'preferredFeatures', ...
    preferredFeatures, 'clusterThresh', 0.2, 'showFigures', true, 'jmiFormulation', false), ...
    'numBootStraps', 1, 'numFeatures', 6);
fl(resfcr.FeatIdx(:, 1))

resfcr2 = BootstrapFsEval(subData, FeatSelClusterRankM([], 'mi', 0, [], 'preferredFeatures', ...
    preferredFeatures, 'clusterThresh', 0.175, 'showFigures', true, 'jmiFormulation', true), ...
    'numBootStraps', 1, 'numFeatures', 6);

resfcr3 = BootstrapFsEval(subData, FeatSelClusterRankM([], naivebc, 0, [], 'preferredFeatures', ...
    preferredFeatures, 'clusterThresh', 0.2, 'showFigures', true, 'jmiFormulation', true), ...
    'numBootStraps', 1, 'numFeatures', 6);

resfcr4 = BootstrapFsEval(subData, FeatSelClusterRankM([], 'mi', 0, [], 'preferredFeatures', ...
    preferredFeatures, 'clusterThresh', 0.7, 'showFigures', true, 'jmiFormulation', true, 'useCorrelation', false), ...
    'numBootStraps', 1, 'numFeatures', 6);

fl(resfcr.FeatIdx(:, 1))

% resj = BootstrapFsEval(subData, FeatSelFeastM([], 'jmi', 0), 'numBootStraps', 10, 'numFeatures', 6);
fl(resfcr2.FeatIdx(:, 1))
origFeats = [9 15 20 23 7 6]; % use in thesis
fl(origFeats)

% NOTE:
% - using nmi / mi criterion does not give same feats as orig formulation
% - using naivebc as criterion but with jmi formulation also does not give
% the orig feats
% - bootstrapping the subData does not give the original features
% - the clfr accuracy is MUCH better with criterion = naivebc and orig formulation than with
% criterion = nmi & JMI.  we should check if this applies to other data...

%% get the old classifier working first
% this loads subData and feats
load 'F:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_OpenCvDTree_clfr_3class.mat'
featsOrig = feats;


% 			float priors[3] = {0.25, 0.5, 0.25};
% 			CvDTreeParams dtreeParams;
% 			dtreeParams.max_depth = 12;
% 			dtreeParams.use_1se_rule = false;
% 			dtreeParams.use_surrogates = false;
% 			dtreeParams.cv_folds = 5;
% 			dtreeParams.truncate_pruned_tree = true;
% 			dtreeParams.priors = priors;
% 			dtreeParams.min_sample_count = 20;


%something up with the mapping so retrain
feats = resfcr2.FeatIdx(:, 1);
p = [1 1 1];
w = subData(:, feats)*opencvdtreec([], 12, {'Priors', [0.25 0.5 0.25], 'MaxDepth', 12, 'Use1seRule', false, ...
        'UseSurrogates', false, 'CVFolds', 5, 'TruncatePrunedTree', true, 'MinSampleCount', min(classsizes(subData))/100});


%% train a classifier (as in ThesisClassification)
% choose a modified selection that is as close to the prev one as possible
% (using clusters)

cdata = subData;
p = [1 1 1];
feats = resfcr2.FeatIdx(:, 1);
feats = origFeats;

w = subData(:, feats)*opencvrtreec([], [], {'Priors', [1 2 1]/4, ...
      'MaxNumOfTreesInTheForest', 5, 'NActiveVars', 4, 'CalcVarImportance', false, 'MaxDepth', 10, 'ForestAccuracy', 0.025})

p = [1 1 1];
w = subData(:, feats)*(scalem([], 'variance')*opencvsvc([], [],...
    {'SVMType', 'C_SVC', 'KernelType', 'RBF', 'Gamma', 25, 'C', .1, 'ClassWeights', double(p)./sum(p)}))

w = subData(:, feats)*opencvdtreec([], 12, {'Priors', p./sum(p), 'MaxDepth', 12, 'Use1seRule', false, ...
        'UseSurrogates', false, 'CVFolds', 5, 'TruncatePrunedTree', true, 'MinSampleCount', min(classsizes(subData))/100});

% [err, cerr, nlabOut] = prcrossval(cdata(:, feats), opencvdtreec([], 12, {'Priors', p./sum(p), 'MaxDepth', 12, 'Use1seRule', false, ...
%         'UseSurrogates', false, 'CVFolds', 5, 'TruncatePrunedTree', true, 'MinSampleCount', min(classsizes(subData))/100}), 10);
% toc
% [cc ccr] = ClfrPerfMeas(cdata, nlabOut);


c = confmat(subData(:, feats)*w);
cn = c./repmat(sum(c,2),1,3)


%% setup file name variables for validating on jv gt
close all hidden; 

imageFileNames = {...
    'F:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0415_rgbn_XCALIB.tif';...
    'F:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0419_rgbn_XCALIB.tif';...
    'F:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321D_319_04_0121_rgbn_XCALIB.tif';...
    'F:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3322C_322_02_0056_rgbn_XCALIB.tif';...
    };

gtFileNames = {...
    'C:\Data\Development\Projects\MSc GeoInformatics\Docs\Trimble\GroundTruthCombinedCorrected_TM21_WGS84.shp';...
    'C:\Data\Development\Projects\MSc GeoInformatics\Docs\Trimble\GroundTruthCombinedCorrected_TM21_WGS84.shp';...
    'C:\Data\Development\Projects\MSc GeoInformatics\Docs\Trimble\GroundTruthCombinedCorrected_TM21_WGS84.shp';...
    'C:\Data\Development\Projects\MSc GeoInformatics\Docs\Trimble\GroundTruthCombinedCorrected_TM23_WGS84.shp';...
    };

% clfrFileNames = {... %empty to use workspace vars i.e. feats and w
%     'F:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_OpenCvDTree_clfr_3class.mat';...
%     'F:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_OpenCvDTree_clfr_3class.mat';...
%     'F:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_OpenCvDTree_clfr_3class.mat';...
%     'F:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_OpenCvDTree_clfr_3class.mat';...
%     };
clfrFileNames = {'','','',''};
global w feats

%% Validate Clfr Against Jv GroundTruth

close all hidden;
resC = {};
res = [];
for g = 1:length(gtFileNames)
    resC{g} = ValidateClfrAgainstJvGroundTruth(imageFileNames{g}, gtFileNames{g}, clfrFileNames{g});
end

res = [resC{:}];

for i = 1:length(res)
      fprintf('%s %f %f\n', res(i).name, res(i).gtCover, res(i).clfCover);
end

fprintf('Canopy cover error Mean(abs) Std(abs): %f (%f)\n', mean(abs([res.gtCover]-[res.clfCover])), std(abs([res.gtCover]-[res.clfCover])));
fprintf('Canopy cover error Median Mad: %f (%f)\n', median([res.gtCover]-[res.clfCover]), mad([res.gtCover]-[res.clfCover]));

%%

