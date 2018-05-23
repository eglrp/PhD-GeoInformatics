function res = BootstrapFsEval(data, fs, varargin)
%BootstrapFsEval evaluates feature selection accuracy and stability
%   Evaluates feature selection accuracy and stability for std prtools
%   approaches featsel* 

numBootStraps = 10;
numFeatures = 6;  % limit the number of features for Tanimoto dist and Clfr accuracy, =0 will use the num of features as resturned by the featsel
ModifyDefaultArgs(varargin)
rng('default') % so we get the same results every time and the same crossval folds.  
            % note this doesn't apply to opencv clfr training which
            % apparently uses its own randomiser (I guess unsuprisingly).
            % There doesn't seem to be any mexopencv access to the c++
            % random number generator
randreset;
g = RandStream('twister'); % this should make it work inside parfor
RandStream.setGlobalStream(g);
s = RandStream('twister', 'Seed', 1); % NB make a local randstream that should make same bootstraps each time

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% scale data for any MI type measure which needs to find histograms
data = 10*(data*scalem(data, 'domain'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w = {};
res.FeatIdx = {};
res.FeatScore = zeros(size(data, 2), numBootStraps); 
res.FeatRank = zeros(size(data, 2), numBootStraps); 
for i = 1:numBootStraps
    fprintf('Boot strap %i of %i\n', i, numBootStraps);
    if true
        [subData1, subData2] = gendat(data, 0.5);  % TO DO: seed the sampling so that the same bootstraps can be repeated for different clfrs
        %[subData1, subData2] = gendat2(data, 0.5, i);  % give same random seeds on each call so that we get identical bootstraps for different calls and threahds
        %randi(10,1,10)
    else  % classify on feat sel data
          % TO DO: seed the sampling so that the same bootstraps can be repeated for different clfrs
        if i == 1
            subData1 = data;
        else
            [subData1] = gendat(data);
        end
        subData2 = subData1;
    end
    subData1 = setprior(subData1, 0);
    subData2 = setprior(subData2, 0);
    tStart = tic;
    [w, r] = subData1*fs;    % TO DO: fs should be spec'd to cross validate
    res.FsDuration(i) = toc(tStart);
    res.FeatIdx{i} = +w;         % TO DO: come up with a consistent way of getting all features ranks etc for fs and be and then limiting the features to the best N

    %note: the below will not work with BE or branch and bound I dont think - these
    %should be importance scores per feature.  BE gives the importance of
    %the model excluding that feature and I dont think branch and bound
    %gives importance.  but we could look more deeply into these things
    %note also: the feature importance in FS is not the importance for that
    %feature but for the model including that feature.  This is not the
    %same thing and importances should be per feature importances for 
    %the correlation stability measure to work.  The models including that
    %feature will be different for different runs so it is not a real
    %indication of that features importance at all.  
    %Is there some way we can extract per feature performance from FS and
    %BE model importances?  We could say that the amount that the model
    %accuracy inc or dec is the feature importance but this is in the 
    %context of the model which can change between runs
    %

    %covert model importances to feature importances
    switch struct(w).name
        case 'Forward FeatSel',
            fScore = diff([0; r(:, 2)]);
            res.FeatScore(abs(r(:, 3)), i) = fScore;
            res.FeatRank(abs(r(:, 3)), i) = r(:, 1);
        case 'Backward FeatSel',
            fScore = diff(r(:, 2));
            res.FeatScore(abs(r(2:end, 3)), i) = fScore;
            res.FeatRank(abs(r(2:end, 3)), i) = r(1:end-1, 1);  % rank features as the order in which they are eliminated
            res.FeatRank(setdiff(1:size(r, 1), abs(r(:,3))), i) = 1;    % rank the final remaining feature as 1
            if numFeatures > 0
                res.FeatIdx{i} = 1:size(r, 1); % a hack for BE - FeatIdx = +w always return the optimal set but we want the full set to select numFeatures from
            end
        case 'B&B FeatSel', %do nothing - we dont have scores 
        case 'Feature Clustering and Ranking'
            resFcr = struct(w).user.user;
            [clustAcc_ clustIdx] = sort(resFcr.ClustAcc);
            res.FeatScore(:, i) = resFcr.FeatClustScore(:);
            res.FeatRank(:, i) = resFcr.FeatClustRank(:);
            res.FeatClustIdx(:, i) = resFcr.FeatClustIdx(:); % hack in extra field for FCR
            res.ClustFeatNLab{i} = resFcr.ClustFeatNLab; % hack in extra field for FCR
            %res.FeatIdx{i} = +w;         % index features by their clusters rather for purposes of stability calcs
            %res.FeatIdx{i} = resFcr.FeatClustIdx(+w);

        otherwise, %assume importances are already per individual feature
            res.FeatScore(abs(r(:, 3)), i) = r(:, 2);
            res.FeatRank(abs(r(:, 3)), i) = r(:, 1);
    end

    % evaluate classifier accuracy
    % Note that opencvrtreec and opencvdtreec use their own random seeding
    % apparently so initialising it outside doesn't help
    % NB Note: these clfr params are optimised for spekboom data (changed to more generic ones now)
            % opencvrtreec([], [], {'Priors', ones(1, getsize(subData2, 3))/getsize(subData2, 3), 'MaxNumOfTreesInTheForest', 5, 'NActiveVars', 4, 'CalcVarImportance', true, 'MaxDepth', 10, 'ForestAccuracy', 0.025}), ...
            % opencvdtreec([], 12, {'Priors', ones(1, getsize(subData2, 3))/getsize(subData2, 3), 'MaxDepth', 12, 'Use1seRule', false, 'UseSurrogates', false, 'CVFolds', 5, 'TruncatePrunedTree', true, 'MinSampleCount', min(classsizes(subData2))/100}), ...
%     clfrs = {opencvsvc([], [], {'SVMType', 'C_SVC', 'KernelType', 'RBF', 'Gamma', 2, 'C', 1, 'ClassWeights', ones(1, getsize(subData2, 3))}), ...
%             opencvrtreec([], [], {'Priors', ones(1, getsize(subData2, 3))/getsize(subData2, 3), 'MaxNumOfTreesInTheForest', 50, 'NActiveVars', 0, 'CalcVarImportance', true, 'MaxDepth', 10, 'ForestAccuracy', 0.025}), ...
%             opencvdtreec([], 12, {'Priors', ones(1, getsize(subData2, 3))/getsize(subData2, 3), 'MaxDepth', 12, 'Use1seRule', false, 'UseSurrogates', false, 'CVFolds', 5, 'TruncatePrunedTree', true, 'MinSampleCount', 10}), ...
%             naivebc, ...
%             opencvknnc([], 3)
%             };
    clfrs = {
%             opencvsvc([], [], {'SVMType', 'C_SVC', 'KernelType', 'RBF', 'Gamma', 2, 'C', 1, 'ClassWeights', ones(1, getsize(subData2, 3))}), ...
%             opencvrtreec([], [], {'Priors', ones(1, getsize(subData2, 3))/getsize(subData2, 3), 'MaxNumOfTreesInTheForest', 50, 'NActiveVars', 0, 'CalcVarImportance', true, 'MaxDepth', 10, 'ForestAccuracy', 0.025}), ...
%             opencvdtreec([], 12, {'Priors', ones(1, getsize(subData2, 3))/getsize(subData2, 3), 'MaxDepth', 12, 'Use1seRule', false, 'UseSurrogates', false, 'CVFolds', 5, 'TruncatePrunedTree', true, 'MinSampleCount', 10}), ...
%             naivebc, ...
            opencvknnc([], 3)
            };

    if numFeatures > 0
        fr = res.FeatRank(res.FeatIdx{i}, i);
        [fr_ fidx] = sort(fr);
        %fprintf('res.FeatIdx{i}: %d, numFeatures: %d\n', length(res.FeatIdx{i}), numFeatures);
        if numFeatures > length(fidx)  % we didn't get as many features as we wanted - something is wrong
            fprintf('ERROR: not enough features - length(fidx): %d, numFeatures: %d\n', length(fidx), numFeatures);
        end
        res.FeatIdx{i} = res.FeatIdx{i}(fidx(1:min(numFeatures, length(fidx))));
    end
    for ci = 1:length(clfrs)
        fprintf('    Clfr eval %i of %i\n', ci, length(clfrs));
        % limit the number of features for fs and be if spec'd
        % features are not necessarily in order of rank, so make sure we
        % choose the best ones 
        % NB: all clfrs get data scaled to unit variance
        [err, cerr, nlabOut] = prcrossval(subData2(:, res.FeatIdx{i}), scalem([], 'variance')*clfrs{ci}, 5, 1);
        c = confmat(getnlab(subData2), nlabOut);
        res.ClfCn{i, ci} = c./repmat(sum(c, 2), 1, size(c, 2));
        res.ClfAcc(i, ci) = mean(diag(res.ClfCn{i, ci}));
        res.ClfName{ci} = struct(clfrs{ci}).name;
    end
end

%use cluster indices rather than feature indices for FCR stability eval
% the problem with this is that the clusters and their numbers may change
% from one bootstrap to the next, so the indexing is not consistent, if we
% wanted to go this way, we would have to do the clustering once and then
% bootstrap the ranking.  which is a little suspect...  so then lets
% renumber clusters from all bootstraps i.e. if clusters are identical,
% they get the same numbers, if not they get new numbers

if strcmpi(struct(w).name, 'Feature Clustering and Ranking')
    res = RenumClustAcrossBootstraps(res);

    if false
        featIdx = res.FeatIdx;
        for i = 1:length(res.FeatIdx)
            res.FeatIdx{i} = res.FeatClustIdx(res.FeatIdx{i}, i);
        end
        res.FeatIdxOrig = featIdx;
    end    
end

res.ClfMeanAcc = mean(res.ClfAcc, 1);
res.ClfStdAcc = std(res.ClfAcc, 1);
res.FsMeanDuration = mean(res.FsDuration);
res = FsStabilityEval(res);
res.FsMapping = fs;
