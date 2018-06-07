%% use clusttools to compare the prototype stability of various clustering
% algorithms

load 'D:\Data\Development\Projects\PhD GeoInformatics\Data\Feature Selection\CompareFsMethodsHs4.mat'
close all
clear my*
clear i m tr ts subData uci* hs* data feats cs idx fl res* cres
randreset
clustMethods = {dclustm, dcluste, dclusth([],[],'a'), dclustk([],[],'kmedoids')};  %dclustf([], 18), 
numBootStraps = 10;
nClust = [12, 7, 10, 10, 20, 20];  %obtained from crude examination of clust output w/o spec'ing num clusters
for mi = 1:length(clustMethods)
    disp(clustMethods{mi});
    for di = 1:length(cdata)
%         clustMethods = {dclustm([], nClust(di)), dcluste([], nClust(di)), dclusth([], nClust(di), 'a'), dclustk([], nClust(di), 'kmedoids')};  %dclustf([], 18), 
        disp(cdataNames(di));
        for bi = 1:numBootStraps        %bootstraps
            data = gendat(cdata{di});
            data = data*scalem(data, 'variance');
            c = 1-abs(corr(+data));
            lab = c*clustMethods{mi};
            nFeats = [];
            for ci = 1:size(lab, 2)  %choose the clustering that is closest to nClust
                nFeats(ci) = length(unique(lab(:, ci)));
            end
            [tmp, idx] = sort(abs(nFeats - nClust(di)));            
            cres{mi, di}.Lab{bi} = lab(:, idx(1));
            cres{mi, di}.FeatIdx{bi} = unique(lab(:, idx(1)))';
        end
        cres{mi, di}.N = size(cdata{di}, 2);
        cres{mi, di} = FsStabilityEval(cres{mi, di});
    end
end

%%


%%
methodStability = [];
methodConsistency = [];
for di = 1:length(cdata)
    for mi = 1:length(clustMethods)
        methodStability(mi, di) = cres{mi,di}.TanimotoStability;
        methodConsistency(mi, di) = cres{mi,di}.Consitency;
    end
end

%%
for di = 1:length(cdata)
    for mi = 1:length(clustMethods)
        numFeats = [];
        for bi = 1:numBootStraps
            for ci = 1:size(cres{mi,di}.Lab{bi}, 2)
                numFeats(bi,ci) = length(unique(cres{mi,di}.Lab{bi}(:, ci)));
            end
        end
        fprintf('Method %s, Data %s\n', struct(clustMethods{mi}).name, cdataNames{di});
        disp(numFeats)
    end
end

% mean(methodStability, 2)