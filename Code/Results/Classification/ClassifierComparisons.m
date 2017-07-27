%%
close all; clear all;
prload('G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\DataAllWin5NoBorder2.mat');
dataAll = changelablist(dataAll, 'Default');
dataAll = setprior(dataAll, 0);
feats = [9 15 20 23 7 6]; %ranked cluster

resAll.label = 'All';
[resAll.err, resAll.cerr, resAll.nlabOut, resAll.tclassf, resAll.tress] = prcrossval(dataAll, opencvdtreec([], 12), 5, 1);

resAll.c = confmat(getnlab(dataAll), resAll.nlabOut);
resAll.cn = resAll.c./repmat(sum(resAll.c, 2), 1, size(resAll.c, 2));
resAll.cnErr = 1 - mean(diag(resAll.cn));

[resAll.crn] = ReduceConfMat(resAll.cn, {[1 3], [2]}, true);
resAll.crnErr = 1 - mean(diag(resAll.crn));


%% Compare separation by XX
res = [];
% compareBy = 'Slope';
compareBy = 'Habitat';

dataAll = changelablist(dataAll, 'Default');
classes = cellstr(getlablist(dataAll));
dataAll = changelablist(dataAll, compareBy);
sepFields = cellstr(getlablist(dataAll));
% sepFields(3) = []; %hack out misspelled valleyt

figure;
scatterd2(dataAll(1:10:end, (1:2)), 'legend')
title('All Classes')

dataAll = changelablist(dataAll, 'Default');
figure;
scatterd2(dataAll(1:10:end, (1:2)), 'legend')
title(['All ' compareBy 's'])
% 
% figure;
% for i = 1:length(classes)
%     dataAll = changelablist(dataAll, compareBy);
%     classIdx = (getnlab(dataAll) == getclassi2(dataAll, classes{i}));
%     subData = dataAll(classIdx,:);
%     subData = changelablist(subData, 'Default');
%     h(i) = subplot(2,2,i);
%     scatterd2(subData(1:5:end, feats(1:2)), 'legend')
%     title(classes{i})
% end
% linkaxes(h, 'xy');

err = [];
figure;
for i = 1:length(sepFields)
    dataAll = changelablist(dataAll, compareBy);
    classIdx = (getnlab(dataAll) == getclassi2(dataAll, sepFields{i}));
    subData = dataAll(classIdx,:);
    subData = changelablist(subData, 'Default');
    h(i) = subplot(1,length(sepFields),i);
    scatterd2(subData(1:5:end, feats(1:2)), 'legend');
    title(sepFields{i})
    
    subData = setprior(subData, 0);
    [tr ts] = gendat(subData(:, [2 5]), 0.5);
    if getsize(tr, 3) ==3
%         [err, cerr, nlabOut, tclassf, tress] = prcrossval(subData, scalem([], 'variance')*opencvsvc([], [], ...
%             {'SVMType', 'C_SVC', 'KernelType', 'RBF', 'Gamma', 5, 'C', 10, 'ClassWeights', double([1; 1; 1])}), 5, 1);
        res(i).label = sepFields{i};
        [res(i).err, res(i).cerr, res(i).nlabOut, res(i).tclassf, res(i).tress] = prcrossval(subData, opencvdtreec([], 12), 5, 1);

        res(i).c = confmat(getnlab(subData), res(i).nlabOut);
        res(i).cn = res(i).c./repmat(sum(res(i).c, 2), 1, size(res(i).c, 2));
        res(i).cnErr = 1 - mean(diag(res(i).cn));

        [res(i).crn] = ReduceConfMat(res(i).cn, {[1 3], [2]}, true);
        res(i).crnErr = 1 - mean(diag(res(i).crn));
        
%         w = tr*qdc;
%         c = confmat(ts*w);
%         cn = c./repmat(sum(c, 2), 1, size(c, 2));
%         err = [err 1-mean(diag(cn))];
    end
end



%% Compare separation by custom match list on XX
res = [];
% compareList = {'Valley', 'Arid'};
% compareBy = 'Habitat';
compareList = {'GroenFontein', 'MatjiesVlei', 'RooiBerg', 'GrootKop'};
compareBy = 'Area';

dataAll = changelablist(dataAll, 'Default');
classes = cellstr(getlablist(dataAll));
dataAll = changelablist(dataAll, compareBy);
% sepFields = cellstr(getlablist(dataAll));
sepLabels = cellstr(getlabels(dataAll));
% sepFields(3) = []; %hack out misspelled valleyt

figure;
scatterd2(dataAll(1:10:end, (1:2)), 'legend')
title('All Classes')

dataAll = changelablist(dataAll, 'Default');
figure;
scatterd2(dataAll(1:10:end, (1:2)), 'legend')
title(['All ' compareBy 's'])

figure;
for i = 1:length(compareList)
    dataAll = changelablist(dataAll, compareBy);
%     classIdx = (getnlab(dataAll) == getclassi2(dataAll, classes{i}));
    classIdx = strmatch(lower(compareList{i}), lower(sepLabels));
    subData = dataAll(classIdx,:);
    subData = changelablist(subData, 'Default');
    h(i) = subplot(1,length(compareList),i);
    scatterd2(subData(1:5:end, feats(1:2)), 'legend')
    title(compareList{i})
end
linkaxes(h, 'xy');

err = [];
figure;
for i = 1:length(compareList)
    dataAll = changelablist(dataAll, compareBy);
%     classIdx = (getnlab(dataAll) == getclassi2(dataAll, sepFields{i}));
    classIdx = strmatch(lower(compareList{i}), lower(sepLabels));
    subData = dataAll(classIdx,:);
    subData = changelablist(subData, 'Default');
    h(i) = subplot(1,length(compareList),i);
    scatterd2(subData(1:5:end, feats(1:2)), 'legend');
    title(compareList{i})
    
    subData = setprior(subData, 0);
    [tr ts] = gendat(subData(:, [2 5]), 0.5);
    if getsize(tr, 3) ==3
%         [err, cerr, nlabOut, tclassf, tress] = prcrossval(subData, scalem([], 'variance')*opencvsvc([], [], ...
%             {'SVMType', 'C_SVC', 'KernelType', 'RBF', 'Gamma', 5, 'C', 10, 'ClassWeights', double([1; 1; 1])}), 5, 1);
        res(i).label = sepFields{i};
        [res(i).err, res(i).cerr, res(i).nlabOut, res(i).tclassf, res(i).tress] = prcrossval(subData, opencvdtreec([], 12), 5, 1);

        res(i).c = confmat(getnlab(subData), res(i).nlabOut);
        res(i).cn = res(i).c./repmat(sum(res(i).c, 2), 1, size(res(i).c, 2));
        res(i).cnErr = 1 - mean(diag(res(i).cn));

        [res(i).crn] = ReduceConfMat(res(i).cn, {[1 3], [2]}, true);
        res(i).crnErr = 1 - mean(diag(res(i).crn));
        
%         w = tr*qdc;
%         c = confmat(ts*w);
%         cn = c./repmat(sum(c, 2), 1, size(c, 2));
%         err = [err 1-mean(diag(cn))];
    end
end

fprintf('Results for all data:\n')
resAll.cn
resAll.crnErr
for i = 1:length(compareList)
    fprintf('Results for %s data:\n', compareList{i});
    res(i).cn
    res(i).crnErr
end
fprintf('Mean over %:\n', compareBy);
mean(cat(3, res(1:length(compareList)).cn), 3)
mean([res(1:length(compareList)).crnErr])

% linkaxes(h, 'xy');
% fprintf('Average error for separation by sepFields: %f\n', mean(err));

% 
% 
% [err, cerr, nlabOut, tclassf, tress] = prcrossval(subData, scalem([], 'variance')*opencvsvc([], [], ...
%     {'SVMType', 'C_SVC', 'KernelType', 'RBF', 'Gamma', 5, 'C', 10, 'ClassWeights', double([1; 1; 1])}), 5, 1);
% 
% c = confmat(getnlab(subData), nlabOut);
% cn = c./repmat(sum(c, 2), 1, size(c, 2))
% 1 - mean(diag(cn))
% 
% [CRn] = ReduceConfMat(cn, {[1 3], [2]}, true)
% 1 - mean(diag(CRn))

%sep by valley and arid
%make sep by area fit into above code


% linkaxes(h, 'xy');
% fprintf('Average error for separation by sepFields: %f\n', mean(err));

% 
% 
% [err, cerr, nlabOut, tclassf, tress] = prcrossval(subData, scalem([], 'variance')*opencvsvc([], [], ...
%     {'SVMType', 'C_SVC', 'KernelType', 'RBF', 'Gamma', 5, 'C', 10, 'ClassWeights', double([1; 1; 1])}), 5, 1);
% 
% c = confmat(getnlab(subData), nlabOut);
% cn = c./repmat(sum(c, 2), 1, size(c, 2))
% 1 - mean(diag(cn))
% 
% [CRn] = ReduceConfMat(cn, {[1 3], [2]}, true)
% 1 - mean(diag(CRn))

%sep by valley and arid
%make sep by area fit into above code
