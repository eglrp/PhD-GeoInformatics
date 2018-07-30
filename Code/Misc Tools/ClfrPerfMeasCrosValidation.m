function [C Cr] = ClfrPerfMeasCrosValidation(data, w, doCv)
if nargin < 3
    doCv = true;
end
nFolds = 10;

randreset;
r = prcrossval(data, [], nFolds, 0);
fprintf('\n')
for f = 1:nFolds
    fprintf('.')
    tr = data(r~=f, :);
    ts = data(r==f, :);
    wtr = tr*w;
    out = ts*wtr;
    nLabIn = getnlab(ts);
    nLabOut = out*nlabeld;
    c(:,:,f) = confmat(nLabIn, nLabOut);
    cn(:,:,f) = c(:,:,f)./repmat(sum(c(:,:,f), 2), 1, size(c(:,:,f), 2));
    [cRn_] = ReduceConfMat(cn, {[1 3], [2]}, false);
    % 1 - mean(diag(CRn))

    [cR(:,:,f)] = ReduceConfMat(c(:,:,f), {[1 3], [2]}, false);
    cRn(:,:,f) = cR(:,:,f)./repmat(sum(cR(:,:,f), 2), 1, size(cR(:,:,f), 2));
end
fprintf('\n')

ll = getlablist(data);
ll = cellstr(ll);

if false
    nLabIn = getnlab(data);

    c = confmat(nLabIn, nLabOut);
    cn = c./repmat(sum(c, 2), 1, size(c, 2))
    % 1 - mean(diag(cn))
    % 1 - sum(diag(c))/sum(c(:))

    [cRn_] = ReduceConfMat(cn, {[1 3], [2]}, false);
    % 1 - mean(diag(CRn))

    [cR] = ReduceConfMat(c, {[1 3], [2]}, false);
    cRn = cR./repmat(sum(cR, 2), 1, size(cR, 2))
    % assert(all(cRn == cRn_))

    C = DispConfMat(c, ll);
    Cr = DispConfMat(cR, ll([1 2]));
else
%     round(cnAvg,3)==round(cnAvg_,3)
    C = DispConfMatCrossValidation(c, ll);
    Cr = DispConfMatCrossValidation(cR, ll([1 2]));

    if false
        cAvg = sum(c, 3);
        cnAvg_ = mean(cn, 3);
        cnAvg = cAvg./repmat(sum(cAvg, 2), 1, size(cAvg, 2));

        cRAvg = sum(cR, 3);
        cRnAvg_ = mean(cRn, 3);
        cRnAvg = cRAvg./repmat(sum(cRAvg, 2), 1, size(cRAvg, 2));

        C = DispConfMatCrossValidation(cAvg, ll);
        Cr = DispConfMatCrossValidation(cRAvg, ll([1 2]));
    end

%     C2 = DispConfMat(cAvg, ll);
%     Cr2 = DispConfMat(cRAvg, ll([1 2]));
end

end

function cc = DispConfMatCrossValidation(cV, ll)
    cV_ = cV;
    prodAccV = zeros(size(cV, 3), size(cV, 2)+1);
    consAccV = zeros(size(cV, 3), size(cV, 1)+1);
    overallErrV = zeros(size(cV, 3), 1);
    overallWErrV = zeros(size(cV, 3), 1);
    kappaV = zeros(size(cV, 3), 1);
    for f = 1:size(cV, 3)
        c = cV(:,:,f);
        c_ = c;
        c(:, end+1) = sum(c, 2);
        c(end+1, :) = sum(c, 1);
        consAccV(f,:) = diag(c(:,:))./c(end, :)';
        prodAccV(f,:) = diag(c(:,:))./c(:, end);
%         c(:, end+1) = 100*prodAcc;
%         c(end+1, :) = 100*[consAcc' nan];

        overallErrV(f) = sum(diag(c_))/sum(sum(c_)); %prior dep
        overallWErrV(f) = mean(prodAccV(f, 1:end-1)); %prior ind / = prior

        c_n = (c_')./sum(sum(c_));
        po = sum(diag(c_n)); %==overallErr
    %     assert(po==overallErr);

        c_n(:, end+1) = sum(c_n, 2);
        c_n(end+1, :) = sum(c_n, 1);
        pc = sum(c_n(end, 1:end-1).*(c_n(1:end-1, end)'));
        pc = c_n(end, 1:end-1)*c_n(1:end-1, end);

        kappaV(f) = (po-pc)/(1-pc);
    end
    cAvg = sum(cV, 3);
    cAvg_ = cAvg;
    
    cnAvg = cAvg./repmat(sum(cAvg, 2), 1, size(cAvg, 2));
%     cRAvg = ReduceConfMat(cAvg, {[1 3], [2]}, false);
%     cRnAvg = cRAvg./repmat(sum(cRAvg, 2), 1, size(cRAvg, 2));
    
    cAvg(:, end+1) = sum(cAvg, 2);
    cAvg(end+1, :) = sum(cAvg, 1);
    consAcc = diag(cAvg)./cAvg(end, :)';
    prodAcc = diag(cAvg)./cAvg(:, end);
    cAvg(:, end+1) = 100*prodAcc;
    cAvg(end+1, :) = 100*[consAcc' nan];
    
    overallErr = sum(diag(cAvg_))/sum(cAvg_(:)); %prior dep
    overallWErr = mean(prodAcc(1:end-1)); %prior ind / = prior
    
    cAvg_n = (cAvg_')./sum(cAvg_(:));
    po = sum(diag(cAvg_n)); %==overallErr
%     assert(po==overallErr);

    cAvg_n(:, end+1) = sum(cAvg_n, 2);
    cAvg_n(end+1, :) = sum(cAvg_n, 1);
    pc = sum(cAvg_n(end, 1:end-1).*(cAvg_n(1:end-1, end)'));
    pc = cAvg_n(end, 1:end-1)*cAvg_n(1:end-1, end);

    kappa = (po-pc)/(1-pc);

    cc = {'', ll{:}, 'Total', 'Prod Acc'};
    cc(2:length(ll)+1, 1) = {(ll{:})}; 
    cc(end+1:end+2, 1) = {'Total'; 'Cons Acc'};

    cc(2:size(cAvg,1)+1, 2:size(cAvg,2)+1) = mat2cell(cAvg, ones(1, size(cAvg,1)), ones(1, size(cAvg,2)));
    
    cc(1, end+1) = {'PA Std %'};
    cc(end+1, 1) = {'CA Std %'}; 
    cc(2:size(prodAccV, 2)+1, end) =  mat2cell(100*std(prodAccV, [], 1)', ones(size(prodAccV, 2), 1), 1);
    cc(end, 2:size(consAccV, 2)+1) =  mat2cell(100*std(consAccV, [], 1), 1, ones(1, size(consAccV, 2)));
    cc(2:size(cAvg,1)+1, 2:size(cAvg,2)+1) = mat2cell(cAvg, ones(1, size(cAvg,1)), ones(1, size(cAvg,2)));
    
    cc(end+1,1:3) = {'Kappa (Std)', kappa, std(kappaV)};
    cc(end+1,1:3) = {'Overall acc (Std%)', 100*(1-overallErr), 100*std(overallErrV)};
    cc(end+1,1:3) = {'Overall =prior acc (Std%)', 100*(1-overallWErr), 100*std(overallWErrV)};
    %estimate the canopy cover error from the conf mat
    cr = cAvg_./repmat(sum(cAvg_,2), 1, size(cAvg_,2));
    pActual = sum(cAvg_, 1)./sum(cAvg_(:)); %actual priors
    
    lims = ones(size(cAvg_,1), 2);
    lims(:,1) = -inf;
    lims(:,2) = inf;
    priorGrid = PriorGrid(lims, 0.1);
    sbIdx = find(strcmpi(ll, 'Spekboom'));
    for i = 1:size(priorGrid,1)
        pHat(i,:) = sum(repmat(priorGrid(i,:)', 1, size(cAvg_, 2)) .* cr, 1);
        ccErr(i) = abs(pHat(i, sbIdx)-priorGrid(i, sbIdx));
    end

    cc(end+1,1:4) = {'CC abs err (all prior combos)', 100*mean(ccErr), 100*std(ccErr), 100*max(ccErr)};

    disp(cc);
%     fprintf('Kappa: %f\n', kappa);
%     fprintf('Overall acc: %f\n', overallErr);
%     fprintf('Overall =prior acc: %f\n\n', overallWErr);
    
end

function cc = DispConfMat(c, ll)
    c_ = c;

    
    c(:, end+1) = sum(c, 2);
    c(end+1, :) = sum(c, 1);
    consAcc = diag(c)./c(end, :)';
    prodAcc = diag(c)./c(:, end);
    c(:, end+1) = 100*prodAcc;
    c(end+1, :) = 100*[consAcc' nan];
    
    overallErr = sum(diag(c_))/sum(c_(:)); %prior dep
    overallWErr = mean(prodAcc(1:end-1)); %prior ind / = prior
    
    c_n = (c_')./sum(c_(:));
    po = sum(diag(c_n)); %==overallErr
%     assert(po==overallErr);

    c_n(:, end+1) = sum(c_n, 2);
    c_n(end+1, :) = sum(c_n, 1);
    pc =  sum(c_n(end, 1:end-1).*(c_n(1:end-1, end)'));
    pc =  c_n(end, 1:end-1)*c_n(1:end-1, end);

    kappa = (po-pc)/(1-pc);

    cc = {'', ll{:}, 'Total', 'Prod Acc'};
    cc(2:length(ll)+1, 1) = {(ll{:})}; 
    cc(end+1:end+2, 1) = {'Total'; 'Cons Acc'};

    cc(2:size(c,1)+1, 2:size(c,2)+1) = mat2cell(c, ones(1, size(c,1)), ones(1, size(c,2)));
    cc(end+1,1:2) = {'Kappa', kappa};
    cc(end+1,1:2) = {'Overall acc', 100*(1-overallErr)};
    cc(end+1,1:2) = {'Overall =prior acc', 100*(1-overallWErr)};
    %estimate the canopy cover error from the conf mat
    cr = c_./repmat(sum(c_,2), 1, size(c_,2));
    pActual = sum(c_, 1)./sum(c_(:)); %actual priors
    
    lims = ones(size(c_,1), 2);
    lims(:,1) = -inf;
    lims(:,2) = inf;
    priorGrid = PriorGrid(lims, 0.1);
    sbIdx = find(strcmpi(ll, 'Spekboom'));
    for i = 1:size(priorGrid,1)
        pHat(i,:) = sum(repmat(priorGrid(i,:)', 1, size(c_, 2)) .* cr, 1);
        ccErr(i) = abs(pHat(i, sbIdx)-priorGrid(i, sbIdx));
    end

    cc(end+1,1:4) = {'CC abs error (over all prior combos)', 100*mean(ccErr), 100*std(ccErr), 100*max(ccErr)};

    disp(cc);
%     fprintf('Kappa: %f\n', kappa);
%     fprintf('Overall acc: %f\n', overallErr);
%     fprintf('Overall =prior acc: %f\n\n', overallWErr);
    
end