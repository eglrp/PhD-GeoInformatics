function [C Cr] = ClfrPerfMeas(data, nLabOut)
nLabIn = getnlab(data);
ll = getlablist(data);
ll = cellstr(ll);
    
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