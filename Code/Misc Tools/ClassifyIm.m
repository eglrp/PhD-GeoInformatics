function outIm = ClassifyIm(imstruct, w, feats)
%     global feats;
%     if (nargin==1 || isempty(featsArg))
%         featsArg = feats;
%     end
    featIm = ExtractFeaturesIm2(double(imstruct.data), 'feats', feats);
    out = prdataset(im2feat(featIm(:,:,feats)))*w*classc;
%     outIm = reshape(+out(:,2), [size(imstruct.data, 1) size(imstruct.data, 2)]);    
    outIm = im2uint8(reshape(+out, [size(imstruct.data, 1) size(imstruct.data, 2) size(w, 2)]));
    if (size(outIm, 3) == 2)
        outIm(:,:,3) = 0;
    end
    fprintf('.');
end