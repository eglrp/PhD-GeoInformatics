function [colourBalImage whiteBal] = ColourBalImage(refIm, srcIm, varargin)

noDataVal = 0;
whiteBal = [];
ModifyDefaultArgs(varargin);

colourBalImage = srcIm;
for i = 1:size(refIm, 3)
    refBand = refIm(:,:,i);
%     refBand = refBand(refBand~=noDataVal);

    srcBand = srcIm(:,:,i);
%     srcBand = srcBand(srcBand~=noDataVal);
    if (isempty(whiteBal))
        whiteBal_(i) = double(mean(refBand(:)))./double(mean(srcBand(:)));
    else
        whiteBal_(i) = whiteBal(i);
    end
    colourBalImage(:,:,i) = srcBand.*whiteBal_(i);
end
whiteBal = whiteBal_;
end

