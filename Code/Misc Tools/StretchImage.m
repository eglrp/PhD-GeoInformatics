function [strecthIm stretchLim] = StretchImage(im, varargin)

maxVal = (2^15-1); %MODIS
noDataVal = 0;
stretchLim = [];

ModifyDefaultArgs(varargin);

strecthIm = zeros(size(im), 'double');

noDataMask = any(im == noDataVal, 3);
noDataMask = noDataMask | any(im == maxVal, 3);
for i = 1:size(im, 3)
    data = im(:,:,i);
    nonodata = data(~noDataMask);
%     nonodata = nonodata(nonodata~=noDataVal);
    nonodata = double(nonodata)./maxVal;
    data(data==maxVal)=0;
    data(data==noDataVal)=0;
    data = double(data)./maxVal;
    if (isempty(stretchLim))
        stretchLim_{i} = stretchlim(nonodata,[0.01 .99]);
%         stretchLim_{i}(1) = 0;
    else
        stretchLim_{i} = stretchLim{i};
    end
    strecthIm(:,:,i) = imadjust(data, stretchLim_{i}, [0 1], 1);    
end
stretchLim = stretchLim_;
end

