function [lutIm] = ApplyIntergraphLut(im, lut)

lutIm = im;
for i = 2:size(lut, 2)
    lut16 = zeros(2^16, 1, 'uint16');
    lut16(1:length(lut(:, i))) = lut(:, i);
    lutIm(:,:,i-1) = intlut(im(:,:,i-1), lut16);
end

end

