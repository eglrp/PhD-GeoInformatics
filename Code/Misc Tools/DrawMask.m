function mask = DrawMask(image, maskFileName, dispMask)

if (nargin == 1)
    maskFileName = 'c:\ngiMask.mat';
end

mask = false(size(image(:,:,1)));

if (exist(maskFileName, 'file'))
    load(maskFileName, 'mask');
    image(repmat(mask, [1 1 3])) = 1;
end

if (nargin == 3)
    dispMask = bwperim(dispMask);
    dispMask(:,:,2:3) = 0;
    image(dispMask) = 1;
end

figure;
imshow(image)

i = 1;
cont = '';
while (~strcmpi(cont,'y'))
    h_ = imfreehand(gca);
    if (~isempty(h_))
        h(i) = h_;
        wait(h(i));
        mask(h(i).createMask) = true;
        cont = input('Press y to stop or anything else to continue\n','s');
        i = i+1;
        save(maskFileName, 'mask');
    end
end

end

