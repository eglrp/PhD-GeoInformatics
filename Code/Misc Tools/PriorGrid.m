function Priors = PriorGrid(PriorRange, StepSize)
%PriorRange is the min and max for each class use +-inf for no limit

if (size(PriorRange,2)~=2)
    error('size(PriorRange,2)~=2');
end

numClasses = size(PriorRange,1);
if (numClasses > 5)
    error('numClasses > 5');
end

switch numClasses
    case 2,
        [p1 p2] = ndgrid(0:StepSize:1, 0:StepSize:1);
        p = [p1(:) p2(:)];
    case 3,
        [p1 p2 p3] = ndgrid(0:StepSize:1, 0:StepSize:1, 0:StepSize:1);
        p = [p1(:) p2(:) p3(:)];
    case 4,
        [p1 p2 p3 p4] = ndgrid(0:StepSize:1, 0:StepSize:1, 0:StepSize:1, 0:StepSize:1);
        p = [p1(:) p2(:) p3(:) p4(:)];
    case 5,
        [p1 p2 p3 p4 p5] = ndgrid(0:StepSize:1, 0:StepSize:1, 0:StepSize:1, 0:StepSize:1, 0:StepSize:1);
        p = [p1(:) p2(:) p3(:) p4(:) p5(:)];
end

idx = sum(p, 2) == 1;
for i = 1:size(PriorRange,1)
    idx = idx & ((p(:,i) >= PriorRange(i,1)) & (p(:,i) <= PriorRange(i,2)));
end
Priors = p(idx,:);
