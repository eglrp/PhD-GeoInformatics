function [outputFeatures, featureJMIMatrix] = JMIC(data, K, varargin)
% JMI_C JMI feature selection with optional computational cost
% A Matlab implementation of the FEAST JMI function with additional logic
% to allow preference of specific features.
% 
% INPUT
%   data    Training dataset
%   K    Number of features to select (default: sort all features)
preferredFeatures = [];
similarityThresh = 0.05;  
ModifyDefaultArgs(varargin)

% fl = cellstr(getfeatlab(data));
noOfFeatures = size(data, 2);
selectedFeatures = false(noOfFeatures, 1);
outputFeatures = -1*ones(K, 1);
classColumn = getnlab(data);
featureJMIMatrix = -1*ones(K, noOfFeatures);
featureMIMatrix = -1*ones(K, noOfFeatures);

classMI = -1*ones(noOfFeatures, 1);
for i = 1:noOfFeatures
    classMI(i) = mi(+data(:, i), classColumn);
end

[maxClassMi, maxClassMiIdx] = max(classMI);

if ~isempty(preferredFeatures)
    possibleFeatures = find(classMI >= (maxClassMi - similarityThresh));
    if length(possibleFeatures) > 1
        % fl(possibleFeatures)
        weightedMi = classMI .* preferredFeatures;
        [m, idx] = min(weightedMi(possibleFeatures));
        fprintf('Choosing %i with MI %f instead of %i with MI %f\n', ...
            possibleFeatures(idx), classMI(possibleFeatures(idx)), ...
            maxClassMiIdx, classMI(maxClassMiIdx));
        maxClassMiIdx = possibleFeatures(idx);
    end
end
outputFeatures(1) = maxClassMiIdx;
selectedFeatures(maxClassMiIdx) = true;

for i = 2:K
%     score = 0.0;
%     currentHighestFeature = 0;
    currentScore = zeros(noOfFeatures, 1);

    for j = 1:noOfFeatures
        if ~selectedFeatures(j)
            currentScore(j) = 0.0;
            for x = 1:i-1
                if (featureJMIMatrix(x, j) == -1)
                    mergedVector = joint(+data(:, [outputFeatures(x) j]), []);
                    featureJMIMatrix(x, j) = mi(mergedVector, classColumn);
%                     featureMIMatrix(x, j) = mi(+data(:, outputFeatures(x)), ...
%                         +data(:, j));
                end
                currentScore(j) = currentScore(j) + featureJMIMatrix(x, j);
            end
%             if (currentScore(j) > score)
%                 score = currentScore;
%                 currentHighestFeature = j;
%             end
        end
    end
    [score, currentHighestFeature] = max(currentScore);
    if ~isempty(preferredFeatures)
        possibleFeatures = find(currentScore >= (score - similarityThresh));
        if length(possibleFeatures) > 1
            weightedScore = score .* preferredFeatures;
            [m, idx] = min(weightedScore(possibleFeatures));
            % fl(possibleFeatures)
            fprintf('Choosing %i with JMI %f instead of %i with JMI %f\n', ...
                possibleFeatures(idx), currentScore(possibleFeatures(idx)), ...
                currentHighestFeature, currentScore(currentHighestFeature));

            currentHighestFeature = possibleFeatures(idx);
        end
    end
    selectedFeatures(currentHighestFeature) = 1;
    outputFeatures(i) = currentHighestFeature;
end

