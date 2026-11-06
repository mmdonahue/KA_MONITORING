function [TP, FP, FN, Precision, Recall, F1, F05] = PrecSummary_v2(spkInds, labInds)
% function [Precision, Recall, F1, F05] = PrecSummary(spkInds, labInds)
%
% PURPOSE:
%   This function is modified from the "PercSummary" function used in Yi
%   et al. J. Neuro. Original notes from JY:
%       classify each detection & ROI according to the ground-truth as follows:
%       for a given ROI, if there is at least one detection inside, call one the
%       TP. All others inside are FP. Similarly, any detections outside an ROI is
%       a FP. Finally, ROIs which lack detection are FN. We do not explicitly
%       consider TNs.
%
% INPUT:
%   spkInds = detected events
%   labInds = labeled events
%       NOTE: You could use spkTms and labTmes, as long as they are
%       measured in the same unit.
%
% OUTPUT:
%   TP = true positives
%   FP = false positives
%   FN = false negatives
%   Precision = TP / (TP + FP)
%   Recall = TP / (TP + FN)
%   F1 = (1+1^2) * ((Precision * Recall) /(1^2 x Preceision) + Recall)
%   F05 = (1+0.5^2) * ((Precision * Recall) /(0.5^2 x Preceision) + Recall)
%
% Modified MMD from code by JDY
% Ewell Lab
% 07/2026


TP = 0;
FP = 0;
FN = 0;

for i = 1:length(labInds)
    % how many detections within this ROI?
    num_dets = length(find((spkInds>=labInds(i,1) & spkInds<labInds(i,2))));

    if num_dets == 0
        FN = FN + 1;
    elseif num_dets == 1
        TP = TP + 1;
    elseif num_dets > 1
        TP = TP + 1;
        FP = FP + num_dets-1;
    end
end

% how many detections outside an ROI?
num_dets = length(find((spkInds<labInds(1,1))));
FP = FP + num_dets;

for i = 1:length(labInds)-1
    % how many detections between end of ROI(i) and start of ROI(i+1)?
    num_dets=length(find((spkInds>=labInds(i,2) & spkInds<labInds(i+1,1))));

    FP = FP + num_dets;
end

num_dets = length(find((spkInds>=labInds(end,2))));
FP = FP + num_dets;

Precision = TP/(TP + FP);
Recall = TP/(TP+FN);
F1 = (2*Precision*Recall)/(Precision+Recall);
F05 = ((1+0.5^2)*Precision*Recall)/((0.5^2)*Precision+Recall);
end