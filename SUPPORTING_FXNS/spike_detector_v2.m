function [spkTms, spkInds] = spike_detector_v2(sig, ts, fs, f_low, f_high, minProm, varargin)
% [spkTms, spkInds] = spike_detector_v2(sig, ts, fs, f_low, f_high, minProm, maxWidth, skewDir)
%
% PURPOSE:
%   This function is modified from the "spike_detector" function used in Yi
%   et al. J. Neuro. 2025. These modifications allow spike indicies to be
%   output, add a maximum width criterion to prevent erroneous spikes
%   from being detected (an issue Meg ran into when first using this
%   function), ...
%
% INPUTS:
%   sig = lfp signal
%   ts = time series for the lfp
%   f_low = freq for low pass filter
%   f_high = frequency for high pass filter
%   minProm = minimum peak prominance for detected events
%   maxWidth = (optional) maximum peak width, in s
%       default = 0.10 s
%   skewDir = direction of IS in the singnal (1 or -1)
%       default = sign(skewness(signal))
%
% OUTPUT:
%   spkTms = times of interictal spikes
%   spkInds = indices of interictal spikes from the lfp array
%
% Modified MMD from code by JDY
% Ewell Lab
% 06/2026

if nargin > 6
    maxWidth = varargin{1};
else
    maxWidth = 0.1; %100 ms
end %check max width
if nargin > 7
    skew_direction = varargin{2};
else
    skew_direction = sign(skewness(sig)); %whether skews neg or pos overall
end %check skewness

maxWidthBins = maxWidth * fs; %convert based on sampling frequency

sig_bp = bandpass(sig, [f_low,f_high], fs);

[~, spkInds] = findpeaks(skew_direction*sig_bp, 'MinPeakProminence', minProm, 'MaxPeakWidth', maxWidthBins);

spkTms = ts(spkInds);

end %function