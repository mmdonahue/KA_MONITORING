function [ts_spike, i_pks] = spike_detector(sig,ts,fs,f_low,f_high,theta)
    sig_bp = bandpass(sig,[f_low,f_high],fs);

    % find NEGATIVE-going peaks
    skew_direction = sign(skewness(sig));               
    [~,i_pks]=findpeaks(skew_direction*sig_bp,"MinPeakProminence",theta);

    ts_spike = ts(i_pks);
end
 
% rms-based
%sig_rms = sqrt(movmean(sig_bp.^2,[15,15]));
%[~,i_pks]=findpeaks(sig_rms,"MinPeakProminence",theta,"MinPeakWidth",5);