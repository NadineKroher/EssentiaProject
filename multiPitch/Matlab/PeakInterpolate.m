function [peak, peakAmp] = PeakInterpolate(ObSpecMag, peak)
% Quadratic interpolation
%
% Input
%   - ObSpecMag         : amplitude spectrum
%   - peak              : frequency bins of detected peaks, column vector 
%                         (peak=1 means that a peak appears at the first positive frequency bin)
% Output
%   - peak              : frequency bins of detected peaks, column vector 
%                         (peak=1 means that a peak appears at the first positive frequency bin)
%   - peakAmp           : amplitude in dB of peaks
% 
% Author: Zhiyao Duan
% Created: 2007

if ~isempty(peak)
    peakAmp = ObSpecMag(peak) - 1/8 *((ObSpecMag(peak+1)-ObSpecMag(peak-1)).^2) ./ (ObSpecMag(peak-1)+ObSpecMag(peak+1)-2*ObSpecMag(peak));
    peak = peak - 1/2*(ObSpecMag(peak+1)-ObSpecMag(peak-1)) ./ (ObSpecMag(peak-1)+ObSpecMag(peak+1)-2*ObSpecMag(peak));
else
    peak = [];
    peakAmp = [];
end
