function [peak, peakNum, ObSpecMagSmoothed] = PeakExtract(ObSpecMag, para, IfPlot)
% Peak detection from amplitude spectrum
% Input
%   - ObSpecMag         : magnitude spectrum
%   - para 
%       - PeakTh        : the absolute threshold (peaks above the maximum amplitude minus this value are considered)
%       - peakTh_rel    : the relative threshold (peaks above the smoothed envelope plus this value are considered)
%       - peakTh_freq   : the frequency limit (peaks below this frequency in bins are considered)
%       - movL          : the moving average threshold
%       - localRange    : the radius of local range to get the local maximum
%   - IfPlot            : plot figures or not
% Output
%   - peak              : frequency bins of detected peaks, column vector (peak=1 means that a peak appears at the first positive frequency bin)
%   - peakNum           : number of detected peaks
%   - ObSpecMagSmoothed : smoothed magnitude spectrum
%
% Author: Zhiyao Duan
% Last modified: 9/4/2009

peakTh = para.peakTh;                                                       % global threshold
peakTh_rel = para.peakTh_rel;                                               % local threshold
peakTh_freq = para.peakTh_freq;                                             % frequency threshold (frequency bins)
movL = para.movL;                                                           % moving average width used to calcualte the smoothed log amplitude spectrum
localRange = para.localRange;                                               % local range width to decide local maximum

peak = [];
peakNum = 0;

if sum(ObSpecMag == 0) == length(ObSpecMag)                                 % if all the elements are zero, then return
    return;
end

SpecLength = size(ObSpecMag, 1);
th = max(ObSpecMag)-peakTh;                                                 % the global threshold
if(IfPlot)
    figure; plot(ObSpecMag);
    hold on; plot(ones(1,SpecLength)*th, 'k'); 
end

% calculate smoothed log amplitude spectrum by moving average
movL = movL-1+mod(movL,2);                                                  % force movL to be odd
ObSpecMagSmoothed = MovingAve(ObSpecMag, movL);
% ObSpecMagSmoothed = smooth(ObSpecMag, movL, 'moving');                    % "smooth" function is in Curve Fitting toolbox. 
if(IfPlot) 
    hold on; plot(ObSpecMagSmoothed, 'r');
    hold on; plot(ObSpecMagSmoothed + peakTh_rel, 'k');
end

ObSpecMagDiff = ObSpecMag - ObSpecMagSmoothed;                              % calculate the relative log amplitude spectrum
if(IfPlot) 
    figure; 
    plot(ObSpecMagDiff);
    close;
end

for i = 2:peakTh_freq
    if ObSpecMag(i) < th || ObSpecMagDiff(i) < peakTh_rel                   % a peak should satisfy the global threshold and local threshold
        continue;
    end
    if ObSpecMag(i) < ObSpecMag(i+1) || ObSpecMag(i) < ObSpecMag(i-1)       % a peak should be a local maximum (in localRange)
        continue;
    end
    if ObSpecMag(i) == max(ObSpecMag(max(i-localRange,1):min(i+localRange,SpecLength)))
        peakNum=peakNum+1;
        peak(peakNum)=i;
    end
end
peak = peak';                                                               % a column vector, peak=1 means that a peak appears at the first positive frequency bin
if IfPlot
    hold on; plot(peak, peakAmp, 'ro');
end