function f0HarmInPeakArea = CheckF0HarmInPeakArea(mpeak, mFharm, peakSize)
% check if the f0 and its harmonics are in the peak area
% return a row vector whose elements are 1 (in) or 0 (not in)
%
% Input:
%   - mpeak             : peak frequency (midi number), a column vector
%   - mFharm            : harmonic frequency (midi number), a column vector
%   - peakSize          : peak region size (default: radius of 0.5 semitone)
% Output:
%   - f0HarmInPeakArea  : a boolean vector
%
% Author: Zhiyao Duan
% Created: 06/04/2012
% Last modified: 06/04/2012

if nargin < 3
    peakSize = 0.5;
end

span_mharm = repmat(mFharm, length(mpeak), 1);
span_mpeak = repmat(mpeak, 1, length(mFharm));
span_err = abs(span_mharm - span_mpeak);

min_span_err = min(span_err, [], 1);
f0HarmInPeakArea = (min_span_err < peakSize);
