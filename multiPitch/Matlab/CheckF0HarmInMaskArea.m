function f0HarmInMaskArea = CheckF0HarmInMaskArea(fpeak, f0, maxHarm)
% check if the f0 and its harmonics are in the masking area
% return a row vector whose elements are 1 (in) or 0 (not in)
% Input
%   fpeak           - peak frequency (Hz), a column vector
%   f0              - F0 (Hz)
%   maxHarm         - upper bound of the number of harmonics
% Output
%   f0HarmInMaskArea    - a boolean vector
%
% Author: Zhiyao Duan
% Created: 06/04/2012
% Last modified: 06/04/2012

span_fharm = repmat(f0 * (1:maxHarm), length(fpeak), 1);
span_fpeak = repmat(fpeak, 1, maxHarm);
span_err = abs(span_fharm - span_fpeak);

min_span_err = min(span_err, [], 1);
f0HarmInMaskArea = (min_span_err < 30);
