function f0HarmExist = CheckF0HarmExist(f0, maxHarm, fs)
% Check if f0 and its harmonics exist in the spectrum
% return a row vector whose elements are 1 (exist) or 0 (not exist)
%
% Input:
%   - f0            : F0 (Hz)
%   - maxHarm       : upper bound of the harmonic number
%   - fs            : sampling frequency (Hz)
% Output:
%   - f0HarmExist   : a boolean vector
%
% Author: Zhiyao Duan
% Created: 06/04/2012
% Last modified: 06/04/2012

f0HarmExist = zeros(1, maxHarm);
tmpf = f0 * (1:maxHarm);
f0HarmExist(tmpf < fs/2) = 1;
