function SpecData = CalculateSTFT(wavData, para)
% Calculate STFT
% Input:
%   - wavData           : wave form mono
%   - para
%       - frameLen      : frame length in points
%       - hop           : hop size in points
%       - window        : a window function vector of frameLength
%       - fftLen        : fft length in points
% Output:
%   - SpecData          : complex spectrogram, dimension = [fftLength/2+1, total frame number]
%
% Author: Zhiyao Duan
% Last modified: 5/24/2013

frameLen = para.frameLen;
hop = para.hop;
win = para.win;
fftLen = para.fftLen;

L = length(wavData);
TotalFrameNum = floor((L-frameLen)/hop + 1);            % the total number of frames

SpecData = zeros(fftLen/2+1, TotalFrameNum);            % to store the spectrum

for frameNum = 1:TotalFrameNum
    startp = hop*(frameNum-1) + 1;
    endp = startp + frameLen - 1;
    CurFrame = wavData(startp:endp) .* win;
    
    CurSpec = fft(CurFrame, fftLen);
    SpecData(:, frameNum) = CurSpec(1:fftLen/2+1);      % take only the first half spectrum
end