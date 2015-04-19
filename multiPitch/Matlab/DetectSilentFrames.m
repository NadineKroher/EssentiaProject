function FrameIndex = DetectSilentFrames(wavData, para)
% Detect silent frames
% Input:
%   - wavData           : wave form mono
%   - para
%       - frameLen      : frame length in points
%       - hop           : hop size in points
%       - window        : a window function vector of frameLength
%       - RMS_th        : RMS threshold to decide if a frame is silent or not
% Output:
%   - FrameIndex        : indicate of which frame is not a noise frame
%
% Author: Zhiyao Duan
% Created: 7/18/2012
% Last modified: 7/18/2012

frameLen = para.frameLen;
hop = para.hop;
win = para.win;
RMS_th = para.RMS_th;

L = length(wavData);
TotalFrameNum = floor((L-frameLen)/hop + 1);            % the total number of frames

FrameIndex = zeros(1, TotalFrameNum);                   % to store the indicative variable that if the frame has enough energy

for frameNum = 1:TotalFrameNum
    startp = hop*(frameNum-1) + 1;
    endp = startp + frameLen - 1;
    CurFrame = wavData(startp:endp) .* win;
    
    if( sqrt(mean(CurFrame.^2)) > RMS_th )              % if the RMS is large enough, then mark it
        FrameIndex(frameNum) = 1;
    end
end