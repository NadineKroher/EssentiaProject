function EstF0s = Notes2F0TimeFreqMat(EstNote, para)
% function EstF0s = Notes2F0TimeFreqMat(EstNote, para)
% Convert the note list representation to the F0 time-frequency matrix
% representation
%
% Input
%   - EstNote       : the note list representation, a column cell with each
%                       cell being a 2*n matrix. The first row is frame number, the second row
%                       is frequency in MIDI number
%   - para
%       - FrameIndex: which frames are silent, used to know the number of frames here
%       - trackNum  : a possible F0 number upper bound in each frame
% Output
%   - EstF0s        : the time-frequency matrix representation of F0s
%
% Author: Zhiyao Duan
% Created: 7/25/2012
% Last modified: 7/25/2012

FrameIndex = para.FrameIndex;
trackNum = para.trackNum;

EstF0s = zeros(trackNum*2, length(FrameIndex));
EstF0Num = zeros(1, length(FrameIndex));

for i = 1:length(EstNote)
    for j = 1:size(EstNote{i}, 2)
        time = EstNote{i}(1, j);
        freq = EstNote{i}(2, j);
        EstF0Num(time) = EstF0Num(time) + 1;
        EstF0s(EstF0Num(time), time) = freq;
    end
end

realMaxF0Num = min(trackNum, length(sum(EstF0s, 2)~=0));
EstF0s = EstF0s(1:realMaxF0Num, :);
