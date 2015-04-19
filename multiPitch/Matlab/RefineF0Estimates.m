function [RefinedF0, RefinedF0Num] = RefineF0Estimates(EstF0, EstF0Num, para)
% Refine the estiamted F0s by majority vote
%
% Input:
%   - EstF0       	: Estimated all F0 matrix
%   - EstF0Num    	: Estimated polyphony for each frame
%   - para        	
%       - midiMin   : the minimum possible F0
%       - midiMax   : the maximum possible F0
%       - binSize   : the size of each frequency bin in semitone in the histogram (binSize=1 for music)
%       - neigSize  : the neighborhood size (#frames, radius not including itself) in the histogram
% Output:
%   - RefinedF0   	: Refined all F0 results
%   - RefinedF0Num	: Refined polyphony
%
% Author: Zhiyao Duan
% Last modified: 6/6/2012

midiMin = para.midiMin;
midiMax = para.midiMax;
binSize = para.binSize;
neigSize = para.neigSize;                                   

[maxPoly, FrameNum] = size(EstF0);
RefinedF0 = zeros(maxPoly, FrameNum);                       % make RefinedF0 the same size as EstF0
RefinedF0Num = zeros(1, FrameNum);                          % initialize RefinedF0Num
hstSize = ceil((midiMax-midiMin)/binSize + 1);              % the size of the histogram

for i=1:FrameNum
    startFrame = max(i-neigSize, 1);                        % the first frame of the neighbourhood
    endFrame = min(i+neigSize, FrameNum);                   % the last frame of the neighbourhood
    
    weightHst = zeros(hstSize, 1);
    weightPitch = zeros(hstSize, 1);
    
    total_weight = 0;
    for j = startFrame:endFrame
        tmp_weight = 1 - abs(j-i)/(neigSize+1);             % weight that this frame contributes
        total_weight = total_weight + tmp_weight;           % total weights of these frames
        RefinedF0Num(i) = RefinedF0Num(i) + tmp_weight * EstF0Num(j);   % weighted sum of polyphony
        
        for k = 1:EstF0Num(j)
            tmp_bin = round((EstF0(k, j)-midiMin)/binSize + 1);     % the bin that EstF0(k,j) falls in
            weightHst(tmp_bin) = weightHst(tmp_bin) + tmp_weight;   % calculate the weight histogram
            weightPitch(tmp_bin) = weightPitch(tmp_bin) + EstF0(k,j) * tmp_weight;  % weighted sum of the pitches in each bin of the histogram
        end
    end
    
    % refine polyphony
    RefinedF0Num(i) = round(RefinedF0Num(i) / total_weight);
    
    % refine pitches
    weightPitch = weightPitch ./ weightHst;
    weightPitch(weightHst == 0) = 0;                            % eliminate NaN values
    [~, robustBins] = sort(weightHst, 'descend');               % sort the weight histogram by descending order
    RefinedF0(1:maxPoly, i) = weightPitch(robustBins(1:maxPoly));
    % Replace the refined pitch with the one estimated in the current frame, if they are similar
    for j = 1:maxPoly
        for k = 1:EstF0Num(i)
            currbin = round((EstF0(k,i)-midiMin)/binSize + 1);  % the bin that EstF0(k,i) falls in
            if robustBins(j) == currbin
                RefinedF0(j, i) = EstF0(k, i);
                break;
            end
        end
    end
end