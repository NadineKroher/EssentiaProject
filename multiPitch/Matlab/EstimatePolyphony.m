function EstF0Num = EstimatePolyphony(EstAllF0MLogLike, para)
% Estimate the polyphony from the likelihood vector.
% The log likelihood generally increases with the number of estimated F0s. 
% The polyphony is the first number that the log likelihood increase exceeds
% the para_poly * total likelihood increase from 1 to maxF0Num
%
% Input:
%   - EstAllF0MLogLike    	: the likelihood vector, a column vector, where each element corresponds to the minus log likelihood given previous F0s
%   - para                	
%       - para_poly         : the threshold to determine the polyphony
%       - maxF0Num          : the upper bound of polyphony
%       - FrameIndex        : indicator vector of which frame is loud enough (value=1)
% Output:
%   - EstF0Num            	: the estimated polyphony
%
% Author: Zhiyao Duan
% Created: 6/6/2012
% Last modified: 6/6/2012

para_poly = para.para_poly;
maxF0Num = para.maxF0Num;
FrameIndex = para.FrameIndex;

EstF0Num = zeros(1, length(FrameIndex));

sals = -EstAllF0MLogLike;
diffsals = EstAllF0MLogLike(1, :) - EstAllF0MLogLike(maxF0Num, :);
diffsals(diffsals < 0) = 0;
S = sals(1, :) + para_poly * diffsals;                                  % the likelihood threshold

for i=1:size(sals, 2)
    if FrameIndex(i) == 0 || (min(EstAllF0MLogLike(:,i))==0 && max(EstAllF0MLogLike(:,i))==0)
        EstF0Num(i) = 0;
    else
        idx = find(sals(:, i) >= S(:, i));                              % find the numbers of F0s that get a likelihood increase larger than the threshold
        if isempty(idx)                                                 % if the likelihood does not increase at all, return polyphony=1
            idx = 1;
        end
        EstF0Num(i) = idx(1);                                           % find the smalles one
    end
end
