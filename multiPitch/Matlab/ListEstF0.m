function F0Info = ListEstF0(EstF0)
% function F0Info = ListEstF0(EstF0)
% Convert estimated F0s matrix to a list
%
% Input
%   - EstF0     : estimated F0s matrix, each column is a frame
% Output
%   - F0Info    : estimated F0s list, each row is [frame number, F0 in midi]
%
% Author: Zhiyao Duan
% Created: 7/25/2012
% Last modified: 7/25/2012

F0Info = zeros(length(EstF0(EstF0~=0)), 2);

[m, n] = size(EstF0);
k = 0;
for i = 1:n
    for j = 1:m
        if EstF0(j, i) ~= 0
            k = k + 1;
            F0Info(k, :) = [i, EstF0(j, i)];
        end
    end
end