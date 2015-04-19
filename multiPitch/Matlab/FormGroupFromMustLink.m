function Group = FormGroupFromMustLink(MustLink)
% Form vertex groups from must-links by breath-first search
% Input
%   - MustLink      : must-link lists for each vertex
% Output
%   - Group         : vertex groups, each row is a group, containing the
%                     indicies of the pitches, in an ascending order
%                     (therefore ascending order of time also)
% 
% Author: Zhiyao Duan
% Created: 5/28/2009
% Last modified: 6/15/2012

Group = cell(1000, 1);
groupNum = 0;
feaNum = length(MustLink);

bUsedMustLink = zeros(feaNum, 1);                   % mark the used must-links

for i = 1:length(MustLink)
    if bUsedMustLink(i) == 0
        % start a new group
        groupNum = groupNum + 1;
        Group{groupNum} = zeros(1, feaNum);
        
        % the first point
        pointNum = 1;                               % the number of points in this group
        Group{groupNum}(pointNum) = i;
        
        % breath-first search to add must-link points into this group
        j = 1;
        while j <= pointNum
            for k = MustLink{Group{groupNum}(j)}
                if ~any(Group{groupNum}==k)         % prevent duplicate     
                    pointNum = pointNum + 1;        % enlarge the current group
                    Group{groupNum}(pointNum) = k;
                end
            end
            bUsedMustLink(Group{groupNum}(j)) = 1;  % this must-link has been used
            j = j + 1;
        end
        
        Group{groupNum} = Group{groupNum}(1:pointNum);
    end
end
Group = Group(1:groupNum);
