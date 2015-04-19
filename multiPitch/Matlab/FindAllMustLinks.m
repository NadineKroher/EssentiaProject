function MustLink = FindAllMustLinks(F0Location, para)
% Form must-links from estimated F0s
% Impose a must-link between two F0s in adjacent frames whose frequencies differ less than Thresh_PeakF0Belong
% Input
%   - F0Location    : Location of each F0, each row corresponds to a F0 
%   - para          : parameters
% Output
%   - MustLink      : Must links
%
% Author: Zhiyao Duan
% Created: 5/28/2009
% Last modified: 6/15/2012

fprintf('Finding all must-links...');

MSL_pd = para.MSL_pd;                                               % pitch difference threshold of must-link
feaNum = size(F0Location, 1);

MustLink = cell(feaNum, 1);

for i = 1:feaNum
    if mod(i, 10) == 0                                              % print progress
        fprintf('.');
    end
    if mod(i, 1000) == 0
        fprintf('\n');
    end
    
    % allocate memory, the number of must-links associated to pitch i will
    % not exceed the total number of pitches
    MustLink{i} = zeros(1, feaNum);                                 
    MustLinkNum = 0;                                                % record the number of must-links associated to pitch i

    Atime = F0Location(i, 1);
    Afreq = F0Location(i, 2);
    
    for j = 1:feaNum
        Btime = F0Location(j, 1);
        Bfreq = F0Location(j, 2);
        
        if (abs(Btime-Atime) == 1) && (abs(Afreq-Bfreq) < MSL_pd)   % must-link criterion
            MustLinkNum = MustLinkNum + 1;
            MustLink{i}(MustLinkNum) = j;
        end
    end
    
    MustLink{i} = MustLink{i}(1:MustLinkNum);
end
fprintf('\n');