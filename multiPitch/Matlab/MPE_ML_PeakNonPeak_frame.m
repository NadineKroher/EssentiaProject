function [FrAllF0, FrAllF0MLogLike] = MPE_ML_PeakNonPeak_frame(peak, peakAmp, peakRelAmp, para)
% MPE estimation using ML+BIC in each frame
%
% Input: 
%   - peak              : detected peak frequencies in MIDI number, a column vector
%   - peakAmp           : detected peak amplitudes in dB, a column vector
%   - peakRelAmp        : detected peak relative log amplitude in each frame, a column vector
% Output: 
%   - FrAllF0           : estimated F0s (in total maxF0Num) in this frame, a column vector
%   - FrAllF0MLogLike   : minus log likelihood obtained by each iteration, a column vector
%
% Author: Zhiyao Duan
% Created: 6/3/2012
% Last modified: 6/4/2012

% parameters
maxF0Num = para.maxF0Num;                                           % the upper bound of polyphony
midiMin = para.midiMin;                                             % the lowest possible frequency of F0 (in midi number)
midiMax = para.midiMax;                                             % the highest possible frequency of F0 (in midi number)
f0step = para.f0step;                                               % F0 search step (in midi number)
Spur_af_Miu = para.Spur_af_Miu;                                     % spurious peak frequency and amplitude mean
Spur_af_Sigma = para.Spur_af_Sigma;                                 % spurious peak frequency and amplitude covariance
dupF0Th = para.dupF0Th;                                             % frequency threshold (in midi number) to decide if two F0 candidates are duplicate
peakNum = length(peak);

% initialization
FrAllF0 = zeros(maxF0Num, 1);
FrAllF0MLogLike = zeros(maxF0Num, 1);

% calculate the spurious peak part likelihood
x = [peakAmp, peak];
Spur_af_Miu = repmat(Spur_af_Miu, size(peakAmp, 1), 1);                      
p_SpurP = 1/(2*pi*(det(Spur_af_Sigma))^0.5)...
    * exp(-1/2 * diag((x-Spur_af_Miu)/Spur_af_Sigma*(x-Spur_af_Miu)'))*1e10;  % model the spurious peak frequencies and amplitudes using a single Gaussian 

% search F0s around the following peaks
snum = round(peakNum/3);                                            % Search space
SearchID = (1:snum)';                                               % the snum lowest frequency peaks
[~, idx] = sort(peakAmp, 'descend');                                % the snum highest absolute amplitude peaks
SearchID = [SearchID; idx(1:snum)];
% [~, idx] = sort(peakRelAmp, 'descend');                             % the snum highest relative amplitude peaks
% SearchID = [SearchID; idx(1:snum)];
SearchID = unique(SearchID);
% form the search space
SearchSpace = [];
for i = 1:length(SearchID)
    peakbin = f0step*round(peak(SearchID(i))/f0step);
    SearchSpace = [SearchSpace, peakbin-1:f0step:peakbin+1, peak(SearchID(i))];
end
SearchSpace = unique(SearchSpace);
idx = SearchSpace>midiMin & SearchSpace<midiMax;
SearchSpace = SearchSpace(idx);

if isempty(SearchSpace)
    return;
end

% start the iterative process of estimating F0s
old_MLogLike = Inf;                                                     % initialize the minus log likelihood as infinity, i.e. no estimate is the the worst estimate
old_p_NormP = -ones(peakNum, 1)*Inf;                                    % initialize the likelihood of each normal peak
old_MLogLike_NP = 0;                                                    % no harmonics yet, therefore the non-peak region minus log likelihood is 0
for m=1:maxF0Num                                                        % iterate from 1 to maxF0Num
    [FrAllF0(m), FrAllF0MLogLike(m), p_NormP, MLogLike_NP]...           % estimate a new F0 by maximizing the increase of the likelihood
        = MPE_ML_PeakNonPeak_frame_iter(SearchSpace, peak, peakAmp, old_MLogLike, old_p_NormP, p_SpurP, old_MLogLike_NP, para);      
    old_MLogLike = FrAllF0MLogLike(m);
    old_p_NormP = p_NormP;
    old_MLogLike_NP = MLogLike_NP;
    if FrAllF0(m) == 0                                                  % no F0 can be found any more to decrease the minus log likelihood
        FrAllF0MLogLike(m:maxF0Num) = FrAllF0MLogLike(m);               % set the rest minus log likelihood to the old one
        break;
    end
    % modify SearchSpace to prevent duplicates in F0 estimates
    [~, idx] = min(abs(peak-FrAllF0(m)));
    idx = abs(SearchSpace-peak(idx))>dupF0Th & abs(SearchSpace-FrAllF0(m))>dupF0Th;
    SearchSpace = SearchSpace(idx);
end
