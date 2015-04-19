function [NewF0, NewMLogLike, p_NormP, MLogLike_NP] = MPE_ML_PeakNonPeak_frame_iter(SearchSpace, peak, peakAmp, old_MLogLike, old_p_NormP, p_SpurP, old_MLogLike_NP, para)
% One greedy search iteration of the Maximum Likelihood-based MPE method.
% It adds one F0 which increases the likelihood most to the set of
% estimated F0s.
%
% Input:
%   - peak        	: detected peak frequencies in MIDI number, a column vector
%   - peakAmp     	: detected peak amplitudes in dB, a column vector
%   - old_MLogLike  : old minus log likelihood, prior calculating the new F0
%   - old_p_NormP   : old normal peak part linear likelihoods of each peak, could be [] or a column vector
%   - P_SpurP       : spurious peak part linear likelihood of each peak, could be [] or a column vector
%   - old_MLogLike_NP   : old non-peak region minus log likelihood
%   - para          : parameters
% Output:
%   - NewF0         : the newly estimated F0
%   - NewMLogLike   : the newly updated minus log likelihood after adding the current F0 to existing F0s
%   - P_NormP       : the newly updated linear normal peak likelihood each peak, a column vector
%   - MLogLike_NP   : the newly updated minus log likelihood of the non-peak region
%
% Author: Zhiyao Duan
% Created: 10/01/2007
% Revised: 02/13/2008 -- add p_NPA (non-peak area probability)
% Revised: 12/15/2008 -- add several switches (bSpuriousPeak, etc.)
% Revised: 06/04/2012 -- change the way of generating F0 hypothesis

% system configurations
bPeakFreqGMM = para.bPeakFreqGMM;
bPeakAmp = para.bPeakAmp;
bSpuriousPeak = para.bSpuriousPeak;
bNonpeak = para.bNonpeak;
bNonpeakProb = para.bNonpeakProb;
bMask = para.bMask;
bHarmonicPrior = para.bHarmonicPrior;

% parameters
fs = para.fs;                                                       % sampling rate
maxHarm = para.maxHarm;                                             % the upper bound of the number of harmonics

% parameters learned from training data
if bSpuriousPeak == 1
    p_s = para.p_s;                                                 % prior probability of observing a spurious peak
else
    p_s = 0;
end
d_mo = para.d_mo;                                                   % the GMM model of p(d) of normal peaks, mo:mean, vo:variance, co: component weights
d_vo = para.d_vo;
d_co = para.d_co;
p_a_fh = para.p_a_fh;                                               % p(a|f,h) of normal peaks
p_e_h = para.p_e_h;                                                 % p(e=1|f0,h), the probability of detecting a harmonic, conditioned on f0 and h
p_e_F0h = para.p_e_F0h;                                             % p(e=1|f0), the probability of detecting a harmonic, only conditioned on f0

% initialization and default return values
MinMLogLike = old_MLogLike;
NewF0 = 0;
NewMLogLike = old_MLogLike;
p_NormP = 0;
MLogLike_NP = 0;

% some calculations
peakNum = length(peak);
FharmDiff = 12 * log2((1:maxHarm));                                 % the relative difference between F0 and its harmonics in midi num
fpeak = midi2hz(peak);                                              % peak frequency in Hz

for mF0 = SearchSpace
    F0 = midi2hz(mF0);

    Ratio = fpeak / F0;                                             % calculate the harmonic number of peak i given F01
    RRatio = round(Ratio);

    % p_NormP_fk: p(fk|F0)=p(dk,h|F0)=p(dk|F0)*p(h|F0), the frequency part of the normal peak likelihood
    err = peak - mF0 - 12*log2(max(RRatio,1));                      % the frequency deviation of all peaks from harmonics of F01
    if bPeakFreqGMM == 1
        p_NormP_f = gmmeval(err', d_mo', d_vo', d_co');             % model the deviation using GMM
        p_NormP_f = p_NormP_f';
    else
        p_NormP_f = 1/sqrt(2*pi)/0.06 * exp(-0.5 * (err/0.06).^2);  % model the deviation using a single gaussian (variance is set to 0.06), learned 0.389 is even worse
    end
    if bHarmonicPrior == 1
        p_NormP_h = 1/maxHarm;                                          % using a uniform distribution gets best results
    else
        p_NormP_h = 1;
    end
    
    % p_NormP_ak: p(ak|fk,F0), the amplitude part of the normal peak likelihood
    if bPeakAmp == 1
        RRatio(RRatio==0) = 1;                                      % peaks that have lower frequencies than F01 are considered the first harmonic
        RRatio(RRatio>maxHarm) = maxHarm;                           % peaks that have higher frequencies than F01 are considered the maxHarm-th harmonic.
        % these operations are to make the following calculations easier
        p_NormP_a = zeros(peakNum,1);
        for k=1:peakNum
            p_NormP_a(k) = p_a_fh(max(round(peakAmp(k)),1), round(peak(k)), RRatio(k));
        end
    else
        p_NormP_a = ones(peakNum, 1);                               % if not using the amplitude information, just set the probability to 1
    end
    
    current_p_NormP = 1e10 * p_NormP_a .* p_NormP_f * p_NormP_h;    % the current likelihood of each peak under F01
                
    idx = (old_p_NormP > current_p_NormP);                          % if a peak was better "explained" by previously estimated F0s rather than the current F0
    current_p_NormP(idx) = old_p_NormP(idx);                        % then do not update the likelihood of the peak
    
    current_MLogLike_P = sum(-log(current_p_NormP * (1-p_s) + p_SpurP * p_s)+10*log(10));   % the current minus log likelihood after adding F01 to previously estimated F0s

    % non-peak part likelihood (p_NPA)
    if bNonpeak == 1
        p_NPA = ones(1, maxHarm);
        % check if a harmonic exist in the spectrum
        f0HarmExist = CheckF0HarmExist(F0, maxHarm, fs);                       
        % check if a harmonic is in the peak area
        mFharm = mF0 + FharmDiff;
        f0HarmInPeakArea = CheckF0HarmInPeakArea(peak, mFharm);                
        
        % if outside the spectrum or in peak area, then the probability of not
        % observing the harmonic in the non-peak region area is 1
        idx = (f0HarmExist==0 | f0HarmInPeakArea==1);                       
        p_NPA(idx) = 1;
        
        % if inside the non-peak region
        if bMask == 1                                             	% consider masking situations
            % check if a harmonic is masked by others
            f0HarmInMaskArea = CheckF0HarmInMaskArea(fpeak, F0, maxHarm);
            % in both non-peak area and masking area
            idx = (f0HarmExist==1 & f0HarmInPeakArea==0) & f0HarmInMaskArea==1; 
            if bNonpeakProb == 1                                    % if use a leaned probability which is conditioned on both F0 and h
                for j = 1:maxHarm
                    if (idx(j) == 1)
                        p_NPA(j) = 1 - 0.5*p_e_F0h(round(mF0), j);   % 1 - 1/2*P(D=1|F0,h) = 1/2*(1 + 1-P(ei=1|F0,h)), two cases either being masked or masking the other
                    end
                end
            else                                                    % if use a learned probability which is only conditioned on h
                p_NPA(idx) = 1 - 0.5*p_e_h(idx);              % 1 - 1/2*P(D=1|h)
            end
            % in non-peak area, but not in masking area
            idx = (f0HarmExist==1 & f0HarmInPeakArea==0) & f0HarmInMaskArea==0;   
            if bNonpeakProb == 1                                    % if use a leaned probability which is conditioned on both F0 and h
                for j = 1:maxHarm
                    if (idx(j) == 1)
                        p_NPA(j) = 1 - p_e_F0h(round(mF0), j);   % P(D=0|F0,h) = 1 - P(D=1|F0,h)
                    end
                end
            else                                                    % if use a learned probability which is only conditioned on h
                p_NPA(idx) = 1 - p_e_h(idx);                  % P(D=0|h) = 1 - P(D=1|h)
            end
        else                                                        % do not consider masking situations
            idx = f0HarmExist==1 & f0HarmInPeakArea==0;             % in non-peak area
            if bNonpeakProb == 1                                    % if use a leaned probability which is conditioned on both F0 and h
                for j = 1:maxHarm
                    if (idx(j) == 1)
                        p_NPA(j) = 1 - p_e_F0h(round(mF0), j);   % 1 - P(D=1|F0,h)
                    end
                end
            else                                                    % if use a learned probability which is only conditioned on h
                p_NPA(idx) = 1 - p_e_h(idx);                  % P(D=0|h) = 1 - P(D=1|h)
            end
        end
        
        current_MLogLike_NP = sum(-log(p_NPA));                     % minus log likelihood of the non-peak area
    else
        current_MLogLike_NP = 0;                                    % if do not consider the non-peak area, then set it to 0 (linear likelihood to 1)
    end
    current_MLogLike_NP = current_MLogLike_NP + old_MLogLike_NP;    % current overall non-peak region likelihood, given all estimated F0s
        
    % whole likelihood
    current_MLogLike = current_MLogLike_P + current_MLogLike_NP;    % current whole likehood
%     current_MLogLike = current_MLogLike_P/peakNum + current_MLogLike_NP/maxHarm;    % current whole likehood
    
    % update return values given the best F0 found in this iteration
    if current_MLogLike < MinMLogLike
        MinMLogLike = current_MLogLike;
        % update return values
        NewF0 = mF0;
        NewMLogLike = current_MLogLike;
        p_NormP = current_p_NormP;
        MLogLike_NP = current_MLogLike_NP;
    end
end
