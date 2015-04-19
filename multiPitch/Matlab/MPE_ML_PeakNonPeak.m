function [EstAllF0, EstAllF0MLogLike] = MPE_ML_PeakNonPeak(PeakData, PeakAmpData, PeakRelAmpData, FrameIndex, para)
% Estimate F0s is all the frame, using ML+BIC
%
% Input:
%   PeakData            - peak frequencies in MIDI number of different frames, 100 * frame no.
%   PeakAmpData         - peak amplitudes in dB of different frames, 100 * frame no.
%   PeakRelAmpData      - peak relative amplitude in dB of different frames, 100 * frame no.
%   para                - parameters
% Output:
%   EstAllF0            - estimated F0s (maxF0Num * L)
%   EstAllF0MLogLike    - the minus log likelihood up to the estiamted F0s (maxF0Num * L)
%
% Author: Zhiyao Duan
% Created: 6/3/2012
% Last modified: 6/3/2012

fprintf('Estimating pitches in each frame...');

maxF0Num = para.maxF0Num;                                         	% upper bound of polyphony
midiMin = para.midiMin;                                            	% lowest possible frequency of F0s

L = size(PeakData, 2);                                          	% frame number

% initialization
EstAllF0 = zeros(maxF0Num, L);
EstAllF0MLogLike = zeros(maxF0Num, L);

for FrameNum = 1:L
%     if FrameNum == 30
%         a=1;
%     end
    fprintf('.');
    if mod(FrameNum, 100) == 0
        fprintf('\n');
    end
    
    if FrameIndex(FrameNum) == 1                                  	% process if the energy of the frame is large enough
        peak = PeakData(:, FrameNum);
        peakAmp = PeakAmpData(:, FrameNum);
        peakRelAmp = PeakRelAmpData(:, FrameNum);
        
%         idx = peak~=0;                                          	% because the peak number in each frame might be less than 100, so there are some 0s
%         peak = peak(idx);
%         peakAmp = peakAmp(idx);
%         peakRelAmp = peakRelAmp(idx);
        
        idx = peak>=midiMin;                                        % only use the peaks whose frequencies are larger than the minimum possible F0
        peak = peak(idx);
        peakAmp = peakAmp(idx);
        peakRelAmp = peakRelAmp(idx);
        
        if ~isempty(peak)
            [FrAllF0, FrAllF0MLogLike] = MPE_ML_PeakNonPeak_frame(peak, peakAmp, peakRelAmp, para); % estiamte F0s for each frame
            EstAllF0(:,FrameNum) = FrAllF0;                       	% gather estimates from different frames
            EstAllF0MLogLike(:,FrameNum) = FrAllF0MLogLike;
        end
    end
end

fprintf('\n');