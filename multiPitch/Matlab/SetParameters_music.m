% Set parameters for music
% "%--" before a parameter indicates that this parameter has been initialized before but used here as well
% "%--" before a comment of a parameter indicates that this parameter is important
%
% Author: Zhiyao Duan
% Created: 6/20/2012
% Last modified: 5/24/2013

% parameters for preprocessing input signal
frameLen = 2048;                            %-- 46ms for fs=44100
zpf = 4;                                    % zero padding factor
win = hamming(frameLen, 'periodic');        % window function
hop = 441;                                  %-- 10ms for fs=44100
RMS_th = sqrt(0.075);                       % frame whose RMS lower than this threshold will be considered silent

% parameters for peak extraction
peakTh = 50;                                % global amplitude threshold (dB)
peakTh_rel = 4;                             % local amplitude threshold (dB)
peakTh_freq = 20000;                        % frequency threshold (Hz)
movL = 400;                                 % moving average width (Hz)
localRange = 0;                             % the range for considering local maximum

% ---------- Estimation ----------------
% parameters for estimating F0s in each single frame
bPeakFreqGMM = 1;                           % 0: modeling frequency deviation using single Gaussian, 1: using GMM
bPeakAmp = 1;                               % 0: not consider peak amplitude, 1: consider peak amplitude
bSpuriousPeak = 1;                          % 0: consider all the peaks as normal peaks, 1: consider spurious peak
bNonpeak = 1;                               % 0: not consider non-peak area, 1: consider non-peak area
bNonpeakProb = 1;                           % 0: constant value, 1: learned value (MAYBE I SHOULD USE 0 HERE)
bMask = 0;                                  % 0: not consider mask are, 1: consider mask area
bHarmonicPrior = 0;                         % 0: not use harmonic prior in peak region likelihood calculation, 1: use harmonic prior. Using it will deemphasize the peak region part.
maxHarm = 50;                               % upper bound of the harmonic number
maxF0Num = 6;                               %-- the number of estimated pitches in each frame in the first place
midiMin = note2midinum('C2');               %-- lowest possible frequency of F0 (midi number)
midiMax = note2midinum('B6');               %-- highest possible frequency of F0 (midi number)
dupF0Th = 0.5;                              % frequency threshold (in midi number) to decide if two F0 candidates are duplicate
f0step = 0.2;                               % F0 search step (midi number)
load(fullfile('model', 'PeakStatValue_music44p1k.mat')); % learned statistics of peak-f0 relations from polyphonic data
load(fullfile('model', 'HarmStatValue_music44p1k.mat')); % learned statistics of harmonic-f0 relations from monophonic data

% parameter for estimating instaneous polyphony
bInstPolyEst = 1;                           % 0: assume each time frame has trackNum of concurrent pitches, 1: estimate instaneous polyphony
para_poly = 0.88;                           % polyphony estimation threshold, percentage of the whole likelihood increase

% parameter for refining the pitches using neighboring frames
bRefine = 1;                                % 0: not refine F0 estiamtes, 1: refine F0 estiamtes using neighbouring frames
binSize = 1;                                % frequency bin size in semitones of the histogram
neigSize = 9;                               %-- the radius of the neighborhood

% parameters for second refinement of the estimated pitches: remove some
% outliers and fill some holes between pitch estimates
bSecondRefine = 0;                          % 0: not refine F0 estimates again, 1: refine F0 estimates again, by removing some outliers and filling some gaps
MSL_pd = 0.3;                               %-- pitch difference threshold of must-link (midi number)
mergeNoteGap = 100;                         %-- threshold for note formation (ms), two notelets with gap less than this threshold can be merged
minNoteLength = 100;                        %-- threshold for note length (ms), notelets shorter than this will be removed