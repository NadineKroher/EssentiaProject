function doMPE(wavFile, outFile, dataType, IfDebug)
% Multiple F0 estimation in each audio frame
%
% Input:
%   - wavFile       : the audio needs to be processed
%   - outFile       : the mpe result file, each row is a frame, values are
%                     in MIDI number
%   - dataType      : 'music' or 'speech'
%   - IfDebug       : 0 - not debug; 1 - debug mode, a number of figures
%                     will be shown
%
% Example usage:
%   doMPE('example\01-AchGottundHerr-short.wav', 'est.f0', 'music', 0);
% or
%   doMPE('example\022_M10_si2182_F02_si803_F10_si2267_cut.wav', 'est.f0', 'speech', 0);
%
% Author: Zhiyao Duan
% Created: 6/3/2012
% Last modified: 5/30/2013

fprintf('\nMulti-pitch estimation on file %s...\n', wavFile);

if nargin < 4
    IfDebug = 0;
end
if IfDebug == 1
    close all;
end

%% Set parameters
% The following scripts set all the parameters for the given datatype. One can
% modify these parameters by editing the corresponding script.
switch dataType
    case 'music'
        SetParameters_music;    % parameters loaded here are good for 44.1kHz music
    case 'speech'
        SetParameters_speech;   % parameters loaded here are good for 16kHz speech
    otherwise
        fprintf('Unknown data type!\n');
        return;
end

%% Multi-pitch Estimation
fprintf('----------Step 1: multi-pitch estimation----------\n');
% Read the wavfiles
[wavData, fs, ~] = ReadWavFile_Mono(wavFile);
if strcmp(dataType, 'speech')
    wavData = resample(wavData,16000,fs);
    fs = 16000;
end
% wavData = wavData(1:fs*5);                                          % TEMP: chunk first 5 seconds for debugging

% Detect silent frames
para = {};
para.frameLen = frameLen;
para.hop = hop;
para.win = win;
para.RMS_th = RMS_th;
FrameIndex = DetectSilentFrames(wavData, para);
if IfDebug == 1
    % plot FrameIndex
    figure; plot(FrameIndex); ylim([0 2]);
end

% STFT
fftLen = 2^nextpow2(frameLen*zpf);                                  % zero padding
para = {};
para.frameLen = frameLen;
para.hop = hop;
para.win = win;
para.fftLen = fftLen;
SpecData = CalculateSTFT(wavData, para);
if IfDebug == 1
    % plot spectrum
    figure;
    imagesc(1:length(FrameIndex), (1:fftLen/2)*fs/fftLen, 20*log10(abs(SpecData(1:fftLen/2, :))));
    set(gca,'YDir','normal');
    title('Spectrogram and ground-truth pitches');
    xlabel('Time (Frame number)');
    ylabel('Frequency (Hz)');
    ylim([0 2000]);
%     % plot the ground-truth pitch
%     % uncomment this block if you have the ground-truth pitches
%     GTF0s = dlmread(strrep(wavFile, '.wav', '.f0'));                % in Hz
%     GTF0s = GTF0s';
%     hold on;
%     plot(1:size(GTF0s,2), GTF0s, '.');
end

% Peak Detection
para = {};
para.peakTh = peakTh;
para.peakTh_rel = peakTh_rel;
para.peakTh_freq = min(fftLen/2, round(peakTh_freq*fftLen/fs));     % frequency threshold (point)
para.movL = round(movL*fftLen/fs);                                  % moving average length (point)
para.localRange = localRange;
[PeakData, PeakAmpData, PeakRelAmpData, PeakNum] = CalculatePeakData(SpecData, FrameIndex, para);
idx = PeakData==0;
PeakData = hz2midi(PeakData*fs/fftLen);                             % change to midi number
PeakData(idx) = 0;
if IfDebug == 1
    % plot peaks
    tempFrame = 200;
    figure; plot((1:fftLen/2)*fs/fftLen, 20*log10(abs(SpecData(1:fftLen/2,tempFrame))));
    hold on; plot(midi2hz(PeakData(1:PeakNum(tempFrame),tempFrame)), PeakAmpData(1:PeakNum(tempFrame),tempFrame), 'ro')
end

% Estimate F0s
para = {};
% system configurations
para.bPeakFreqGMM = bPeakFreqGMM;
para.bPeakAmp = bPeakAmp;
para.bSpuriousPeak = bSpuriousPeak;
para.bNonpeak = bNonpeak;
para.bNonpeakProb = bNonpeakProb;
para.bMask = bMask;
para.bHarmonicPrior = bHarmonicPrior;
% parameters
para.maxHarm = maxHarm;                                                 % harmonic number upper bound
para.f0step = f0step;                                                   % f0 search step: f0step * 100% of the peak frequency
para.maxF0Num = maxF0Num;                                               % polyphony upper bound
para.fs = fs;
para.midiMin = midiMin;                                                 % f0 lower bound (e.g. midinum('C2');
para.midiMax = midiMax;                                                 % fo higher bound (e.g. midinum('B6');
para.dupF0Th = dupF0Th;                                                 % frequency threshold (in midi number) to decide if two F0 candidates are duplicate
% learned parameters
para.p_s = p_s;                                                         % spurious peak probability
para.Spur_af_Miu = Spur_af_Miu;                                         % spurious peak mean vector
para.Spur_af_Sigma = Spur_af_Sigma;                                     % spurious peak covariance matrix
para.d_mo = d_mo;                                                       % normal peak GMM mean vectors
para.d_vo = d_vo;                                                       % normal peak GMM variance vectors
para.d_co = d_co;                                                       % normal peak GMM component weights vectors
para.p_a_fh = p_a_fh;                                                   % p(a|f,h) of normal peaks
% para.p_h_f0 = p_h_f0;                                                 % P(h|F0) of harmonics
para.p_e_h = p_e_h;                                                     % P(D|h) of harmonics
para.p_e_F0h = p_e_F0h;                                                 % P(D|h, f0) of harmonics
[EstAllF0, EstAllF0MLogLike] = MPE_ML_PeakNonPeak(PeakData, PeakAmpData, PeakRelAmpData, FrameIndex, para);
if IfDebug==1
    figure;
    imagesc(1:length(FrameIndex), (1:fftLen/2)*fs/fftLen, 20*log10(abs(SpecData(1:fftLen/2, :))));
    set(gca,'YDir','normal');
    title('Before polyphony estimation');
    xlabel('Time (Frame number)');
    ylabel('Frequency (Hz)');
    ylim([0 2000]);
    hold on;
    plot(1:size(EstAllF0,2), midi2hz(EstAllF0), '.');
end

% Estimate polyphony
if bInstPolyEst == 1
    para = {};
    para.para_poly = para_poly;
    para.maxF0Num = maxF0Num;
    para.FrameIndex = FrameIndex;
    EstF0Num = EstimatePolyphony(EstAllF0MLogLike, para);
    if IfDebug==1
        figure;
        imagesc(1:length(FrameIndex), (1:fftLen/2)*fs/fftLen, 20*log10(abs(SpecData(1:fftLen/2, :))));
        set(gca,'YDir','normal');
        title('After polyphony estimation');
        xlabel('Time (Frame number)');
        ylabel('Frequency (Hz)');
        ylim([0 2000]);
        hold on;
        EstF0 = zeros(maxF0Num, size(EstAllF0,2));
        for i = 1:length(FrameIndex)
            if EstF0Num(i) > 0
                EstF0(1:min(EstF0Num(i),maxF0Num), i) = EstAllF0(1:min(EstF0Num(i),maxF0Num), i);
            end
        end
        plot(1:size(EstF0,2), midi2hz(EstF0), '.');
    end
else
    EstF0Num = min(sum(EstAllF0~=0, 1), maxF0Num*ones(size(FrameIndex)));
end

% Refine F0 estiamtes
if bRefine == 1
    para = {};
    para.midiMin = midiMin;
    para.midiMax = midiMax;
    para.binSize = binSize;
    para.neigSize = neigSize;                                           % neighbour frame radius to refine current F0 estimates
    [RefinedAllF0, RefinedF0Num] = RefineF0Estimates(EstAllF0, EstF0Num, para);
    EstAllF0 = RefinedAllF0;
    EstF0Num = RefinedF0Num;

    if IfDebug==1
        figure;
        imagesc(1:length(FrameIndex), (1:fftLen/2)*fs/fftLen, 20*log10(abs(SpecData(1:fftLen/2, :))));
        set(gca,'YDir','normal');
        title('After refinement');
        xlabel('Time (Frame number)');
        ylabel('Frequency (Hz)');
        ylim([0 2000]);
        EstF0 = zeros(maxF0Num, size(EstAllF0,2));
        for i = 1:length(FrameIndex)
            if EstF0Num(i) > 0
                EstF0(1:min(EstF0Num(i),maxF0Num), i) = EstAllF0(1:min(EstF0Num(i),maxF0Num), i);
            end
        end
        hold on;
        plot(1:size(EstF0,2), midi2hz(EstF0), '.');
    end
end

EstF0 = zeros(maxF0Num, size(EstAllF0,2));
for i = 1:length(FrameIndex)
    if EstF0Num(i) > 0
        EstF0(1:min(EstF0Num(i),maxF0Num), i) = EstAllF0(1:min(EstF0Num(i),maxF0Num), i);
    end
end

% Second refinement of F0 estimates
% Remove outliers, fill gaps
if bSecondRefine == 1
    % List all estimated F0s and their frames
    F0Info = ListEstF0(EstF0);
    
    % Find pairwise must-link constraints and form must-link groups
    para = {};
    para.MSL_pd = MSL_pd;                                               % Frequency deviation threshold to decide must links
    AllMustLink = FindAllMustLinks(F0Info, para);                   % Find all must-links that the domain knowledge can provide
    PreGroup = FormGroupFromMustLink(AllMustLink);                      % Form groups from must-links
    
    % Merge groups into notes and remove small groups
    para = {};
    para.mergeNoteGap = mergeNoteGap*fs/1000/hop;                   % change to frame number
    para.minNoteLength = minNoteLength*fs/1000/hop;                 % change to frame number
    EstNote = FormingNotesForSomeTrack(0, PreGroup, zeros(length(PreGroup),1), F0Info, para);

    % Back to F0 time-frequency representation
    para = {};
    para.FrameIndex = FrameIndex;
    para.trackNum = maxF0Num;
    EstF0 = Notes2F0TimeFreqMat(EstNote, para);
    
    if IfDebug==1
        figure;
        imagesc(1:length(FrameIndex), (1:fftLen/2)*fs/fftLen, 20*log10(abs(SpecData(1:fftLen/2, :))));
        set(gca,'YDir','normal');
        title('After second refinement');
        xlabel('Time (Frame number)');
        ylabel('Frequency (Hz)');
        ylim([0 2000]);
        hold on;
        plot(1:size(EstF0,2), midi2hz(EstF0), '.');
    end    
end

% output multi-pitch estimation result
idx = sum(EstF0~=0, 2)~=0;
EstF0 = EstF0(idx, :);
dlmwrite(outFile, EstF0', '\t');

fprintf('Done!\n');