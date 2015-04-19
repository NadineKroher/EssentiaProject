% Calculate statistics of the existence of harmonics in the non-peak
% regions from monophonic speech recordings.
%
% This code calls several functions (gmmest.m, gmmplot.m, and gmmval.m) in
% the GMM toolbox written by Dan Ellis. These functions can be downloaded
% at http://www.ee.columbia.edu/~dpwe/muscontent/matlab/
%
% Author: Zhiyao Duan
% Created: 05/31/2012
% Last modified: 5/25/2013

%% Parameters
clc; clear; close all;

IfDebug = 0;            % in the debug mode, some figures will be plotted

frameLen = 512;         % 32ms for fs=16000Hz
hop = 160;              % 10ms for fs=16000
zpf = 4;                % zero padding factor
win = hamming(frameLen, 'periodic');    % window function
RMS_th = sqrt(0.075);   % RMS threshold to classify if a frame is silent
peakTh = 50;            % peak detection absolute value threshold (dB)
peakTh_rel = 4;         % peak detection relative value threshold (dB)
peakTh_freq = 16000;    % peak detection frequency upperbound (Hz)
movL = 400;             % peak detection moving average window length (Hz)
localRange = 0;         % peak detection threshold.
maxHarm = 50;           % maximum harmonic number considered
NoteMin = 30;           % minimum F0 considered, 46Hz
NoteMax = 66;           % maximam F0 considered, 370Hz
FharmDiff = 12 * log2((1:maxHarm)); % the relative difference between F0 and its harmonics in midi num

%% Do on each wav file
% This path contains monophonic speech recordings for training. Each
% recording (a .wav file) should have its ground-truth F0 file stored in
% the .f0 file of the same name as the .wav file. Each line of the .f0 file
% stores a single number of the F0 (in Hz) in the corresponding frame.
% Note: the length and hop sie of the time frames in the F0 file should
% match the parameters listed above.
dirPath = 'D:\Work\My_data\MultiPitch_Speech\PTDB-TUG-mixed\Train\Mix_1_0\';

PANum = zeros(1, maxHarm);      % count of harmonics that are in the peak region
TotalNum = zeros(1, maxHarm);   % total count of harmonics (considering harmonics whose frequency is lower than the Nyquist sampling frequency)
HarmDevStat = zeros(1, maxHarm);% frequency deviation from the closest peak

PANum_f0 = zeros(100, maxHarm);
TotalNum_f0 = zeros(100, maxHarm);

wavNo = 0; % the number of wav files
wavAllName = dir(dirPath);
for fnum = 1:size(wavAllName,1) % each file
    if( ~isempty(strfind(wavAllName(fnum).name, '.wav')) ) %find wav files
        wavName = wavAllName(fnum).name;
        wavFile = strcat(dirPath, wavName);
        wavNo = wavNo + 1;
        fprintf(1, 'File %d is going to be processed.\n', wavNo);
        if wavNo < 0    % to skip some wrong files, modify the right-hand-side number
            continue;
        end

        % Read the wavfile, and resample it to 16kHz
        [wavData, fs, ~] = ReadWavFile_Mono(wavFile);
        wavData = resample(wavData, 16000, fs);
        fs = 16000;
        
        % Detect silent frames
        para = {};
        para.frameLen = frameLen;
        para.hop = hop;
        para.win = win;
        para.RMS_th = RMS_th;
        FrameIndex = DetectSilentFrames(wavData, para);
        
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
            xlabel('Time (Frame number)');
            ylabel('Frequency (Hz)');
            % plot FrameIndex
            figure; plot(FrameIndex); ylim([0 2]);
        end        
        
        % Peak Detection
        para = {};
        para.peakTh = peakTh;
        para.peakTh_rel = peakTh_rel;
        para.peakTh_freq = min(fftLen/2, round(peakTh_freq*fftLen/fs));     % frequency threshold (point)
        para.movL = round(movL*fftLen/fs);                                  % moving average length (point)
        para.localRange = localRange;
        [PeakData, PeakAmpData, PeakRelAmpData, PeakNum] = CalculatePeakData(SpecData, FrameIndex, para);
        PeakData = PeakData * fs / fftLen;                                  % change to Hz
        if IfDebug == 1
            % plot peaks
            figure; plot((1:fftLen/2)*fs/fftLen, 20*log10(abs(SpecData(1:fftLen/2,296))));
            hold on; plot(PeakData(:,296), PeakAmpData(:,296), 'ro')
        end

        % Read ground-truth pitches
        GTF0s = dlmread(strrep(wavFile, '.wav', '.f0'));                    % in Hz
        GTF0s = GTF0s';
        
        % Collect peak-f0 relations
        FrameNum = min(size(GTF0s,2), length(FrameIndex));
        idx = (sum(GTF0s(:,1:FrameNum), 1)~=0 & FrameIndex(1:FrameNum)~=0);
        UseFrames = 1:FrameNum;
        UseFrames = UseFrames(idx);
        for i = UseFrames
            tmpPeakData = hz2midi(PeakData(1:PeakNum(i), i));
            tmpf0 = hz2midi(GTF0s(i));
            
            % check if the f0 and its harmonics exist (smaller than the nyguist frequency)
            f0HarmExist = CheckF0HarmExist(midi2hz(tmpf0), maxHarm, fs);        
            TotalNum = TotalNum + f0HarmExist;
            TotalNum_f0(round(tmpf0), :) = TotalNum_f0(round(tmpf0), :) + f0HarmExist;
            
            % check if the f0 and its harmonics are in the peak area
            tmpmFharm = tmpf0 + FharmDiff;
            f0HarmInPeakArea = CheckF0HarmInPeakArea(tmpPeakData, tmpmFharm);
            f0HarmExistInPeakArea = f0HarmInPeakArea .* f0HarmExist;
            PANum = PANum + f0HarmExistInPeakArea;
            PANum_f0(round(tmpf0), :) = PANum_f0(round(tmpf0), :) + f0HarmExistInPeakArea;            
        end
       
        if mod(wavNo, 10) == 0
            save 'model_learning\HarmStat_speech_temp.mat' PANum TotalNum;
        end
    end
end
save 'model_learning\HarmStat_speech16k.mat' PANum TotalNum PANum_f0 TotalNum_f0;   % raw data
delete 'model_learning\HarmStat_speech_temp.mat';

%% Calculate statistics from the raw data
load 'model_learning\HarmStat_speech16k.mat';

% p(e=1 | h): p_e_h
p_e_h = PANum ./ TotalNum;
TotalNum_f0(TotalNum_f0 == 0) = 10000;
% p(e=1 | f0,h): p_e_F0h
p_e_F0h = PANum_f0 ./ TotalNum_f0;
p_e_F0h(NoteMin, :) = 0.9999;
p_e_F0h(p_e_F0h == 1) = 0.9999;
h = fspecial('gaussian',[3 3],1.5); % Gaussian smoothing
p_tmp = conv2(p_e_F0h, h, 'same');
p_e_F0h(NoteMin+1:NoteMax-1, 2:49) = p_tmp(NoteMin+1:NoteMax-1, 2:49);

save 'model_learning\HarmStatValue_speech16k.mat' p_e_h p_e_F0h;

%% plot p(e=1 | f0,h): p_e_F0h
figure; h = axes('FontSize',30);
[X,Y] = meshgrid(1:maxHarm, 36:95);
Z = 1 - p_e_F0h(36:95, :);
Z(60,50) = 1.1;
mesh(X,Y,Z); colormap(gray); 
xlabel('Harmonic number','FontSize',30); 
ylabel('F0 (MIDI number)','FontSize',30); 
zlabel('Probability','FontSize',30);
axis([0 50 36 95 0 1]);
