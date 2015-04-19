% Calculate statistics of the existence of harmonics in the non-peak
% regions from note samples of the Iowa dataset.
%
% This code calls the YIN pitch detection toolbox, written by Alain de
% Cheveigné. This toolbox can be downloaded at
% http://www.ee.columbia.edu/~dpwe/LabROSA/doc/.
% This code also calls several functions (gmmest.m, gmmplot.m, and gmmval.m) in
% the GMM toolbox written by Dan Ellis. These functions can be downloaded
% at http://www.ee.columbia.edu/~dpwe/muscontent/matlab/
%
% Author: Zhiyao Duan
% Last modified: 5/25/2013

%% Parameters
clc; clear; close all;

IfDebug = 0;            % in the debug mode, some figures will be plotted

Num_FileToSkip = 467;   % Number of files to skip from the beginning in training
                        % This is good for avoiding files that has problems
                        % to deal with.
frameLen = 2048;        % Frame length (Points), 46ms for fs=44100Hz
hop = 1024;             % Frame hop (Points), 23ms for fs=44100Hz
zpf = 4;                % zero padding factor
window = hamming(FrameLength, 'periodic');  % window function
RMS_th = 0.075;         % RMS threshold to decide if use this file or not
peakTh = 50;            % peak detection absolute value threshold (dB)
peakTh_rel = 4;         % peak detection relative value threshold (dB)
peakTh_freq = 20000;    % peak detection frequency upperbound (Hz)
movL = 400;             % moving average width (Hz)
localRange = 0;         % the range for considering local maximum
maxHarm = 50;           % maximum harmonic number considered  
NoteMin = 36;           % minimum F0 considered, C2
NoteMax = 95;           % maximum F0 considered, B6
FharmDiff = 12 * log2((1:maxHarm)); % The relative difference between F0 and its harmonics in midi num

% YIN parameters
P.minf0 = 40;                                                       % Hz - minimum expected F0 (default: 30 Hz)
P.maxf0 = 2200;                                                     % Hz - maximum expected F0 (default: SR/(4*dsratio))
P.thresh = 0.1;                                                     % threshold (default: 0.1)
% P.relfag = 1;                                                       % if ~0, thresh is relative to min of difference function (default: 1)
P.hop = 32;                                                         % s - interval between estimates (default: 32)
% P.range = [1 4096];                                                 % samples - range of samples ([start stop]) to process
%   P.bufsize:  samples - size of computation buffer (default: 10000)
P.sr = 44100;                                                       % Hz - sampling rate (usually taken from file header)
% P.wsize = 4096;                                                     % samples - integration window size (defaut: SR/minf0)
%	P.lpf:		Hz - intial low-pass filtering (default: SR/4)
%	P.shift		0: shift symmetric, 1: shift right, -1: shift left
%	(default: 0)

%% Do on each wav file
% This path contains note samples for training. Each note sample is in a
% wave file, whose name is the note name.
dirPath = 'W:\My_data\ChordData\SampleMaterial\Training\';          % The directory of the training samples

TotalNum = zeros(1, maxHarm);                                       % Total number of valid harmonics (smaller than the Nyquist frequency) for each harmonic
PANum = zeros(1, maxHarm);                                          % Total number of valid harmonics that in the peak area for each harmonic
TotalNum_f0 = zeros(100, maxHarm);                                  % Total number of valid harmonics (smaller than the Nyquist frequency) for each harmonic, given F0
PANum_f0 = zeros(100, maxHarm);                                     % Total number of valid harmonics that in the peak area for each harmonic, given F0

load 'HarmStat_mono.mat'                                            % Skip wrong files encountered in YIN

wavNo = 0;% the number of wav files
wavAllName = dir(dirPath);
for fnum = 1:size(wavAllName,1) % each note file
    if( ~isempty(strfind(wavAllName(fnum).name, '.wav')) )          % Find wav files
        wavName = wavAllName(fnum).name;
        wavFile = strcat(dirPath, wavName);

        % Get the tag of f0 to limit the search space
        tmpposarray = strfind(wavName, '_');
        tmppos = max(tmpposarray);
        wavNote = wavName(tmppos+1:tmppos+3);
        if wavNote(3)=='.'
            wavNote(3)='';
        end
        wavMidi = note2midinum(wavNote);

        if (wavMidi < 36 || wavMidi > 95)                           % Out of range, then skip
            continue;
        end

        wavNo = wavNo + 1;
        fprintf(1, 'File %d is going to be processed.\n', wavNo);
        if wavNo <= Num_FileToSkip                                  % To skip some wrong files encountered in YIN
            continue;
        end

        % Read the wavfile, and select the middle portion of it
        [wavData, fs, ~] = ReadWavFile_Mono(wavFile);
        startp = 8193;
        endp = length(wavData)-8192;
        endp = startp + floor((endp-startp)/hop)*hop;
        wavDataSel = wavData(startp:endp+1);                        % Data selected for training
        
        % Detect silent frames
        para = {};
        para.frameLen = frameLen;
        para.hop = hop;
        para.win = win;
        para.RMS_th = RMS_th;
        FrameIndex = DetectSilentFrames(wavDataSel, para);
        
        % STFT
        fftLen = 2^nextpow2(frameLen*zpf);                                  % zero padding
        para = {};
        para.frameLen = frameLen;
        para.hop = hop;
        para.win = win;
        para.fftLen = fftLen;
        SpecData = CalculateSTFT(wavDataSel, para);
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
        para.movL = movL;
        para.localRange = localRange;
        [PeakData, PeakAmpData, PeakRelAmpData, PeakNum] = CalculatePeakData(specData, FrameIndex, para);
        PeakData = PeakData * fs / fftLen;                                  % change to Hz
        if IfDebug == 1
            % plot peaks
            figure; plot((1:fftLen/2)*fs/fftLen, 20*log10(abs(SpecData(1:fftLen/2,296))));
            hold on; plot(PeakData(:,296), PeakAmpData(:,296), 'ro')
        end
        realSegNum = size(PeakData, 2);                             % Number of frames

        % Calculate f0 using YIN
        P.minf0 = midi2hz(wavMidi-1);
        P.maxf0 = midi2hz(wavMidi+1);
        R = yin(wavDataSel, P);         %% calling YIN here. You can download it from YIN's webpage.
        f0 = zeros(1, realSegNum);
        for segnum = 1:realSegNum
            startp = 1 + (segnum-1) * hop/P.hop;
            endp = startp + FrameLength/P.hop - 1;
            tmpf0 = R.f0(startp:endp);                              % F0s in this frame
            [minAp, idx] = min(R.ap(startp:endp) ./ R.pwr(startp:endp));    % Get the best F0 estimate of this frame
            f0(segnum) = 2^tmpf0(idx)*440;                          % Convert to Hz
        end
        f0 = repmat(f0, size(PeakData, 1), 1);

        % Collect peak-f0 relations
        PeakHarm = PeakData ./ f0;                                  % Peak harmonic numbers
        for segnum = 1:realSegNum                                   % For each frame
            tmpPeakData = PeakData(:, segnum);
            idx = (tmpPeakData~=0) & (round(PeakHarm(:, segnum))>=1);   % Valid peaks
            if (idx == 0)
                continue;
            end
            tmpPeakData = hz2midi(tmpPeakData(idx));                % Convert peak frequencies to midi numbers
            tmpf0 = f0(1, segnum);
            
            % Check if the f0 and its harmonics exist (smaller than the nyguist frequency)
            f0HarmExist = CheckF0HarmExist(tmpf0, maxHarm, fs);        
            TotalNum = TotalNum + f0HarmExist;
            TotalNum_f0(round(hz2midi(tmpf0)), :) = TotalNum_f0(round(hz2midi(tmpf0)), :) + f0HarmExist;
            
            % check if the f0 and its harmonics are in the peak area
            tmpmFharm = hz2midi(tmpf0) + FharmDiff;
            f0HarmInPeakArea = CheckF0HarmInPeakArea(tmpPeakData, tmpmFharm);
            f0HarmExistInPeakArea = f0HarmInPeakArea .* f0HarmExist;
            PANum = PANum + f0HarmExistInPeakArea;
            PANum_f0(round(hz2midi(tmpf0)), :) = PANum_f0(round(hz2midi(tmpf0)), :) + f0HarmExistInPeakArea;            

            if round(hz2midi(tmpf0)) == 39
                a = 1;
            end
        end
       
        if mod(wavNo, 10) == 0
            save 'HarmStat_mono.mat' PANum PANum_f0 TotalNum TotalNum_f0;   % Save results every 10 frames
        end
    end
end
% p(e=1 | h): p_e_h
p_e_h = PANum ./ TotalNum;                                    
TotalNum_f0(TotalNum_f0 == 0) = 10000;                              % 10000 is an arbitary number to make sure the demoninator not be zero
% p(e=1 | f0,h): p_e_F0h
PANum_f0(1:35, :) = 10000;
p_e_F0h = PANum_f0 ./ TotalNum_f0;
p_e_F0h(35, :) = 0.9999;
p_e_F0h(p_harmExist_f0 == 1) = 0.9999;
h = fspecial('gaussian',[3 3],1.5);                                 % Gaussian smoothing
p_tmp = conv2(p_e_F0h, h, 'same');
p_e_F0h(36:94, 2:49) = p_tmp(36:94, 2:49);

save 'HarmStatValue_music44p1k.mat' p_e_h p_e_F0h;           % Save the statistics
%% plot P(D|f0)
figure; h = axes('FontSize',30);
[X,Y] = meshgrid(1:maxHarm, 36:95);
Z = p_e_F0h(36:95, :);
Z(60,1) = 1.1;
mesh(X,Y,Z); colormap(gray); 
xlabel('Harmonic number','FontSize',30); 
ylabel('F0 (MIDI number)','FontSize',30); 
zlabel('Probability','FontSize',30);
axis([0 50 36 95 0 1]);
