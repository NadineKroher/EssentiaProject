% Calculate statistics of the frequencies and amplitudes of the peaks of
% music from both monophonic and polyphonic chords, which were randomly
% mixed using Iowa instrument note samples.
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

Num_FileToSkip = 0;     % Number of files to skip from the beginning in training
                        % This is good for avoiding files that has problems
                        % to deal with.
                        
frameLen = 2048;        % Frame length (Points), 46ms for fs=44100Hz
hop = 1024;             % Frame hop (Points), 23ms for fs=44100Hz
zpf = 4;                % zero padding factor
win = hamming(frameLen, 'periodic');  % window function
RMS_th = 0.075;         % RMS threshold to decide if use this file or not
peakTh = 50;            % peak detection absolute value threshold (dB)
peakTh_rel = 12;        % peak detection relative value threshold (dB)
peakTh_freq = 16000;    % peak detection frequency upperbound (Hz)
movL = 400;             % moving average width (Hz)
localRange = 1;         % the range for considering local maximum
maxHarm = 50;           % maximum harmonic number considered  
NoteMin = 36;           % minimum F0 considered, C2
NoteMax = 95;           % maximum F0 considered, B6

% YIN parameters
P.minf0 = 40;                                                           % Hz - minimum expected F0 (default: 30 Hz)
P.maxf0 = 2200;                                                         % Hz - maximum expected F0 (default: SR/(4*dsratio))
P.thresh = 0.1;                                                         % threshold (default: 0.1)
% P.relfag = 1;                                                           % if ~0, thresh is relative to min of difference function (default: 1)
P.hop = 32;                                                             % s - interval between estimates (default: 32)
% P.range = [1 4096];                                                     % samples - range of samples ([start stop]) to process
%   P.bufsize:  samples - size of computation buffer (default: 10000)
P.sr = 44100;                                                           % Hz - sampling rate (usually taken from file header)
% P.wsize = 4096;                                                       % samples - integration window size (defaut: SR/minf0)
%	P.lpf:		Hz - intial low-pass filtering (default: SR/4)
%	P.shift		0: shift symmetric, 1: shift right, -1: shift left
%	(default: 0)

%% Do on each wav file
% The training path contains randomly mixed chords from note samples in the
% SamplePath.
TrainPath = 'D:\Work\My_data\ChordData\Training\RandomChords\';              % The directory of training chords
SamplePath = 'D:\Work\My_data\ChordData\SampleMaterial\Training\';           % The directory of the samples of the training chords

PeakStat = [];                                                          % To store the training statistics

numPeaks = size(PeakStat, 1);                                           % Number of current detected peaks        
sizePeakStat = size(PeakStat, 1);                                       % Size of the PeakStat table (dynamic resizing)

fileAllName = dir(TrainPath);                                           % The training file name list
for i = 1:length(fileAllName)
    fileAllName(i).name = strcat(TrainPath, fileAllName(i).name);
end
fileAllName = fileAllName(3:end);

fnum = 1;                                                               % Index of the open list head
wavNo = 0;                                                              % The number of training wav files
while (fnum <= length(fileAllName))
    if( fileAllName(fnum).isdir )                                       % If it is a folder, enqueue all the files in this folder to the open list
        tmpAllName = dir(fileAllName(fnum).name);
        for i = 1:length(tmpAllName)
            if  strcmp(tmpAllName(i).name, '.') || strcmp(tmpAllName(i).name, '..') % Prevent to add the folder itself to the list
                continue;
            else
                tmpAllName(i).name = strcat(fileAllName(fnum).name, '\', tmpAllName(i).name);
                fileAllName = [fileAllName; tmpAllName(i)];             % Enqueue open list
            end
        end
        fnum = fnum + 1;
        continue;
    elseif( isempty(strfind(fileAllName(fnum).name, '.wav')) )          % If it is a file, but not *.wav
        fnum = fnum + 1;
        continue;
    else                                                                % If is is a *.wav file
        wavMixture = fileAllName(fnum).name;
        fprintf(1, 'wavMixture name: %s.\n', wavMixture);
        
        wavSamples = GetSampleFilesFromMixRef(wavMixture, SamplePath);  % Find all the samples that compose this wavfile
        samplenum = length(wavSamples);                                 % The number of samples composing this chord
                
        if (min([wavSamples.midinum]) < 36 || max([wavSamples.midinum]) > 95) % Check if the samples are qualified
            continue;
        end

        %---- begin process ----
        wavNo = wavNo + 1;                                              % The counter of the training chords increases by 1
        fprintf(1, 'File %d are going to process.\n', wavNo); 
        fnum = fnum + 1;                                                % Dequeue open list
        if wavNo <= Num_FileToSkip                                      % Skip some wrong files encountered in YIN
            continue;
        end
        
        % Read the the training chord and its samples
        [wavDataMixture, fs, ~] = ReadWavFile_Mono(wavMixture);
        wavDataSamples = zeros(size(wavDataMixture, 1), samplenum);     % Allocate the array to store the wavDataSamples
        for i = 1:samplenum
            tmpwavData = ReadWavFile_Mono(wavSamples(i).filename);
            wavDataSamples(1:length(tmpwavData), i) = tmpwavData;
        end
        
        % Take the middle frames (exclude the transients) to learn
        startp = 8193;
        endp = length(wavDataMixture)-8192;
        endp = startp + floor((endp-startp)/hop)*hop;
        wavDataMixtureSel = wavDataMixture(startp:endp+1);
        wavDataSamplesSel = wavDataSamples(startp:endp+1, :);
        
        % Detect silent frames
        para = {};
        para.frameLen = frameLen;
        para.hop = hop;
        para.win = win;
        para.RMS_th = RMS_th;
        FrameIndex = DetectSilentFrames(wavDataMixtureSel, para);
  
        % STFT
        fftLen = 2^nextpow2(frameLen*zpf);                                  % zero padding
        para = {};
        para.frameLen = frameLen;
        para.hop = hop;
        para.win = win;
        para.fftLen = fftLen;
        SpecData = CalculateSTFT(wavDataMixtureSel, para);
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
        para.movL = round(movL*fftLen/fs);                                  % moving average length (poin
        para.localRange = localRange;
        [PeakData, PeakAmpData, PeakRelAmpData, PeakNum] = CalculatePeakData(SpecData, FrameIndex, para);
        PeakData = PeakData * fs / fftLen;                              % Change the unit from Point to Hz
        if IfDebug == 1
            % plot peaks
            figure; plot((1:fftLen/2)*fs/fftLen, 20*log10(abs(SpecData(1:fftLen/2,296))));
            hold on; plot(PeakData(:,296), PeakAmpData(:,296), 'ro')
        end
        realSegNum = size(PeakData, 2);                                 % Number of frames
        
        % Calculate f0s using YIN
        AllF0s = zeros(samplenum, realSegNum);
        for i = 1:samplenum
            P.minf0 = midi2hz(wavSamples(i).midinum - 1);               % Search range
            P.maxf0 = midi2hz(wavSamples(i).midinum + 1);
            % YIN
            R = yin(wavDataSamplesSel(:, i), P);    %% calling YIN here. You can download it from YIN's webpage
            for segnum = 1:realSegNum
                startp = 1 + (segnum-1) * hop/P.hop;
                endp = startp + frameLen/P.hop - 1;
                tmpf0 = R.f0(startp:endp);                              % The estimated F0s in this frame
                [minAp, idx] = min(R.ap(startp:endp) ./ R.pwr(startp:endp));    % Take the best F0 estimates in this frame
                AllF0s(i, segnum) = 2^tmpf0(idx)*440;                   % Change the unit from MIDI number to Hz
            end
        end
        
        % start to stat
        for i = 1:realSegNum
            % Get F0s in this frame
            f0s = AllF0s(:, i)';                                        % A row vector, each one is a F0 in that frame
            f0s = f0s(f0s ~= 0);                                        % Eliminate 0 F0s in each frame
            f0s_num = length(f0s);                                      % Number of F0s in each frame
            
            % Get peaks and their amplitudes in this frame
            peaks = PeakData(:, i);                                     % A column vector, each one is a peak frequency in that frame
            peaks_amp = PeakAmpData(peaks~=0, i);                       % Peak amplitudes
            peaks = peaks(peaks ~= 0);                                  % Eliminate 0 frequency peaks in each frame
            peaks_num = length(peaks);
            
            % Change both F0s and peaks into matrices
            span_f0s = repmat(f0s, peaks_num, 1);                       % Change to matrix, to find the minimum peaks_harmerr
            span_peaks = repmat(peaks, 1, f0s_num);
            
            % Calculate peak frequency deviations
            span_peaks_harm = round(span_peaks ./ span_f0s);
            span_peaks_harm(span_peaks_harm == 0) = 1;                  %%%%% use all the peaks or not, 1: use; 0: not use peaks whose frequency are smaller than F0s
            span_peaks_harmerr = span_peaks - span_f0s .* span_peaks_harm;
            
            if f0s_num == 1
                peaks_harm = span_peaks_harm;
                peaks_f0 = span_f0s;
            else
                [peaks_harmerr, idx] = min(abs(span_peaks_harmerr), [], 2); % Calculate the minimum peak frequency deviation using different F0s
                peaks_harm = zeros(peaks_num, 1);
                peaks_f0 = zeros(peaks_num, 1);
                for j = 1:peaks_num
                    peaks_harm(j) = span_peaks_harm(j, idx(j));         % Determine the F0 and harmonic number of each peak
                    peaks_f0(j) = f0s(idx(j));
                end
            end
            
            % Change peaks, peaks_f0 and peaks_harmerr to midinum scale
            for j = 1:peaks_num
                peaks_f0(j) = hz2midi(peaks_f0(j));
                peaks(j) = hz2midi(peaks(j));
            end
            peaks_harmerr = peaks - peaks_f0 - 12*log2(peaks_harm);
            
%             if max(abs(peaks_harmerr)) > 1
%                 x = 0;
%             end

            % Record peaks in the current frame
            curPeakStat = [peaks_amp, peaks, peaks_f0, peaks_harm, peaks_harmerr];
            curPeakStat = curPeakStat(peaks_harm >= 1, :);
            num_curPeakStat = size(curPeakStat, 1);
            
            while numPeaks + num_curPeakStat > sizePeakStat             % Dynamic resize
                PeakStat = [PeakStat; zeros(sizePeakStat+1, 5)];
                sizePeakStat = sizePeakStat * 2 + 1;
            end
            PeakStat(numPeaks + 1 : numPeaks + num_curPeakStat, :) = [peaks_amp, peaks, peaks_f0, peaks_harm, peaks_harmerr];
            numPeaks = numPeaks + num_curPeakStat;
        end    
        
        if mod(wavNo, 10) == 0                                          % Save results every 10 frames
            PeakStat = PeakStat(1:numPeaks, :);
            save 'model_learning\PeakStat_music_tmp.mat' PeakStat;
        end        
    end
end
save 'model_learning\PeakStat_music44p1k.mat' PeakStat;                          % Save the final results
delete 'model_learning\PeakStat_music_tmp.mat';
%% Calculate statistics from the raw data
load 'model_learning\PeakStat_music44p1k.mat' PeakStat;

PeakStat_Norm = PeakStat(PeakStat(:,4)<=maxHarm & abs(PeakStat(:,5))<=0.5, :); % Normal peaks
corrcoef(PeakStat_Norm)	% Calculate the correlation coefficients

PeakStat_Spur = PeakStat(abs(PeakStat(:, 5))>0.5, :); % Spurious peaks

% p(s=1)
p_s = size(PeakStat_Spur, 1) / size(PeakStat, 1); % The ratio of the spurious peaks
% p(a,f | s=1)
Spur_af_Miu = mean(PeakStat_Spur(:, 1:2));
Spur_af_Sigma = cov(PeakStat_Spur(:, 1:2));
% draw contour of p(a,f | s=1)
ctrs{1}=0:80;
ctrs{2}=0:136;
[ConN, C]=hist3(PeakStat_Spur(:,1:2),ctrs);
h = fspecial('gaussian',[9 11],2); % Gaussian smoothing
ConNN = conv2(ConN, h);
[X, Y] = meshgrid(10:136,0:60);
Z = ConNN(5:65,16:142)/sum(sum(ConNN)); Z(1,1)=-0.0000001;
figure; axes('FontSize',30); meshc(X,Y,Z); colormap('gray');
axis([10 136 0 60 -0.002 0.002]);
ylabel('Amplitude (dB)','FontSize',30); 
xlabel('Frequency (MIDI number)','FontSize',30); 
zlabel('Probability density','FontSize',30);

[Spur_af_mo,Spur_af_vo,Spur_af_co,Spur_af_jx] = gmmest(PeakStat_Spur(:,1:2),2,[],[],100,-1);

%% Calculate the 3-d p(a,f,h)
% count all (a,f,h) triples of normal peaks
Num_afh = zeros(80, 137, 51);
for i = 1:size(PeakStat_Norm, 1)
    Amp = max(round(PeakStat_Norm(i,1)),1);
    f = round(PeakStat_Norm(i, 2));
    h = PeakStat_Norm(i, 4);
    Num_afh(Amp, f, h) = Num_afh(Amp, f, h) + 1;
end       
Num_afh_Norm = smooth3(Num_afh,'gaussian',[11 11 5], 2);
Num_fh = sum(Num_afh_Norm, 1);
Num_fh(Num_fh==0) = 1;
% p(a | f,h,s=0)
p_a_fh = zeros(80, 137, 51);
for i = 1:80
    p_a_fh(i,:,:) = Num_afh_Norm(i,:,:)./Num_fh;
end
p_a_fh(p_a_fh < 1e-5) = 1e-5;
% p(a,f | s=0)
Num_af = sum(Num_afh_Norm, 3);
p_af = Num_af / size(PeakStat_Norm, 1); % p_af
p_af(p_af < 1e-5) = 1e-5;

% draw marginal distributions -----------------------
% draw p(a,f | s=0)
af = zeros(80, 137);
for i=1:80
    for j=1:137
        af(i,j)=sum(Num_afh_Norm(i,j,:));
    end
end
figure; axes('FontSize',30);
[X,Y] = meshgrid(30:137, 1:80);
Z = af(:,30:137)/size(PeakStat_Norm, 1);
mesh(X,Y,Z); colormap(gray); 
ylabel('Amplitude (dB)','FontSize',30); 
xlabel('MIDI number','FontSize',30); 
zlabel('Probability density','FontSize',30);
axis([30 137 0 80 0 0.0015]);
% draw p(a,h | s=0)
ah = zeros(80, 51);
for i=1:80
    for j=1:51
        ah(i,j)=sum(Num_afh_Norm(i,:,j));
    end
end
figure;  axes('FontSize',30);
[X,Y] = meshgrid(1:51, 1:80);
Z = ah/size(PeakStat_Norm, 1);
mesh(X,Y,Z); colormap(gray); 
ylabel('Amplitude (dB)','FontSize',30); 
xlabel('Harmonic number','FontSize',30); 
zlabel('Probability density','FontSize',30);
axis([0 51 0 80 0 0.004]);
% draw p(f,h)
fh = zeros(137, 51);
for i=1:137 
    for j=1:51 
        fh(i,j)=Num_fh(1,i,j);
    end 
end
figure; axes('FontSize',30);
[X,Y] = meshgrid( 1:51,30:137);
Z = fh(30:137,:)/size(PeakStat_Norm, 1);
mesh(X,Y,Z); colormap(gray); 
ylabel('MIDI number','FontSize',30); 
xlabel('Harmonic number','FontSize',30); 
zlabel('Probability density','FontSize',30);
axis([0 51 30 137 0 0.0015]);
% ------------------------------------------------------

%% calculate p(d|f0), p(d|h), calculate gmm model of the deviations of the normal peaks
% p(d|f0)
idx = round(PeakStat_Norm(:,3))<55 & round(PeakStat_Norm(:,3))>=45; N = size(idx(idx==1),1);
[n,xout]=hist(PeakStat_Norm(idx, 5), -0.5:0.001:0.5);
figure; axes('FontSize', 20); 
bar(xout, n/sum(n)/0.001, 'FaceColor', [0.5,0.5,0.5]);
axis([-0.2 0.2 0 40]);
% p(d|h)
for i = 1:2 %maxHarm
    idx = round(PeakStat_Norm(:,4)) == i; N = size(idx(idx==1),1);
    [n,xout]=hist(PeakStat_Norm(idx, 5), -0.5:0.001:0.5);
    figure; axes('FontSize', 20); 
    bar(xout, n/sum(n)/0.001, 'FaceColor', [0.5,0.5,0.5]);
    axis([-0.2 0.2 0 40]);
end
% gmm to get mo, vo and co
[d_mo,d_vo,d_co,jx] = gmmest(PeakStat_Norm(:,5),4,[],[],100,-1);
axis([-0.2 0.2 0 35]);

%% calculate the harmonic existence probability p(h|f0)
% count all (h, f0) pairs of normal peaks
Num_hf0 = zeros(maxHarm, 100);
for i = 1:size(PeakStat_Norm, 1)
    f0 = round(PeakStat_Norm(i, 3));
    h = PeakStat_Norm(i, 4);
    Num_hf0(h, f0) = Num_hf0(h, f0) + 1;
end
h = fspecial('gaussian',[5 5],1.5); % Gaussian smoothing
Num_hf0_S = conv2(Num_hf0, h, 'same');
% p(h | f0)
p_h_f0 = zeros(maxHarm, 100);
for i = 1:100
    if sum(Num_hf0(:, i)) ~= 0
        p_h_f0(:, i) = Num_hf0(:, i) / sum(Num_hf0(:, i));
    end
end
p_h_f0(p_h_f0 < 1e-5) = 1e-5;
% draw p(hi, f0)
figure; h = axes('FontSize',30);
[X,Y] = meshgrid(36:95, 1:maxHarm);
Z = Num_hf0_S(:, 36:95)/size(PeakStat_Norm, 1);
mesh(X,Y,Z); colormap(gray); 
xlabel('F0 (MIDI number)','FontSize',30); 
ylabel('Harmonic number','FontSize',30); 
zlabel('Probability density','FontSize',30);
axis([36 95 0 51 0 0.002]);
%% save
save 'model_learning\PeakStatValue_music44p1k.mat' p_s Spur_af_Miu Spur_af_Sigma p_a_fh d_mo d_vo d_co;
save 'model_learning\PeakStatValue_music44p1k_full.mat' p_s Spur_af_Miu Spur_af_Sigma Spur_af_mo Spur_af_vo Spur_af_co p_a_fh p_h_f0 p_af d_mo d_vo d_co;
