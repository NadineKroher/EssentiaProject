% Calculate statistics of the frequencies and amplitudes of the peaks of
% speech from both monophonic and polyphonic training data.
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

Num_FileToSkip = 0;     % Number of files to skip from the beginning in training
                        % This is good for avoiding files that has problems
                        % to deal with.

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

%% Do on each wav file
% This path contains monophonic speech recordings for training. Each
% recording (a .wav file) should have its ground-truth F0 file stored in
% the .f0 file of the same name as the .wav file. Each line of the .f0 file
% stores a single number of the F0 (in Hz) in the corresponding frame.
% Note: the length and hop sie of the time frames in the F0 file should
% match the parameters listed above.
TrainPath = 'D:\Work\My_data\MultiPitch_Speech\PTDB-TUG-mixed\Train\';

PeakStat = [];
numPeaks = size(PeakStat, 1);       % number of current detected peaks
sizePeakStat = size(PeakStat, 1);   % size of the PeakStat table

% find all subfolders
fileAllName = dir(TrainPath);
for i = 1:length(fileAllName)
    fileAllName(i).name = strcat(TrainPath, fileAllName(i).name);
end
fileAllName = fileAllName(3:end);

fnum = 1; % index of the open list head
wavNo = 0;% the number of training wav files
while (fnum <= length(fileAllName))
    if( fileAllName(fnum).isdir ) % folder
        tmpAllName = dir(fileAllName(fnum).name);
        for i = 1:length(tmpAllName)
            if  strcmp(tmpAllName(i).name, '.') || strcmp(tmpAllName(i).name, '..') % prevent to add the same fold to the list
                continue;
            else
                tmpAllName(i).name = strcat(fileAllName(fnum).name, '\', tmpAllName(i).name);
                fileAllName = [fileAllName; tmpAllName(i)]; % enqueue open list
            end
        end
        fnum = fnum + 1;
        continue;
    elseif( isempty(strfind(fileAllName(fnum).name, '.wav')) ) % is a file, but not *.wav
        fnum = fnum + 1;
        continue;
    else % *.wav file
        wavMixture = fileAllName(fnum).name;
        fprintf(1, 'wavMixture name: %s.\n', wavMixture);

        %---- begin process ----
        wavNo = wavNo + 1;
        fprintf(1, 'File %d is going to be processed.\n', wavNo); 
        fnum = fnum + 1;       % dequeue open list
        if wavNo < Num_FileToSkip    % to skip some wrong files encountered in YIN
            continue;
        end
        
        % Read the wavfile, and resample it to 16kHz
        [wavDataMixture, fs, ~] = ReadWavFile_Mono(wavMixture);
        wavDataMixture = resample(wavDataMixture, 16000, fs);
        fs = 16000;
        
        % Detect silent frames
        para = {};
        para.frameLen = frameLen;
        para.hop = hop;
        para.win = win;
        para.RMS_th = RMS_th;
        FrameIndex = DetectSilentFrames(wavDataMixture, para);
        
        % STFT
        fftLen = 2^nextpow2(frameLen*zpf);                                  % zero padding
        para = {};
        para.frameLen = frameLen;
        para.hop = hop;
        para.win = win;
        para.fftLen = fftLen;
        SpecData = CalculateSTFT(wavDataMixture, para);
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
            figure; plot((1:fftLen/2)*fs/fftLen, 20*log10(abs(SpecData(1:fftLen/2,300))));
            hold on; plot(PeakData(:,300), PeakAmpData(:,300), 'ro')
        end
        
        % Read ground-truth pitches
        GTF0s = dlmread(strrep(wavMixture, '.wav', '.f0'));                 % in Hz
        GTF0s = GTF0s';

        % Collect peak-f0 relations
        FrameNum = min(size(GTF0s,2), length(FrameIndex));
        idx = (sum(GTF0s(:,1:FrameNum), 1)~=0 & FrameIndex(1:FrameNum)~=0);
        UseFrames = 1:FrameNum;
        UseFrames = UseFrames(idx);
        for i = UseFrames
            fprintf(1, 'Frame %d is going to be processed...\n', i);
            f0s = GTF0s(:, i)';                                             % a row vector, each one is a F0 in that frame
            f0s = f0s(f0s ~= 0);
            f0s_num = length(f0s);            
            
            peaks = PeakData(:, i);                                         % a column vector, each one is a peak frequency in that frame
            peaks_amp = PeakAmpData(peaks~=0, i);
            peaks = peaks(peaks ~= 0);
            peaks_num = length(peaks);            
            
            span_f0s = repmat(f0s, peaks_num, 1);                           % change to matrix, to find the minimum peaks_harmerr
            span_peaks = repmat(peaks, 1, f0s_num);
            
            span_peaks_harm = round(span_peaks ./ span_f0s);
            span_peaks_harm(span_peaks_harm == 0) = 1;                      % every peak needs to find a harmonic
            span_peaks_harmerr = span_peaks - span_f0s .* span_peaks_harm;
            
            if f0s_num == 1
                peaks_harm = span_peaks_harm;
                peaks_f0 = span_f0s;
            else                                                            % if more than one F0, use the best one to calculate the frequency deviation of a peak
                [peaks_harmerr, idx] = min(abs(span_peaks_harmerr), [], 2);
                peaks_harm = zeros(peaks_num, 1);
                peaks_f0 = zeros(peaks_num, 1);
                for j = 1:peaks_num
                    peaks_harm(j) = span_peaks_harm(j, idx(j));
                    peaks_f0(j) = f0s(idx(j));
                end
            end
            
            % change peaks, peaks_f0 and peaks_harmerr to midinum scale
            peaks_f0 = hz2midi(peaks_f0);
            peaks = hz2midi(peaks);
            peaks_harmerr = peaks - peaks_f0 - 12*log2(peaks_harm);

            curPeakStat = [peaks_amp, peaks, peaks_f0, peaks_harm, peaks_harmerr];
            curPeakStat = curPeakStat(peaks_harm >= 1, :);
            num_curPeakStat = size(curPeakStat, 1);
            
            while numPeaks + num_curPeakStat > sizePeakStat % resize
                PeakStat = [PeakStat; zeros(sizePeakStat+1, 5)];
                sizePeakStat = sizePeakStat * 2 + 1;
            end
            PeakStat(numPeaks + 1 : numPeaks + num_curPeakStat, :) = [peaks_amp, peaks, peaks_f0, peaks_harm, peaks_harmerr];
            numPeaks = numPeaks + num_curPeakStat;
        end    
        
        if mod(wavNo, 10) == 0 % save
            save 'model_learning\PeakStat_speech_tmp.mat' PeakStat numPeaks;
        end        
    end
end
PeakStat = PeakStat(1:numPeaks, :);
save 'model_learning\PeakStat_speech16k.mat' PeakStat;     % raw data
delete 'model_learning\PeakStat_speech_tmp.mat';

%% Calculate statistics from the raw data
load 'model_learning\PeakStat_speech16k.mat';

% normal peaks, s=0
PeakStat_Norm = PeakStat(PeakStat(:,4)<=maxHarm & abs(PeakStat(:,5))<=0.5, :);
corrcoef(PeakStat_Norm)
% spurious peaks, s=1
PeakStat_Spur = PeakStat(abs(PeakStat(:, 5))>0.5, :);

% p(s=1)
p_s = size(PeakStat_Spur, 1) / size(PeakStat, 1); % the ratio of the noisy peaks
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
figure; axes('FontSize',20); meshc(X,Y,Z); colormap('gray');
axis([10 136 0 60 -0.002 0.003]);
ylabel('Log amplitude (dB)','FontSize',20); 
xlabel('Log frequency (MIDI number)','FontSize',20); 
zlabel('Probability density','FontSize',20);

[Spur_af_mo,Spur_af_vo,Spur_af_co,Spur_af_jx] = gmmest(PeakStat_Spur(:,1:2),2,[],[],100,-1);

%% Calculate the 3-d p(a,f,h)
% count all (a,f,h) triples of normal peaks
Num_afh = zeros(80, 137, 51);
for i = 1:size(PeakStat_Norm, 1)
    amp = max(round(PeakStat_Norm(i,1)),1);
    f = round(PeakStat_Norm(i, 2));
    h = PeakStat_Norm(i, 4);
    Num_afh(amp, f, h) = Num_afh(amp, f, h) + 1;
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
axis([0 51 0 80 0 0.008]);
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
save 'model_learning\PeakStatValue_speech16k.mat' p_s Spur_af_Miu Spur_af_Sigma p_a_fh d_mo d_vo d_co;
save 'model_learning\PeakStatValue_speech16k_full.mat' p_s Spur_af_Miu Spur_af_Sigma Spur_af_mo Spur_af_vo Spur_af_co p_a_fh p_h_f0 p_af d_mo d_vo d_co;
