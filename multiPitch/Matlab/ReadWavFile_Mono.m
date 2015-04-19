function [wavData, fs, nBits] = ReadWavFile_Mono(wavFile, TimeLength)
% Read wavfile into a single track, truncate, zero mean, and normalize RMS.
% Average the two channels if the file is stereo.
% Input:
%   - wavFile   : file name
%   - TimeLength: time length in ms
% Output:
%   - wavData   : mono data, column vector
%   - fs        : sampling rate
%   - nBits     : number of bits
% Author: Zhiyao Duan
% Modified: 5/24/2013

[wavData, fs, nBits] = wavread(wavFile);
wavData = mean(wavData, 2);         % convert stereo to mono
L = length(wavData);

if nargin == 2
    tmpwavLength = floor(TimeLength*fs/1000);
    L = min(L, tmpwavLength);       % truncate wavData if it is longer than input
    wavData = wavData(1:L);
end

wavData = wavData - mean(wavData);  % zero mean
wavData = wavData./sqrt(mean(wavData.^2));  % normalize RMS

% figure; plot(1:L, wavData);
