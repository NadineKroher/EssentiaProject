% Read polyphony and the file names of the samples which compose the
% wavFile (mixture)
% Zhiyao Duan
% 11/14/2008

% Input:
%   -wavFile: the wavFile name (including path) 
%   -samplePath: the path where the samples are stored
% Ouput:
%   -Samples: a struct array
%       -midinum: the midi number of the sample
%       -filename: the sample file name (including path)

function Samples = GetSampleFilesFromMixRef(wavFile, samplePath)
textName = strrep(wavFile,'.wav','.txt');
fid = fopen(textName,'r');
Polyphony = fscanf(fid, '%d', 1);

% Initialize the size of the struct array
Samples(Polyphony).notename = '';
Samples(Polyphony).filename = '';

if Polyphony == 1
    while(fscanf(fid, '%c', 1)~=10)% find the "new line"
    end
    Samples(1).midinum = note2midinum(fscanf(fid, '%3c', 1));
    Samples(1).filename = wavFile;
    fclose(fid);
    return;
end

for i=1:Polyphony
    while(fscanf(fid, '%c', 1)~=10)% find the "new line"
    end
    Samples(i).midinum = note2midinum(fscanf(fid, '%3c', 1));
    fscanf(fid, '%f', 1);
    fscanf(fid, '%f', 1);
    fscanf(fid, '%s', 1);
    fscanf(fid, '%s', 1);
    filename = fscanf(fid, '%s', 1);
    
    
    tmpposarray = strfind(filename, '\');
    tmppos = max(tmpposarray);
    
    Samples(i).filename = strcat(samplePath, filename(tmppos+1 : end));
end
fclose(fid);
