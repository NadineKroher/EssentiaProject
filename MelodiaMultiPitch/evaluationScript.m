%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% EVALUATE RAW PITCH ACCURACY %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

% file paths
filenames=cell(4,1);
filenames{1}='violin.csv';
filenames{2}='voice.csv';
filenames{3}='flute.csv';
filenames{4}='sax.csv';
root='/Users/GinSonic/MTG/EssentiaProject/Data/Monophonic';

% metrics init
numFrames=0;
rpaVamp=0;
rpaYin=0;
rpaEss=0;
rpaPoly=0;

for m=1:length(filenames);
    % load files
    file=filenames{m};
    gtPath=strcat(root,'/f0GroundTruth/', file);
    gt=load(gtPath);
    time=gt(:,1);
    gt=gt(:,2);
    essentiaPath=strcat(root,'/f0MonoMelodiaEssentia/',file);
    essentia=load(essentiaPath);
    essentia=essentia(:,2);
    vampPath=strcat(root,'/f0MonoMelodiaVamp/',file);
    vamp=load(vampPath);
    vamp=vamp(:,2);
    polyPath=strcat(root,'/f0PolyMelodiaEssentia/',file);
    poly=load(polyPath);
    poly=poly(:,2);
    yinPath=strcat(root,'/f0Yin/',file);
    yin=load(yinPath);
    yin=yin(:,2);
    
    % set to same length
    len=length(gt);
    lDiff=len-size(essentia,1);
    if lDiff>0
        essentia=[essentia; zeros(lDiff,1)];
    end
    lDiff=len-size(vamp,1);
    if lDiff>0
        vamp=[vamp; zeros(lDiff,1)];
    end
    lDiff=len-size(poly,1);
    if lDiff>0
        poly=[poly; zeros(lDiff,1)];
    end
    lDiff=len-size(yin,1);
    if lDiff>0
        yin=[yin; zeros(lDiff,1)];
    end
    
    % evaluate
    for i=1:length(gt)
        numFrames=numFrames+1;
       if abs(1200*log2(yin(i)/gt(i)))<75 || (yin(i)==0 && gt(i)==0)
           rpaYin=rpaYin+1;
       end  
       if abs(1200*log2(vamp(i)/gt(i)))<75 || (vamp(i)==0 && gt(i)==0)
           rpaVamp=rpaVamp+1;
       end 
       if abs(1200*log2(essentia(i)/gt(i)))<75 || (essentia(i)==0 && gt(i)==0)
           rpaEss=rpaEss+1;
       end
       if abs(1200*log2(poly(i)/gt(i)))<75 || (poly(i)==0 && gt(i)==0)
           rpaPoly=rpaPoly+1;
       end 
    end
    
    
end

rpaPoly=rpaPoly/numFrames
rpaVamp=rpaVamp/numFrames
rpaEss=rpaEss/numFrames
rpaYin=rpaYin/numFrames


