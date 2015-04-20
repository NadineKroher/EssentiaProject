//
//  main.cpp
//  XCodeExampleProject
//
//  Created by Nadine Kroher on 14/03/15.
//  Copyright (c) 2015 MTG. All rights reserved.
//

#include <iostream>
#include "essentia.h"
#include "taglib.h"
#include "fftw3.h"
#include <essentia/algorithmfactory.h>
#include <essentia/essentiamath.h>
#include <essentia/pool.h>

using namespace essentia;
using namespace std;
using namespace essentia::standard;

///////////////////////
// GLOBAL PARAMETERS //
///////////////////////

// FFT + spectral Peaks
int frameLen=2048;
int zpf=4;
int hop=441;

// F0-estimation
bool bPeakFreqGMM = true; // 0: modeling frequency deviation using single Gaussian, 1: using GMM
bool bPeakAmp = true; // 0: not consider peak amplitude, 1: consider peak amplitude
bool bSpuriousPeak = true; // 0: consider all the peaks as normal peaks, 1: consider spurious peak
bool bNonpeak = true; // 0: not consider non-peak area, 1: consider non-peak area
bool bNonpeakProb = true; // 0: constant value, 1: learned value (MAYBE I SHOULD USE 0 HERE)
bool bMask = false; // 0: not consider mask are, 1: consider mask area
bool bHarmonicPrior = false; // 0: not use harmonic prior in peak region likelihood calculation, 1: use harmonic prior. Using it will deemphasize the peak region part.
int maxHarm = 50; // upper bound of the harmonic number
int maxF0Num = 6; // the number of estimated pitches in each frame in the first place
int midiMin = 48; // lowest possible frequency of F0 (midi number)
int midiMax = 107; // highest possible frequency of F0 (midi number)
float dupF0Th = 0.5; // frequency threshold (in midi number) to decide if two F0 candidates are duplicate
float f0step = 0.2; // F0 search step (midi number)

// polyphony estimation
int bInstPolyEst=1; // 0: assume each time frame has trackNum of concurrent pitches, 1: estimate instaneous polyphony
float para_poly=0.88; // polyphony estimation threshold, percentage of the whole likelihood increase

// post-processing
bool bRefine = true; // 0: not refine F0 estiamtes, 1: refine F0 estiamtes using neighbouring frames
int binSize=1; // frequency bin size in semitones of the histogram
int neigSize=9; // radius of the neighborhood
bool bSecondRefine=false; // 0: not refine F0 estimates again, 1: refine F0 estimates again, by removing some outliers and filling some gaps
float MSL_pd=0.3; // pitch difference threshold of must-link (midi number)
float mergeNoteGap=100.0; // threshold for note formation (ms), two notelets with gap less than this threshold can be merged
float minNoteLength=100.0; // threshold for note length (ms), notelets shorter than this will be removed



int main(int argc, const char * argv[]) {
    
    // ESSENTIA INIT
    essentia::init();
    AlgorithmFactory& factory = standard::AlgorithmFactory::instance();
    
    // LOAD AUDIO
    string audioFilename="/Users/GinSonic/MTG/EssentiaProject/Data/multiPitch/sine.wav";

    Algorithm* audioLoader = factory.create("MonoLoader","filename", audioFilename,"sampleRate", 44100);
    vector<Real> audioSamples;
    audioLoader->output("audio").set(audioSamples);
    audioLoader->compute();
    delete audioLoader;
    
    ///////////////////////////////////////////////////
    // STAGE 1 - 2: FFT and spectral peak extraction //
    ///////////////////////////////////////////////////
    
    // algorithm setup
    vector<Real> frame, windowedFrame, spectrum, fSpectrum, peakBins, peakMags, relMags;
    
    Algorithm* frameCutter=factory.create("FrameCutter", "frameSize", frameLen, "hopSize", hop, "startFromZero", true);
    frameCutter->input("signal").set(audioSamples);
    frameCutter->output("frame").set(frame);
    frameCutter->reset();
    
    Algorithm* window=factory.create("Windowing", "type", "hann", "zeroPadding", (zpf-1)*frameLen);
    window->input("frame").set(frame);
    window->output("frame").set(windowedFrame);
    
    Algorithm* spec=factory.create("Spectrum", "size", zpf*frameLen);
    spec->input("frame").set(windowedFrame);
    spec->output("spectrum").set(spectrum);
    
    Algorithm* mov=factory.create("MovingAverage", "size", 40);
    mov->input("signal").set(spectrum);
    mov->output("signal").set(fSpectrum);
    
    Algorithm* spectralPeaks=factory.create("PeakDetection", "range", zpf*frameLen*0.5, "minPosition", 1, "maxPosition", zpf*frameLen*0.5, "threshold",0.0); // we don't use spectralPeaks, since we need the bins and not the corresponding frequencies
    spectralPeaks->input("array").set(spectrum);
    spectralPeaks->output("positions").set(peakBins);
    spectralPeaks->output("amplitudes").set(peakMags);
    
    // output storage
    vector<vector<Real> > PeakData; // peak frequency bins for each frames
    vector<vector<Real> > PeakAmpData; // peak magnitudes for each frames
    vector<vector<Real> > RelPeakAmpData; // rel. peak magnitudes for each frame
    vector<int> PeakNum; // number of peaks in each frame
    int numFrames=0;
    
    // frame-wise processing
    while (true){
        
        relMags.clear();
        frameCutter->compute(); // get a frame
        
        if (!frame.size()){ // end of track is reached
            break;
        }
        numFrames++;
        
        if (isSilent(frame)){
            peakBins.clear();
            peakMags.clear();
            PeakData.push_back(peakBins);
            PeakAmpData.push_back(peakMags);
            continue;
        }
        
        // spectral peaks
        window->compute();
        spec->compute();
        mov->compute();
        
        for (int i=0; i<spectrum.size(); i++){
            spectrum[i]=20*log10(spectrum[i]*817); // dB magnitude + scaling for Matlab comparison
            fSpectrum[i]=20*log10(fSpectrum[i]*817);
        }
        
        spectralPeaks->compute();
        
        for (int i=0; i<peakMags.size(); i++){ // convert magnitude
            relMags.push_back(spectrum[round(peakBins[i])]-fSpectrum[round(peakBins[i])]);
            //cout << spectrum[round(peakBins[i])] << "   " << fSpectrum[round(peakBins[i])] << endl;
        }
        
        PeakData.push_back(peakBins);
        PeakAmpData.push_back(peakMags);
        RelPeakAmpData.push_back(relMags);
    }
    
    
    return 0;
}
