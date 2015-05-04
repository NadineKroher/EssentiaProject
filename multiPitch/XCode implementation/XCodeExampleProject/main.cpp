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
    vector<Real> frame, windowedFrame, spectrum, fSpectrum, peakBins, peakMags, peakBinsP, peakMagsP, relMagsP;
    
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
    
    Algorithm* mov=factory.create("MovingAverage","size", int(400*zpf*frameLen/44100));
    mov->input("signal").set(spectrum);
    mov->output("signal").set(fSpectrum);
    
    Algorithm* spectralPeaks=factory.create("SpectralPeaks","magnitudeThreshold",-50);
    spectralPeaks->input("spectrum").set(spectrum);
    spectralPeaks->output("frequencies").set(peakBins);
    spectralPeaks->output("magnitudes").set(peakMags);
    
    // output storage
    vector<vector<Real> > PeakData; // peak frequency bins for each frames
    vector<vector<Real> > PeakAmpData; // peak magnitudes for each frames
    vector<vector<Real> > RelPeakAmpData; // rel. peak magnitudes for each frame
    vector<int> PeakNum; // number of peaks in each frame
    int numFrames=0;
    
    // frame-wise processing
    while (true){
        
        peakBins.clear();
        peakBinsP.clear();
        peakMags.clear();
        peakMagsP.clear();
        relMagsP.clear();
        
        frameCutter->compute(); // get a frame
        
        if (!frame.size()){ // end of track is reached
            break;
        }
        numFrames++;
        
        if (isSilent(frame)){
            peakBinsP.clear();
            peakMagsP.clear();
            relMagsP.clear();
            PeakData.push_back(peakBinsP);
            PeakAmpData.push_back(peakMagsP);
            RelPeakAmpData.push_back(relMagsP);
            continue;
        }
        
        // spectral peaks
        window->compute();
        spec->compute();
        float ma=spectrum[argmax(spectrum)];
        // get log spectrum
        for (int ii=0; ii<spectrum.size(); ii++){
            spectrum[ii]=20*log10(spectrum[ii]/ma);
        }
        spectralPeaks->compute();
        mov->compute();
        for (int ii=0; ii<peakBins.size(); ii++){
            int bin=round(zpf*frameLen*peakBins[ii]/44100); // get relative magnitude
            float rMag=peakMags[ii]-fSpectrum[bin];
            float mPeak=69+12*log2(peakBins[ii]/440); // conversion to MIDI number
            if (rMag>5 && mPeak>midiMin){ // select only peaks above the local threshold and above min. considered F0
                peakBinsP.push_back(mPeak);
                peakMagsP.push_back(peakMags[ii]);
                relMagsP.push_back(rMag);
            }
        }

        PeakData.push_back(peakBinsP);
        PeakAmpData.push_back(peakMagsP);
        RelPeakAmpData.push_back(relMagsP);
    }
    
    
    
    /////////////////////////////////////////////////////////
    // STAGE 3-4: F0 candidate selection and f0 estimation //
    /////////////////////////////////////////////////////////
    vector<vector<float> > EstAllF0, EstAllLogLike;
    
    // spurious peak frequency and amplitude mean
    vector<float> Spur_af_Miu;
    Spur_af_Miu.push_back(23.0106);
    Spur_af_Miu.push_back(82.0925);
    
    // spurious peak frequency and amplitude covariance
    vector<vector<float> > Spur_af_Sigma;
    vector<float> s0, s1;
    s0.push_back(86.7795);
    s0.push_back(-89.5364);
    s1.push_back(-89.5364);
    s1.push_back(481.6379);
    Spur_af_Sigma.push_back(s0);
    Spur_af_Sigma.push_back(s1);
    float det=s0[0]*s1[1]-s0[1]*s1[0]; // determant
    vector<float> inv0, inv1; // inverse
    vector<vector<float> >inv;
    inv0.push_back(s1[1]/det);
    inv0.push_back(-s0[1]/det);
    inv1.push_back(-s1[0]/det);
    inv1.push_back(s0[0]/det);
    inv.push_back(inv0);
    inv.push_back(inv1);
    
    for (int i=0; i<PeakData.size(); i++){
        
        vector<float> FrAllF0, FrLogLike;
        if (!PeakData[i].size()){ // no peaks found in this frame
            EstAllF0.push_back(FrAllF0);
            EstAllLogLike.push_back(FrLogLike);
            continue;
        }
        
        vector<float> p_SpurP;
        for (int ii=0; ii<PeakData[i].size(); ii++){
            
            // spurious peak likelihood (single Gaussian model)
            float x1=PeakData[i][ii];
            float x2=PeakAmpData[i][ii];
            float m1=Spur_af_Miu[0];
            float m2=Spur_af_Miu[1];
            float inv11=inv0[0];
            float inv12=inv0[1];
            float inv21=inv1[0];
            float inv22=inv1[1];
            float exponent=((x1-m1)*inv11+(x2-m2)*inv21)*(x1-m1)+((x1-m1)*inv12+(x2-m2)*inv22)*(x2-m2);
            exponent*=0.5;
            float nom=1/sqrt((det*2*3.14)*(det*2*3.14));
            float val=nom*exp(exponent);
            p_SpurP.push_back(val);
            
            int peakNum=PeakData[i].size();
            int snum=int(peakNum/3);
        }
        
    }
    return 0;
}
