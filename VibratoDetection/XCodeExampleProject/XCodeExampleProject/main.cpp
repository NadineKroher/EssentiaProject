//
//  main.cpp
//
//  Script for testing vocal vibrato detection as implemented in Melodia.
//


#include <iostream>

#include "essentia.h"
#include "taglib.h"
#include "fftw3.h"
#include <essentia/algorithmfactory.h>
#include <essentia/essentiamath.h>
#include <essentia/pool.h>

using namespace std;
using namespace essentia;
using namespace essentia::standard;

int main(int argc, const char * argv[]) {
    
    // load f0
    string filename="/Users/GinSonic/MTG/EssentiaProject/Data/vibrato/test.csv";
    ifstream csvfile(filename);
    string line;
    vector<float> time, pitch;
    while (getline(csvfile,line))
    {
        std::istringstream ss(line);
        std::string token;
        getline(ss, token, ',');
        time.push_back(stof(token.c_str()));
        getline(ss, token, ',');
        pitch.push_back(stof(token.c_str()));
    }
    
    // get contour start and end indices
    vector<float> startC, endC;
    if (pitch[0]>0){
        startC.push_back(0);
    }
    for (int i=0; i<pitch.size()-1; i++){
        if (pitch[i+1]>0 && pitch[i]==0){
            startC.push_back(i+1);
        }
        if (pitch[i+1]==0 && pitch[i]>0){
            endC.push_back(i);
        }
    }
    if (endC.size()<startC.size()){
        endC.push_back(pitch.size()-1);
    }
    
    // Parameters
    Real sampleRate=44100;
    int hopSize=128;
    Real vibSampleRate=sampleRate / hopSize;
    int vibFrameSize=int(0.350 * vibSampleRate);
    int vibHopSize=1;
    int vibZeroPaddingFactor=4;
    int vibFFTSize=vibFrameSize * vibZeroPaddingFactor;
    vibFFTSize=pow(2, ceil(log(vibFFTSize)/log(2)));
    Real vibMinFreq=5.0;
    Real vibMaxFreq=8.0;
    Real vibdBDropLobe=15;
    Real vibdBDropSecondPeak=20;
    
    essentia::init();
    AlgorithmFactory& factory = standard::AlgorithmFactory::instance();
    Algorithm* v = factory.create("Vibrato");
    
    vector<float> freq, ext;
    v->input("sampleRate").set(Real(44100/128));
    v->input("pitch").set(pitch);
    v->output("vibratoFrequency").set(freq);
    v->output("vibratoExtend").set(ext);
    
    v->compute();
    /*
    // Algorithm setup
    AlgorithmFactory& factory = standard::AlgorithmFactory::instance();
    Algorithm* frameCutter=factory.create("FrameCutter", "frameSize", vibFrameSize, "hopSize", vibHopSize, "startFromZero", true);
    Algorithm* spectrum=factory.create("Spectrum", "size", vibFFTSize);
    Algorithm* window=factory.create("Windowing", "type", "hann");
    Algorithm* spectralPeaks=factory.create("SpectralPeaks", "sampleRate", vibSampleRate, "maxPeaks", 3, "orderBy", "magnitude");
    
    vector<Real> vibFreq;
    vibFreq.resize(pitch.size());
    for (int i=0; i<pitch.size(); i++){
        vibFreq[i]=0.0;
    }
    
    // iterate over contour segments
    for (int i=0; i<startC.size(); i++){
        
        // get a segment in cents
        vector<Real> contour;
        for (int ii=startC[i]; ii<=endC[i]; ii++){
            contour.push_back(1200*log2(pitch[ii]/55.0));
        }
        
        // setup algorithm I/O
        vector<Real> frame;
        frameCutter->input("signal").set(contour);
        frameCutter->output("frame").set(frame);
        vector<Real> windowedFrame;
        window->input("frame").set(frame);
        window->output("frame").set(windowedFrame);
        vector<Real> vibSpectrum;
        spectrum->input("frame").set(windowedFrame);
        spectrum->output("spectrum").set(vibSpectrum);
        vector<Real> peakFrequencies, peakMagnitudes;
        spectralPeaks->input("spectrum").set(vibSpectrum);
        spectralPeaks->output("frequencies").set(peakFrequencies);
        spectralPeaks->output("magnitudes").set(peakMagnitudes);
        frameCutter->reset();
        
        int frameNo=0;
        int minFrameNo=0;
        
        // frame-wise processing
        while (true){
            bool vibrato=true;
            
            //get a frame
            frameCutter->compute();
            frameNo++;
            
            if(!frame.size()){
                break;
            }
            
            // subtract mean pitch from frame
            Real m=mean(frame, 0, frame.size()-1);
            for (int ii=0; ii<frame.size(); ii++){
                frame[ii]-=m;
            }
            
            // spectral peaks
            window->compute();
            spectrum->compute();
            spectralPeaks->compute();
            
            int numberPeaks = peakFrequencies.size();
            if (!numberPeaks) {
                vibrato=false;
            }
            
            if (peakFrequencies[0] < vibMinFreq || peakFrequencies[0] > vibMaxFreq) {
                vibrato=false;
            }
            
            if (numberPeaks > 1) {  // there is at least one extra peak
                if (peakFrequencies[1] <= vibMaxFreq) {
                    vibrato=false;
                }
                if (20 * log10(peakMagnitudes[0]/peakMagnitudes[1]) < vibdBDropLobe) {
                    vibrato=false;
                }
            }
            
            if (numberPeaks > 2) {  // there is a second extra peak
                if (peakFrequencies[2] <= vibMaxFreq) {
                    vibrato=false;
                }
                if (20 * log10(peakMagnitudes[0]/peakMagnitudes[2]) < vibdBDropSecondPeak) {
                    vibrato=false;
                }
            }
            
            if ((frame[argmax(frame)]-frame[argmin(frame)])>250){
                vibrato=false;
            }
            
            if(vibrato){
                for (int ii=startC[i]+frameNo; ii<startC[i]+frameNo+vibFrameSize; ii++){
                vibFreq[ii]=1;
                }
            }
        }
    }
    
    string csvFilenameMono="/Users/GinSonic/MTG/EssentiaProject/Data/vibrato/testV.csv";
    ofstream outFileMono;
    outFileMono.open(csvFilenameMono);
    for (int ii=0; ii<vibFreq.size(); ii++){
        outFileMono << float(ii)*hopSize/sampleRate << ", " << vibFreq[ii] << endl;
    }
    
    // clean up
    delete frameCutter;
    delete spectralPeaks;
    delete spectrum;
    delete window;
    
    */
    
    string csvFilenameMono="/Users/GinSonic/MTG/EssentiaProject/Data/vibrato/testV.csv";
    ofstream outFileMono;
    outFileMono.open(csvFilenameMono);
    for (int ii=0; ii<freq.size(); ii++){
        outFileMono << float(ii)*hopSize/sampleRate << ", " << freq[ii] << ", " << ext[ii] << endl;
    }
    
    essentia::shutdown();
    
    return 0;
}
