//////////////////////////////////////////////////////////////////////
//  Extracts the pitch using the monophonic melodia implementation. //
//////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include "essentia.h"
#include "taglib.h"
#include "fftw3.h"
#include "avcodec.h"
#include <essentia/algorithmfactory.h>
#include <essentia/essentiamath.h>
#include <essentia/pool.h>

using namespace std;
using namespace essentia;
using namespace essentia::standard;

int main(int argc, const char * argv[]) {
    
    // file path
    string audioFilename="/Users/GinSonic/MTG/EssentiaProject/Data/monoExamples/bassLoop.wav";
    
    // parameters
    int sampleRate = 44100;
    int frameSize = 2048;
    int hopSize = 128;
    Real voicingTolerance=0.2;
    
    // algorithm setup
    essentia::init();
    AlgorithmFactory& factory = standard::AlgorithmFactory::instance();
    
    // MELODIA
    Algorithm* audio = factory.create("MonoLoader","filename", audioFilename,"sampleRate", sampleRate);
    Algorithm* predmelMono = factory.create("PitchMelodia", "sampleRate", sampleRate, "frameSize", frameSize, "hopSize", hopSize);
    Algorithm* el = factory.create("EqualLoudness","sampleRate", sampleRate);
    //Algorithm* predmel = factory.create("PredominantMelody", "sampleRate", sampleRate, "frameSize", frameSize, "hopSize", hopSize, "voicingTolerance", voicingTolerance);
    
    // Algorithm I/O setup
    vector<Real> audioVec, audioEQ, pitch, pitchConfidence, pitchMono, pitchConfidenceMono;
    audio->output("audio").set(audioVec);
    el->input("signal").set(audioVec);
    el->output("signal").set(audioEQ);
    //predmel->input("signal").set(audioVec);
    //predmel->output("pitch").set(pitch);
    //predmel->output("pitchConfidence").set(pitchConfidence);
    predmelMono->input("signal").set(audioEQ);
    predmelMono->output("pitch").set(pitchMono);
    predmelMono->output("pitchConfidence").set(pitchConfidenceMono);
    
    // Compute
    audio->compute();
    el->compute();
    //predmel->compute();
    predmelMono->compute();

    /*
    // YIN FFT
    Algorithm* frameCutter = factory.create("FrameCutter",
                                            "frameSize", frameSize,
                                            "hopSize", hopSize,
                                            "startFromZero", false);
    
    Algorithm* window = factory.create("Windowing",
                                       "type", "hann",
                                       "zeroPadding", 0);
    
    Algorithm* spectrum = factory.create("Spectrum",
                                         "size", frameSize);
    
    Algorithm* pitchDetect = factory.create("PitchYinFFT",
                                            "frameSize", frameSize,
                                            "sampleRate", sampleRate);
    
    // set frameCutter:
    vector<Real> frame;
    frameCutter->input("signal").set(audioVec);
    frameCutter->output("frame").set(frame);
    
    // set windowing:
    vector<Real> windowedframe;
    window->input("frame").set(frame);
    window->output("frame").set(windowedframe);
    
    // set spectrum:
    vector<Real> spec;
    spectrum->input("frame").set(windowedframe);
    spectrum->output("spectrum").set(spec);
    
    // set pitch extraction:
    Real thisPitch = 0., thisConf = 0;
    pitchDetect->input("spectrum").set(spec);
    pitchDetect->output("pitch").set(thisPitch);
    pitchDetect->output("pitchConfidence").set(thisConf);
    
    vector<Real> pitchYin;
    // process:
    while (true) {
        frameCutter->compute();
        
        if (!frame.size())
            break;
        
        
        if (isSilent(frame)){
            pitchYin.push_back(0);
            continue;
        }
        
        window->compute();
        spectrum->compute();
        pitchDetect->compute();
        pitchYin.push_back(thisPitch);
    }

    */
    // write to file
    /*
    string csvFilename="/Users/GinSonic/MTG/EssentiaProject/Data/Monophonic/f0PolyMelodiaEssentia/sax.csv";
    ofstream outFile;
    outFile.open(csvFilename);
    for (int ii=0; ii<pitch.size(); ii++){
        outFile << float(ii)*hopSize/sampleRate << ", " << pitch[ii] << endl;
    }
    */
    string csvFilenameMono="/Users/GinSonic/MTG/EssentiaProject/Data/monoExamples/bassLoop.csv";
    ofstream outFileMono;
    outFileMono.open(csvFilenameMono);
    for (int ii=0; ii<pitchMono.size(); ii++){
        outFileMono << float(ii)*hopSize/sampleRate << ", " << pitchMono[ii] << endl;
    }
    
    /*
    string csvFilenameYin="/Users/GinSonic/MTG/EssentiaProject/Data/Monophonic/f0Yin/sax.csv";
    ofstream outFileYin;
    outFileYin.open(csvFilenameYin);
    for (int ii=0; ii<pitchYin.size(); ii++){
        outFileYin << float(ii)*hopSize/sampleRate << ", " << pitchYin[ii] << endl;
    }
    */
    
    // clean up
    delete audio;
    //delete predmel;
    delete predmelMono;
    //delete frameCutter;
    //delete pitchDetect;
    //delete window;
    //delete spectrum;
    
    essentia::shutdown();
    
    return 0;
}
