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
    string audioFilename="/Users/GinSonic/MTG/EssentiaProject/Data/audioMono/track1.wav";
    
    // parameters
    int sampleRate = 44100;
    int frameSize = 2048;
    int hopSize = 128;
    Real voicingTolerance=0.2;
    
    // algorithm setup
    essentia::init();
    AlgorithmFactory& factory = standard::AlgorithmFactory::instance();
    Algorithm* audio = factory.create("MonoLoader","filename", audioFilename,"sampleRate", sampleRate);
    Algorithm* predmel = factory.create("MelodiaMonophonic", "sampleRate", sampleRate, "frameSize", frameSize, "hopSize", hopSize, "voicingTolerance", voicingTolerance);
    
    // Algorithm I/O setup
    vector<Real> audioVec, pitch, pitchConfidence;
    audio->output("audio").set(audioVec);
    predmel->input("signal").set(audioVec);
    predmel->output("pitch").set(pitch);
    predmel->output("pitchConfidence").set(pitchConfidence);
    
    // Compute
    audio->compute();
    predmel->compute();

    // write to file
    string csvFilename="/Users/GinSonic/MTG/EssentiaProject/Data/audioMono/track1.csv";
    ofstream outFile;
    outFile.open(csvFilename);
    for (int ii=0; ii<pitch.size(); ii++){
        outFile << float(ii)*hopSize/sampleRate << ", " << pitch[ii] << endl;
    }
    // clean up
    essentia::shutdown();
    
    return 0;
}
