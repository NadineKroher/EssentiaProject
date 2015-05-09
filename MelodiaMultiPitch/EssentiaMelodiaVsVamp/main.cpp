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
    string audioFilename="/Users/GinSonic/MTG/EssentiaProject/MelodiaMultiPitch/EssentiaMelodiaVsVamp/test2.wav";
    
    // parameters
    int sampleRate = 44100;
    int frameSize = 2048;
    int hopSize = 128;
    
    // algorithm setup
    essentia::init();
    AlgorithmFactory& factory = standard::AlgorithmFactory::instance();
    
    // MELODIA
    Algorithm* audio = factory.create("MonoLoader","filename", audioFilename,"sampleRate", sampleRate);
    Algorithm* predmelMulti = factory.create("MultiPitchMelodia", "sampleRate", sampleRate, "frameSize", frameSize, "hopSize", hopSize, "minFrequency", 120);
    Algorithm* el = factory.create("EqualLoudness","sampleRate", sampleRate);
    
    // Algorithm I/O setup
    vector<Real> audioVec, audioEQ;
    vector<vector<Real> >pitchMulti;
    audio->output("audio").set(audioVec);
    el->input("signal").set(audioVec);
    el->output("signal").set(audioEQ);
    predmelMulti->input("signal").set(audioEQ);
    predmelMulti->output("pitch").set(pitchMulti);
    
    // Compute
    audio->compute();
    el->compute();
    predmelMulti->compute();
    
    // write to file
    string csvFilenameMulti="/Users/GinSonic/MTG/EssentiaProject/MelodiaMultiPitch/EssentiaMelodiaVsVamp/test2.csv";
    ofstream outFileMulti;
    outFileMulti.open(csvFilenameMulti);
    for (int ii=0; ii<pitchMulti.size(); ii++){
        for (int jj=0; jj<pitchMulti[ii].size(); jj++){
            outFileMulti <<  float(ii)*hopSize/sampleRate << ", " << pitchMulti[ii][jj] << endl;
        }
    }
    
    // clean up
    delete audio;
    delete predmelMulti;

    
    essentia::shutdown();
    
    return 0;
}
