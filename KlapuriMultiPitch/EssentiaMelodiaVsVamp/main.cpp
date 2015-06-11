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
    string audioFilename="/Users/GinSonic/MTG/EssentiaProject/KlapuriMultiPitch/EssentiaMelodiaVsVamp/test.wav";
    
    // parameters
    int sampleRate = 44100;
    int frameSize = 2048;
    int hopSize = 128;
    
    // algorithm setup
    essentia::init();
    AlgorithmFactory& factory = standard::AlgorithmFactory::instance();
    
    Algorithm* audio = factory.create("MonoLoader","filename", audioFilename,"sampleRate", sampleRate);
    Algorithm* multiPitch = factory.create("MultiPitchKlapuri", "sampleRate", sampleRate);
    vector<Real> audioVec;
    vector<vector<Real> >pitchMulti;
    
    audio->output("audio").set(audioVec);
    multiPitch->input("signal").set(audioVec);
    multiPitch->output("pitch").set(pitchMulti);
    
    // Compute
    audio->compute();
    multiPitch->compute();
    
    // write to file
    string csvFilenameMulti="/Users/GinSonic/MTG/EssentiaProject/KlapuriMultiPitch/EssentiaMelodiaVsVamp/test.csv";
    ofstream outFileMulti;
    outFileMulti.open(csvFilenameMulti);
    for (int ii=0; ii<pitchMulti.size(); ii++){
        cout << pitchMulti[ii] << endl;
        for (int jj=0; jj<pitchMulti[ii].size(); jj++){
            
            outFileMulti <<  float(ii)*hopSize/sampleRate << ", " << pitchMulti[ii][jj] << endl;
            cout << float(ii)*hopSize/sampleRate << ", " << pitchMulti[ii][jj] << endl;
        }
    }
    
    // clean up
    delete audio;

    
    essentia::shutdown();
    
    return 0;
}
