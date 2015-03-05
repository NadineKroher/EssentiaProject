///////////////////////////////////////////////////////////////////////////////
//  Extracts the predominant melody using Essentia for comparison with Vamp. //
///////////////////////////////////////////////////////////////////////////////

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
    
    string root="/Users/GinSonic/MTG/EssentiaProject/Data/audio/Track0";
    string rootOut="/Users/GinSonic/MTG/EssentiaProject/Data/PM_Essentia/Track0";
    
    // register the algorithms.
    essentia::init();
    
    for (int i=1; i<6; i++){
        
        // path
        string audioFilename=root+to_string(i)+".wav";
        cout << audioFilename << endl;
        
        // parameters (same as in Vamp, unknown set to default)
        int sampleRate = 44100;
        int frameSize = 2048;
        int hopSize = 128;
        Real voicingTolerance=0.2;
        
        // Algorithm parameter setup
        AlgorithmFactory& factory = standard::AlgorithmFactory::instance();
        Algorithm* audio = factory.create("MonoLoader","filename", audioFilename,"sampleRate", sampleRate);
        Algorithm* predmel = factory.create("PredominantMelody", "sampleRate", sampleRate, "frameSize", frameSize, "hopSize", hopSize, "voicingTolerance", voicingTolerance);
        
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
        string csvFilename=rootOut+to_string(i)+".csv";
        ofstream outFile;
        outFile.open(csvFilename);
        for (int ii=0; ii<pitch.size(); ii++){
            outFile << float(ii)*hopSize/sampleRate << ", " << pitch[ii] << endl;
        }
        
        // clean up
        delete audio;
        delete predmel;
    
    }
    
    return 0;
}
