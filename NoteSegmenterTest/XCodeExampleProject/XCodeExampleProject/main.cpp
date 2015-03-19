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

int main(int argc, const char * argv[]) {
    
    
    string audioFilename="/Users/GinSonic/MTG/EssentiaProject/Data/monoExamples/voice.wav";
    
    // register the algorithms.
    essentia::init();
    
    // parameters (same as in Vamp, unknown set to default)
    int sampleRate = 44100;
    int frameSize = 2048;
    int hopSize = 128;

    // Algorithm parameter setup
    AlgorithmFactory& factory = standard::AlgorithmFactory::instance();
    Algorithm* audio = factory.create("MonoLoader","filename", audioFilename,"sampleRate", sampleRate);
    Algorithm* melExt = factory.create("PitchMelodia", "sampleRate", sampleRate, "frameSize", frameSize, "hopSize", hopSize);
    Algorithm* noteSeg = factory.create("PitchContourSegmentation");
    
    // Algorithm I/O setup
    vector<Real> audioVec, pitch, pitchConfidence, onset, dur;
    vector<int> qPitch;
    audio->output("audio").set(audioVec);
    melExt->input("signal").set(audioVec);
    melExt->output("pitch").set(pitch);
    melExt->output("pitchConfidence").set(pitchConfidence);
    noteSeg->input("pitch").set(pitch);
    noteSeg->input("signal").set(audioVec);
    noteSeg->output("onset").set(onset);
    noteSeg->output("duration").set(dur);
    noteSeg->output("MIDIpitch").set(qPitch);

    // Compute
    audio->compute();
    melExt->compute();
    noteSeg->compute();
    
    // write f0 to file
    string csvFilenameMono="/Users/GinSonic/MTG/EssentiaProject/Data/monoExamples/voice.csv";
    ofstream outFileMono;
    outFileMono.open(csvFilenameMono);
    for (int ii=0; ii<pitch.size(); ii++){
        outFileMono << float(ii)*hopSize/sampleRate << ", " << pitch[ii] << endl;
    }
    ofstream outFileNotes;
    outFileNotes.open("/Users/GinSonic/MTG/EssentiaProject/Data/monoExamples/voice.notes");
    for (int ii=0; ii<onset.size(); ii++){
         outFileNotes << onset[ii] << ", " << dur[ii] << ", " << qPitch[ii] << endl;
    }
    return 0;
}
