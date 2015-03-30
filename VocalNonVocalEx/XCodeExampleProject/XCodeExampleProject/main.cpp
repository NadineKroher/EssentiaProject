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
#include <string.h>

using namespace std;
using namespace essentia;
using namespace essentia::standard;
vector <vector<float> > features;

int main(int argc, const char * argv[]) {
    
    string root="/Users/GinSonic/MTG/EssentiaProject/Data/vocalSegmentation/";
    
    // register the algorithms.
    essentia::init();
    
    // parameters
    int sampleRate = 44100;
    int frameSize = 1024;
    
    // Algorithm parameter setup
    AlgorithmFactory& factory = standard::AlgorithmFactory::instance();
    
    for (int i=0; i<10; i++){
        
        string filenameA=root+to_string(i)+".wav";
        string filenameC=root+to_string(i)+"GT.csv";
        
        // load audio
        Algorithm* audio = factory.create("MonoLoader","filename", filenameA,"sampleRate", sampleRate);
        vector<float> audioStream;
        audio->output("audio").set(audioStream);
        audio->compute();
        
        // load csv
        ifstream in(filenameC.c_str());
        string line;
        vector<float> timeInd;
        while (getline(in,line))
        {
            string val=line.substr(0, line.size()-10);
            timeInd.push_back(stof(val));
        }
        
        // get voicing ground truth
        vector<float> startC, endC;
        for (int ii=0; ii<timeInd.size()-1; ii+=2){
            startC.push_back(timeInd[ii]);
            endC.push_back(timeInd[ii+1]);
        }
        
        int numFrames=floor(float(audioStream.size())/float(frameSize));
        float voicing[numFrames];
        memset(voicing, 0, numFrames*sizeof(float));
        for (int ii=0; ii<startC.size(); ii++){
            int startS=round(startC[ii]*float(sampleRate)/float(frameSize));
            int endS=round(endC[ii]*float(sampleRate)/float(frameSize));
            for (int jj=startS; jj<=endS; jj++){
                voicing[jj]=1;
            }
        }
        
        string csvFilenameMono=root+to_string(i)+"GTF.csv";
        ofstream outFileMono;
        outFileMono.open(csvFilenameMono);
        for (int ii=0; ii<numFrames; ii++){
            outFileMono << float(ii)*frameSize/sampleRate << ", " << voicing[ii] << endl;
        }
        
        // extract features
        vector<Real> frame, windowedFrame;
        vector<Real> spectrum, mfccCoeffs, mfccBands;
        
        Algorithm* fc    = factory.create("FrameCutter",
                                          "frameSize", frameSize,
                                          "hopSize", frameSize);
        
        Algorithm* w     = factory.create("Windowing",
                                          "type", "blackmanharris62");
        
        Algorithm* spec  = factory.create("Spectrum");
        Algorithm* mfcc  = factory.create("MFCC");
        
        fc->input("signal").set(audioStream);
        fc->output("frame").set(frame);
        w->input("frame").set(frame);
        w->output("frame").set(windowedFrame);
        spec->input("frame").set(windowedFrame);
        spec->output("spectrum").set(spectrum);
        mfcc->input("spectrum").set(spectrum);
        mfcc->output("bands").set(mfccBands);
        mfcc->output("mfcc").set(mfccCoeffs);
        
        int frameNo=-1;
        vector<float> mv;
        
        for (int ii=0; ii<13; ii++){
            mv.push_back(0.0);
        }
        vector<vector<float> >featuresLoc;
        while (true) {
            
            vector<float> feature;
            // compute a frame
            fc->compute();
            frameNo++;
            // if it was the last one (ie: it was empty), then we're done.
            if (!frame.size()) {
                break;
            }
            
            // if the frame is silent, just drop it and go on processing
            if (isSilent(frame)) continue;
            
            w->compute();
            spec->compute();
            mfcc->compute();
            for (int ii=0; ii<mfccCoeffs.size(); ii++){
                feature.push_back(mfccCoeffs[ii]);
                mv[ii]+=mfccCoeffs[i];
            }
            feature.push_back(voicing[frameNo]);
            features.push_back(feature);
        }
        
        // clean up
        delete fc;
        delete w;
        delete spec;
        delete mfcc;
        delete audio;
    }
    
   
    
    // write csv file
    string csvFilenameMono="/Users/GinSonic/MTG/EssentiaProject/Data/vocalSegmentation/features.csv";
    ofstream outFileMono;
    outFileMono.open(csvFilenameMono);
    outFileMono << "MFCC1, MFCC2, MFCC3, MFCC4, MFCC5, MFCC6, MFCC7, MFCC8, MFCC9, MFCC10, MFCC11, MFCC12, MFCC13, class; " << endl;
    for (int ii=0; ii<features.size(); ii++){
        for (int jj=0; jj<features[ii].size()-1; jj++){
            outFileMono << features[ii][jj] << ", ";
        }
        if (features[ii][features[ii].size()-1]==0){
            outFileMono << "vocals" << "; " << endl;
        }else{
            outFileMono << "instr" << "; " << endl;
        }
    }
    
    essentia::shutdown();
    
    return 0;
}
