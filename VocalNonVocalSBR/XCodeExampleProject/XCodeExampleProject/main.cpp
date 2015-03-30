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

using namespace std;
using namespace essentia;
using namespace essentia::standard;

int main(int argc, const char * argv[]) {
    
    // SETUP //
    essentia::init();
    AlgorithmFactory& factory = AlgorithmFactory::instance();
    
    string root="/Users/GinSonic/MTG/EssentiaProject/Data/vocalSegmentation/";
    
    // parameters
    int sampleRate = 44100;
    int frameSize = 2048;

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
        
        // extract features
        Algorithm* frameCutterL = factory.create("FrameCutter", "frameSize", frameSize, "hopSize", frameSize);
        Algorithm* wL     = factory.create("Windowing","type", "blackmanharris62");
        Algorithm* specL  = factory.create("Spectrum");
        
        // spectral band ratio parameters
        float f1Low=80.0;
        float f1High=400.0;
        float f2Low=500.0;
        float f2High=6000.0;
        int b1Low=round((0.5*frameSize*f1Low)/(0.5*sampleRate));
        int b1High=round((0.5*frameSize*f1High)/(0.5*sampleRate));
        int b2Low=round((0.5*frameSize*f2Low)/(0.5*sampleRate));
        int b2High=round((0.5*frameSize*f2High)/(0.5*sampleRate));
        
        // extract spectral band ratio
        vector<float> windowedFrame, spectrum, frame, sbr;
        
        // connect algorithms
        frameCutterL->input("signal").set(audioStream);
        frameCutterL->output("frame").set(frame);
        wL->input("frame").set(frame);
        wL->output("frame").set(windowedFrame);
        specL->input("frame").set(windowedFrame);
        specL->output("spectrum").set(spectrum);
        
        while(true){
            frameCutterL->compute();
            if (!frame.size()){
                break;
            }
            // sepctral band ratio
            wL->compute();
            specL->compute();
            float highMag, lowMag;
            lowMag=0;
            highMag=0;
            for (int i=b1Low; i<=b1High; i++){
                lowMag+=spectrum[i];
            }
            for (int i=b2Low; i<=b2High; i++){
                highMag+=spectrum[i];
            }
            if (lowMag!=0){
                sbr.push_back(20*log10(highMag/lowMag));
            }else{
                sbr.push_back(0);
            }
        }

        delete wL;
        delete specL,
        delete frameCutterL;
        delete audio;
        
        // write csv file
        string csvFilenameMono=root+to_string(i)+"SBR.csv";
        ofstream outFileMono;
        outFileMono.open(csvFilenameMono);
        for (int ii=0; ii<sbr.size(); ii++){
            outFileMono << float(ii)*frameSize/sampleRate << ", " << sbr[ii] << ", " << voicing[ii] << endl;
        }

    }
    
    return 0;
}
