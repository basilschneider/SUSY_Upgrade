/*
 * SUSY_Upgrade_Plotter.h
 *
 *  Created on: 24 Aug 2016
 *      Author: jkiesele
 */

#ifndef SUSY_Upgrade_Plotter_H_
#define SUSY_Upgrade_Plotter_H_

#include "interface/basicAnalyzer.h"
#include "interface/sampleCollection.h"
#include "interface/stackPlotter.h"
#include "classes/DelphesClasses.h"

typedef std::tuple<int,std::string,double> cutflowtuple;

class SUSY_Upgrade_Plotter: public d_ana::basicAnalyzer{
    public:
        SUSY_Upgrade_Plotter():d_ana::basicAnalyzer(){}
        ~SUSY_Upgrade_Plotter(){}


    private:
        void analyze(size_t id);
        bool isBasicSelection();

        void postProcess();

        std::vector<Electron>* skimelecs;
        std::vector<MissingET>* skimmet;

        // Vector of 3-tuple storing cutflow: Cut number, Cut name, Number of events
        std::vector<cutflowtuple> cutflow;

};





#endif /* SUSY_Upgrade_Plotter */
