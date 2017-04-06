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
#include "classes/DelphesClasses.h"


class SUSY_Upgrade_Plotter: public d_ana::basicAnalyzer{
    public:
        SUSY_Upgrade_Plotter():d_ana::basicAnalyzer(){}
        ~SUSY_Upgrade_Plotter(){}


    private:
        void analyze(size_t id);

        void postProcess();
};





#endif /* SUSY_Upgrade_Plotter */
