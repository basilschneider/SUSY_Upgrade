/*
 * SUSY_Upgrade_Plotter.cpp
 *
 *  Created on: 24 Aug 2016
 *      Author: jkiesele
 */


#include <iostream>
#include "interface/SUSY_Upgrade_Plotter.h"

int main(int argc, char* argv[]){

    if(argc!=2){
        std::cout << "SUSY_Upgrade_Plotter: need exactly one input file" <<std::endl;
        exit (-1);
    }


    std::string inputfile=argv[1];

    SUSY_Upgrade_Plotter analyser;


    analyser.readConfigFile(inputfile);

    analyser.start();

    return 1;
}
