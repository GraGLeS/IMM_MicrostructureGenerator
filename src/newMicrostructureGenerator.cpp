/*
 * newMicrostructurGenerator.cpp
 *
 *  Created on: May 2, 2016
 *      Author: jn538375
 */

#include <iostream>
#include "microStructureHdl.h"
#include "Settings.h"
using namespace std;

int main(int argc, char *argv[]) {
	try {
		if (argc > 1)
			Settings::readXML(argv[1]);
		else
			Settings::readXML();
	} catch (exception& e) {
		cout << "Unable to parse parameters file! Details:\n";
		cout << e.what() << endl;
		cout << "Simulation will now halt." << endl;
		return 0;
	}

	microStructureHdl myHdl;

	//myHdl.testprng( 100000, 1.0e+14, 5.0e+13 );
	myHdl.initEnvironment();
	myHdl.ReadAdditionalInputFiles();
	myHdl.GeneratePolycrystallineStructureOfGrains();
	myHdl.DistributeGrainOriAndSEE();

	myHdl.GenerateSubgrainStructureInEachGrain();
	myHdl.DistributeSubgrainOrientations();
	myHdl.DistributeSubgrainSEE();

	myHdl.SaveData();
	myHdl.SaveDetailedDiagnostics();
	myHdl.SaveParenthood();
	myHdl.PlotIPF2DSection();
	myHdl.ReportProfile();
	cout << "Program finished successfully!";
	return 0;
}

