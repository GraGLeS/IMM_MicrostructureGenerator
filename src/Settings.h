/*
 * Settings.h
 *
 *  Created on: Apr 25, 2016
 *      Author: Jonathan Nierengarten
 */
//*********************************************************************************************************************************************************************************************
#ifndef SETTINGS_H_
#define SETTINGS_H_
using namespace std;
#include <string>

#define MINIMUM_DISCRETIZATION		10

//*********************************************************************************************************************************************************************************************

enum E_TEXTURE {
	E_DEFAULT_TEXTURE,
	E_USE_PREFERENCEORI
};

enum E_MICROGENMODE {
	E_VORONOI,
	E_PREDEFINEDSTRUCTURE,
	E_DEFAULT_GENERATOR
};

enum E_CRYSTAL_STRUCTURE { //Corresponds to the CRYSTAL STRUCTURE of the sample the user wants to work with
	E_HCP, //Hexagonal compact
	E_FCC, //Face-centered cubic
	E_BCC, //Body-centered cubic
	E_DEFAULT_STRUCTURE
};

enum E_PLOT_DIMENSION {
	E_3D,
	E_2D,
	E_DEFAULT_DIMENSION
};

//*********************************************************************************************************************************************************************************************

class Settings {
public:

	static unsigned int NumberOfGrains; //Number of wished grains in the sample
	static unsigned int NumberOfSubgrains;
	static unsigned int NumberOfGridpoints;
	static unsigned long NumberOfPointsPerSubGrain; //Number of interpolating points for each grain
	static unsigned long MaximumNumberOfThreads;
	//choose Texture Gen Mode
	static double SubgrainOriScatter;
	static double StoredElasticEnergyMax; //Energy of a given dislocation in the sample
	static double StoredElasticEnergyMin;
	static double StoredElasticScatterGrain;
	static double StoredElasticScatterSubgrain;
	static bool ExecuteInParallel; //Activate the execution in parallel
	static bool PlotIPF2DSection;
	static double PlotWindowXMin;
	static double PlotWindowXMax;
	static double PlotWindowYMin;
	static double PlotWindowYMax;

	//static bool UseOrientationSpace;
	static E_MICROGENMODE MicroGenMode;
	static E_CRYSTAL_STRUCTURE CrystalStructure;
	static E_TEXTURE TextureGEN;
	static E_PLOT_DIMENSION PlotDimension;

	static string ReadFromFilename;
	static string AdditionalFilename;

	static void readXML(string filename = "");

};

#endif /* SETTINGS_H_ */
