/*
 * microStructureHdl.h
 *
 *  Created on: Mai 2, 2016
 *      Author: Jonathan Nierengarten
 */

#ifndef MICROSTRUCTUREHDL_H_
#define MICROSTRUCTUREHDL_H_

#include <string>
#include "dimensionalBuffer.h"
#include "Eigen/Dense"
#include "myQuaternion.h"
#include "utilities.h"
#include "Settings.h"

using namespace std;
using namespace Eigen;
class Grains;
class SubGrain;
class IterativeGrainScheduler;
class mathMethods;
class myQuaternion;
struct myPreferenceOri;
class randomClass;


#define REDCHAN						0
#define GREENCHAN					1
#define BLUECHAN					2
#define ALPHACHAN					3
#define UCHAR_RANGE_MIN				0
#define UCHAR_RANGE_MAX				255
#define IPFZMAPPING_MAXSIZE			10000


struct RGB
{
	unsigned char R;
	unsigned char G;
	unsigned char B;
	unsigned char ALPHA;
    RGB() : R(0), G(0), B(0), ALPHA(255) {} //white as the default
};


class timeLogger {

public:
	timeLogger();
	~timeLogger();

	void logev(const string title, double time);
	unsigned int get_entries( void );
	vector <double> times;
	vector <string> titles;

private:
	unsigned int entries;
};


/*!
 * \class microStructureHdl
 * \brief Class encapsulating a Level Set Box.
 *
 * LSbox class contains a container of type unsigned it:
 * for each voxel the container stores a grain ID <br>
 *
 * The class is used to organize all routines performed on the hosted grains <br>
 *
 *
 *
 */
class microStructureHdl {

public:
	microStructureHdl();
	virtual ~microStructureHdl();
	void ReadAdditionalInputFiles();
	void readPreferenceOrientationFromFile();
	void GeneratePolycrystallineStructureOfGrains();
	void GenerateSubgrainStructureInEachGrain();
	void VoroGEN();
	void readParticleFile();
	void VoroGenPseudoPeriodic();
	void initializeGrains(vector<vector<Eigen::Vector3d>> hulls,
			vector<double> grainVolume);
	void find_neighbors();
	void updateGlobalVoxelContainer();
	void Execute_SubgrainConstruction();

	void DistributeGrainOriAndSEE();
	void DistributeSubgrainOrientations();
	void readTextureFileVPSC();
	myPreferenceOri findNextPreferenceOrientation(myQuaternion ori);
	void DistributeSubgrainSEE();
	void SaveData();
	void SaveDetailedDiagnostics();
	void SaveParenthood();
	void CreateColormap( struct RGB* thecolormap );
	unsigned int CountNumberOfSubgrains();
	void PlotIPF2DSection();
	void ReportProfile();

	void initEnvironment();
	void initNUMABindings();
	void saveTexture();    //Plots density of orientation
	void plotGrains();
	void copyContainer();

	void testprng( unsigned int n, double mu, double sigma );

private:
	DimensionalBuffer<unsigned int>* m_container;
	vector<Grains*> m_grains;
	double m_h;
	IterativeGrainScheduler* m_grainScheduler;
	vector<myQuaternion>* m_OrientationSpace;
	vector<Vector3d>* m_ParticlesPositions;
	vector<myPreferenceOri>* m_PreferenceOrientations;
	randomClass* m_seqRND;
	timeLogger myprofiler;
};

#endif /* MICROSTRUCTUREHDL_H_ */
