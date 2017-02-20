/*
 * microStructureHdl.cpp
 *
 *  Created on: Mai 2, 2016
 *      Author: Jonathan Nierengarten
 */

#include "microStructureHdl.h"
#include "Grains.h"
#include "voro++/include/voro++/voro++.hh"
#include "utilities.h"
#include <stdexcept>
#include "IterativeGrainScheduler.h"
#include <fstream>
#include <iostream>
#include <string>
#include <omp.h>
#include <numa.h>
#include "mymath.h"
#include "RTree.h"
#include "SubGrain.h"
#include "random.h"
#include <iomanip>
#include "lodepng.h"

using namespace voro;
using namespace std;
using namespace Eigen;

//Definition of the needed functions **********************************************************************************************************************************************************

timeLogger::timeLogger() {
	entries = 0;
}

timeLogger::~timeLogger() {
	titles.clear();
	times.clear();
}

void timeLogger::logev(const string title, double time) {
	titles.push_back(title);
	times.push_back(time);
	entries++;
}

unsigned int timeLogger::get_entries(void) {
	return entries;
}

microStructureHdl::microStructureHdl() {

	m_container = NULL;
	m_grainScheduler = NULL;
	m_OrientationSpace = NULL;
	m_ParticlesPositions = NULL;
	m_seqRND = NULL;

	int subgrains = Settings::NumberOfSubgrains;
	if (subgrains <= 0)
		subgrains = 1;

	if (Settings::PlotDimension != E_2D) {

		if (Settings::MicroGenMode == E_PREDEFINEDSTRUCTURE) {
			Settings::NumberOfGridpoints = pow(
					(subgrains * Settings::NumberOfGrains * 27), (1 / 3.0))
					* pow(4 / 3 * PI, 1. / 3) / 2
					* Settings::NumberOfPointsPerSubGrain;
		} else {
			Settings::NumberOfGridpoints = pow(
					(subgrains * Settings::NumberOfGrains), (1 / 3.0))
					* pow(4 / 3 * PI, 1. / 3) / 2
					* Settings::NumberOfPointsPerSubGrain;
			cout << "Number of Gridpoints: " << Settings::NumberOfGridpoints
					<< endl;
		}
		m_container = new DimensionalBuffer<unsigned int>(0, 0, 0,
				Settings::NumberOfGridpoints, Settings::NumberOfGridpoints,
				Settings::NumberOfGridpoints);

	} else {
		if (Settings::MicroGenMode == E_PREDEFINEDSTRUCTURE) {
			Settings::NumberOfGridpoints = pow(
					(subgrains * Settings::NumberOfGrains * 9), (1 / 2.0))
					* sqrt(PI) / 2 * Settings::NumberOfPointsPerSubGrain;
		} else {
			Settings::NumberOfGridpoints = pow(
					(subgrains * Settings::NumberOfGrains), (1 / 2.0))
					* sqrt(PI) / 2 * Settings::NumberOfPointsPerSubGrain;
			cout << "Number of Gridpoints: " << Settings::NumberOfGridpoints
					<< endl;
		}
		m_container = new DimensionalBuffer<unsigned int>(0, 0, 0,
				Settings::NumberOfGridpoints, Settings::NumberOfGridpoints, 1);

	}
	m_h = 1. / Settings::NumberOfGridpoints;
	m_seqRND = new randomClass(-3000);
	//vector<randomClass*>* m_threadlocalRND = new vector<randomClass*> [omp_get_max_threads() - 1];
}

microStructureHdl::~microStructureHdl() {
	delete m_container;
	delete m_grainScheduler;
	delete m_OrientationSpace;
	delete m_ParticlesPositions;
	delete m_seqRND;
}

void microStructureHdl::GeneratePolycrystallineStructureOfGrains() {
	double gtime = omp_get_wtime();
	switch (Settings::MicroGenMode) {
	case E_VORONOI: {
		VoroGEN();
		break;
	}
	case E_PREDEFINEDSTRUCTURE: {
		VoroGenPseudoPeriodic();
		break;
	}
	case E_VPSC: {
		VoroGEN();
		break;
	}
	default: {
		break;
	}

	}
	copyContainer();

	myprofiler.logev("GenerateGrainHull", (omp_get_wtime() - gtime));

//plotGrains();
}

void microStructureHdl::GenerateSubgrainStructureInEachGrain() {
	double gtime = omp_get_wtime();

	Execute_SubgrainConstruction();

	myprofiler.logev("GenerateSubgrainHull", (omp_get_wtime() - gtime));
}

void microStructureHdl::readTextureFileVPSC() {
	FILE * OriFromFile;
	OriFromFile = fopen(Settings::ReadFromFilename.c_str(), "r");
	int id, N = 0;
	char c;
	// count number of orientations
	do {
		c = fgetc(OriFromFile);
		if (c == '\n')
			N++;
	} while (c != EOF);
	rewind(OriFromFile);
	N -= 4;
	// read over header
	int i = 1;
	do {
		c = fgetc(OriFromFile);
		if (c == '\n')
			i++;
	} while (i <= 4);
	double rho, phi1, PHI, phi2;
	m_OrientationSpace = new vector<myQuaternion>;
	m_ParticlesPositions = new vector<Vector3d>;
	for (int i = 0; i < N; i++) {
		fscanf(OriFromFile, "%lf \t %lf \t %lf \t %lf \n", &phi1, &PHI, &phi2,
				&rho);
		m_OrientationSpace->push_back(myQuaternion(phi1, PHI, phi2));
	}
	fclose(OriFromFile);

	cout
			<< "ParticleFile VPSC read successfully with m_OrientationSpace.size = "
			<< m_OrientationSpace->size() << endl;
}

void microStructureHdl::readParticleFile() {
	FILE * OriFromFile;
	OriFromFile = fopen(Settings::ReadFromFilename.c_str(), "r");
	int id, N = 0;
	char c;
// count number of orientations
	do {
		c = fgetc(OriFromFile);
		if (c == '\n')
			N++;
	} while (c != EOF);
	rewind(OriFromFile);

// read over header
	double vol, x, y, z, boxSize;
	double q[4];
	m_OrientationSpace = new vector<myQuaternion>;
	m_ParticlesPositions = new vector<Vector3d>;
	for (int i = 0; i < N - 1; i++) {
		if (i == 0) {
			fscanf(OriFromFile, "%lf \n", &boxSize);
			do {
				c = fgetc(OriFromFile);
			} while (c != '\n');
			continue;
		}
		if (i == 1) {
			do {
				c = fgetc(OriFromFile);
			} while (c != '\n');
			continue;
		}
		//else {
		fscanf(OriFromFile,
				"%d \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \n",
				&id, &x, &y, &z, &vol, &q[0], &q[1], &q[2], &q[3]);
		m_OrientationSpace->push_back(myQuaternion(q[0], q[1], q[2], q[3]));
		m_ParticlesPositions->push_back(
				Vector3d(x / boxSize, y / boxSize, z / boxSize));
		//}
	}
	fclose(OriFromFile);

	cout << "ParticleFile read successfully with m_OrientationSpace.size = "
			<< m_OrientationSpace->size() << endl;
//for ( unsigned int j = 0; j < m_OrientationSpace->size(); j++ ) cout << (*m_OrientationSpace)[j].get_q0() << ";" <<  (*m_OrientationSpace)[j].get_q1() << endl;
}

void microStructureHdl::VoroGenPseudoPeriodic() {
	cout << "Started Voro Gen Pseudo Periodic" << endl;
	m_grains.resize((Settings::NumberOfGrains * 27) + 1);
	bool randbedingung = false; // bei false ist der container halb offen?! d.h. gitterwert mit 1 werden keinem partikel zugeordnet
	voronoicell_neighbor c;
	int blocks = (int) (pow((Settings::NumberOfGrains * 27 / 8), (1 / 3.)) + 1);
	if (blocks < 1)
		blocks = 1;
	voro::container con(0, 1, 0, 1, 0, 1, blocks, blocks, blocks, randbedingung,
			randbedingung, randbedingung, 8);

	c_loop_all vl(con);

	/**********************************************************/
	int k_limit = 3;
	if (Settings::PlotDimension == E_2D)
		k_limit = 1;
	int region = 0;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < k_limit; k++) {
				for (int id = 0; id < Settings::NumberOfGrains; id++) {
					double x = (*m_ParticlesPositions)[id][0] / 3.
							+ (i * 1. / 3);
					double y = (*m_ParticlesPositions)[id][1] / 3.
							+ (j * 1. / 3);
					double z = (*m_ParticlesPositions)[id][2] / 3.
							+ (k * 1. / 3);
					if (Settings::PlotDimension == E_2D)
						z = 0.;
					int mid = id + (region * Settings::NumberOfGrains);
					con.put(mid, x, y, z);
				}
				region++;
			}
		}
	}

	/**********************************************************/
	Settings::NumberOfGrains *= 27;
	/**********************************************************/
	vector<vector<Eigen::Vector3d> > initialHulls;
	vector<double> grainVolume;
	vector<double> cellCoordinates;
	if (vl.start()) {
		initialHulls.resize(Settings::NumberOfGrains + 1);
		grainVolume.resize(Settings::NumberOfGrains + 1);
		do {
			double cur_x, cur_y, cur_z;
			con.compute_cell(c, vl);
			//new: get the grain_id
			unsigned int box_id = vl.pid() + 1;
			vl.pos(cur_x, cur_y, cur_z);
			grainVolume[box_id] = c.volume();
			c.vertices(cur_x, cur_y, cur_z, cellCoordinates); //Erstellung des Polyheders
			for (unsigned int i = 0; i < cellCoordinates.size() / 3; i++) {
				initialHulls.at(box_id).push_back(
						Eigen::Vector3d(cellCoordinates.at(3 * i + 1),
								cellCoordinates.at(3 * i),
								cellCoordinates.at(3 * i + 2)));
			}
		} while (vl.inc());

		con.draw_particles("VoronoyP.gnu");
		con.draw_cells_gnuplot("VoronoyC.gnu");
		double x, y, z, rx, ry, rz;
		int cell_id;
		int z_limit = Settings::NumberOfGridpoints;
		if (Settings::PlotDimension == E_2D)
			z_limit = 1;
		for (int k = 0; k < z_limit; k++) {
			for (int j = 0; j < Settings::NumberOfGridpoints; j++) {
				for (int i = 0; i < Settings::NumberOfGridpoints; i++) {
					y = j * m_h;
					x = i * m_h;
					z = k * m_h;
					if (con.find_voronoi_cell(x, y, z, rx, ry, rz, cell_id)) {
						unsigned int box_id = cell_id + 1;
						m_container->setValueAt(j, i, k, box_id);
					} else {
						m_container->setValueAt(j, i, k, 0);
					}
				}
			}
		}
	} else {
		throw runtime_error("Voronoy container error at start() method!");
	}

//generate grainboxes by evaluating the initial hulls objects
	initializeGrains(initialHulls, grainVolume);

}

void microStructureHdl::VoroGEN() {
	cout << "Started Voro Gen" << endl;
	m_grains.resize(Settings::NumberOfGrains + 1);
	bool randbedingung = false; // bei false ist der container halb offen?! d.h. gitterwert mit 1 werden keinem partikel zugeordnet
//	if (randbedingung == false)
//		Settings::NumberOfGridpoints = 1;

	voronoicell_neighbor c; //Voronoi Zelle die ihre Nachbaren kennt
//Erstellung eines Containers mit gewissen Dimensionen und Randbedingungen
	int blocks = (int) (pow((Settings::NumberOfGrains / 8), (1 / 3.)) + 1);
	if (blocks < 1)
		blocks = 1;
	voro::container con(0, 1, 0, 1, 0, 1, blocks, blocks, blocks, randbedingung,
			randbedingung, randbedingung, 8);

	c_loop_all vl(con); //Schleife (Iteration)

	/**********************************************************/
// Randomly add particles into the container
	for (int i = 0; i < Settings::NumberOfGrains; i++) {
		double x = m_seqRND->MersenneTwister(); //->parkMiller();
		double y = m_seqRND->MersenneTwister(); //parkMiller();
		double z = m_seqRND->MersenneTwister(); //parkMiller();
		if (Settings::PlotDimension == E_2D)
			z = 0.0;
		con.put(i, x, y, z); //Jedem Korn i werden die Koordinaten x,y,z zugeordnet
	}

	/**********************************************************/

	vector<vector<Eigen::Vector3d> > initialHulls;
	vector<double> grainVolume;
	vector<double> cellCoordinates;
	if (vl.start()) {
		initialHulls.resize(Settings::NumberOfGrains + 1);
		grainVolume.resize(Settings::NumberOfGrains + 1);
		do {
			double cur_x, cur_y, cur_z;
			con.compute_cell(c, vl);
			//new: get the grain_id
			unsigned int box_id = vl.pid() + 1;
			vl.pos(cur_x, cur_y, cur_z);
			grainVolume[box_id] = c.volume();
			c.vertices(cur_x, cur_y, cur_z, cellCoordinates); //Erstellung des Polyheders
			for (unsigned int i = 0; i < cellCoordinates.size() / 3; i++) {
				initialHulls.at(box_id).push_back(
						Eigen::Vector3d(cellCoordinates.at(3 * i + 1),
								cellCoordinates.at(3 * i),
								cellCoordinates.at(3 * i + 2)));
			}
		} while (vl.inc());

		double x, y, z, rx, ry, rz;
		int cell_id;
		int z_limit = Settings::NumberOfGridpoints;
		if (Settings::PlotDimension == E_2D)
			z_limit = 1;
		for (int k = 0; k < z_limit; k++) {
			for (int j = 0; j < Settings::NumberOfGridpoints; j++) {
				for (int i = 0; i < Settings::NumberOfGridpoints; i++) {
					y = j * m_h;
					x = i * m_h;
					z = k * m_h;
					if (con.find_voronoi_cell(x, y, z, rx, ry, rz, cell_id)) {
						unsigned int box_id = cell_id + 1;
						m_container->setValueAt(j, i, k, box_id);
					} else {
						m_container->setValueAt(j, i, k, 0);
					}
				}
			}
		}
	} else {
		throw runtime_error("Voronoy container error at start() method!");
	}

//generate grainboxes by evaluating the initial hulls objects
	initializeGrains(initialHulls, grainVolume);
}

void microStructureHdl::Execute_SubgrainConstruction() {
	cout << "Start subgrain construction " << endl;
#pragma omp parallel
	{
		randomClass* threadlocalRNG = new randomClass();
		//threadlocalRNG->initPM(-1 * (omp_get_thread_num()) - 1);
		uint32_t overflowint = (uint32_t) pow(2.0, 31);
		threadlocalRNG->initMT(overflowint - omp_get_thread_num() - 1);
		vector<unsigned int>& workload = m_grainScheduler->getThreadWorkload(
				omp_get_thread_num());
//		cout << "Thread: " << omp_get_thread_num() << " has received workload" << endl;
		for (auto id : workload) {
//			cout << id << endl;
			if (id <= Settings::NumberOfGrains)
				m_grains[id]->SubGrainConstructor(*threadlocalRNG);
		}
		delete threadlocalRNG;
	}
}

void microStructureHdl::copyContainer() {
	for (auto it : m_grains)
		if (it != NULL)
			it->copyContainerToGrain(m_container);
}

void microStructureHdl::find_neighbors() {
	RTree<unsigned int, int, 3, float> tree;
	int min[3], max[3];
	for (unsigned int i = 1; i <= Settings::NumberOfGrains; i++) {
		if (m_grains[i] == NULL)
			continue;
		min[0] = m_grains[i]->getMinX();
		min[1] = m_grains[i]->getMinY();
		min[2] = m_grains[i]->getMinZ();
		max[0] = m_grains[i]->getMaxX();
		max[1] = m_grains[i]->getMaxY();
		max[2] = m_grains[i]->getMaxZ();
		tree.Insert(min, max, i);
	}
	for (unsigned int id = 1; id <= Settings::NumberOfGrains; id++) {
		m_grains[id]->computeDirectNeighbours(tree);
	}
}

void microStructureHdl::DistributeGrainOriAndSEE() {
	double gtime = omp_get_wtime();

	//###MK::readPreferenceOrientationFromFile();
	cout << "Sample grain orientations" << endl;
	vector<Grains*>::iterator it;
	for (it = ++m_grains.begin(); it != m_grains.end(); it++) {
		myQuaternion ori;
		if (Settings::MicroGenMode == E_VORONOI
				|| Settings::MicroGenMode == E_VPSC) {
			unsigned int mOrientations = m_OrientationSpace->size();
// 			unsigned int randomOri = m_OrientationSpace->size();
// 			while (randomOri >= mOrientations) {
// 				randomOri = m_seqRND->MersenneTwister()
// 						* m_OrientationSpace->size(); //parkMiller()
// 			}
		    ori = (*m_OrientationSpace)[(*it)->get_ID()];
// 			ori = (*m_OrientationSpace)[randomOri]; //pick randomly from list of predefined orientations

//			cout << "Grain " << (*it)->get_ID()
//					<< " becomes mapped on UserDefinedOrientation " << randomOri
//					<< endl;

			//ori.randomOriShoemakeQuat(*m_seqRND);
		} else if (Settings::MicroGenMode == E_PREDEFINEDSTRUCTURE) { //from predefined orientations
			ori = (*m_OrientationSpace)[((*it)->get_ID() - 1) % 27];
		} else {
			ori.randomOriShoemakeQuat(*m_seqRND); //random orientations
		}

		(*it)->set_Orientation(ori);
		(*it)->set_PreforiProperties(findNextPreferenceOrientation(ori));
		double gr_see_mu = (*it)->get_SEEFromPrefori();
		double gr_see_sig = (*it)->get_SEEGrainScatterFromPrefOri();
//cout << "Grain = " << (*it)->get_ID() << " prefori mu/sig/orisig = " << gr_see_mu << ";" << gr_see_sig << ";" << (*it)->get_OriScatterFromPrefOri() << "----prefori q0 " <<  (*it)->get_PrefOriQuatQ0() << ";" << (*it)->get_PrefOriQuatQ1() << ";" << (*it)->get_PrefOriQuatQ2() << ";" << (*it)->get_PrefOriQuatQ3() << endl;
		(*it)->set_SEE(m_seqRND->r4_nor(gr_see_mu, gr_see_sig)); //set stored elastic energy of grain to a specific value from a normal distribution
	}
	cout << "Sample sub-grain orientation from reference orientation of grain"
			<< endl;

	myprofiler.logev("Distribute grain preference oris",
			(omp_get_wtime() - gtime));
}

void microStructureHdl::DistributeSubgrainOrientations() {
	double gtime = omp_get_wtime();

#pragma omp parallel
	{
		randomClass* threadlocalRNG = new randomClass();
		//threadlocalRNG->initPM(-1 * (omp_get_thread_num()) - 1);
		uint32_t overflowint = (uint32_t) pow(2.0, 31);
		threadlocalRNG->initMT(overflowint - omp_get_thread_num() - 1);

		vector<unsigned int>& workload = m_grainScheduler->getThreadWorkload(
				omp_get_thread_num());
		for (auto id : workload) {
			if (id <= Settings::NumberOfGrains) {
				m_grains[id]->generateSubGrainOri(*threadlocalRNG);
			}
		}
		delete threadlocalRNG;
	}

	myprofiler.logev("Distribute sub-grain orientations",
			(omp_get_wtime() - gtime));
}

void microStructureHdl::ReadAdditionalInputFiles() {
	double gtime = omp_get_wtime();

	switch (Settings::MicroGenMode) {
	case E_VPSC: {
		readTextureFileVPSC();
		break;
	}
	case E_PREDEFINEDSTRUCTURE: {
		readParticleFile();
		break;
	}
	default: {
		readParticleFile();
		break;
	}
	}
	readPreferenceOrientationFromFile();

	myprofiler.logev("ReadAdditionalInput", (omp_get_wtime() - gtime));
}

void microStructureHdl::readPreferenceOrientationFromFile() {
	FILE* file;
	file = fopen(Settings::AdditionalFilename.c_str(), "r");
	int N = 0;
	char c;
// count number of orientations
	do {
		c = fgetc(file);
		if (c == '\n')
			N++;
	} while (c != EOF);
	cout << " Found " << N << " Preference Oris in File" << endl;
	rewind(file);
	m_PreferenceOrientations = new vector<myPreferenceOri>;
	for (int i = 0; i < N; i++) {
		double phi1, PHI, phi2, SEE, oriScatter, SEEGrainScatter,
				SEESubgrainScatter, RelSubgrainSizeScaler;
		fscanf(file, "%lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf ",
				&phi1, &PHI, &phi2, &oriScatter, &SEE, &SEEGrainScatter,
				&SEESubgrainScatter, &RelSubgrainSizeScaler); //SEE reads as mu and SEEScatter as sigma
		char c;
		//read comments
//		do {
//			c = fgetc(file);
//		} while (c != '\n');
		m_PreferenceOrientations->push_back(
				myPreferenceOri(phi1 * PI / 180.0, PHI * PI / 180.0,
						phi2 * PI / 180.0, oriScatter * PI / 180.0, SEE,
						SEEGrainScatter, SEESubgrainScatter,
						RelSubgrainSizeScaler)); //define a preference orientation with specific instead of default properties in terms of mu, sigma and oriscatter
//		cout << "PrefRefs = " << phi1 << ";" << PHI << ";" << phi2 << ";"
//				<< oriScatter << ";" << SEE << ";" << SEEGrainScatter << ";"
//				<< SEESubgrainScatter << ";" << RelSubgrainSizeScaler << endl;
	}
	fclose(file);
	cout
			<< "PrefRefDefault = "
			<< Settings::SubgrainOriScatter
			<< ";"
			<< (Settings::StoredElasticEnergyMax
					+ Settings::StoredElasticEnergyMin) * 0.5 << ";"
			<< Settings::StoredElasticScatterGrain << ";" << 1.0 << endl;
}

myPreferenceOri microStructureHdl::findNextPreferenceOrientation(
		myQuaternion ori) {
	if (Settings::TextureGEN == E_USE_PREFERENCEORI) {
		double min = 15.0 * PI / 180.0;
		myPreferenceOri* i;
		unsigned int idx = 0;
		for (unsigned int p = 0; p < m_PreferenceOrientations->size(); p++) {
			double misori = ori.misorientationCubicQxQ(
					&((*m_PreferenceOrientations)[p].ori));
			if (misori <= min) {
				min = misori;
				idx = p;
			}
		}
		if (min <= 10.0 * PI / 180.0) {
			//return *i;
//cout << "Finding preference orientation with min/idx = " << min << ";" << idx << endl;
			return myPreferenceOri(ori,
					(*m_PreferenceOrientations)[idx].subgrainsScatterOri,
					(*m_PreferenceOrientations)[idx].SEE,
					(*m_PreferenceOrientations)[idx].SEEGrainScatter,
					(*m_PreferenceOrientations)[idx].SEESubgrainScatter,
					(*m_PreferenceOrientations)[idx].RelSubgrainSizeScaling);

		} else {
//cout << "Check for preference orientation but too far off = " << min << endl;
			return myPreferenceOri(ori); //in case no close enough preference orientation was found, thus the ori itself is returned as a PreferenceOrientation
		}
	} else {
//cout << "Simply returning preference orientation " << endl;
		return myPreferenceOri(ori);
	}
}

void microStructureHdl::DistributeSubgrainSEE() {
	double gtime = omp_get_wtime();
	cout << "Sample Stored Elastic Energy" << endl;
#pragma omp parallel
	{
		randomClass* threadlocalRNG = new randomClass();
		//threadlocalRNG->initPM(-1 * (omp_get_thread_num()) - 1);
		uint32_t overflowint = (uint32_t) pow(2.0, 31);
		threadlocalRNG->initSHR3(overflowint - omp_get_thread_num() - 1);
		threadlocalRNG->initR4Uni(overflowint - omp_get_thread_num() - 1);
		threadlocalRNG->r4_nor_setup();

		vector<unsigned int>& workload = m_grainScheduler->getThreadWorkload(
				omp_get_thread_num());
		for (auto id : workload) {
			if (id <= Settings::NumberOfGrains) {
				m_grains[id]->generateSubGrainSEE(*threadlocalRNG);
			}
		}
		delete threadlocalRNG;
	}

	myprofiler.logev("Distribute stored elastic energy",
			(omp_get_wtime() - gtime));
}

void microStructureHdl::SaveData() {
	double asciiftime = omp_get_wtime();
	cout << "Pipe data into file" << endl;
	string filename = string("Microstructure") + string(".uds");
	ofstream file;
	file.open(filename.c_str());
	int offset = 0;
	vector<Grains*>::iterator itG;
	double NrGrains = 0;
	if (Settings::NumberOfSubgrains != 0)
		for (itG = ++m_grains.begin(); itG != m_grains.end(); itG++) {
			if ((*itG)->m_SubGrains.size() > 1)
				NrGrains += ((*itG)->m_SubGrains.size() - 1);
			else
				NrGrains++;
		}
	else
		NrGrains = Settings::NumberOfGrains;
	file << "|| Settings || *(sf) || Parameter, Value" << endl;
	file << "FirstID" << "\t" << 1 << endl;
	file << "DX" << "\t" << Settings::NumberOfPointsPerSubGrain << endl;
	file << "DY" << "\t" << Settings::NumberOfPointsPerSubGrain << endl;
	if (Settings::PlotDimension != E_2D)
		file << "DZ" << "\t" << Settings::NumberOfPointsPerSubGrain << endl;
	file << "NX" << "\t" << Settings::NumberOfGridpoints << endl;
	file << "NY" << "\t" << Settings::NumberOfGridpoints << endl;
	if (Settings::PlotDimension != E_2D)
		file << "NZ" << "\t" << Settings::NumberOfGridpoints << endl;
	file << "NGrains" << "\t" << NrGrains << endl;

	if (Settings::PlotDimension != E_2D) {
		file
				<< "|| Points || *(iffffffiiiiiiff) || ID, x, y, z, bunge1, bunge2, bunge3, xmin, xmax, ymin, ymax, zmin, zmax, vol, stored"
				<< endl;
	} else {
		file
				<< "|| Points || *(ifffffiiiiff) || ID, x, y, bunge1, bunge2, bunge3, xmin, xmax, ymin, ymax, vol, stored"
				<< endl;
	}

	int newoffset = 0;
	for (itG = ++m_grains.begin(); itG != m_grains.end(); itG++) {
		if ((Settings::NumberOfSubgrains != 0)
				&& (*itG)->m_SubGrains.size() > 1) {

			newoffset = (*itG)->copySubgrainsToGlobalContainer(m_container,
					offset); // implicit rehashing of all ID's
			offset = --newoffset;
			vector<SubGrain*>::iterator it;
			for (it = ++((*itG)->m_SubGrains.begin());
					it != (*itG)->m_SubGrains.end(); it++) {
				// ID, x, y, z, bunge1, bunge2, bunge3, xmin, xmax, ymin, ymax, zmin, zmax, volume, stored
				myQuaternion* ori = (*it)->get_Orientation();
				double* bunge = ori->Quaternion2Euler();
				file
						<< (*it)->get_ID()
						<< "\t"
						<< (((*it)->getMaxX() - (*it)->getMinX()) / 2
								+ (*it)->getMinX())
						<< "\t"
						<< (((*it)->getMaxY() - (*it)->getMinY()) / 2
								+ (*it)->getMinY()) << "\t";
				if (Settings::PlotDimension != E_2D) {
					file
							<< (((*it)->getMaxZ() - (*it)->getMinZ()) / 2
									+ (*it)->getMinZ()) << "\t";
				}
				file << bunge[0] << "\t" << bunge[1] << "\t" << bunge[2] << "\t"
						<< (*it)->getMinX() << "\t" << (*it)->getMaxX() << "\t"
						<< (*it)->getMinY() << "\t" << (*it)->getMaxY() << "\t";
				if (Settings::PlotDimension != E_2D) {
					file << (*it)->getMinZ() << "\t" << (*it)->getMaxZ()
							<< "\t";
				}
				file << (*it)->get_Volume() << "\t" << (*it)->get_SEE() << endl;
				delete[] bunge;
			}
		} else {
			// ID, x, y, z, bunge1, bunge2, bunge3, xmin, xmax, ymin, ymax, zmin, zmax, volume, stored
			myQuaternion ori = (*itG)->getOri();
			double* bunge = ori.Quaternion2Euler();
			newoffset = (*itG)->copySubgrainsToGlobalContainer(m_container,
					offset);
			offset = newoffset;
			file
					<< (*itG)->get_ID()
					<< "\t"
					<< (((*itG)->getMaxX() - (*itG)->getMinX()) / 2
							+ (*itG)->getMinX())
					<< "\t"
					<< (((*itG)->getMaxY() - (*itG)->getMinY()) / 2
							+ (*itG)->getMinY()) << "\t";
			if (Settings::PlotDimension != E_2D) {
				file
						<< (((*itG)->getMaxZ() - (*itG)->getMinZ()) / 2
								+ (*itG)->getMinZ()) << "\t";
			}
			file << bunge[0] << "\t" << bunge[1] << "\t" << bunge[2] << "\t"
					<< (*itG)->getMinX() << "\t" << (*itG)->getMaxX() << "\t"
					<< (*itG)->getMinY() << "\t" << (*itG)->getMaxY() << "\t";
			if (Settings::PlotDimension != E_2D) {
				file << (*itG)->getMinZ() << "\t" << (*itG)->getMaxZ() << "\t";
			}
			file << int((*itG)->get_Volume() / m_h / m_h + 0.5) << "\t"
					<< (*itG)->getSEE() << endl;
			delete[] bunge;

		}
	}
//write grain infos

	file.close();

	myprofiler.logev("Writing ASCII UDS", (omp_get_wtime() - asciiftime));
	double binftime = omp_get_wtime();

	cout << "created and piped: " << offset
			<< " subgrains to binary file. Move on ... " << endl;
	stringstream filename2;
	filename2.str("");
	filename2 << "Container" << ".raw";
	FILE* binaryFile;

	int size;
	if (Settings::PlotDimension != E_2D)
		size = pow(Settings::NumberOfGridpoints, 3);
	else
		size = pow(Settings::NumberOfGridpoints, 2);
	binaryFile = fopen(filename2.str().c_str(), "w");
	fwrite(m_container->getRawData(), sizeof(unsigned int), size, binaryFile);
	fclose(binaryFile);

	myprofiler.logev("Writing BINARY Container", (omp_get_wtime() - binftime));
}

void microStructureHdl::SaveDetailedDiagnostics() {
	double asciidiagn = omp_get_wtime();
	cout << "Pipe diagnostics into file" << endl;
	string filename = string("MicrostructureDiagnostics") + string(".uds");
	ofstream file;
	file.open(filename.c_str());
	int offset = 0;
	vector<Grains*>::iterator itG;
	unsigned int NrGrains = 0;
	if (Settings::NumberOfSubgrains != 0)
		for (itG = ++m_grains.begin(); itG != m_grains.end(); itG++) {
			if ((*itG)->m_SubGrains.size() > 1)
				NrGrains += ((*itG)->m_SubGrains.size() - 1);
			else
				NrGrains++;
		}
	else
		NrGrains = Settings::NumberOfGrains;
	file << "|| Settings || *(sf) || Parameter, Value" << endl;
	file << "FirstID" << "\t" << 1 << endl;
	file << "DX" << "\t" << Settings::NumberOfPointsPerSubGrain << endl;
	file << "DY" << "\t" << Settings::NumberOfPointsPerSubGrain << endl;
	if (Settings::PlotDimension != E_2D)
		file << "DZ" << "\t" << Settings::NumberOfPointsPerSubGrain << endl;
	file << "NX" << "\t" << Settings::NumberOfGridpoints << endl;
	file << "NY" << "\t" << Settings::NumberOfGridpoints << endl;
	if (Settings::PlotDimension != E_2D)
		file << "NZ" << "\t" << Settings::NumberOfGridpoints << endl;
	file << "NGrains" << "\t" << NrGrains << endl;

	if (Settings::PlotDimension != E_2D) {
		file
				<< "|| Subgrains || *(iffffffiiiiiiffifffffff) || ID, x, y, z, bunge1, bunge2, bunge3, xmin, xmax, ymin, ymax, zmin, zmax, vol, stored, parent ID, disori to parent meas, disori to parent target, parent bunge1, parent bunge2, parent bunge3, parent stored, parent stored grain scatter, parent stored subgrain scatter, parent size scaler"
				<< endl;
	} else {
		file
				<< "|| Subgrains || *(ifffffiiiiffifffffff) || ID, x, y, bunge1, bunge2, bunge3, xmin, xmax, ymin, ymax, vol, stored, parent ID, disori to parent meas, disori to parent target, parent bunge1, parent bunge2, parent bunge3, parent stored, parent stored grain scatter, parent subgrain scatter, parent size scaler"
				<< endl;
	}

	for (itG = ++m_grains.begin(); itG != m_grains.end(); itG++) {
		if ((Settings::NumberOfSubgrains != 0)
				&& (*itG)->m_SubGrains.size() > 1) {

			vector<SubGrain*>::iterator it;

			myQuaternion parentori = (*itG)->getOri();
			double* parentbunge = parentori.Quaternion2Euler();
			double oriscatter = (*itG)->get_OriScatterFromPrefOri();
			double see = (*itG)->getSEE();
			double seegrscatter = (*itG)->get_SEEGrainScatterFromPrefOri();
			double seesgrscatter = (*itG)->get_SEESubgrainScatterFromPrefOri();
			double sizescaler = (*itG)->get_RelSizeScalingFromPrefori();

			for (it = ++((*itG)->m_SubGrains.begin());
					it != (*itG)->m_SubGrains.end(); it++) {
				// ID, x, y, z, bunge1, bunge2, bunge3, xmin, xmax, ymin, ymax, zmin, zmax, volume, stored
				myQuaternion* ori = (*it)->get_Orientation();
				double* bunge = ori->Quaternion2Euler();
				file
						<< (*it)->get_ID()
						<< "\t"
						<< (((*it)->getMaxX() - (*it)->getMinX()) / 2
								+ (*it)->getMinX())
						<< "\t"
						<< (((*it)->getMaxY() - (*it)->getMinY()) / 2
								+ (*it)->getMinY()) << "\t";
				if (Settings::PlotDimension != E_2D) {
					file
							<< (((*it)->getMaxZ() - (*it)->getMinZ()) / 2
									+ (*it)->getMinZ()) << "\t";
				}
				file << bunge[0] << "\t" << bunge[1] << "\t" << bunge[2] << "\t"
						<< (*it)->getMinX() << "\t" << (*it)->getMaxX() << "\t"
						<< (*it)->getMinY() << "\t" << (*it)->getMaxY() << "\t";
				if (Settings::PlotDimension != E_2D) {
					file << (*it)->getMinZ() << "\t" << (*it)->getMaxZ()
							<< "\t";
				}
				file << (*it)->get_Volume() << "\t" << (*it)->get_SEE() << "\t"
						<< (*itG)->get_oldID() << "\t"
						<< (double) parentori.misorientationCubicQxQ(ori)
						<< "\t" << oriscatter << "\t" << parentbunge[0] << "\t"
						<< parentbunge[1] << "\t" << parentbunge[2] << "\t"
						<< see << "\t" << seegrscatter << "\t" << seesgrscatter
						<< "\t" << sizescaler << endl;
				delete[] bunge;
			}

			delete[] parentbunge;

		} else {
			// ID, x, y, z, bunge1, bunge2, bunge3, xmin, xmax, ymin, ymax, zmin, zmax, volume, stored
			myQuaternion ori = (*itG)->getOri();
			double* bunge = ori.Quaternion2Euler();

			file
					<< (*itG)->get_ID()
					<< "\t"
					<< (((*itG)->getMaxX() - (*itG)->getMinX()) / 2
							+ (*itG)->getMinX())
					<< "\t"
					<< (((*itG)->getMaxY() - (*itG)->getMinY()) / 2
							+ (*itG)->getMinY()) << "\t";
			if (Settings::PlotDimension != E_2D) {
				file
						<< (((*itG)->getMaxZ() - (*itG)->getMinZ()) / 2
								+ (*itG)->getMinZ()) << "\t";
			}
			file << bunge[0] << "\t" << bunge[1] << "\t" << bunge[2] << "\t"
					<< (*itG)->getMinX() << "\t" << (*itG)->getMaxX() << "\t"
					<< (*itG)->getMinY() << "\t" << (*itG)->getMaxY() << "\t";
			if (Settings::PlotDimension != E_2D) {
				file << (*itG)->getMinZ() << "\t" << (*itG)->getMaxZ() << "\t";
			}
			file << int((*itG)->get_Volume() / m_h / m_h + 0.5) << "\t"
					<< (*itG)->getSEE() << "\t" << (*itG)->get_oldID() << "\t"
					<< 0.0 << "\t" << (*itG)->get_OriScatterFromPrefOri()
					<< "\t" << bunge[0] << "\t" << bunge[1] << "\t" << bunge[2]
					<< "\t" << (*itG)->getSEE() << "\t"
					<< (*itG)->get_SEEGrainScatterFromPrefOri() << "\t"
					<< (*itG)->get_SEESubgrainScatterFromPrefOri() << "\t"
					<< (*itG)->get_RelSizeScalingFromPrefori() << endl;
			delete[] bunge;
		}
	}

	file.close();

	myprofiler.logev("Writing ASCII Diagnostics",
			(omp_get_wtime() - asciidiagn));
}

void microStructureHdl::SaveParenthood() {
	//stores all subgrain IDs and their parent grain IDs to become in a postprocessing of simulation data capable of
	//identifying all sub-grains which are adjacent to a grain-grain boundary even though no neighbors were kept track of
	double ptime = omp_get_wtime();

	//count grains
	vector<Grains*>::iterator itG;
	vector<SubGrain*>::iterator it;
	unsigned int NrGrains = 0;
	if (Settings::NumberOfSubgrains != 0)
		for (itG = ++m_grains.begin(); itG != m_grains.end(); itG++) {
			if ((*itG)->m_SubGrains.size() > 1)
				NrGrains += ((*itG)->m_SubGrains.size() - 1);
			else
				NrGrains++;
		}
	else
		NrGrains = Settings::NumberOfGrains;

	unsigned int* rawdata = NULL;
	rawdata = new unsigned int[2 * NrGrains]; //SubgrainID1, ParentID1, SubgrainID2, ParentID2, ..., SubgrainIDN, ParentIDN

	//populate
	unsigned int i = 0;
	for (itG = ++m_grains.begin(); itG != m_grains.end(); itG++) {
		if ((Settings::NumberOfSubgrains != 0)
				&& (*itG)->m_SubGrains.size() > 1) {
			for (it = ++((*itG)->m_SubGrains.begin());
					it != (*itG)->m_SubGrains.end(); it++) {
				rawdata[i] = (*it)->get_ID();
				rawdata[i + 1] = (*itG)->get_oldID();
				i += 2;
			}
		} else {
			rawdata[i] = (*itG)->get_ID();
			rawdata[i + 1] = (*itG)->get_oldID();
			i += 2;
		}
	}
	if (i != 2 * NrGrains) {
		cout << "Determination of parenthood resulted in inconsistencies!"
				<< endl;
		delete[] rawdata;
		return;
	}

//##MK::Test, can be deleted for ( i = 0; i < 2*NrGrains; i += 2 ) { cout << "SGID/GID = " << rawdata[i] << "\t\t" << rawdata[i+1] << endl; }

	//store in binary file
	stringstream parentfn;
	parentfn << "Parenthood.bin";

	FILE* binaryFile;
	binaryFile = fopen(parentfn.str().c_str(), "w");
	fwrite(rawdata, sizeof(unsigned int), 2 * NrGrains, binaryFile);
	fclose(binaryFile);

	delete[] rawdata,

	myprofiler.logev("Writing BINARY Parenthood", (omp_get_wtime() - ptime));
	cout << "Storing of parenthood in binary file was successful." << endl;
}

//void microStructureHdl::CreateColormap(struct RGB* thecolormap) {
//	double gtime = omp_get_wtime();
//
//#pragma omp parallel
//	{
//		vector<Grains*>::iterator itG;
//		int gid = 0;
//		//	for (itG = ++m_grains.begin(); itG != m_grains.end(); itG++) {
//		//	if ( gid % omp_get_num_threads() == omp_get_thread_num() ) {
//		vector<unsigned int>& workload = m_grainScheduler->getThreadWorkload(
//				omp_get_thread_num());
//		for (auto id : workload) {
//			if ((Settings::NumberOfSubgrains != 0)
//					&& m_grains[id]->m_SubGrains.size() > 1) {
//				vector<SubGrain*>::iterator it;
//				for (it = ++(m_grains[id]->m_SubGrains.begin());
//						it != m_grains[id]->m_SubGrains.end(); it++) {
//					unsigned int id = (*it)->get_ID() - 1;
//					myQuaternion* ori = (*it)->get_Orientation();
//
//					unsigned char argb[3] = { UCHAR_RANGE_MIN, UCHAR_RANGE_MIN,
//							UCHAR_RANGE_MIN };
//					ori->quat2ipfz(argb);
////#pragma omp critical
////{
////	std::cout << "thread/sid/R/G/B/colormap[id].ALPHA = " << omp_get_thread_num() << ";" << "\t\t" << id << "\t\t" << (int) argb[0] << ";" << (int) argb[1] << ";" << (int) argb[2] << "\t\t\t" << ori->get_x() << ";" << ori->get_y() << std::endl;
////}
//
//					thecolormap[id].R = argb[REDCHAN];
//					thecolormap[id].G = argb[GREENCHAN];
//					thecolormap[id].B = argb[BLUECHAN];
//					thecolormap[id].ALPHA = UCHAR_RANGE_MAX; //=255;
//
//				} //for all subgrains
//			} else {
//				unsigned int id = m_grains[id]->get_ID() - 1;
//				myQuaternion ori = m_grains[id]->getOri();
//				unsigned char argb[3] = { UCHAR_RANGE_MIN, UCHAR_RANGE_MIN,
//						UCHAR_RANGE_MIN };
//
//				ori.quat2ipfz(argb);
//
//				thecolormap[id].R = argb[REDCHAN];
//				thecolormap[id].G = argb[GREENCHAN];
//				thecolormap[id].B = argb[BLUECHAN];
//				thecolormap[id].ALPHA = UCHAR_RANGE_MAX; //=255;
////#pragma omp critical
////{
////std::cout << "thread/gid/R/G/B/colormap[id].ALPHA = " << omp_get_thread_num() << ";" << "\t\t" << id << "\t\t" << (int) argb[0] << ";" << (int) argb[1] << ";" << (int) argb[2] << "\t\t\t" << ori.get_x() << ";" << ori.get_y() << std::endl;
////}
//			}
//		}
//	} //end of parallel region
//
//	myprofiler.logev("Parallel IPF color mapping", (omp_get_wtime() - gtime));
//}
void microStructureHdl::CreateColormap( struct RGB* thecolormap )
{
        double gtime = omp_get_wtime();

#pragma omp parallel
        {
                vector<Grains*>::iterator itG;
                int gid = 0;
                for (itG = ++m_grains.begin(); itG != m_grains.end(); itG++) {
                        if ( gid % omp_get_num_threads() == omp_get_thread_num() ) {
                                if ((Settings::NumberOfSubgrains != 0) && (*itG)->m_SubGrains.size() > 1) {
                                        vector<SubGrain*>::iterator it;
                                        for (it = ++((*itG)->m_SubGrains.begin()); it != (*itG)->m_SubGrains.end(); it++) {
                                                unsigned int id = (*it)->get_ID() - 1;
                                                myQuaternion* ori = (*it)->get_Orientation();

                                                unsigned char argb[3] = { UCHAR_RANGE_MIN, UCHAR_RANGE_MIN, UCHAR_RANGE_MIN };
                                                ori->quat2ipfz( argb );
//#pragma omp critical
//{
//      std::cout << "thread/sid/R/G/B/colormap[id].ALPHA = " << omp_get_thread_num() << ";" << "\t\t" << id << "\t\t" << (int) argb[0] << ";" << (int) argb[1] << ";" << (int) argb[2] << "\t\t\t" << ori->get_x() << ";" << ori->get_y() << std::endl;
//}

                                                thecolormap[id].R = argb[REDCHAN];
                                                thecolormap[id].G = argb[GREENCHAN];
                                                thecolormap[id].B = argb[BLUECHAN];
                                                thecolormap[id].ALPHA = UCHAR_RANGE_MAX; //=255;


                                        } //for all subgrains
                                }
                                else {
                                        unsigned int id = (*itG)->get_ID()-1;
                                        myQuaternion ori = (*itG)->getOri();
                                        unsigned char argb[3] = { UCHAR_RANGE_MIN, UCHAR_RANGE_MIN, UCHAR_RANGE_MIN };

                                        ori.quat2ipfz( argb );

                                        thecolormap[id].R = argb[REDCHAN];
                                        thecolormap[id].G = argb[GREENCHAN];
                                        thecolormap[id].B = argb[BLUECHAN];
                                        thecolormap[id].ALPHA = UCHAR_RANGE_MAX; //=255;
//#pragma omp critical
//{
//std::cout << "thread/gid/R/G/B/colormap[id].ALPHA = " << omp_get_thread_num() << ";" << "\t\t" << id << "\t\t" << (int) argb[0] << ";" << (int) argb[1] << ";" << (int) argb[2] << "\t\t\t" << ori.get_x() << ";" << ori.get_y() << std::endl;
//}
                                }
                        } //end the work I do
                        gid++;
                }
        } //end of parallel region

        myprofiler.logev("Parallel IPF color mapping", (omp_get_wtime() - gtime));
}

unsigned int microStructureHdl::CountNumberOfSubgrains() {
	unsigned int NumberOfSubgrains = 0;
	vector<Grains*>::iterator itG;
	for (itG = ++m_grains.begin(); itG != m_grains.end(); itG++) {
		if ((Settings::NumberOfSubgrains != 0)
				&& (*itG)->m_SubGrains.size() > 1) {
			vector<SubGrain*>::iterator it;
			for (it = ++((*itG)->m_SubGrains.begin());
					it != (*itG)->m_SubGrains.end(); it++)
				NumberOfSubgrains++;

		} else {
			NumberOfSubgrains++;
		}
	}
	cout << "Total number of subgrains = " << NumberOfSubgrains << endl;
	return NumberOfSubgrains;
}

void microStructureHdl::PlotIPF2DSection() {
//MK::requires to become executed after an ID rehashing for the subgrains has been performed!

	if (Settings::PlotIPF2DSection == true) {
		double gtime = omp_get_wtime();

		unsigned int xmi = Settings::PlotWindowXMin
				* (double) Settings::NumberOfGridpoints;
		unsigned int xmx = Settings::PlotWindowXMax
				* (double) Settings::NumberOfGridpoints;
		unsigned int ymi = Settings::PlotWindowYMin
				* (double) Settings::NumberOfGridpoints;
		unsigned int ymx = Settings::PlotWindowYMax
				* (double) Settings::NumberOfGridpoints;
		if (xmx >= m_container->getMaxX())
			xmx--; //bounds check
		if (ymx >= m_container->getMaxY())
			ymx--;
		if ((xmx - xmi + 1) > IPFZMAPPING_MAXSIZE) {
			xmi = 0;
			xmx = IPFZMAPPING_MAXSIZE;
			ymi = 0;
			ymx = IPFZMAPPING_MAXSIZE;
		}
		unsigned int z = 0;
		unsigned int xlen = xmx - xmi + 1;

		cout << "Generating a IPF mapping of the domain on " << xmi << ";"
				<< xmx << ";" << ymi << ";" << ymx << " width "
				<< (xmx - xmi + 1) << " x " << (ymx - ymi + 1) << endl;

		unsigned int SubgrainsInTotal = CountNumberOfSubgrains();
		struct RGB* colormap = new struct RGB[SubgrainsInTotal];

//for ( unsigned int id = 0; id < SubgrainsInTotal; id++ ) {cout << "thread/0/R/G/B/colormap[id].ALPHA = " << omp_get_thread_num() << ";" << id << ";" << (int) colormap[id].R << ";" << (int) colormap[id].G << ";" << (int) colormap[id].B << "--" << (int) colormap[id].ALPHA << endl;}cout << endl << endl;

		CreateColormap(colormap);

		cout << "Colormapping of orientations was successful" << endl;

//cout << endl << endl; for ( unsigned int id = 0; id < SubgrainsInTotal; id++ ) { cout << "thread/sid/R/G/B/colormap[id].ALPHA = " << omp_get_thread_num() << ";" << id << ";" << (int) colormap[id].R << ";" << (int) colormap[id].G << ";" << (int) colormap[id].B << "--" << (int) colormap[id].ALPHA << endl;}

		//implicit barrier after pragma omp parallel assures colormap is filled before filling the IPF mapping with cached color values
		gtime = omp_get_wtime();

		unsigned char* ipf = new unsigned char[4 * (xmx - xmi + 1)
				* (ymx - ymi + 1)];
		unsigned int pxid = 0;
		unsigned int pxc = 0;
		//cout << m_container->getMinX() << ";" << m_container->getMaxX() << "\t\t" << m_container->getMinY() << ";" << m_container->getMaxY() << "\t\t" << m_container->getMinZ() << ";" << m_container->getMaxZ() << endl;
		for (unsigned int y = ymi; y <= ymx; y++) {
			for (unsigned int x = xmi; x <= xmx; x++) {
				//row the y, column the x, and depth the z coordinate of the element
				pxid = (m_container->getValueAt(y, x, z)) - 1;
				pxc = 4 * ((x - xmi) + ((y - ymi) * xlen));
				ipf[pxc + REDCHAN] = colormap[pxid].R;
				ipf[pxc + GREENCHAN] = colormap[pxid].G;
				ipf[pxc + BLUECHAN] = colormap[pxid].B;
				ipf[pxc + ALPHACHAN] = 255;
			}
		}

		myprofiler.logev("Mapping of IDs to colors", (omp_get_wtime() - gtime));

		cout
				<< "Mapping of sub-grain IDs to colors was successful, writing image..."
				<< endl;
		gtime = omp_get_wtime();

		ostringstream fname;
		fname << "Microstructure.IPFZ.png";

		lodepng::encode(fname.str().c_str(), ipf, (xmx - xmi + 1),
				(ymx - ymi + 1)); //utilize lodePNG to plot

		delete[] ipf;
		delete[] colormap;

		myprofiler.logev("Plotting the structure", (omp_get_wtime() - gtime));
	}
}

void microStructureHdl::ReportProfile() {
	cout << endl << "Approximate OpenMP omp_get_wtime report (seconds)" << endl;
	for (unsigned int e = 0; e < myprofiler.get_entries(); e++) {
		cout << myprofiler.titles[e] << "\t\t" << setprecision(6)
				<< myprofiler.times[e] << endl;
	}
}

struct NUMANode {
	int num_cpus;
	int numa_cpus[64];
};

unsigned int my_numa_bitmask_weight(const struct bitmask *mask) {
	unsigned int weight = 0;
	for (unsigned int j = 0; j < mask->size; j++) {
		if (numa_bitmask_isbitset(mask, j)) {
			weight++;
		}
	}
	return weight;
}

void microStructureHdl::initEnvironment() {

//Set up correct Maximum Number of threads
	if (Settings::ExecuteInParallel) {
		Settings::MaximumNumberOfThreads = omp_get_max_threads();

	} else {
		Settings::MaximumNumberOfThreads = 1;
		omp_set_num_threads(Settings::MaximumNumberOfThreads);
	}

//These lines might need to be moved if spatial distribution of grains is utilized
//At best the grain scheduler should be configurable through the parameters file

//m_grainScheduler = new IterativeGrainScheduler(Settings::MaximumNumberOfThreads, Settings::NumberOfParticles);

//choose grain scheduler:
//	if (Settings::GrainScheduler == E_ITERATIVE) {
//		m_grainScheduler = new IterativeGrainScheduler(
//				Settings::MaximumNumberOfThreads, Settings::NumberOfParticles);
//	} else if (Settings::GrainScheduler == E_SQUARES) {
//		m_grainScheduler = new SquaresGrainScheduler(
//				Settings::MaximumNumberOfThreads, Settings::NumberOfParticles);
//	} else if (Settings::GrainScheduler == E_DEFAULT_SCHEDULER) {
//		m_grainScheduler = new IterativeGrainScheduler(
//				Settings::MaximumNumberOfThreads, Settings::NumberOfParticles);
//	}
	m_grainScheduler = new IterativeGrainScheduler(
			Settings::MaximumNumberOfThreads, Settings::NumberOfGrains);
	initNUMABindings();
}

void microStructureHdl::initNUMABindings() {
	vector<NUMANode> nodes;
	nodes.reserve(16);
	numa_available();
// returns a mask of CPUs on which the current task is allowed to run.
	bitmask* mask = numa_get_run_node_mask();
	bitmask* cpus = numa_allocate_cpumask();
	for (unsigned int j = 0; j < mask->size; j++) {
		if (numa_bitmask_isbitset(mask, j)) {
			printf("We are allowed to used node %d\n", j);
			NUMANode node;
			memset(&node, 0xFF, sizeof(node));
			//converts a node number to a bitmask of CPUs.
			//The user must pass a bitmask structure with a mask buffer long enough to represent all possible cpu's
			numa_node_to_cpus(j, cpus);
			node.num_cpus = my_numa_bitmask_weight(cpus);
			int cpuCounter = 0;
			for (unsigned int i = 0; i < cpus->size; i++) {

				if (numa_bitmask_isbitset(cpus, i)
						&& numa_bitmask_isbitset(numa_all_cpus_ptr, i)) {
					node.numa_cpus[cpuCounter] = i;
					cpuCounter++;
				}
			}
			nodes.push_back(node);
		}
	}
	numa_free_cpumask(cpus);
#pragma omp parallel
	{
		int threadID = omp_get_thread_num();
		for (unsigned int i = 0; i < nodes.size(); i++) {
			if (threadID < nodes.at(i).num_cpus) {
#pragma omp critical
				{
					printf("Will bind thread %d to cpu %d\n",
							omp_get_thread_num(),
							nodes.at(i).numa_cpus[threadID]);
					cpu_set_t set;
					CPU_ZERO(&set);
					CPU_SET(nodes.at(i).numa_cpus[threadID], &set);
					int res = sched_setaffinity(0, sizeof(set), &set);
					printf(res == 0 ? "Successful\n" : "Failed\n");
				}
				break;
			}
			threadID -= nodes.at(i).num_cpus;
		}
	}
}

void microStructureHdl::initializeGrains(vector<vector<Eigen::Vector3d>> hulls,
		vector<double> grainVolume) {
	cout << "Started Initializing Grains " << endl;
	m_grainScheduler->buildThreadWorkloads(hulls, Settings::NumberOfGridpoints);

	bool exceptionHappened = false;
	string error_message;
	m_grains.resize(Settings::NumberOfGrains + 1);

#pragma omp parallel
	{
		vector<unsigned int>& workload = m_grainScheduler->getThreadWorkload(
				omp_get_thread_num());
		for (auto id : workload) {
			if (id <= Settings::NumberOfGrains)
				try {
					Grains* newGrain = new Grains(id, hulls[id], m_container,
							grainVolume[id]);
					m_grains[id] = newGrain;

				} catch (exception& e) {
#pragma omp critical
					{
						exceptionHappened = true;
						error_message += string("Grain ")
								+ to_string((unsigned long long) id)
								+ " in its constructor! Reason : " + e.what()
								+ string("\n");
					}
				}

			if (exceptionHappened) {
				throw runtime_error(error_message);
			}
		}
	}
}

//void microStructureHdl::updateGlobalVoxelContainer() {
//#pragma omp parallel
//	{
//		vector<unsigned int>& workload = m_grainScheduler->getThreadWorkload(
//				omp_get_thread_num());
//
//		//man darf hier parallel ausführen, aber ist das performant???
//		for (auto id : workload) {
//			if (id <= Settings::NumberOfGrains) {
//				m_grains[id]->copySubgrainsToGlobalContainer(m_container);
//			}
//		}
//
//	}
//}

void microStructureHdl::saveTexture() {

	string filename = string("Texture") + string(".txt");

	FILE* output = fopen(filename.c_str(), "wt");

	for (auto it : m_grains) {
		for (auto sub : it->get_SubGrains()) {

			myQuaternion* ori = sub->get_Orientation();
			double* bunge = ori->Quaternion2Euler();
			fprintf(output, "%lf \t%lf \t%lf\t %lf \n", bunge[0], bunge[1],
					bunge[2], sub->get_Volume());
			delete[] bunge;
		}
	}

	fclose(output);
}

void microStructureHdl::plotGrains() {
	vector<Grains*>::iterator it;
	for (it = ++m_grains.begin(); it != m_grains.end(); it++) {
		(*it)->plotLocalContainer();
	}
}

void microStructureHdl::testprng(unsigned int n, double mu, double sigma) {
//get some output from Marsaglia's SHR3 powered Ziggurat for standardnormal distributed
	/*for ( int tid = 0; tid < 128; tid++ ) {
	 m_seqRND->initSHR3( (uint32_t) pow(2.0, 31) - tid - 1 );
	 m_seqRND->initR4Uni( (uint32_t) pow(2.0, 31) - tid - 1 );
	 m_seqRND->r4_nor_setup();
	 for ( unsigned int i = 0; i < n; i++ ) {
	 //rescale
	 cout << setprecision(6) << (m_seqRND->r4_nor() * sigma) + mu << endl;
	 }
	 }*/

	for (int tid = 0; tid < 128; tid++) {
		randomClass agen;
		agen.initMT((uint32_t) pow(2.0, 31) - tid - 1);
		for (unsigned int i = 0; i < n; i++) {
			cout << setprecision(6) << agen.MersenneTwister() << endl;
		}
	}
}

