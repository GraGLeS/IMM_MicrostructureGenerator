/*
 * Grains.h
 *
 *  Created on: May 4, 2016
 *      Author: jn538375
 */

#ifndef GRAINS_H_
#define GRAINS_H_

#include "RTree.h"
#include <string>
#include "dimensionalBuffer.h"
#include "myQuaternion.h"
#include "Eigen/Dense"
#include "utilities.h"
#include "Settings.h"

class mathMethods;
class container;
class SubGrain;
struct myPreferenceOri {
	myQuaternion ori;
	double subgrainsScatterOri;
	double SEE;
	double SEEGrainScatter;
	double SEESubgrainScatter;
	double RelSubgrainSizeScaling; //specifies the inverse of how strongly the average sub-grain diameter should be changed
	myPreferenceOri(myQuaternion _ori) { //preference orientation identifying default properties for all grains far away from all distinguished components
		ori = _ori;
		SEE = 0.5 * (Settings::StoredElasticEnergyMax + Settings::StoredElasticEnergyMin); //interpreted as mu of the distribution
		;
		SEEGrainScatter = Settings::StoredElasticScatterGrain; //interpreted as the sigma of the distribution
		SEESubgrainScatter = Settings::StoredElasticScatterSubgrain;
		subgrainsScatterOri = Settings::SubgrainOriScatter*PI/180.0; //interpreted as the sigma of a rAyleigh distribution
		RelSubgrainSizeScaling = 1.0;
	}
	myPreferenceOri(double phi1, double PHI, double phi2, double _oriScatter,
			double _SEE, double _SEEGrainScatter, double _SEESubgrainScatter, double _RelSizeScaler) :
			SEE(_SEE), SEEGrainScatter(_SEEGrainScatter), SEESubgrainScatter(_SEESubgrainScatter), subgrainsScatterOri(_oriScatter), RelSubgrainSizeScaling(_RelSizeScaler) {
		double euler[3] = { phi1, PHI, phi2 };
		ori.euler2Quaternion(euler);
	}
	myPreferenceOri(myQuaternion _ori, double _oriScatter, double _SEE, double _SEEGrainScatter, double _SEESubgrainScatter, double _RelSizeScaler) {
		ori = _ori;
		subgrainsScatterOri = _oriScatter;
		SEE = _SEE;
		SEEGrainScatter = _SEEGrainScatter;
		SEESubgrainScatter = _SEESubgrainScatter;
		RelSubgrainSizeScaling = _RelSizeScaler;
	}
	myPreferenceOri& operator=(const myPreferenceOri &rhs) {
		if (this != &rhs) //oder if (*this != rhs)
				{
			ori = rhs.ori;
			SEE = rhs.SEE;
			SEEGrainScatter = rhs.SEEGrainScatter;
			SEESubgrainScatter = rhs.SEESubgrainScatter;
			subgrainsScatterOri = rhs.subgrainsScatterOri;
			RelSubgrainSizeScaling = rhs.RelSubgrainSizeScaling;
		}
		return *this;
	}
};

using namespace std;

class Grains {
public:
	friend class SubGrain;
	//Constructors
	Grains();
	Grains(myQuaternion ori, double SEE);
	Grains(int id, vector<Eigen::Vector3d>& hull,
			DimensionalBuffer<unsigned int>* container, double volume);
	Grains(int id, int minX, int maxX, int minY, int maxY, int minZ, int maxZ,
			DimensionalBuffer<unsigned int>* container);

	//Destructor
	virtual ~Grains();

	void SubGrainConstructor(randomClass& r);
	//set functions
	void set_Orientation(myQuaternion ori);

	void computeDirectNeighbours(
			const RTree<unsigned int, int, 3, float>& tree);
	void copyVoxelData(DimensionalBuffer<unsigned int> *container);
	int copySubgrainsToGlobalContainer(
			DimensionalBuffer<unsigned int> *container, int offset);

	void plotDensityOfOrientation();
	void plotLocalContainer();
	void copyContainerToGrain(DimensionalBuffer<unsigned int>* container);

	void generateSubGrainOri(randomClass& r);
	void generateSubGrainSEE(randomClass& r);
	int rehashAllSubGrains(int offset);
	//get functions
	inline myQuaternion getOri() {
		return *m_orientation;
	}
	inline double get_Volume() {
		return m_Volume;
	}
	inline double getSEE() {
		return m_SEE;
	}
	inline void set_SEE(double see) {
		m_SEE = see;
	}
	inline int getMinX() const {
		return m_localContainer->getMinX();
	}
	inline int getMaxX() const {
		return m_localContainer->getMaxX();
	}
	inline int getMinY() const {
		return m_localContainer->getMinY();
	}
	inline int getMaxY() const {
		return m_localContainer->getMaxY();
	}
	inline int getMinZ() const {
		return m_localContainer->getMinZ();
	}
	inline int getMaxZ() const {
		return m_localContainer->getMaxZ();
	}
	inline vector<SubGrain*>& get_SubGrains() {
		return m_SubGrains;
	}
	inline unsigned int get_ID() {
		return (unsigned int)m_ID;
	}
	inline unsigned int get_oldID() {
		return (unsigned int)m_oldID;
	}
	inline void set_PreforiProperties(myPreferenceOri PrefOri){
		m_PrefOri = new myPreferenceOri(PrefOri);
	}
	inline double get_SEEFromPrefori(){
		return m_PrefOri->SEE;
	}
	inline double get_SEEGrainScatterFromPrefOri(){
		return m_PrefOri->SEEGrainScatter;
	}
	inline double get_SEESubgrainScatterFromPrefOri(){
		return m_PrefOri->SEESubgrainScatter;
	}
	inline double get_OriScatterFromPrefOri(){
		return m_PrefOri->subgrainsScatterOri;
	}
	inline double get_PrefOriQuatQ0() {
		return m_PrefOri->ori.get_q0();
	}
	inline double get_PrefOriQuatQ1() {
		return m_PrefOri->ori.get_q1();
	}
	inline double get_PrefOriQuatQ2() {
		return m_PrefOri->ori.get_q2();
	}
	inline double get_PrefOriQuatQ3() {
		return m_PrefOri->ori.get_q3();
	}
	inline double get_RelSizeScalingFromPrefori() {
		return m_PrefOri->RelSubgrainSizeScaling;
	}
	vector<SubGrain*> m_SubGrains;
private:
	myQuaternion* m_orientation;
	DimensionalBuffer<unsigned int>* m_localContainer;
	vector<unsigned int> m_NeighborCandidates;
	double m_SEE; //stored elastic energy
	double m_Volume;
	int m_ID;
	int m_oldID;		//mean to store clear and unambiguous the id that the parent grain had before the global ids were rehashed
	myPreferenceOri* m_PrefOri;
};

#endif /* GRAINS_H_ */
