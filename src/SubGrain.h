/*
 * SubGrain.h
 *
 *  Created on: 07.06.2016
 *      Author: jn538375
 */

#ifndef SUBGRAIN_H_
#define SUBGRAIN_H_

#include <string>
#include "dimensionalBuffer.h"
#include "myQuaternion.h"
#include "RTree.h"
class container;
class Grains;

using namespace std;

class SubGrain {
public:
	//Constructors
	SubGrain();
	SubGrain(myQuaternion ori, double SEE, double volume);
	SubGrain(int id, vector<double> cellcoordinates, double volume, Grains* owner);

	//Destructor
	virtual ~SubGrain();

	//Set Functions
	void setID(int id);
	inline void set_Orientation(myQuaternion ori) {
		if (m_orientation == NULL)
			m_orientation = new myQuaternion(); //##MK::potential memory leak, because specificallyDisorient has already allocated heap object
		*m_orientation = ori;
	}
	inline myQuaternion* get_Orientation() {
		return m_orientation;
	}
	inline double get_Volume() {
		return m_Volume;
	}
	inline int get_ID() {
		return m_ID;
	}
	inline double get_SEE() {
		return m_SEE;
	}
	inline void set_SEE(double see) {
		m_SEE = see;
	}
	inline int getMinX() const {
		return m_xmin;
	}
	inline int getMaxX() const {
		return m_xmax;
	}
	inline int getMinY() const {
		return m_ymin;
	}
	inline int getMaxY() const {
		return m_ymax;
	}
	inline int getMinZ() const {
		return m_zmin;
	}
	inline int getMaxZ() const {
		return m_zmax;
	}


private:
	myQuaternion* m_orientation;
	double m_SEE;
	double m_Volume;
	int m_ID;
	int m_xmin;
	int m_xmax;
	int m_ymin;
	int m_ymax;
	int m_zmin;
	int m_zmax;
	Grains* m_owner; //MK::do not delete only backreference

};

#endif /* SUBGRAIN_H_ */
