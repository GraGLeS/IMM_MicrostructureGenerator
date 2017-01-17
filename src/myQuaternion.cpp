/*
 * myQuaternion.cpp
 *
 *  Created on: 17.09.2015
 *      Author: miessen
 */
#include "myQuaternion.h"
#include "mymath.h"
#include "utilities.h"

myQuaternion::myQuaternion(double phi1, double PHI, double phi2){
	double euler[3]={phi1,PHI,phi2};
	euler2Quaternion(euler);
}

myQuaternion::myQuaternion(void) :
		q0(1), q1(0), q2(0), q3(0), x(FAIL_MYMATH_NUMERIC), y(FAIL_MYMATH_NUMERIC) {
}

myQuaternion::myQuaternion(double q0, double q1, double q2, double q3) :
		q0(q0), q1(q1), q2(q2), q3(q3), x(FAIL_MYMATH_NUMERIC), y(FAIL_MYMATH_NUMERIC) {
}

myQuaternion& myQuaternion::operator=(const myQuaternion & rhs) {
	if (this != &rhs) //oder if (*this != rhs)
			{
		q0 = rhs.q0;
		q1 = rhs.q1;
		q2 = rhs.q2;
		q3 = rhs.q3; //Copy-Konstruktor
		x = rhs.x;
		y = rhs.y;
	}
	return *this; //return reference on the object itself
}

myQuaternion::myQuaternion(double alpha, double x, double y, double z,
		bool radiant) {

	//##MK::URGENT CHECK NECESSARY IT IS VERY LIKELY INCONSISTENT WITH OUR BUNGE INTERPRETATION AND THAT WE USE "ZXZ" convention!

	double _normVector = 1.0 / sqrt(SQR(x) + SQR(y) + SQR(z));
	if (!radiant)
		alpha = alpha * PI / 180;

	// alpha in radiant
	// Representing Attitude: Euler Angles, Unit Quaternions, and Rotation Vectors,	James Diebel, 2006
	q0 = cos(alpha / 2);
	q1 = x * sin(alpha / 2) * _normVector;
	q2 = y * sin(alpha / 2) * _normVector;
	q3 = z * sin(alpha / 2) * _normVector;

	x = FAIL_MYMATH_NUMERIC;
	y = FAIL_MYMATH_NUMERIC;
}


myQuaternion::~myQuaternion(void) {
}

double myQuaternion::getNorm(void) {
	double qNorm = sqrt(SQR(q1) + SQR(q2) + SQR(q3) + SQR(q0));
	return qNorm;
}
void myQuaternion::Normalize() {
	double qNorm = getNorm();
	q0 /= qNorm;
	q1 /= qNorm;
	q2 /= qNorm;
	q3 /= qNorm;
}


void myQuaternion::Invert() {
	q1 *= -1.0;
	q2 *= -1.0;
	q3 *= -1.0;
}

myQuaternion myQuaternion::Inverse() {
	myQuaternion r(q0, -q1, -q2, -q3);
	return r;
}

myQuaternion myQuaternion::operator*(const myQuaternion & rhs) {
	//MK::ok, mathematically multiplies quaternions q and p, active or passive rotation convention does not matter
	//verified via http://mathworld.wolfram.com/Quaternion.html as well as http://www.mathworks.de/de/help/aeroblks/quaternionmultiplication.html
	//equivalent to the multiplication given in Grimmer, H, 1974 Acta Cryst A30, 685-688
	//mind vector cross product notation qp = q0p0 - qbar*pbar + q0 pbar + p0 qbar + qbar cross pbar with bar vector quantities parts of the quaternion
	myQuaternion r;
	r.q0 = +q0 * rhs.q0 - q1 * rhs.q1 - q2 * rhs.q2 - q3 * rhs.q3;
	r.q1 = +q1 * rhs.q0 + q0 * rhs.q1 - q3 * rhs.q2 + q2 * rhs.q3;
	r.q2 = +q2 * rhs.q0 + q3 * rhs.q1 + q0 * rhs.q2 - q1 * rhs.q3;
	r.q3 = +q3 * rhs.q0 - q2 * rhs.q1 + q1 * rhs.q2 + q0 * rhs.q3;
	return r;
}

void myQuaternion::Sort() {
	double arr[4] = { q0, q1, q2, q3 };
	int last = 4 - 2;
	int isChanged = 1;

	while (last >= 0 && isChanged) {
		isChanged = 0;
		for (int k = 0; k <= last; k++)
			if (arr[k] > arr[k + 1]) {
				SWAP(arr[k], arr[k + 1]);
				isChanged = 1;
			}
		last--;
	}
	q3 = arr[0];
	q2 = arr[1];
	q1 = arr[2];
	q0 = arr[3];
}

myQuaternion myQuaternion::misorientationQuaternionCubic(myQuaternion* p) {
	//MK::ok
	myQuaternion qm1; //Inverse of quaternion q

	//Inverse of quaternion q is the same like the conjugate for unit quaternions
	qm1 = Inverse();

	myQuaternion result; //Resulting misorientation quaternion, m = pq-1

	result = qm1 * (*p); //MK:: was p post applied to qm1 but that was not consistent with MTex, in particular not with -(q^-1) = -q0,q1,q2,q3

	//Now, we have to determine the smallest angle, following Grimmer H, Acta Cryst 1974, A30 pp685-688
	double r0[6][4]; //There are 12 possible angles

	//Note: The notation r0 is due to the definition of the quaternion which lie
	//in the fundamental zone, this vector possesses the smallest angle, in such a way
	//that r0 is actually the scalar part of this quaternion

	double a = result.q0;
	double b = result.q1;
	double c = result.q2;
	double d = result.q3;

	double fac = 1.0 / sqrt(2.0); //0.70710678;

	//six fundamental quaternions
	r0[0][0] = (a + b) * fac;
	r0[0][1] = (a - b) * fac;
	r0[0][2] = (c + d) * fac;
	r0[0][3] = (c - d) * fac;
	r0[1][0] = (a + c) * fac;
	r0[1][1] = (a - c) * fac;
	r0[1][2] = (b + d) * fac;
	r0[1][3] = (b - d) * fac;
	r0[2][0] = (a + d) * fac;
	r0[2][1] = (a - d) * fac;
	r0[2][2] = (b + c) * fac;
	r0[2][3] = (b - c) * fac;
	r0[3][0] = (a + b + c + d) * 0.5;
	r0[3][1] = (a + b - c - d) * 0.5;
	r0[3][2] = (a - b + c - d) * 0.5;
	r0[3][3] = (a - b - c + d) * 0.5;
	r0[4][0] = (a + b + c - d) * 0.5;
	r0[4][1] = (a + b - c + d) * 0.5;
	r0[4][2] = (a - b + c + d) * 0.5;
	r0[4][3] = (a - b - c - d) * 0.5;
	r0[5][0] = a;
	r0[5][1] = b;
	r0[5][2] = c;
	r0[5][3] = d;

	int mi = 0;
	double max = 0.0;

	for (int i = 0; i < 6; i++) //Determing the quaternion with the maximal component and the component itself
		for (int j = 0; j < 4; j++) {
			if (fabs(r0[i][j]) > max) {
				max = fabs(r0[i][j]);
				mi = i;
			}
		}

	result.q0 = fabs(r0[mi][0]); //Disorientation requires all components positive
	result.q1 = fabs(r0[mi][1]);
	result.q2 = fabs(r0[mi][2]);
	result.q3 = fabs(r0[mi][3]);

	result.Sort(); //Sorting into ascending order, because a desorientation in the SST
	//requires a quaternion with q0>=q1>=q2>=q3 which represents a minimal
	return result;
	//additionally it is required that rq3 >= sum of all others and rq3 * (sqrt(2)-1) >= rq2
}

void myQuaternion::euler2Quaternion(double *euler) {
	/*20130326MK convention: Bunge ZXZ, which represents the (3,1,3) case analyzed in: Diebel 2006
	 Representing Attitude: Euler Angles, Unit Quaternions, and Rotation Vectors, James Diebel, Stanford University, Stanford, California 94301-9010, Email: diebel@stanford.edu
	 //cos(a+b) = c(a+b) = cacb-sasb
	 //cos(a-b) = c(a-b) = cacb+sasb
	 //sin(a+b) = s(a+b) = sacb+casb
	 //sin(a-b) = s(a-b) = sacb-casb
	 */

	double p1 = euler[0];
	double t = euler[1];
	double p2 = euler[2];

	Real co1 = cos(t / 2);
	Real s1 = sin(t / 2);

	//double test[4]={0};
	//quaternion2Euler( p,test);

	q0 = co1 * cos((p1 + p2) / 2);
	q1 = s1 * cos((p1 - p2) / 2);
	q2 = s1 * sin((p1 - p2) / 2);
	q3 = co1 * sin((p1 + p2) / 2);
}

double* myQuaternion::Quaternion2Euler(void) {
	//convention: Bunge, ZXZ, equal to case (3,1,3) as
	//  analyzed in Diebel, James, 2006:
	//  Representing Attitude: Euler Angles, Unit Quaternions and Rotation Vectors
	//Gimbal lock situation analyzed following the line of Melcher et. al., Conversion of EBSD data by a quaternion based algorithm....
	//TECHNISCHE MECHANIK, 30, 4, (2010), 401--413
	//dont forget to define QUAT2EUL_ETA 1e-20

	Real PHI, sP, phi1, phi2;

	Real cosPHI = SQR(q3) - SQR(q2) - SQR(q1) + SQR(q0);

	Real y0 = 2 * q1 * q3 - 2 * q0 * q2;
	Real x0 = 2 * q2 * q3 + 2 * q0 * q1;
	Real y1 = 2 * q1 * q3 + 2 * q0 * q2;
	Real x1 = -2 * q2 * q3 + 2 * q0 * q1;

	if (cosPHI > 1.)
		cosPHI = 1.;

	if (SQR(1. - cosPHI) <= QUAT2EUL_ETA)
		PHI = 0.;
	else
		PHI = acos(cosPHI);

	sP = sin(PHI); //handle the gimbal lock situation that arouses as a quarternion does not define a Bunge Euler angle uniquely

	if (sP != 0) {
		phi2 = atan2(y0 / sP, x0 / sP);
		phi1 = atan2(y1 / sP, x1 / sP);
	} else {
		phi1 = atan2((2 * q1 * q2 + 2 * q0 * q3),
		SQR(q0) + SQR(q1) - SQR(q2) - SQR(q3));
		phi2 = 0.;
	}

	//without additional sample and crystal symmetry the Euler space is symmetric to 0 <= phi1 <= 2*_PI_ , 0 <= PHI <= _PI_, 0 <= phi2 <= 2*_PI_
	if (phi1 < 0.0)
		phi1 += 2 * PI;
	if (phi2 < 0.0)
		phi2 += 2 * PI;

	double *euler = new double[3];
	euler[2] = phi2; //following the notation order used in Diebel, James 2006
	euler[1] = PHI;
	euler[0] = phi1;
	return euler;
}

double* myQuaternion::Quaternion2EulerConst(void) const {
	//convention: Bunge, ZXZ, equal to case (3,1,3) as
	//  analyzed in Diebel, James, 2006:
	//  Representing Attitude: Euler Angles, Unit Quaternions and Rotation Vectors
	//Gimbal lock situation analyzed following the line of Melcher et. al., Conversion of EBSD data by a quaternion based algorithm....
	//TECHNISCHE MECHANIK, 30, 4, (2010), 401--413
	//dont forget to define QUAT2EUL_ETA 1e-20

	Real PHI, sP, phi1, phi2;

	Real cosPHI = SQR(q3) - SQR(q2) - SQR(q1) + SQR(q0);

	Real y0 = 2 * q1 * q3 - 2 * q0 * q2;
	Real x0 = 2 * q2 * q3 + 2 * q0 * q1;
	Real y1 = 2 * q1 * q3 + 2 * q0 * q2;
	Real x1 = -2 * q2 * q3 + 2 * q0 * q1;

	if (cosPHI > 1.)
		cosPHI = 1.;

	if (SQR(1. - cosPHI) <= QUAT2EUL_ETA)
		PHI = 0.;
	else
		PHI = acos(cosPHI);

	sP = sin(PHI); //handle the gimbal lock situation that arouses as a quarternion does not define a Bunge Euler angle uniquely

	if (sP != 0) {
		phi2 = atan2(y0 / sP, x0 / sP);
		phi1 = atan2(y1 / sP, x1 / sP);
	} else {
		phi1 = atan2((2 * q1 * q2 + 2 * q0 * q3),
		SQR(q0) + SQR(q1) - SQR(q2) - SQR(q3));
		phi2 = 0.;
	}

	//without additional sample and crystal symmetry the Euler space is symmetric to 0 <= phi1 <= 2*_PI_ , 0 <= PHI <= _PI_, 0 <= phi2 <= 2*_PI_
	if (phi1 < 0.0)
		phi1 += 2 * PI;
	if (phi2 < 0.0)
		phi2 += 2 * PI;

	double *euler = new double[3];
	euler[2] = phi2; //following the notation order used in Diebel, James 2006
	euler[1] = PHI;
	euler[0] = phi1;
	return euler;
}


void myQuaternion::randomOriShoemakeQuat(randomClass& r) {
	//K. Shoemake, Graphic Gems III (editor D. Kirk) CalTech pp124-134
	double X0 = r.MersenneTwister(); //parkMiller();
	double X1 = r.MersenneTwister(); //parkMiller();
	double X2 = r.MersenneTwister(); //parkMiller();

	double r1 = sqrt(1 - X0);
	double r2 = sqrt(X0);
	double theta1 = 2 * _PI_ * X1;
	double theta2 = 2 * _PI_ * X2;

	q0 = r1 * sin(theta1); //w
	q1 = r1 * cos(theta1); //v
	q2 = r2 * sin(theta2); //x
	q3 = r2 * cos(theta2); //z

	Normalize();
}


myQuaternion* myQuaternion::randomMisorientationShoemake(double theta,
		randomClass& thernd) {
	double q[4] = { 0, 0, 0, 0 };
	double qcrit = cos(0.5 * theta);

	if (theta < MINIMUM_ANGULAR_SPREAD) {
		qcrit = cos(0.5 * MINIMUM_ANGULAR_SPREAD);
	}

	while (q[0] < qcrit) {

		double X0 = thernd.MersenneTwister(); //parkMiller();
		double X1 = thernd.MersenneTwister(); //parkMiller();
		double X2 = thernd.MersenneTwister(); //parkMiller();

		double r1 = sqrt(1 - X0);
		double r2 = sqrt(X0);
		double theta1 = 2 * _PI_ * X1;
		double theta2 = 2 * _PI_ * X2;

		q[0] = r1 * sin(theta1); //w
		q[1] = r1 * cos(theta1); //v
		q[2] = r2 * sin(theta2); //x
		q[3] = r2 * cos(theta2); //z

		double qnorm = sqrt( SQR(q[0]) + SQR(q[1]) + SQR(q[2]) + SQR(q[3]));

		//normalize
		q[0] = q[0] / qnorm;
		q[1] = q[1] / qnorm;
		q[2] = q[2] / qnorm;
		q[3] = q[3] / qnorm;
	}

	myQuaternion* qr = new myQuaternion(q[0], q[1], q[2], q[3]);
	return qr;
}


myQuaternion* myQuaternion::randomRotationSO3Shoemake( randomClass& thernd )
{
	double rot[4] = { 0, 0, 0, 0 };

	//K. Shoemake, Graphic Gems III (editor D. Kirk) CalTech pp124-134
	double X0 = thernd.MersenneTwister(); //parkMiller();
	double X1 = thernd.MersenneTwister(); //parkMiller();
	double X2 = thernd.MersenneTwister(); //parkMiller();

	double r1 = sqrt(1 - X0);
	double r2 = sqrt(X0);
	double theta1 = 2 * _PI_ * X1;
	double theta2 = 2 * _PI_ * X2;

	rot[0] = r1 * sin(theta1); //w
	rot[1] = r1 * cos(theta1); //v
	rot[2] = r2 * sin(theta2); //x
	rot[3] = r2 * cos(theta2); //z

	double rotnorm = sqrt( SQR(rot[0]) + SQR(rot[1]) + SQR(rot[2]) + SQR(rot[3]));

	myQuaternion* qr = new myQuaternion(rot[0], rot[1], rot[2], rot[3]);
	return qr;
}


double myQuaternion::misorientationCubicQxQ(myQuaternion* p) {
	int i;
	myQuaternion qm1 = (*this).Inverse(); //Inverse of quaternion this
	myQuaternion r = qm1 * (*p); //Resulting quaternion, rotation of the two previous quaternions pq-1

	//Now, we have to determine the smallest angle.

	Real r0[6][4]; //There are 12 possible angles

	//Note: The notation r0 is due to the definition of the quaternion which lie
	//in the fundamental zone, this vector possesses the smallest angle, in such a way
	//that r0 is actually the scalar part of this quaternion

	double a = r.q0;
	double b = r.q1;
	double c = r.q2;
	double d = r.q3;

	Real fac = 1.0 / sqrt(2.0); //0.70710678;

	//six fundamental quaternions
	r0[0][0] = (a + b) * fac;
	r0[0][1] = (a - b) * fac;
	r0[0][2] = (c + d) * fac;
	r0[0][3] = (c - d) * fac;
	r0[1][0] = (a + c) * fac;
	r0[1][1] = (a - c) * fac;
	r0[1][2] = (b + d) * fac;
	r0[1][3] = (b - d) * fac;
	r0[2][0] = (a + d) * fac;
	r0[2][1] = (a - d) * fac;
	r0[2][2] = (b + c) * fac;
	r0[2][3] = (b - c) * fac;
	r0[3][0] = (a + b + c + d) * 0.5;
	r0[3][1] = (a + b - c - d) * 0.5;
	r0[3][2] = (a - b + c - d) * 0.5;
	r0[3][3] = (a - b - c + d) * 0.5;
	r0[4][0] = (a + b + c - d) * 0.5;
	r0[4][1] = (a + b - c + d) * 0.5;
	r0[4][2] = (a - b + c + d) * 0.5;
	r0[4][3] = (a - b - c - d) * 0.5;
	r0[5][0] = a;
	r0[5][1] = b;
	r0[5][2] = c;
	r0[5][3] = d;

	Real omega = 0.0;

	for (i = 0; i < 6; i++)
		for (int j = 0; j < 4; j++)
			if (fabs(r0[i][j]) > omega)
				omega = fabs(r0[i][j]);

	if (omega > 1.0) //avoid singularity of acos function
		omega = (Real) (int) omega;

	omega = 2 * acos(omega);
	return omega;
}

myQuaternion* myQuaternion::specificallyDisorientednewOriFromReference(
		double sigma_rayl_dis2bunge, randomClass& r) {
	myQuaternion* cand = NULL;

	double theta;
	double prob;
	double coinFlip;
	bool foundvalid = false;
	//acceptance rejection assuming the subgrains to follow a Rayleigh disorientation distribution

	while (foundvalid == false) {

		cand = this->randomRotationSO3Shoemake( r ); //MK::allocates heap and returns pointer to candidate...

		myQuaternion newori = (*cand) * (*this);

		theta = this->misorientationCubicQxQ( &newori );
		prob = 1.0 - exp(-0.5 * SQR(theta / sigma_rayl_dis2bunge));

		//acceptance rejectance scheme
		coinFlip = r.MersenneTwister();

		if (prob <= coinFlip) {
			foundvalid = true;
			break;
		}

		//cand was unsuccessful, so we need to find another one
		delete cand;
		cand = NULL;
	}

	//cand points to a successful rotor to get a properly disoriented orientation to the parent grain
	myQuaternion* quat = new myQuaternion();
	*quat = (*cand) * (*this);

	delete cand;
	return quat; //return pointer to properly allocated quaternion on the heap...
}


/*MK::OBSOLETE!
 * myQuaternion* myQuaternion::newOrientationFromReference( randomClass& r) {
	myQuaternion* newori = new myQuaternion();
	myQuaternion* ori = NULL; //##MK::memory leak

	ori = randomRotationSO3Shoemake( r );

	*newori = (*ori) * (*this);
	return newori;
}*/


void myQuaternion::QuatOnVector3D( double * q, double* v, double* r )
{
	//in accordance with Spieß2009 and Morawiec Pospiech 1989

	//get rotation matrix components
	double r11 = SQR(q0) + SQR(q1) - SQR(q2) - SQR(q3);
	double r12 = 2*(q1*q2 - q0*q3);
	double r13 = 2*(q1*q3 + q0*q2);

	double r21 = 2*(q1*q2 + q0*q3);
	double r22 = SQR(q0) - SQR(q1) + SQR(q2) - SQR(q3);
	double r23 = 2*(q2*q3 - q0*q1);

	double r31 = 2*(q1*q3 - q0*q2);
	double r32 = 2*(q2*q3 + q0*q1);
	double r33 = SQR(q0) - SQR(q1) - SQR(q2) + SQR(q3);

	//x,y,z vector3d v
	r[0] = (r11 * v[0]) + (r12 * v[1]) + (r13 * v[2]);
	r[1] = (r21 * v[0]) + (r22 * v[1]) + (r23 * v[2]);
	r[2] = (r31 * v[0]) + (r32 * v[1]) + (r33 * v[2]);
}


void myQuaternion::multiplyQuaternions( double *q, double* p, double* r )
{
	//MK::ok, mathematically multiplies quaternions q and p, active or passive rotation convention does not matter
	//verified via http://mathworld.wolfram.com/Quaternion.html as well as http://www.mathworks.de/de/help/aeroblks/quaternionmultiplication.html
	//equivalent to the multiplication given in Grimmer, H, 1974 Acta Cryst A30, 685-688
	//mind vector cross product notation qp = q0p0 - qbar*pbar + q0 pbar + p0 qbar + qbar cross pbar with bar vector quantities parts of the quaternion
	r[0] = + q[0] *	p[0]	- q[1] * p[1]	- q[2] *	p[2]	- q[3] *	p[3];
	r[1] = + q[1] *	p[0]	+ q[0] * p[1]	- q[3] *	p[2]	+ q[2] *	p[3];
	r[2] = + q[2] *	p[0]	+ q[3] * p[1]	+ q[0] *	p[2]	- q[1] *	p[3];
	r[3] = + q[3] *	p[0]	- q[2] * p[1]	+ q[1] *	p[2]	+ q[0] *	p[3];
}


void myQuaternion::project2fundamentalregion_ipfz( void )
{
	//projects itself into fundamental region and stores location x,y in standard triangle

	//MK::ok, define crystallographic m-3m symmetry operators in quaternion representation
	double qsymm[SYMMETRIES_IN_CUBIC][4];

	double _sqrt2 = 1.0 / pow ( 2.0, 0.5);
	double half = 0.5;

	//<100>90*degrees four-fold symmetries
	qsymm[0][0] = 1.0;				qsymm[0][1] = 0.0;				qsymm[0][2] = 0.0;				qsymm[0][3] = 0.0;			/*identity*/
	qsymm[1][0] = _sqrt2;			qsymm[1][1] = _sqrt2;			qsymm[1][2] = 0.0;				qsymm[1][3] = 0.0;			/*[100]90*/
	qsymm[2][0] = 0.0;				qsymm[2][1] = 1.0;				qsymm[2][2] = 0.0;				qsymm[2][3] = 0.0;			/*[100]180*/
	qsymm[3][0] = -1.0 * _sqrt2;	qsymm[3][1] = _sqrt2;			qsymm[3][2] = 0.0;				qsymm[3][3] = 0.0;			/*[100]270*/
	qsymm[4][0] = _sqrt2;			qsymm[4][1] = 0.0;				qsymm[4][2] = _sqrt2;			qsymm[4][3] = 0.0;			/*[010]90*/ //MTex -sq2 0 -sq2 0
	qsymm[5][0] = 0.0;				qsymm[5][1] = 0.0;				qsymm[5][2] = 1.0;				qsymm[5][3] = 0.0;			/*[010]180*/ //##MK::be careful, this is consistent with Pomana it is inconsistent with Mtex stating 0,0,-1,0 but q = -q are the same quaternions
	qsymm[6][0] = -1.0 * _sqrt2;	qsymm[6][1] = 0.0;				qsymm[6][2] = _sqrt2;			qsymm[6][3] = 0.0;			/*[010]270*/
	qsymm[7][0] = _sqrt2;			qsymm[7][1] = 0.0;				qsymm[7][2] = 0.0;				qsymm[7][3] = _sqrt2;		/*[001]90*/
	qsymm[8][0] = 0.0;				qsymm[8][1] = 0.0;				qsymm[8][2] = 0.0;				qsymm[8][3] = 1.0;			/*[001]180*/
	qsymm[9][0] = -1.0 * _sqrt2;	qsymm[9][1] = 0.0;				qsymm[9][2] = 0.0;				qsymm[9][3] = _sqrt2;		/*[001]270*/

	//<110>180*degrees two-fold symmetries
	qsymm[10][0] = 0.0;				qsymm[10][1] = _sqrt2;			qsymm[10][2] = _sqrt2;			qsymm[10][3] = 0.0;			/*[110]180*/
	qsymm[11][0] = 0.0;				qsymm[11][1] = _sqrt2;			qsymm[11][2] = -1.0 * _sqrt2;	qsymm[11][3] = 0.0;			/*[1-10]180*/
	qsymm[12][0] = 0.0;				qsymm[12][1] = _sqrt2;			qsymm[12][2] = 0.0;				qsymm[12][3] = _sqrt2;		/*[101]180*/
	qsymm[13][0] = 0.0;				qsymm[13][1] = -1.0 * _sqrt2;	qsymm[13][2] = 0.0;				qsymm[13][3] = _sqrt2;		/*[-101]180*/ //mtex332 0 sq2 0 -sq2
	qsymm[14][0] = 0.0;				qsymm[14][1] = 0.0;				qsymm[14][2] = _sqrt2;			qsymm[14][3] = _sqrt2;		/*[011]180*/ //mtex 0 0 -sq -sq
	qsymm[15][0] = 0.0;				qsymm[15][1] = 0.0;				qsymm[15][2] = -1.0 * _sqrt2;	qsymm[15][3] = _sqrt2;		/*[110]180*/

	//<111>120*degrees, three-fold symmetries
	qsymm[16][0] = half;			qsymm[16][1] = half;			qsymm[16][2] = half;			qsymm[16][3] = half;		/*[111]120*/
	qsymm[17][0] = -1.0 * half;		qsymm[17][1] = half;			qsymm[17][2] = half;			qsymm[17][3] = half;		/*[111]240*/
	qsymm[18][0] = half;			qsymm[18][1] = half;			qsymm[18][2] = -1.0 * half;		qsymm[18][3] = half;		/*[1-11]240*/
	qsymm[19][0] = -1.0 * half; 	qsymm[19][1] = half;			qsymm[19][2] = -1.0 * half;		qsymm[19][3] = half;		/*[1-11]240*/
	qsymm[20][0] = half;			qsymm[20][1] = -1.0 * half;		qsymm[20][2] = half;			qsymm[20][3] = half;		/*[-111]120*/
	qsymm[21][0] = -1.0 * half;		qsymm[21][1] = -1.0 * half;		qsymm[21][2] = half;			qsymm[21][3] = half;		/*[-111]240*/ //mtex h h -h -h
	qsymm[22][0] = half;			qsymm[22][1] = -1.0 * half;		qsymm[22][2] = -1.0 * half;		qsymm[22][3] = half;		/*[-1-11]120*/ //Mtex332-h h h -h
	qsymm[23][0] = -1.0 * half;		qsymm[23][1] = -1.0 * half;		qsymm[23][2] = -1.0 * half;		qsymm[23][3] = half;		/*[-1-11]240*/
//cout << "All m-3m fcc crystal symmetry quaternion operators loaded successfully." << endl;

	double nd[3] = {0.0, 0.0, 1.0}; //z-direction normal vector

	//project to standard triangle by calculating symmetric variants, first ndp
	double triposp[SYMMETRIES_IN_CUBIC][2];

	for (int s = 0; s < SYMMETRIES_IN_CUBIC; s++) {

		myQuaternion qq = myQuaternion( qsymm[s][0], qsymm[s][1], qsymm[s][2], qsymm[s][3] );
		myQuaternion qtestqsym = (*this) * qq;

		//myQuaternion qm1 = (*this).Inverse(); //Inverse of quaternion this
		//myQuaternion r = qm1 * (*p); //Resulting quaternion, rotation of the two previous quaternions pq-1

		qtestqsym.Invert();

		double hbar[3];
       	//in accordance with Spieß2009 and Morawiec Pospiech 1989
		double q0 = qtestqsym.q0;
		double q1 = qtestqsym.q1;
		double q2 = qtestqsym.q2;
		double q3 = qtestqsym.q3;

       	//x,y,z vector3d v
       	hbar[0] = ((SQR(q0) + SQR(q1) - SQR(q2) - SQR(q3))	 * nd[0]) 	+ ( 2*(q1*q2 - q0*q3) 						* nd[1]) 	+ (2*(q1*q3 + q0*q2)	 * nd[2]);
       	hbar[1] = (2*(q1*q2 + q0*q3)						 * nd[0])	+ ((SQR(q0) - SQR(q1) + SQR(q2) - SQR(q3)) 	* nd[1]) 	+ (2*(q2*q3 - q0*q1) * nd[2]);
       	hbar[2] = (2*(q1*q3 - q0*q2)						 * nd[0])	+ (2*(q2*q3 + q0*q1) 						* nd[1]) 	+ (SQR(q0) - SQR(q1) - SQR(q2) + SQR(q3) * nd[2]);

       	//##MK::antipodal consistency with MTex3.3.2
		if ( hbar[2] < 0.0 ) {
			hbar[0] *= -1.0;
			hbar[1] *= -1.0;
			hbar[2] *= -1.0;
		}

		//now get the azimuth angle theta and assure no throw out of range exceptions by the acos function, but before range limiting
		if( hbar[2] > (1.0 - DOUBLE_ACCURACY) ) hbar[2]= 1.0;
		if( hbar[2] < (-1.0 + DOUBLE_ACCURACY) ) hbar[2] = -1.0;

		double theta = acos( hbar[2] ); //result is [0.0 <= theta <= pi]

		//the fact that theta goes to pi can cause tan singularities if theta --> pi/2 periodicities
		double tt = tan(theta/2);
		if ( tt > (1.0 - DOUBLE_ACCURACY) ) { tt = 1.0; }

		double psi = 0.0; //capture discontinuity of the atan2 at 0,0
		if ( fabs(hbar[0]) > DOUBLE_ACCURACY && fabs(hbar[1]) > DOUBLE_ACCURACY ) {
			psi = atan2( hbar[1], hbar[0] );
		}

		//psi is limited against pi
		double epps = 1e-10;
		if ( psi > ( _PI_ - epps ) ) { psi = _PI_; }
		if ( psi < ( (-1.0 * _PI_) + epps ) ) { psi = -1.0 * _PI_; }


		triposp[s][0] = tt * cos(psi);
		triposp[s][1] = tt * sin(psi);
	}

	//check which are in the standard triangle and closest to the center in the standard triangle x >= y > 0
	double minnorm = 10.0;
	double norm = 10.0;
	double xmin = FAIL_MYMATH_NUMERIC;
	double ymin = FAIL_MYMATH_NUMERIC;

	//double sstr = pow( 2, 0.5 ) - 1.0;
	//sstr = pow ( sstr + EPS_ENVIRONMENT, 2); //radius of the sst circle

	for	( int s = 0; s < SYMMETRIES_IN_CUBIC; s++) {
		double x2y2 = SQR(triposp[s][0]) + SQR(triposp[s][1]);
		norm = pow( x2y2, 0.5 );

		if ( norm <= minnorm && triposp[s][0] >= 0.0 && triposp[s][1] >= 0.0 && triposp[s][0] >= triposp[s][1] ) { //&& x2y2 <= sstr ) {
			xmin = triposp[s][0];
			ymin = triposp[s][1];
			minnorm = norm;
		}
	}

	this->x = xmin;
	this->y = ymin;
}


void myQuaternion::quat2ipfz( unsigned char *rgb )
{
	project2fundamentalregion_ipfz();

	double position[2] = {this->x, this->y };

//cout << "Position=" << position[0] << ";" << position[1] << endl;

	if ( position[0] == FAIL_MYMATH_NUMERIC || position[1] == FAIL_MYMATH_NUMERIC ) {
		//color in black, a color otherwise not utilized any mismatch and failure in the function!
		rgb[0] = 0;
		rgb[1] = 0;
		rgb[2] = 0;
		return;
	}

	//heuristic approach to catch numeric cases too close to the IPFZ coloring triangle vertices
	//"red"
	double xr = 0.0;							double yr = 0.0; //red
	double xg = pow(2.0, 0.5) - 1.0;			double yg = 0.0; //green
	double xb = (0.5 * pow(3.0, 0.5)) - 0.5;	double yb = xb; //blue

	double xo = position[0];					double yo = position[1];

	//heuristic catch to avoid running in on-the-vertex-location singularities
	/*guideline implementation of the scaling approach Molodov
	scaling suggestion from K. Molodov HKL and TSL implement different white
	points and different strategies with respect to how to rescale the RGB
	mixing, source unknown so heuristic decision, no problem: because
	an IPF colorcoding for the m-3m, 1 symmetry is possible purpose is to show
	qualitative location, coding necessarily is unpractical for sharp ipfz
	fiber textures in particular if close to the vertices, then a power-law
	rescale might exaggerate contrast the stronger the closer to the vertices
	but the less distinct for orientations which are located near the white point*/

	if ( fabs(xo-xr) < EPS_PROXIMITY_IPFZ && fabs(yo-yr) < EPS_PROXIMITY_IPFZ ) {
		rgb[0] = 255; //assign pure RED
		rgb[1] = 0;
		rgb[2] = 0;
		return;
	}
	if ( fabs(xo-xg) < EPS_PROXIMITY_IPFZ && fabs(yo-yg) < EPS_PROXIMITY_IPFZ ) {
		rgb[0] = 0;
		rgb[1] = 255; //assign pure GREEN
		rgb[2] = 0;
		return;
	}
	if ( fabs(xo-xb) < EPS_PROXIMITY_IPFZ && fabs(yo-yb) < EPS_PROXIMITY_IPFZ ) {
		rgb[0] = 0;
		rgb[1] = 0;
		rgb[2] = 255; //assign pure BLUE
		return;
	}

	//identify RGB color in stereographic triangle
	double k1g = (xg-xo)*(yr-yb) - (yg-yo)*(xr-xb);
	//k1g can only be come 0 if: then the lines are parallel or coincident  xg-xo and yg-yo are either both zero or both terms the same because or naturally yr-yb = xr-xb = const  != 0catch k1g <= DOUBLE_ACCURACY)
	//which based on the shape of the triangle can be safely assumed

	double xgbar = ( (xg*yo - yg*xo)*(xr-xb) - (xg-xo)*(xr*yb - yr*xb) ) / k1g;
	double ygbar = ( (xg*yo - yg*xo)*(yr-yb) - (yg-yo)*(xr*yb - yr*xb) ) / k1g;
	double absggbar = fabs(pow( (SQR(xg-xgbar) + SQR(yg-ygbar)) , 0.5));
	double absogbar = fabs(pow( (SQR(xo-xgbar) + SQR(yo-ygbar)) , 0.5));

	double k1b = (xb-xo)*(yr-yg) - (yb-yo)*(xr-xg);
	double xbbar =( (xb*yo - yb*xo)*(xr-xg) - (xb-xo)*(xr*yg - yr*xg) ) / k1b;
	double ybbar =( (xb*yo - yb*xo)*(yr-yg) - (yb-yo)*(xr*yg - yr*xg) ) / k1b;
	double absbbbar = fabs(pow( (SQR(xb-xbbar) + SQR(yb-ybbar)), 0.5));
	double absobbar = fabs(pow( (SQR(xo-xbbar) + SQR(yo-ybbar)) , 0.5));

	double k1r = pow( (yo/xo), 2.0 );
	//x0 == zero already excluded by R vertex proximity check and the coordiante system choice x>=y
	double xrbar = pow( (2+k1r), 0.5);
	xrbar = xrbar - 1;
	xrbar /= (k1r + 1.0);
	double yrbar = yo/xo * xrbar;
	double absrrbar = fabs(pow( (SQR(xr-xrbar) + SQR(yr-yrbar)), 0.5));
	double absorbar = fabs(pow( (SQR(xo-xrbar) + SQR(yo-yrbar)), 0.5));

	//IPF color stretch approach
	double rrggbb[3];
	rrggbb[0] = pow( (absorbar/absrrbar), IPF_COLOR_STRETCH_R );
	rrggbb[1] = pow( (absogbar/absggbar), IPF_COLOR_STRETCH_G );
	rrggbb[2] = pow( (absobbar/absbbbar), IPF_COLOR_STRETCH_B );

	double maxx = rrggbb[0];
	if ( rrggbb[1] > maxx ) maxx = rrggbb[1];
	if ( rrggbb[2] > maxx ) maxx = rrggbb[2];

	//K. Molodov got a better agrreement with the HKL colo though by anisotropically IPF_COLOR_STRETCHING via reverse engineering
	int rr = rrggbb[0] * (1.0 / maxx) * 255;
	int gg = rrggbb[1] * (1.0 / maxx) * 255;
	int bb = rrggbb[2] * (1.0 / maxx) * 255;

	//"RRGGBB=" << rr << ";" << gg << ";" << bb << endl;

	rgb[0] = rr;
	rgb[1] = gg;
	rgb[2] = bb;
}
