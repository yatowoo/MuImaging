// Class - MuEvent Implementation

#include "stdafx.h"
#include "MuEvent.h"

MuEvent::MuEvent() :
	fEventID(-1),fMomentum(0.0*MeV),fYr(1.),
	fSAngle(0.0*rad),fPoCA(NULL),fPc(NULL),fQc(NULL),
	fData(NULL),fDxz(NULL),fDyz(NULL),fCovar(NULL)
{
	// empty constructor
}

MuEvent::MuEvent(_ID ev_id, _Point3 hits[], 
	double ev_p, int nHits): 
		fEventID(ev_id), fMomentum(ev_p),fYr(1.),
		fSAngle(0.0*rad),fPoCA(NULL),fPc(NULL),fQc(NULL),
		fData(NULL),fDxz(NULL),fDyz(NULL),fCovar(NULL)
{	
	if(fEventID < 0 || fMomentum < SYS::MOMENTUM_MIN || hits == NULL
		|| nHits != SYS::NUM_DETECTOR )
		fEventID = -1;		// Invalid Event
	else{
		for(int idx = 0 ; idx < SYS::NUM_DETECTOR ; idx++)
			fHits[idx].SetXYZ(hits[idx].X(),hits[idx].Y(),hits[idx].Z());
		
		fData = new TMatrixD(2,1);
		fDxz = new TMatrixD(2,1);
		fDyz = new TMatrixD(2,1);
		this->CalcPoCA();

		fCovar = new TMatrixD(2,2);
		this->fCovXZ = new TMatrixD(2,2);
		this->fCovYZ = new TMatrixD(2,2);
	}

}

MuEvent::~MuEvent(){

	if(fPoCA != NULL){
		delete fPoCA; fPoCA = NULL;
	}
	if(fPc != NULL){
		delete fPc; fPc = NULL;
	}
	if(fQc != NULL){
		delete fQc; fQc = NULL;
	}
	if(fData != NULL){
		delete fData; fData = NULL;
	}
	if(fCovar != NULL){
		delete fCovar; fCovar = NULL;
	}
}

bool MuEvent::CalcPoCA(){
		
	if(fEventID == -1)
		return false;
	
	_Point3& p0 = fHits[0];
	_Point3& p1 = fHits[1];
	_Point3& q0 = fHits[2];
	_Point3& q1 = fHits[3];

	_Vector3 u = p1 - p0;
	_Vector3 v = q1 - q0;
	_Vector3 w = p1 - q0;

	// Length of Flight
	this->fLength[0] = u.Mag();
	this->fLength[1] = w.Mag();
	this->fLength[2] = v.Mag();

	// Generate Scattering Info
		// Scattering Angle
	fSAngle = u.Angle(v);
	(*fData)[0][0] = fSAngle;
		// Transport Deflection
	_Point3 exitRef(p1.X(),p1.Y(),q0.Z());

	_Point3 exitOrigin = p1 + u*((SYS::Z_POS[2] - SYS::Z_POS[1])/u.Z());
	exitOrigin = exitOrigin - exitRef;
	
	_Vector3 exitReal = q0 - exitRef;
	
	_Vector3 deflect = exitReal - exitOrigin;
	(*fData)[1][0] = deflect.Mag();
		// Projection Scattering
		// Oxz
	(*fDxz)[0][0] = std::atan(v.X()/v.Z())-std::atan(u.X()/u.Z());
	(*fDxz)[1][0] = deflect.X() ;//* std::cos((*fDxz)[0][0]) 
									/// std::cos(vecOutProj.Theta());
		// Oyz
	(*fDyz)[0][0] = std::atan(v.Y()/v.Z())-std::atan(u.Y()/u.Z());
	(*fDyz)[1][0] = deflect.Y() ;//* std::cos((*fDyz)[0][0]) 
									/// std::cos(vecOutProj.Theta());
	
	// PoCA Reconstrution
	if(std::fabs(fSAngle) > SYS::ANGLE_CUT){
		this->fEventID = -1;
		return false;
	}

	if( std::fabs(fSAngle) < SYS::ANGLE_TOLERANCE)
		return false;

	double a = u * u;
	double b = u * v;
	double c = v * v;
	double d = u * w;
	double e = v * w;
	
	_Point3 pc = p1 + ((b*e-c*d)  / (a*c-b*b))  * u;
	_Point3 qc = q0 + ((a*e-b*d)/(a*c-b*b)) * v;
	
	fPc = new _Point3(pc);
	fQc = new _Point3(qc);
	fPoCA = new _Point3(0.5*(pc+qc));

	return true;
}//CalcPoCA