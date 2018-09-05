// YatoTypes.h : include Yato's usual type definition
//

#pragma once

// ROOT Library
#include "TVirtualFitter.h"
#include "TFitter.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TF2.h"
#include "TVector3.h"
#include "TMatrixD.h"

// YATO'S INCLUDE
#include "PhysicalSymbol.h"

class MuVoxel;
class MuEvent;
class MuImage;

// YATO'S TYPEDEF
typedef long _ID;

typedef double _Len;
typedef double _Angle;
typedef double _P;
typedef double _SDensity;

typedef TVector3 _Index3;
typedef TVector3 _Point3;
typedef TVector3 _Vector3;

typedef std::list<long> _IDList;
typedef std::map<_Len, _Point3*> _PMap;
typedef std::map<_ID, MuVoxel*> _VMap;
typedef std::map<_ID, MuEvent*> _EVMap;

// SYSTEM PARAMETERS
namespace SYS {
	// System Set-up
	static const int NUM_DETECTOR = 4;
	/* PSMTS / CRIPT
	static const _Len X_START = -50*cm, X_END = 50*cm;
	static const _Len Y_START = -50*cm, Y_END = 50*cm;
	static const _Len Z_POS[NUM_DETECTOR] = {-120*cm, -70*cm, 30*cm, 80*cm};
	*/
	//*/ MRPC CRIS
	static const _Len X_START = -20*cm, X_END = 20*cm;
	static const _Len Y_START = -20*cm, Y_END = 20*cm;
	static const _Len Z_POS[NUM_DETECTOR] = {-90*cm, -20*cm, 20*cm, 90*cm};	
	//*/
	/* CTX
	static const _Len X_START = 0*cm, X_END = 90*cm;
	static const _Len Y_START = 0*cm, Y_END = 25.8*cm;
	static const _Len Z_POS[NUM_DETECTOR] = {-175.9*cm, -105.2*cm, -70.7*cm, 0*cm};	
	*/
	// Tolerance & Default Parameter
	static const _Angle ANGLE_TOLERANCE = 1.*mrad;
	static const _Angle ANGLE_CUT = 200*mrad;
	static const _Len PATH_TOLERANCE = 1*mm;
	static const _P MOMENTUM_P0 = 3.0*GeV/c;
	static const _P MOMENTUM_MIN = 100.0*MeV/c;
	static const _Len ATOMIC_VOXEL = 1.0*mm;
	static const _Len PRE_VOXEL = 5*cm;
	// ML-EM Parameter
	static const long MAX_ITERATION = 20;
	static const double LOGP_TOLERANCE = 1.0e-6;
	// Radiation Length
	static const _Len L0_AIR = 3.039e+4*cm;
	static const _Len L0_U235 = 0.3166*cm;
	static const _Len L0_LEAD = 0.5612*cm;
	static const _Len L0_IRON = 1.757*cm; 
	static const _Len L0_CONCRETE = 11.55*cm;
	// Scattering Density (V.L.Highly)
	static const _SDensity LAMBDA_AIR = 22.09*mrad*mrad / L0_AIR;
};