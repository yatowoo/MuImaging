/* Cosmic-Ray Muon Imaging Program
 *
 * Class : MuVoxel
 *
 * Author: Yitao Wu   2016/04/01
*/

/*************************************************************************
 * Copyright (C) 2014-2016, Yitao Wu, wu.topgun@gmail.com                *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see GPLv3.									 *
 *************************************************************************/

#ifndef YATO_MuVoxel
#define YATO_MuVoxel

#include "stdafx.h"

#include "MuEvent.h"

typedef std::map<_ID, TMatrixD*> _WMap;

class MuImage;

class MuVoxel {

private:
	_ID fVoxelID;	// Generate by MuImage
	_Point3* fPi;
	_Point3* fPf;	// Voxel Border
	_SDensity fLambda;
	_WMap fWeights;
	_WMap fWxz;
	_WMap fWyz;

public:
	MuVoxel();
	MuVoxel(_ID, _Point3&, _Point3&);
	virtual ~MuVoxel();

	_ID GetVID(){ return fVoxelID;}
	bool GetBorder(_Point3& pi, _Point3& pf);
	_SDensity GetSDensity(){ return fLambda;}
	void SetSDensity(_SDensity lambda){ fLambda = lambda;}
	bool GetWeight(_ID ev_id, TMatrixD& mW);
	bool GetWxz(_ID ev_id, TMatrixD& mW);
	bool GetWyz(_ID ev_id, TMatrixD& mW);
	_WMap& GetWeights(){return this->fWeights;}
	_WMap& GetWxz(){return this->fWxz;}
	_WMap& GetWyz(){return this->fWyz;}

	TMatrixD* GenWeightMatrix(_Len,_Len);
	bool AddEvent(_ID, _Len L, _Len T, _Vector3& vecIN, _Point3& ptIN, _Point3& ptOUT);
	bool AddPoCA(_ID, _Angle);
	bool AddPoCA(_ID, _SDensity,_Len,_Len);
	bool AddPoCA(MuEvent*);
	void CalcSDensity(int);
};

#endif