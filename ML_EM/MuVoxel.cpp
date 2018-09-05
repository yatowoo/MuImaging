// Class - MuVoxel Implementation

#include "stdafx.h"
#include "MuVoxel.h"

MuVoxel::MuVoxel() : fVoxelID(-1)
{
	fPi = NULL;
	fPf = NULL;
	fLambda = SYS::LAMBDA_AIR;
}

MuVoxel::MuVoxel(_ID vid, _Point3 &pi, _Point3 &pf)
	: fVoxelID(vid)
{
	fPi = NULL;
	fPf = NULL;
	fLambda = SYS::LAMBDA_AIR;
	if( fVoxelID >= 0 &&
		pi[0] < pf[0] && pi[1] < pf[1] && pi[2] < pf[2]){
			fPi = new _Point3(pi[0],pi[1],pi[2]);
			fPf = new _Point3(pf[0],pf[1],pf[2]);
	}
	else
		fVoxelID = -1;
}

MuVoxel::~MuVoxel(){
	delete fPi;
	delete fPf;
}

bool MuVoxel::GetBorder(_Point3& pi, _Point3& pf){
	if(fPi == NULL || fPf == NULL)
		return false;
	pi.SetXYZ(fPi->X(),fPi->Y(),fPi->Z());
	pf.SetXYZ(fPf->X(),fPf->Y(),fPf->Z());
	return true;
}

TMatrixD* MuVoxel::GenWeightMatrix(_Len L, _Len T){
	TMatrixD* mTmp = new TMatrixD(2,2);

	(*mTmp)[0][0] = L;
	(*mTmp)[0][1] = L*L/2.0+ L*T;
	(*mTmp)[1][0] = L*L/2.0+ L*T;
	(*mTmp)[1][1] = L*L*L/3.0+L*L*T+L*T*T;

	return mTmp;
}

bool MuVoxel::AddEvent(_ID ev_id, _Len pathL, _Len pathT,
						_Vector3& vecIN, _Point3& ptIN, _Point3& ptOUT){
	if(fWeights[ev_id] != NULL)
		return false;
	
	//TMatrixD* mTmp = new TMatrixD(2,2);
	// Total
	/*
	(*mTmp)[0][0] = pathL;
	(*mTmp)[0][1] = pathL*pathL/2.0+ pathL*pathT;
	(*mTmp)[1][0] = pathL*pathL/2.0+ pathL*pathT;
	(*mTmp)[1][1] = std::pow(pathL,3.0)/3.0
						+pathL*pathL*pathT+pathL*pathT*pathT;
	*/
	this->fWeights[ev_id] = GenWeightMatrix(pathL,pathT);
	
	// Oxz
	//mTmp = new TMatrixD(2,2);
	_Vector3& vProj = vecIN.Unit();
	vProj.SetY(0.0);
	this->fWxz[ev_id] = GenWeightMatrix(pathL*vProj.Mag(),pathT*vProj.Mag());
	//Oyz
	vProj = vecIN.Unit();
	vProj.SetX(0.0);
	this->fWyz[ev_id] = GenWeightMatrix(pathL*vProj.Mag(),pathT*vProj.Mag());
	//mTmp = NULL;

	return true;
}
bool MuVoxel::AddPoCA(_ID ev_id, _Angle dTheta, _Len LOFxz, _Len LOFyz){
	if(fWeights[ev_id] == NULL)
		return false;
	/*
	_Len pathL = (*fWeights[ev_id])[0][0];
	if(pathL < SYS::PATH_TOLERANCE)
		pathL = SYS::PATH_TOLERANCE;

	this->fLambda += dTheta * dTheta / pathL ;
	*/
	_Len Lxz = (*fWxz[ev_id])[0][0];
	_Len Lyz = (*fWyz[ev_id])[0][0];
	// NEED UPGRADE
	
	if (Lxz < fPf->Z() - fPi->Z())
		Lxz = SYS::PRE_VOXEL;
	if (Lyz < fPf->Z() - fPi->Z())
		Lyz = SYS::PRE_VOXEL;

	this->fLambda += 0.5* dTheta * dTheta *(Lxz/(LOFxz*LOFxz)+Lyz/(LOFyz*LOFyz));

	return true;
}
bool MuVoxel::AddPoCA(_ID ev_id, _Angle dTheta){
	if(fWeights[ev_id] == NULL){
		std::cout<<" !!!WARNNING - Muon_"<<ev_id<<" did not get through voxel - "
			<< this->fVoxelID<<std::endl;
		return false;
	}
	/*
	_Len pathL = (*fWeights[ev_id])[0][0];
	if(pathL < SYS::PATH_TOLERANCE)
		pathL = SYS::PATH_TOLERANCE;

	this->fLambda += dTheta * dTheta / pathL ;
	*/
	_Len Lxz = (*fWxz[ev_id])[0][0];
	_Len Lyz = (*fWyz[ev_id])[0][0];
	// NEED UPGRADE
	
	if (Lxz < fPf->Z() - fPi->Z())
		Lxz = SYS::PRE_VOXEL;
	if (Lyz < fPf->Z() - fPi->Z())
		Lyz = SYS::PRE_VOXEL;

	this->fLambda += 0.5* dTheta * dTheta *(1./Lxz+1./Lyz);

	return true;
}

bool MuVoxel::AddPoCA(MuEvent* ev){
	int ev_id = ev->GetEventID();
	if(fWeights[ev_id] == NULL){
		std::cout<<" !!!WARNNING - Muon_"<<ev_id<<" did not get through voxel - "
			<< this->fVoxelID<<std::endl;
		return false;
	}
	/*
	_Len pathL = (*fWeights[ev_id])[0][0];
	if(pathL < SYS::PATH_TOLERANCE)
		pathL = SYS::PATH_TOLERANCE;

	this->fLambda += dTheta * dTheta / pathL ;
	*/
	TMatrixD data(2,1);
	ev->GetDxz(data);
	_Angle dThetaXZ = data[0][0]/ev->GetYr();
	ev->GetDyz(data);
	_Angle dThetaYZ = data[0][0]/ev->GetYr();

	_Len Lxz = (*fWxz[ev_id])[0][0];
	_Len Lyz = (*fWyz[ev_id])[0][0];

	// NEED UPGRADE
	if (Lxz < fPf->Z() - fPi->Z())
		Lxz = SYS::PRE_VOXEL;
	if (Lyz < fPf->Z() - fPi->Z())
		Lyz = SYS::PRE_VOXEL;

	this->fLambda += 0.5* (dThetaXZ*dThetaXZ/Lxz+dThetaYZ*dThetaYZ/Lyz);

	return true;
}

void MuVoxel::CalcSDensity(int entriesLimit){
	if(!fWeights.empty() && fWeights.size() > entriesLimit)
		fLambda /= fWeights.size();
	else
		fLambda = SYS::LAMBDA_AIR;
}

bool MuVoxel::GetWeight(_ID ev_id, TMatrixD& mW){
	if(fWeights[ev_id]){
		mW = (*(fWeights[ev_id]));
		return true;
	}
	else{
		std::cout<<" [+] Error - Weight[" 
					<< ev_id <<"]["<<this->fVoxelID<<"]"
					<< " Not Exist!" << std::endl;
		return false;
	}
}
bool MuVoxel::GetWxz(_ID ev_id, TMatrixD& mW){
	if(fWxz[ev_id]){
		mW = (*(fWxz[ev_id]));
		return true;
	}
	else{
		std::cout<<" [+] Error - Weight[" 
					<< ev_id <<"]["<<this->fVoxelID<<"]"
					<< " Not Exist!" << std::endl;
		return false;
	}
}
bool MuVoxel::GetWyz(_ID ev_id, TMatrixD& mW){
	if(fWyz[ev_id]){
		mW = (*(fWyz[ev_id]));
		return true;
	}
	else{
		std::cout<<" [+] Error - Weight[" 
					<< ev_id <<"]["<<this->fVoxelID<<"]"
					<< " Not Exist!" << std::endl;
		return false;
	}
}