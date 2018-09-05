// Class - MuImage Implementation

#include "stdafx.h"
#include "MuImage.h"

MuImage::MuImage() :
	fNIteration(-1), fLnP(0.0),
	fVoxelSize(SYS::PRE_VOXEL){

	fPi = new _Point3(SYS::X_START,SYS::Y_START,SYS::Z_POS[1]);
	fPf = new _Point3(SYS::X_END,SYS::Y_END,SYS::Z_POS[2]);
	fMatError = new TMatrixD(2,2);
	(*fMatError)[0][0] = 14.*14.;
	(*fMatError)[0][1] = 20.;
	(*fMatError)[1][0] = 20.;
	(*fMatError)[1][1] = 10.*10;
	this->InitialVoxels();
}

MuImage::~MuImage(){
	delete fPi; fPi = NULL;
	delete fPf; fPf = NULL;
	delete fMatError; fMatError = NULL;
}

_ID MuImage::GetVoxelID(_Point3& pos,_Index3& idx){
	if(pos[0] < fPi->X() || pos[1] < fPi->Y() || pos[2] < fPi->Z())
		return -1;
	if(pos[0] >= fPf->X() || pos[1] >= fPf->Y() || pos[2] > fPf->Z())
		return -1;

	long idxX = (long)((pos.X() - fPi->X())/fVoxelSize);
	long idxY = (long)((pos.Y() - fPi->Y())/fVoxelSize);
	long idxZ = (long)((pos.Z() - fPi->Z())/fVoxelSize);
	idx.SetXYZ(idxX,idxY,idxZ);

	static long nBinX = (fPf->X() - fPi->X())/fVoxelSize;
	static long nBinY = (fPf->Y() - fPi->Y())/fVoxelSize;
	static long nBinZ = (fPf->Z() - fPi->Z())/fVoxelSize;

	if(idxZ == nBinZ)	idxZ -=1;	// exit out

	return (idxX + idxY * nBinX + idxZ * nBinX * nBinY);
}

_ID MuImage::GetVoxelID(_Point3& pi){
	_Index3 idx;
	return this->GetVoxelID(pi,idx);
}

bool MuImage::InitialVoxels(){
	MuVoxel* voxelTmp = NULL;
	_Point3 pi, pf;
	_ID vid = 0;

	for(_Len zTmp = fPi->Z(); zTmp < fPf->Z() ; zTmp += SYS::PRE_VOXEL){
		for(_Len yTmp = fPi->Y(); yTmp < fPf->Y() ; yTmp += SYS::PRE_VOXEL){
			for(_Len xTmp = fPi->X(); xTmp < fPf->X() ; xTmp += SYS::PRE_VOXEL){
				pi.SetXYZ(xTmp,yTmp,zTmp);
				pf.SetXYZ(xTmp+SYS::PRE_VOXEL,yTmp+SYS::PRE_VOXEL,zTmp+SYS::PRE_VOXEL);
				vid = this->GetVoxelID(pi);
				voxelTmp = new MuVoxel(vid,pi,pf);
				fMapVoxels[vid] = voxelTmp;
			}
		}
	}

	voxelTmp = NULL;
	return false;
}
//
// Scan-Line Method
//
// DO NOT insert pi into ptMap
void MuImage::CrossBorder(_Point3& pi, _Point3& pf, _PMap& ptMap){

	_Index3 idx_pi, idx_pf;
	GetVoxelID( pi, idx_pi);
	GetVoxelID( pf, idx_pf);

	_Vector3 lineVec = pf - pi;

	ptMap[pf.Z()] = new _Point3(pf);
	
	// Locate all the point which line crossing voxel boundary
	_ID idxNow = 0;
	_Point3* pNow = NULL;
	for(int idxNow = 1 + std::min(idx_pi.X(),idx_pf.X());
		idxNow <= std::max(idx_pi.X(),idx_pf.X());idxNow++){
		_Len xNow = fPi->X() + idxNow * fVoxelSize;
		pNow = new _Point3;
		*pNow = (pi + (xNow - pi.X())/lineVec.X() * lineVec);
		ptMap[pNow->Z()] = pNow;
	}
	for(int idxNow = 1 + std::min(idx_pi.Y(),idx_pf.Y());
				idxNow <= std::max(idx_pi.Y(),idx_pf.Y());idxNow++){
		_Len yNow = fPi->Y() + idxNow * fVoxelSize;
		pNow = new _Point3;
		*pNow = (pi + (yNow - pi.Y())/lineVec.Y()* lineVec);
		ptMap[pNow->Z()] = pNow;
	}
	for(int idxNow = 1 + std::min(idx_pi.Z(),idx_pf.Z());
				idxNow <= std::max(idx_pi.Z(),idx_pf.Z());idxNow++){
		_Len zNow = fPi->Z() + idxNow * fVoxelSize;
		pNow = new _Point3;
		*pNow = (pi + (zNow - pi.Z())/lineVec.Z()* lineVec);
		ptMap[pNow->Z()] = pNow;
	}
}

bool MuImage::AddEvent(MuEvent* ev){

	if(ev == NULL || ev->GetEventID() == -1)
		return false;
	// build incident vector
	_Point3 muIN;ev->GetHit(0,&muIN);
	_Point3 muOUT;ev->GetHit(1,&muOUT);
	_Vector3 vecIN = muOUT - muIN;
	// get in&out point
	ev->GetHit(1,&muIN);
	ev->GetHit(2,&muOUT);
	if(GetVoxelID(muIN) == -1 || GetVoxelID(muOUT) == -1)
		return false;

	//YYYYYYYYYYYYYYYYYYYY
	//
	//	Key Algorithm 1st - Positioning
	//		
	//	Detail : Calculate all the points that muon cross with borders of each voxel
	//
	//YYYYYYYYYYYYYYYYYYYY 
	_PMap ptsBorder;
	_Point3 pc,qc, poca;
	_ID idPoCA = -1;		// Locate PoCA & act a flag
	// Pattern Recognition
	if(ev->HasPoCA()){
		// Event with processible scatering angle

		ev->GetPoCA(&poca);
		ev->GetPc(&pc);
		ev->GetQc(&qc);

		idPoCA = this->GetVoxelID(poca);	
		if(idPoCA != -1){
			// PoCA in Volume
			// 1. PcQc in Volume
			if( GetVoxelID(pc) == -1 || GetVoxelID(qc) == -1){
				return false;
				//idPoCA = 0;
				//this->CrossBorder(muIN,poca,ptsBorder);
				//this->CrossBorder(poca, muOUT, ptsBorder);				
			}
			// 2. Pc or Qc out of Volume
			else{
				this->CrossBorder(muIN,pc,ptsBorder);
				this->CrossBorder(pc,qc,ptsBorder);
				this->CrossBorder(qc, muOUT, ptsBorder);			
			}
		}else{
			// 3. PoCA out of Volume
			return false;
			//idPoCA = 0;
			//this->CrossBorder(muIN, muOUT, ptsBorder);
		}
	}
	// scattering angle too small
	else
		this->CrossBorder(muIN, muOUT, ptsBorder);

	this->fMapEvents[ev->GetEventID()] = ev;
	

	//YYYYYYYYYYYYYYYYYYYY
	//
	//	Key Algorithm 2rd - Mark
	//		
	//	Detail : mark all the voxel the muon get through with weight matrix
	//		-->Weight Matrix : calculate path segment and build L-T
	//
	//YYYYYYYYYYYYYYYYYYYY 
	_Len pathL = 0.0, pathT = 0.0;
	_ID vidNow, vidNext;
	_PMap::iterator itr = ptsBorder.begin();
	_Point3 *ptLast = &muIN, *ptNow= itr->second, *ptNext=NULL;

	// calculate flight length in volume
	/*
	_Len LOF = 0.0;
	
	for(itr = ptsBorder.begin() ; itr != ptsBorder.end() ; itr++){
		ptNext = itr->second;
		_Vector3 segment = *ptLast - *ptNext;
		ptLast = ptNext;
		LOF += segment.Mag();
	}
	ev->SetLOF(LOF);
	*/
	// Start 
	itr = ptsBorder.begin();
	ptLast = &muIN, ptNow= itr->second, ptNext=NULL;
	for( itr++; itr != ptsBorder.end() ; itr++){
		ptNext = itr->second;

		pathL += (ptNow->Z() - ptLast->Z())/vecIN.Z() * vecIN.Mag();
		
		vidNow = GetVoxelID((*ptLast + *ptNow)*0.5);
		vidNext = GetVoxelID((*ptNow + *ptNext)*0.5);
		if( vidNow != vidNext){
			pathT = (fPf->Z() - ptNow->Z())/vecIN.Z() * vecIN.Mag();
			// Insert Matrix element into voxel
				// NEED UPGRADE
			fMapVoxels[vidNow]->AddEvent(ev->GetEventID(), pathL, pathT,vecIN,*ptLast,*ptNow);
			// Insert Voxel index into event
			//LOF += pathL;
			ev->AddVoxel(vidNow);
			//
			pathL = 0.0;
		}
		if(fPf->Z() - ptNext->Z() < SYS::PATH_TOLERANCE){
			// Volume Bottom
			pathL += (ptNext->Z() - ptNow->Z())/vecIN.Z() * vecIN.Mag();
			pathT = 0.0;
			// Insert Matrix element into voxel
				// NEED UPGRAD
			fMapVoxels[vidNext]->AddEvent(ev->GetEventID(), pathL, pathT,vecIN,*ptLast,*ptNow);
			// Insert Voxel index into event
			//LOF += pathL;
			ev->AddVoxel(vidNext);
		}
		ptLast = ptNow; 
		ptNow = ptNext;
	}
	
	// delete - all elem(new _Point3) in ptsBorder
	ptLast = NULL;
	ptNow = NULL;
	ptNext = NULL;
	for( itr = ptsBorder.begin(); itr != ptsBorder.end() ; itr++){
		ptNow = itr->second;
		ptsBorder[itr->first] = NULL;
		delete ptNow;
	}
	ptNow = NULL;

	//YYYYYYYYYYYYYYYYYYYY
	//
	//	Key Algorithm 3rd - Allocation
	//		
	//	Detail : fill each voxel that scattering happend in patterns
	//
	//YYYYYYYYYYYYYYYYYYYY
	if(idPoCA > 0)
		//fMapVoxels[idPoCA]->AddPoCA(ev->GetEventID(), ev->GetSAngle()/ev->GetYr());
		fMapVoxels[idPoCA]->AddPoCA(ev);
	/* NEED UPGRADE
	if(idPoCA == 0){
		_IDList* vids = ev->GetVoxelList();
		_Vector3 flight = muOUT - muIN;
		_Len Lxz = sqrt(flight.Mag2()-flight.Y()*flight.Y());
		_Len Lyz = sqrt(flight.Mag2()-flight.X()*flight.X());
		for(_IDList::iterator itr = vids->begin(); itr != vids->end() ; itr++){
			_ID evid = ev->GetEventID();
			fMapVoxels[*itr]->AddPoCA(evid, ev->GetSAngle(),Lxz,Lyz);
		}
	}
	*/

	return true;
}

// Calc Lambda /= EventCount
void MuImage::PreImaging(int entriesLimit){
	this->fNIteration = 0;
	_VMap::iterator itr;
	for(itr = fMapVoxels.begin() ; itr != fMapVoxels.end() ; itr ++)
		itr->second->CalcSDensity(entriesLimit);
}

void MuImage::PrintSDTable(ostream& outSD){

	_ID vID = 0;
	static long nBinX = (fPf->X() - fPi->X())/fVoxelSize;
	static long nBinY = (fPf->Y() - fPi->Y())/fVoxelSize;
	static long nBinZ = (fPf->Z() - fPi->Z())/fVoxelSize;

	std::cout<<" [-] Output Scatter Density Table.(Oxy)  ！ START!!!"<<std::endl;
	for(int iz = 0 ; iz < nBinZ ; iz ++){
		for(int iy = 0 ; iy < nBinY ; iy ++){
			for(int ix = 0 ; ix < nBinX ; ix ++){
				vID = ix + iy * nBinX + iz * nBinX * nBinY;
				outSD<<std::setw(6)<<(int)(fMapVoxels[vID]->GetSDensity()/(mrad*mrad/cm))<<"\t";
			}// for x index
			outSD<<"\n";
		}//for y index
	}//for z index
	std::cout<<" [-] Output Scatter Density Table.  ！ OVER!!!"<<std::endl;
}

void MuImage::PrintEntriesTable(ostream& outENT){

	_ID vID = 0;
	static long nBinX = (fPf->X() - fPi->X())/fVoxelSize;
	static long nBinY = (fPf->Y() - fPi->Y())/fVoxelSize;
	static long nBinZ = (fPf->Z() - fPi->Z())/fVoxelSize;

	std::cout<<" [-] Output PoCA Entry Table.(Oxy) ！ START!!!"<<std::endl;
		for(int iz = 0 ; iz < nBinZ ; iz ++){
	for(int iy = 0 ; iy < nBinY ; iy ++){
			for(int ix = 0 ; ix < nBinX ; ix ++){
				vID = ix + iy * nBinX + iz * nBinX * nBinY;
				outENT<<std::setw(6)<<fMapVoxels[vID]->GetWeights().size()<<"\t";
			}// for x index
			outENT<<"\n";
		}//for y index
	}//for z index
	std::cout<<" [-] Output PoCA Entry Table. ！ OVER!!!"<<std::endl;
}

// Calc Covariance Matrix for each event
	// Call all Voxels passed an get SD & Weight
void MuImage::GenerateCovar(){

	//TMatrixD sigma(2,2);
	TMatrixD Sxz(2,2);
	TMatrixD Syz(2,2);
	//TMatrixD weight(2,2);
	TMatrixD Wxz(2,2);
	TMatrixD Wyz(2,2);
	_EVMap::iterator evItr;
	_IDList* vids = NULL;
	MuVoxel* voxel = NULL;
	for(evItr = fMapEvents.begin() ; evItr != fMapEvents.end() ; evItr ++){
		vids = evItr->second->GetVoxelList();
		for(_IDList::iterator voItr = vids->begin() ; voItr != vids->end() ; voItr++){
			voxel = fMapVoxels[*voItr];
			
			//voxel->GetWeight(evItr->first,weight);
			//sigma += voxel->GetSDensity() * weight;

			voxel->GetWxz(evItr->first,Wxz);
			Sxz += voxel->GetSDensity() * Wxz;

			voxel->GetWyz(evItr->first,Wyz);
			Syz += voxel->GetSDensity() * Wyz;
		}
		//evItr->second->SetInverseCovar(sigma);
		Sxz = pow(evItr->second->GetYr(),2)*Sxz + *fMatError;
		evItr->second->SetCovXZ(Sxz);
		Syz = pow(evItr->second->GetYr(),2)*Syz + *fMatError;
		evItr->second->SetCovYZ(Syz);
	}
}

// Update SD Table by each voxel
void MuImage::UpdateVoxel(){

	_ID idEV = 0;
	//TMatrixD data(2,1);
	TMatrixD Dxz(2,1);
	TMatrixD Dyz(2,1);
	//TMatrixD dataT(1,2);
	TMatrixD DxzT(1,2);
	TMatrixD DyzT(1,2);
	//TMatrixD inSigma(2,2);
	TMatrixD Sxz(2,2);
	TMatrixD Syz(2,2);
	//TMatrixD weight(2,2);
	TMatrixD Wxz(2,2);
	TMatrixD Wyz(2,2);

	TMatrixD tmp;
	_SDensity lambda = 0.0;

	_VMap::iterator vItr;
	_WMap::iterator wItr;
	for(vItr = fMapVoxels.begin() ; vItr != fMapVoxels.end() ; vItr ++){
		lambda = vItr->second->GetSDensity();
		double sum = 0.0;

		// Proj Oxz
		_WMap& w = vItr->second->GetWxz();
		if(w.size() < 1) continue;
		for(wItr = w.begin() ; wItr != w.end() ; wItr ++){
			idEV = wItr->first;
			Wxz = *(wItr->second);
			this->fMapEvents[idEV]->GetDxz(Dxz);
			DxzT.Transpose(Dxz);
			this->fMapEvents[idEV]->GetCovXZ(Sxz);
			tmp.ResizeTo(1,1);
			tmp = (DxzT*Sxz * Wxz*Sxz*Dxz);
			sum += tmp[0][0]*pow(fMapEvents[idEV]->GetYr(),2);
			tmp.ResizeTo(2,2);
			tmp = Sxz * Wxz;
			sum -= (tmp[0][0] + tmp[1][1])*pow(fMapEvents[idEV]->GetYr(),2);
		}//for each event through this voxel

		//Proj Oyz
		w = vItr->second->GetWyz();
		if(w.size() < 1) continue;
		for(wItr = w.begin() ; wItr != w.end() ; wItr ++){
			idEV = wItr->first;
			Wyz = *(wItr->second);
			this->fMapEvents[idEV]->GetDyz(Dyz);
			DyzT.Transpose(Dyz);
			this->fMapEvents[idEV]->GetCovYZ(Syz);
			tmp.ResizeTo(1,1);
			tmp = (DyzT*Syz * Wyz*Syz*Dyz);
			sum += tmp[0][0]*pow(fMapEvents[idEV]->GetYr(),2);
			tmp.ResizeTo(2,2);
			tmp = Syz * Wyz;
			sum -= (tmp[0][0] + tmp[1][1])*pow(fMapEvents[idEV]->GetYr(),2);
		}//for each event through this voxel

		sum *= lambda/w.size();
		sum += 1;
		if(sum*lambda < SYS::LAMBDA_AIR)
			sum = SYS::LAMBDA_AIR / lambda;
		this->fLnP += w.size() * (std::log(lambda)+sum);
		vItr->second->SetSDensity(sum * lambda);
	}//for each voxel
}

void MuImage::MLEMIterate(int nItr){
	if(this->fNIteration < 0)
		this->PreImaging();
	for(int i = 0 ; i < nItr ; i++){
		std::cout << "------>Iterations " << i+1 <<" Start." <<std::endl;
		this->fLnP = 0.0;
		this->GenerateCovar();
		std::cout << " [-] Covariance Generated. " <<std::endl;
		this->UpdateVoxel();
		std::cout << " [-] AllVoxels Update. " <<std::endl;
		std::cout << "<------Iterations " << i <<",\t"
					<< "LogP : " << fLnP << std::endl;
		this->fNIteration ++;
	}
}

_SDensity MuImage::GetLambdaMax(){
	_SDensity lambdaMax = 0.,lambdaTmp;
	_VMap::iterator vItr;
	for(vItr = fMapVoxels.begin() ; vItr != fMapVoxels.end() ; vItr ++){
		lambdaTmp = vItr->second->GetSDensity();
		if(lambdaTmp > lambdaMax)
			lambdaMax = lambdaTmp;
	}
	return lambdaMax;
}