// SPoCA.cpp : Code for ROOT
//

#include "stdafx.h"
#include "TApplication.h"
#include "TSystem.h"

#include<iostream>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<string>
#include<cmath>
#include<map>

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TH3.h"
#include "TH2.h"

#define MAX(x,y) (x>y?x:y)
#define MIN(x,y) (x<y?x:y)

using namespace std;

const static double PI = 3.14159265358;

const Float_t ANGLE_LIMITED = 0.001;	// rad
const Int_t BIN_X = 30, BIN_Y = 10, BIN_Z = 10;
const Float_t X_START = -20., X_END= 20.0;	// cm
const Float_t Y_START = -20., Y_END= 20.;	// cm
const Float_t Z_START = -20., Z_END= 20.;	// cm
const Float_t Z_POS[4] = {-90., -20., 20.,90.};

typedef TVector3 Index3;
typedef TVector3 Point3;
typedef TVector3 Vector3;

// Muon Tomography Image Struct
struct MuImage {
	Float_t xi;
	Float_t xf;
	Float_t yi;
	Float_t yf;
	Float_t zi;
	Float_t zf;
	Int_t binX;
	Int_t binY;
	Int_t binZ;
	Long_t size;
	Float_t* sDenstiy;
	Long_t* countVoxel;
};

// Locate Point in Volume & Get the voxel id/index
Long_t GetVoxelID( Index3& idx, Point3& pNow, MuImage& image){
	// Check
	if (pNow.X() < image.xi || pNow.X() > image.xf )
		return -1;
	if (pNow.Y() < image.yi || pNow.Y() > image.yf )
		return -1;
	if (pNow.Z() < image.zi || pNow.Z() > image.zf )
		return -1;

	Long_t vId;
	Int_t idxX, idxY, idxZ;

	idxX = (pNow.X() - image.xi) * image.binX / (image.xf - image.xi);
	idxY = (pNow.Y() - image.yi) * image.binY / (image.yf - image.yi);
	idxZ = (pNow.Z() - image.zi) * image.binZ / (image.zf - image.zi);
	if(pNow.Z() == image.zf)
		idxZ--;

	idx.SetXYZ(idxX, idxY, idxZ);
	vId = idxX + idxY * image.binX + idxZ * image.binX * image.binY;
	return vId;
}
Long_t GetVoxelID(Point3& pNow, MuImage& image){
	Index3 idxTmp;
	return GetVoxelID(idxTmp, pNow, image);
}

// Calculate the path muon through each voxel in an event
// by Scan-Line xyz
Bool_t CalcPath(Point3& pi, Point3& pf, Float_t* ev_path, MuImage& image){
	
	// Avoid the boundary
	if(pf.Z() == image.zf)
		pf.SetZ(pf.Z() - 0.001);

	Index3 idx_pi, idx_pf;
	Long_t voxelID = 0;
	voxelID = GetVoxelID( idx_pi, pi, image);
	if(voxelID < 0 )
		return false;
	voxelID =GetVoxelID( idx_pf, pf, image);
	if(voxelID < 0 )
		return false;

	TVector3 lineVec = pf - pi;
	lineVec = lineVec.Unit();

	Float_t dX = (image.xf - image.xi)/image.binX ;
	Float_t dY = (image.yf - image.yi)/image.binY ;
	Float_t dZ = (image.zf - image.zi)/image.binZ ;

	map<Float_t, Point3*> mapNode;	// key - Point3.Z()
	mapNode[pi.Z()] = &pi;
	mapNode[pf.Z()] = &pf;
	
	// Locate all the point which line crossing voxel boundary
	Int_t idxNow = 0;
	Point3* pNow = NULL;
	for(Int_t idxNow = 1 + MIN(idx_pi.X(),idx_pf.X());
				idxNow <= MAX(idx_pi.X(),idx_pf.X());idxNow++){
		Float_t xNow = image.xi + idxNow * dX;
		pNow = new Point3;
		*pNow = (pi + (xNow - pi.X())/lineVec.X() * lineVec);
		mapNode[pNow->Z()] = pNow;
	}
	for(Int_t idxNow = 1 + MIN(idx_pi.Y(),idx_pf.Y());
				idxNow <= MAX(idx_pi.Y(),idx_pf.Y());idxNow++){
		Float_t yNow = image.yi + idxNow * dY;
		pNow = new Point3;
		*pNow = (pi + (yNow - pi.Y())/lineVec.Y()* lineVec);
		mapNode[pNow->Z()] = pNow;
	}
	for(Int_t idxNow = 1 + MIN(idx_pi.Z(),idx_pf.Z());
				idxNow <= MAX(idx_pi.Z(),idx_pf.Z());idxNow++){
		Float_t zNow = image.zi + idxNow * dZ;
		pNow = new Point3;
		*pNow = (pi + (zNow - pi.Z())/lineVec.Z()* lineVec);
		mapNode[pNow->Z()] = pNow;
	}

	// Calc Path Segement
	map<Float_t, Point3*>::iterator itr;

	Point3* pNext = NULL;
	for(itr = mapNode.begin() ; itr != mapNode.end() && pNext != &pf ; itr++){
		if(itr == mapNode.begin()){
			pNow = &pi;
			continue;
		}

		pNext = itr->second;
		voxelID = GetVoxelID(((*pNow)+(*pNext))*0.5, image);
		if(voxelID < 0 )
			return false;
		ev_path[voxelID] += (pNext->Z()-pNow->Z()) * lineVec.Z();

		pNow = pNext;		
		
	}

	return true;
}

// Reconstruct the closeset approching point of scattering
Double_t CalcPoCA(TVector3& pPoCA, TVector3* vPos, MuImage& image){
	TVector3& p0 = vPos[0];
	TVector3& p1 = vPos[1];	
	TVector3& q0 = vPos[2];
	TVector3& q1 = vPos[3];
	
	TVector3 u = p1 - p0;
	TVector3 v = q1 - q0;
	TVector3 w = p0 - q0;

	static Float_t* path = NULL;
	if(path == NULL)
		path = (Float_t*)malloc(image.size * sizeof(Float_t));
	for(Int_t i = 0 ; i < image.size ; i ++)
		path[i] = 0.0;	
	
	Float_t angle = u.Angle(v);
	angle = std::atan(v.X()/v.Z())-std::atan(u.X()/u.Z());
	if( fabs(angle) > ANGLE_LIMITED){

		Float_t a = u * u;
		Float_t b = u * v;
		Float_t c = v * v;
		Float_t d = u * w;
		Float_t e = v * w;
		
		TVector3 pc = p0 + ((b*e-c*d)  / (a*c-b*b))  * u;
		TVector3 qc = q0 + ((a*e-b*d)/(a*c-b*b)) * v;
	
		pPoCA.SetX( pc.x() / 2.0 + qc.x() / 2.0);
		pPoCA.SetY( pc.y() / 2.0 + qc.y() / 2.0);
		pPoCA.SetZ ( pc.z() / 2.0 + qc.z() / 2.0);
	
		Long_t voxelID = GetVoxelID(pPoCA,image);
		if(voxelID < 0 || GetVoxelID(pc,image) < 0  || GetVoxelID(qc,image) < 0)
			return -angle;

		CalcPath(p1, pc, path, image);
		CalcPath(pc, qc, path, image);
		CalcPath(qc, q0, path, image);

		Float_t dL = (Z_END - Z_START)/BIN_Z;//path[voxelID];
		if(dL < 0.1)
			dL = 0.1;
		image.sDenstiy[voxelID] += pow(10.0,6.0)*angle * angle / dL;
	}
	else
		CalcPath(p1, q0, path, image);

	for(Int_t i = 0 ; i < image.size ; i++)
		if(path[i] > 0.1)
			image.countVoxel[i]++;

	return angle;
}

int main(int argc, char* argv[])
{

	MuImage theImage = {X_START, X_END, Y_START, Y_END, Z_START, Z_END, \
							BIN_X, BIN_Y, BIN_Z};
	theImage.size = theImage.binX * theImage.binY * theImage.binZ ;
	theImage.countVoxel = (Long_t*)malloc(theImage.size * sizeof(Long_t));
	theImage.sDenstiy = (Float_t*)malloc(theImage.size * sizeof(Float_t));
	for(Int_t i = 0 ; i < theImage.size ; i ++){
		theImage.countVoxel[i] = 0;
		theImage.sDenstiy[i] = 0.0;
	}
	
	TApplication rootApp("rootApp",&argc, argv);
	
	TFile *fileData = new TFile("../Data/CRYTracks_RealEmpty.root");
	TTree *pTree = (TTree*)fileData->Get("CRY");
	/*
	ifstream fin("E:\\Users\\torre\\OneDrive\\CODE\\ROOT\\DATA\\MRPC_Record.dat");
	if(!fin.is_open()){
		cout<<"!!! File Not Found !!!" <<endl;
		return 1;
	}*/
	ofstream fout("../Data/CRIS_Empty_SPoCA_SDTable_Oxz.txt");

	TCanvas* mu_canvas = new TCanvas("mu_c","Canvas for CRIS",800,600);
	TH1* hAng = new TH1D("angle","Scattering Angle of Muon",1500,-50*PI/2,50*PI/2);
	hAng->Draw();

	// Set Address for Tree Branch
	Double_t x[4], y[4], t[4];
	Double_t angle, ke;
	
	stringstream ssBrName;
	string strBrName;
	for( Int_t i = 0 ; i < 4 ; i++){
		ssBrName<<"X"<<i;
		ssBrName>>strBrName;
		ssBrName.clear();
		const char* tmpBrXName = strBrName.c_str();
		pTree->SetBranchAddress(tmpBrXName, x+i);
		
		ssBrName<<"Y"<<i;
		ssBrName>>strBrName;
		ssBrName.clear();
		const char* tmpBrYName = strBrName.c_str();
		pTree->SetBranchAddress(tmpBrYName, y+i);

		ssBrName<<"T"<<i;
		ssBrName>>strBrName;
		ssBrName.clear();
		const char* tmpBrTName = strBrName.c_str();
		pTree->SetBranchAddress(tmpBrYName, t+i);
		
	}
	pTree->SetBranchAddress("KE", &ke);
	

	Point3 pPos[4];
	Point3 pPoCA;
	Long_t countPoCA = 0, countTotal = 0;

	for(Long_t i = 0 ; i < 150000 && i < pTree->GetEntries(); i++){
		countTotal ++;

		pTree->GetEntry(i);
		//for(int j = 0 ; j < 4 ; j++)
			//fin>>x[j]>>y[j]>>ke;

		for(Int_t j = 0 ; j < 4 ; j++){

			pPos[j].SetXYZ(x[j]/10.0,y[j]/10.0,Z_POS[j]);

			
		}

		angle = CalcPoCA(pPoCA, pPos, theImage);
		hAng->Fill((1000*angle));
		//hAng->Fill(x[2]);

		if(angle > 0)
			countPoCA++;

		if(countTotal % 1000 == 0){

			mu_canvas->Modified();
			mu_canvas->Update();

			cout<<"PoCA Event / Total Entries : "
				<<countPoCA<<" / "<<countTotal
				<<endl;
		}
		if( gSystem->ProcessEvents())
			return -1;
	}

	cout<<"PoCA Event / Total Entries : "<<countPoCA<<" / "<<countTotal<<endl;
	cout<<"SAngle FWHM : "<< hAng->GetBinCenter(hAng->FindLastBinAbove(hAng->GetMaximum()/2.))
		-
		hAng->GetBinCenter(hAng->FindFirstBinAbove(hAng->GetMaximum()/2.))
		<<endl;

	for(Long_t i = 0 ; i < theImage.size ; i ++){
		if(theImage.countVoxel[i] > 0)
			theImage.sDenstiy[i] /= theImage.countVoxel[i];
		else
			theImage.sDenstiy[i] = 0.0;
	}
	Long_t voxelID;
	//Int_t idxY = 0;
	for(Int_t idxZ = 0 ; idxZ < theImage.binZ ; idxZ++){
		for(Int_t idxY = 0 ; idxY < theImage.binY ; idxY++){
			for(Int_t idxX = 0 ; idxX < theImage.binX ; idxX++){
				voxelID = idxX + idxY * theImage.binX + idxZ * theImage.binX * theImage.binY;
				fout<<setw(6)<<(int)(theImage.sDenstiy[voxelID])<<"\t";
				//fout<<setw(6)<<(int)theImage.countVoxel[voxelID]<<"\t";
			}
			fout<<endl;
		}//forY
	}//forZ

	//fileData->Close();
	//fin.close();
	fout.close();
	std::cout<<"ALL Over!"<<std::endl;

	rootApp.Run();
	return 0;
}

