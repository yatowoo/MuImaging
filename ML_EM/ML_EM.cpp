// ML_EM.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include "TCanvas.h"

#include "MuEvent.h"
#include "MuVoxel.h"
#include "MuImage.h"

//#define QUICK_CHECK
#define IMAGING
//#define POSITIONING
//#define QA		// Quality Assessment

// Generate Momentum Info with TOF
double MomentumK(double t[SYS::NUM_DETECTOR], MuEvent* ev){
	_Len L[3];
	L[0] = ev->L(0);
	L[1] = ev->LOF();
	L[2] = ev->L(2);
	double beta = 0.1*((t[2]-t[0])/(L[0]+L[1])
						+(t[3]-t[1])/(L[1]+L[2])
						+(t[3]-t[0])/(L[0]+L[1]+L[2]));
	if(beta<1)
		return -1;
	//return beta;
	return 1/sqrt(beta*beta-1);
}

int main(int argc, char* argv[])
{

	// Output file
	//ofstream fout("../Data/Image/Entry_EmptyNew_200mrad_5cm_10min.txt");
	ofstream fout("_tmp.txt");

	// ROOT System
	TApplication rootApp("rootApp",&argc, argv);
	// ROOT Tree & Histogram
	TFile *fileData = new TFile("../Data/Tracks/CRYTracks_Empty_New.root");
	TTree *pTree = (TTree*)fileData->Get("CRY");
	TCanvas *muC = new TCanvas("c_mu","Canvas for muon imaging",1200,600);
	
	TH1 *hGamma = new TH1D("Yc","Muon Momentum (pc/Eu)",100,0,50);
	TH1 *hYato = new TH1D("err","Momentum Estimation (dp/p)",50,-1,4); 
	TH1 *hAngle = new TH1D("s","Scattering Angle - Oxz",200,-200,200); 
	TH1 *hDis = new TH1D("d","Deflection - Oxz",100,-50,50);
	TH1 *hCov = new TH1D("cov","Angle x Deflection - Oxz",200,-500,1000);
	TH2 *hMul = new TH2D("mul","Joint Distribution - Oxz",200,-1000,1000,200,-200,200);

	muC->Divide(3,2);
	muC->cd(1);hGamma->Draw();
	muC->cd(2);hAngle->Draw();
	muC->cd(3);hCov->Draw();
	muC->cd(4);hYato->Draw();
	muC->cd(5);hDis->Draw();
	muC->cd(6);hMul->Draw();

	// Set Address for Tree Branch
	double x[SYS::NUM_DETECTOR], y[SYS::NUM_DETECTOR],t[SYS::NUM_DETECTOR];
	double ke;
	
	std::stringstream ssBrName;
	std::string strBrName;
	for( Int_t i = 0 ; i < SYS::NUM_DETECTOR ; i++){
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
		pTree->SetBranchAddress(tmpBrTName, t+i);
		
	}
	pTree->SetBranchAddress("KE", &ke);
	
	//
	// Read Entry from TTree & Add Event
#ifdef QUICK_CHECK
	for(int nRun = 0 ; nRun < 100 ;nRun ++){
#endif
	MuImage theImage;	
	_Point3 ptPos[SYS::NUM_DETECTOR];
	MuEvent* evTmp = NULL;
	Long_t nEff = 0, nTotal = 0;
	std::cout << "------> Event Read Start"<< std::endl;
#ifdef QUICK_CHECK
	for(Long_t i = nRun*1500 ; i < (nRun+1)*1500 && i< pTree->GetEntries() ; i++){
#endif
#ifndef QUICK_CHECK
	for(Long_t i = 0 ; i < 1000000 && i< pTree->GetEntries() ; i++){
#endif
		nTotal ++;
		pTree->GetEntry(i);
		for(Int_t j = 0 ; j < SYS::NUM_DETECTOR ; j++)
			ptPos[j].SetXYZ(x[j]*mm,y[j]*mm,SYS::Z_POS[j]);

		// Momentum
		_Vector3 muonIN = ptPos[1] - ptPos[0];
		static double Eu = 105.6587; // MeV
		double gamma = (ke*1000/Eu + 1);
		//fout<<gamma<<"\t"<<eta<<"\t"<<data[0][0]<<"\n";

		evTmp = new MuEvent(i,ptPos);
		if(evTmp->GetEventID() != -1 ){


			// Fill Histogram with Scattering Data in Oxz
			static TMatrixD data(2,1);
			evTmp->GetData(data);

			// Generate Momentum Info
			//hGamma->Fill(sqrt(gamma*gamma-1));
			if(data[0][0]>30*mrad){
				hGamma->Fill(sqrt(gamma*gamma-1));
				double Kp = MomentumK(t,evTmp);
				if(Kp > 0 && Kp < 30.)
					hYato->Fill(Kp/sqrt(gamma*gamma-1)-1);
				if(Kp < 30.)
					evTmp->SetYr(Kp*Kp);
			}

			//hYato->Fill(pow(gamma/(gamma*gamma-1),2)/(std::cos(data[0][0])));
			evTmp->GetDxz(data);

			//hGamma->Fill(gamma/(gamma*gamma-1)/sqrt(cos(data[0][0])));

			hAngle->Fill(data[0][0]/mrad);
			hDis->Fill(data[1][0]);
			hCov->Fill(data[0][0]/mrad*data[1][0]/mm);
			hMul->Fill(data[0][0]/mrad,data[1][0]);

#ifdef QA
			if(evTmp->HasPoCA()){
#else
			if(theImage.AddEvent(evTmp)){
#endif
				nEff++;
			}//addEvent
		}//no empty event

		// Dynamic Showing
		if(nTotal % 1000 == 0 ){
			muC->cd(1);hGamma->Draw();
			muC->cd(2);hAngle->Draw();
			muC->cd(3);hCov->Draw();
			muC->cd(4);hYato->Draw();
			muC->cd(5);hDis->Draw();
			muC->cd(6);hMul->Draw();
			muC->Modified();
			muC->Update();
			std::cout <<"Eff / Total : "<< nEff << " / " << nTotal <<std::endl;
		}
		if(gSystem->ProcessEvents())
			break;
	}//ForeachEntryOfTree
	// Fit & Draw
	muC->cd(1);hGamma->Draw();
	muC->cd(2);hAngle->Draw();
	muC->cd(3);hCov->Draw();
	muC->cd(4);hYato->Fit("gaus","","",-0.5,0.5);hYato->Draw();
	muC->cd(5);hDis->Draw();
	muC->cd(6);hMul->Draw();
	muC->Modified();
	muC->Update();
	muC->SaveAs("../Data/Utils/EmptyNew.root","");
	if(gSystem->ProcessEvents())
		return 1;
	std::cout <<"Eff / Total : "<< nEff << " / " << nTotal <<std::endl
				<<"<------ Event Read Over "<< std::endl;

#ifdef IMAGING
	// Imaging
	theImage.PreImaging(30);
	//theImage.MLEMIterate(10);
#endif
#ifdef QUICK_CHECK
	theImage.PreImaging(30);
	std::cout<<theImage.GetLambdaMax()/(mrad*mrad/cm)<<"\n";
	fout<<theImage.GetLambdaMax()/(mrad*mrad/cm)<<"\n";
#endif
#ifdef IMAGING
	//theImage.PrintSDTable(fout);
#endif
	theImage.PrintEntriesTable(fout);
	//
#ifdef QUICK_CHECK
}// end of nRun
#endif

	fileData->Close();
	fout.close();
	rootApp.Run();

	return 0;
}

