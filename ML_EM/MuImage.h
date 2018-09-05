/* Cosmic-Ray Muon Imaging Program
 *
 * Class : MuImage
 * 
 * Author: Yitao Wu   2016/04/01
 */

/*************************************************************************
 * Copyright (C) 2014-2016, Yitao Wu, wu.topgun@gmail.com                *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see GPLv3.									 *
 *************************************************************************/

#ifndef YATO_MuImage
#define YATO_MuImage

#include "stdafx.h"

#include "MuEvent.h"
#include "MuVoxel.h"

class MuImage {

private:
	_Point3* fPi;			//
	_Point3* fPf;			// Border of Imaging Volume
	TMatrixD* fMatError;		// System Error Matrix(2,2)
	_Len fVoxelSize;
	_VMap fMapVoxels;		// Map of Voxel in the image
							// Key by VoxelID
	_EVMap fMapEvents;		// Map of events in the image
							// Key by EventID
	int fNIteration;
	double fLnP;

public:
	MuImage();
	~MuImage();

	bool AddEvent(MuEvent*);
	void PreImaging(int entriesLimit = 3);
	void MLEMIterate(int NItr = 1);

	void PrintSDTable(ostream& outSD = std::cout);
	void PrintEntriesTable(ostream& outENT= std::cout);

	_SDensity GetLambdaMax();

private:
	_ID GetVoxelID(_Point3&);
	_ID GetVoxelID(_Point3&,_Index3&);
	bool InitialVoxels();
	void CrossBorder(_Point3&, _Point3&, _PMap&);
	// Parallel
	void GenerateCovar();
	void UpdateVoxel();


};

#endif