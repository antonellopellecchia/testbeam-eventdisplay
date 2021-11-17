//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file visualization/standalone/src/StandaloneVisAction.cc
/// \brief Implementation of the StandaloneVisAction class
//
//

#include <map>
#include <array>
#include <iostream>
#include <math.h>

#include "StandaloneVisAction.hh"

#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"
#include "G4Polyhedron.hh"
#include "G4Box.hh"
#include "G4Trap.hh"
#include "G4SubtractionSolid.hh"
#include "G4SystemOfUnits.hh"
#include "G4Polyline.hh"
#include "G4Point3D.hh"
#include "G4Circle.hh"

#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TGraphErrors.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void StandaloneVisAction::ReadEvent() {
  std::cout << "Reading event from ROOT file..." << std::endl;
  if (fEventRead) {
    std::cout << "Event already read." << std::endl;
    return;
  }

  TFile rechitFile(fInputPath.c_str(), "READ");
  TTree *rechitTree = (TTree *) rechitFile.Get("rechitTree");

  // rechits 2d variables
  int nrechits2d;
  std::vector<int> *vecRechit2DChamber = new std::vector<int>();
  std::vector<double> *vecRechit2D_X_Center = new std::vector<double>();
  std::vector<double> *vecRechit2D_Y_Center = new std::vector<double>();
  std::vector<double> *vecRechit2D_X_Error = new std::vector<double>();
  std::vector<double> *vecRechit2D_Y_Error = new std::vector<double>();
  std::vector<double> *vecRechit2D_X_ClusterSize = new std::vector<double>();
  std::vector<double> *vecRechit2D_Y_ClusterSize = new std::vector<double>();
  
  rechitTree->SetBranchAddress("nrechits2d", &nrechits2d);
  rechitTree->SetBranchAddress("rechit2DChamber", &vecRechit2DChamber);
  rechitTree->SetBranchAddress("rechit2D_X_center", &vecRechit2D_X_Center);
  rechitTree->SetBranchAddress("rechit2D_Y_center", &vecRechit2D_Y_Center);
  rechitTree->SetBranchAddress("rechit2D_X_error", &vecRechit2D_X_Error);
  rechitTree->SetBranchAddress("rechit2D_Y_error", &vecRechit2D_Y_Error);
  rechitTree->SetBranchAddress("rechit2D_X_clusterSize", &vecRechit2D_X_ClusterSize);
  rechitTree->SetBranchAddress("rechit2D_Y_clusterSize", &vecRechit2D_Y_ClusterSize);

  rechitTree->GetEntry(fEventNumber);

  // fit track
  int zStart = -(697+254+294+200)*mm;
  int zEnd = (170+697+200)*mm;
  TF1 trackX("fTrackX", "[0]+[1]*x", zStart, zEnd);
  TF1 trackY("fTrackY", "[0]+[1]*x", zStart, zEnd);
  trackX.SetParameters(0., 0.);
  trackY.SetParameters(0., 0.);

  TGraphErrors graphX;
  TGraphErrors graphY;
  graphX.SetName("gTrackX");
  graphY.SetName("gTrackY");

  for (int irechit=0; irechit<nrechits2d; irechit++) {
    // add to graph for fitting
    int chamber = vecRechit2DChamber->at(irechit);
    graphX.SetPoint(irechit, coordinateMap[chamberNames[chamber]][2], vecRechit2D_X_Center->at(irechit));
    graphX.SetPointError(irechit, 10, vecRechit2D_X_Error->at(irechit));
    graphY.SetPoint(irechit, coordinateMap[chamberNames[chamber]][2], vecRechit2D_Y_Center->at(irechit));
    graphY.SetPointError(irechit, 10, vecRechit2D_Y_Error->at(irechit));

    // create hit markers
    fHits.push_back(G4Circle(G4Point3D(
      vecRechit2D_X_Center->at(irechit)*mm,
      vecRechit2D_Y_Center->at(irechit)*mm,
      coordinateMap[chamberNames[chamber]][2]
    )));

    std::cout << 
    chamberNames[chamber] << " " <<
    vecRechit2D_X_Center->at(irechit)*mm << " " << 
    vecRechit2D_Y_Center->at(irechit)*mm << " " << 
    coordinateMap[chamberNames[chamber]][2] << std::endl;
  }

  // fit track and propagate to external points
  graphX.Fit(&trackX, "SQ");
  graphY.Fit(&trackY, "SQ");

  fTrack.push_back(G4Point3D(trackX.Eval(zStart), trackY.Eval(zStart), zStart));
  fTrack.push_back(G4Point3D(trackX.Eval(zEnd), trackY.Eval(zEnd), zEnd));

  std::cout << "Event read." << std::endl;
  fEventRead = true;
}

void StandaloneVisAction::Draw() {
  G4VVisManager* pVisManager = G4VVisManager::GetConcreteInstance();
  if (pVisManager) {

    this->ReadEvent();

    for (auto chamber:coordinateMap) {
      const std::string chamberName = chamber.first;
      const std::array<double, 3> chamberCoordinates = chamber.second;
      pVisManager->Draw(
        G4Box(
          chamberName,
          sizeMap[chamberName][0], sizeMap[chamberName][1], sizeMap[chamberName][2]
        ),
        G4VisAttributes(G4Colour(1,1,0)),
        G4Translate3D(chamberCoordinates[0], chamberCoordinates[1], chamberCoordinates[2])
      );
    }

    //G4Polyline track;
    G4Colour red(1.0, 0.0, 0.0);
    G4VisAttributes trackAttributes(red);
    fTrack.SetVisAttributes(trackAttributes);
    pVisManager->Draw(fTrack);

    G4VisAttributes hitAttributes(G4Colour(0.,0.,1.));
    for (G4Circle hit:fHits) {
      hit.SetVisAttributes(hitAttributes);
      pVisManager->Draw(hit);
    }

    //G4Trap ge21Chamber("ge21Chamber", ge21SizeMap[3], ge21SizeMap[2], ge21SizeMap[1], ge21SizeMap[0]);
    G4Trap ge21Chamber("ge21Chamber",
      0.5*ge21SizeMap[0], 0.5*ge21SizeMap[1], 0.5*ge21SizeMap[3], 0.5*ge21SizeMap[3],
      0.5*ge21SizeMap[2]
    );
    //G4RotationMatrix xRotation;
    //xRotation.rotateX(-M_PI_2);
    pVisManager->Draw(
      ge21Chamber,
      G4VisAttributes(G4Colour(0,1,0)),
      G4RotateX3D(-M_PI_2)
      //G4Translate3D(ge21CoordinateMap[0], ge21CoordinateMap[1], ge21CoordinateMap[2]),
    );
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
