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
/// \file visualization/standalone/include/StandaloneVisAction.hh
/// \brief Definition of the StandaloneVisAction class
//
//

#ifndef STANDALONEVISACTION_HH
#define STANDALONEVISACTION_HH

#include<vector>

#include "G4SystemOfUnits.hh"
#include "G4VUserVisAction.hh"
#include "G4Polyline.hh"
#include "G4Circle.hh"

class StandaloneVisAction: public G4VUserVisAction {
    
  virtual void Draw();

  public:
    StandaloneVisAction(std::string inputPath, int eventNumber) {
      fInputPath = inputPath;
      fEventNumber = eventNumber;
    }

    void ReadEvent();

  private:
    std::string fInputPath;
    int fEventNumber;

    bool fEventRead = false;
    G4Polyline fTrack;
    std::vector<G4Circle> fHits;

    std::vector<std::string> chamberNames {"BARI-001", "BARI-002", "BARI-003", "BARI-004"};
  
    std::map<std::string, std::array<double, 3>> coordinateMap {
      {"BARI-001", std::array<double, 3>{0, 0, -(697+254+294)*mm}},
      {"BARI-002", std::array<double, 3>{0, 0, -(254+294)*mm}},
      {"BARI-003", std::array<double, 3>{0, 0, 170*mm}},
      {"BARI-004", std::array<double, 3>{0, 0, (170+697)*mm}}
    };
    std::map<std::string, std::array<double, 3>> sizeMap {
      {"BARI-001", std::array<double, 3>{10*cm, 10*cm, 1*cm}},
      {"BARI-002", std::array<double, 3>{10*cm, 10*cm, 1*cm}},
      {"BARI-003", std::array<double, 3>{10*cm, 10*cm, 1*cm}},
      {"BARI-004", std::array<double, 3>{10*cm, 10*cm, 1*cm}}
    };

    std::array<double, 3> ge21CoordinateMap {0., 0., 0.};
    std::array<double, 4> ge21SizeMap {488.92*mm, 643.43*mm, 431.5*mm, 1*cm};
};

#endif

