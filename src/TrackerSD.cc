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
//
/// \file B2/B2a/src/TrackerSD.cc
/// \brief Implementation of the B2::TrackerSD class

#include "TrackerSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4AnalysisManager.hh"

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackerSD::TrackerSD(const G4String& name,
                     const G4String& hitsCollectionName)
 : G4VSensitiveDetector(name)
{
  collectionName.insert(hitsCollectionName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackerSD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection

  fHitsCollection
    = new TrackerHitsCollection(SensitiveDetectorName, collectionName[0]);

  // Add this collection in hce

  G4int hcID
    = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection( hcID, fHitsCollection );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool TrackerSD::ProcessHits(G4Step* aStep,
                                     G4TouchableHistory*)
{
  // energy deposit
  G4double edep = aStep->GetTotalEnergyDeposit();

  if (edep==0.) return false;

  auto newHit = new TrackerHit();

  newHit->SetTrackID  (aStep->GetTrack()->GetTrackID());
  G4int number=-1;
  if(collectionName[0] == "DetCollection1") number=1;
  if(collectionName[0] == "DetCollection2") number=0;
  newHit->SetChamberNb(number);
  newHit->SetEdep(edep);
  newHit->SetPos (aStep->GetPostStepPoint()->GetPosition());

  fHitsCollection->insert( newHit );

  //newHit->Print();

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackerSD::EndOfEvent(G4HCofThisEvent*)
{
  auto analysisManager = G4AnalysisManager::Instance();
  
  
  
  G4int nofHits = fHitsCollection->entries();
  G4double deposition = 0;
  G4int chamberno = -1;
  G4int copynumber = -1;
  
  for( G4int jcount = 0; jcount < nofHits; jcount++)
  {
    deposition += (*fHitsCollection)[jcount]->GetEdep();
    G4cout << "deposition:" << deposition << G4endl;
    G4int chamberno = (*fHitsCollection)[jcount]->GetChamberNb(); //bu önceden newchamberno idi
    /*
    if(newchamberno != chamberno && chamberno >-1) 
    {
    	G4cout << jcount << " chamberno: " << chamberno << "->" << newchamberno << G4endl;
    chamberno = newchamberno;
    }
    */
    
    G4int newcopynumber = (*fHitsCollection)[jcount]->GetCopyNumber(); //bu önceden newcopynumber idi
    /*
    if(newcopynumber != copynumber && copynumber >-1) 
    {
    G4cout << jcount << " chamberno: " << chamberno << "->" << newchamberno << G4endl;
    copynumber = newcopynumber;
    }
    */
      
  }
  
  if(deposition != 0.)
  {
    analysisManager->FillNtupleDColumn(0, 0, deposition);
    analysisManager->FillNtupleDColumn(0, 1, nofHits);
    analysisManager->FillNtupleDColumn(0, 2, chamberno);
    analysisManager->FillNtupleDColumn(0, 3, copynumber);
    analysisManager->AddNtupleRow(0);
  }

  
  //G4ThreeVector position = (*fHitsCollection)->GetPos();
  //analysisManager->FillNtupleDColumn(0, 1, position.getX());
  //analysisManager->FillNtupleDColumn(0, 2, position.getY());
  //analysisManager->FillNtupleDColumn(0, 3, position.getZ());

  
  ///////////////////////////////
    if ( verboseLevel>1 ) {
  //G4int nofHits = fHitsCollection->entries();
     G4cout << G4endl
            << "-------->Hits Collection: in this event they are " << nofHits
            << " hits in the tracker chambers: " << G4endl;
     for ( G4int icount=0; icount<nofHits; icount++ ){ (*fHitsCollection)[icount]->Print();
     G4cout << (*fHitsCollection)[icount] << G4endl;}
  }
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}

