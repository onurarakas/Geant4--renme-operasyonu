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
#include "G4RootAnalysisManager.hh"


//#include "TStyle.h"
//#include "TCanvas.h"
//#include "TH1F.h"
//#include "TFile.h"
//#include "TLegend.h"
//#include "TGraph.h"

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

  /*TCanvas* c1 = new TCanvas("c1", "Sensitive Detectors", 800, 600);
   
  c1->Divide(2, 1);
   
  h1->SetFillColorAlpha(kRed, 0.04);
  c1->cd(1);
  h2->Draw();
  c1->cd(2);
  h1->Draw();
  */
  
  // energy deposit
  G4double edep = aStep->GetTotalEnergyDeposit();

  if (edep==0.) return false;

  auto newHit = new TrackerHit();

  newHit->SetTrackID  (aStep->GetTrack()->GetTrackID());
  G4int number=-1;
  G4int copyno=-1;
  if(collectionName[0] == "DetCollection1") number=1;
  if(collectionName[0] == "DetCollection2") number=0;
  
  auto analysisManager = G4AnalysisManager::Instance();
    
  analysisManager->CreateH1("1","SD1",100,0,100);
  
  analysisManager->CreateH1("2","SD2",100,0,100);
  
  copyno = aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber();

  newHit->SetChamberNb(number);
  newHit->SetEdep(edep);
  newHit->SetCopyNumber(copyno);
  newHit->SetPos (aStep->GetPostStepPoint()->GetPosition());

  fHitsCollection->insert( newHit );

  //newHit->Print();

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackerSD::EndOfEvent(G4HCofThisEvent*)
{
  auto analysisManager = G4AnalysisManager::Instance();
  
  
  // Create a canvas to display the histograms
  TCanvas* canvas = new TCanvas("canvas", "Energy Deposition", 800, 600);

  // Divide the canvas into two parts for the histograms
  canvas->Divide(1, 2);
  
  TH1F* hist1 = new TH1F("hist1", "Energy Deposition - SD1", 100, 0, 100);
  TH1F* hist2 = new TH1F("hist2", "Energy Deposition - SD2", 100, 0, 100);
  

  G4int nofHits = fHitsCollection->entries();
  G4double deposition = 0;
  G4double edepo = 0;
  G4int chamberno = -1;
  G4int copynumber = -1;
  
  for( G4int jcount = 0; jcount < nofHits; jcount++)
  {
    deposition += (*fHitsCollection)[jcount]->GetEdep();
    G4cout << "deposition:" << deposition << G4endl;
    G4int newchamberno = (*fHitsCollection)[jcount]->GetChamberNb(); //bu önceden newchamberno idi
    
    if(newchamberno != chamberno && chamberno >-1) 
    {
    	G4cout << jcount << " chamberno: " << chamberno << "->" << newchamberno << G4endl;
    
    }
    chamberno = newchamberno;
    
    G4int newcopynumber = (*fHitsCollection)[jcount]->GetCopyNumber(); //bu önceden newcopynumber idi
    
    if(newcopynumber != copynumber && copynumber >-1) 
    {
    G4cout << jcount << " chamberno: " << chamberno << "->" << newchamberno << G4endl;

    }
    copynumber = newcopynumber;
      
  }
  
  if(deposition != 0.)
  {
    analysisManager->FillNtupleDColumn(0, 0, deposition);
    analysisManager->FillNtupleDColumn(0, 1, nofHits);
    analysisManager->FillNtupleDColumn(0, 2, chamberno);
    analysisManager->FillNtupleDColumn(0, 3, copynumber);
    analysisManager->AddNtupleRow(0);


    if(chamberno == 1)
    {
      //analysisManager->FillH1(1, deposition);
      
      hist1->Fill(deposition);
      
    }
    
    else if(chamberno == 2)
    {
      hist2->Fill(deposition);
    }
    
  }
  
  
  // Set the histogram titles and axis labels
  hist1->SetTitle("Sensitive Detector 1");
  hist1->GetXaxis()->SetTitle("Energy Deposition");
  hist1->GetYaxis()->SetTitle("Entries");

  hist2->SetTitle("Sensitive Detector 2");
  hist2->GetXaxis()->SetTitle("Energy Deposition");
  hist2->GetYaxis()->SetTitle("Entries");

  // Draw the histograms
  canvas->cd(1);
  hist1->Draw();
  canvas->Update();

  canvas->cd(2);
  hist2->Draw();
  canvas->Update();

  // Save the canvas to a ROOT file
  canvas->SaveAs("energy_deposition.root");

  // Clean up memory
  delete hist1;
  delete hist2;
  delete canvas;
  //analysisManager->Write();
  //analysisManager->CloseFile();

  
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

