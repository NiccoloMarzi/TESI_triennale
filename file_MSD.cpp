#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TLegend.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TRandom.h"
#include "TStyle.h"
#include "TTree.h"
#include "TCanvas.h"
#include "THStack.h"
#include <vector>
#include <TChain.h>

void loop(){
    // /////apertura del file newgeom///////////////////////////
    TFile *file1 = new TFile("tree4306_newgeom_MAR2022.root"); //apro il file ROOT
    //file1->ls(); //listing del file
    TTree *Tout1 = (TTree*)file1->Get("Tree;5");
    //Tout1->Print(); //listing del contenuto del tree
    Int_t entries1 = Tout1->GetEntries(); //numero di entries nel tree
    
    
    // /////apertura del file pileup////////////////////////////
    TFile *file2 = new TFile("tree4306_pileup_MAR2022.root"); //apro il file ROOT
    //file2->ls(); //listing del file
    TTree *Tout2 = (TTree*)file2->Get("Tree;3");
    //Tout2->Print(); //listing del contenuto del tree
    Int_t entries2 = Tout2->GetEntries(); //numero di entries nel tree
    
    
    /*TChain *chain = new TChain("tree");
    chain->Add("tree4306_newgeom_MAR2022.root");
    chain->Add("tree4307_newgeom_MAR2022.root");
    TTree *tree = chain;*/
    
    // /////riempimento istogrammi//////////////////////////////
    Int_t TWPoints;
    Tout1->SetBranchAddress("TWPoints",&TWPoints);
    
    vector<int> *MSDPoints = 0;
    TBranch *bMSDPoints = 0;
    Tout1->SetBranchAddress("MSDPoints", &MSDPoints, &bMSDPoints);
    
    Bool_t Frag;
    Tout1->SetBranchAddress("Frag", &Frag);
    
    Bool_t SCPileup;
    Tout2->SetBranchAddress("SCPileup",&SCPileup);
    
    Double_t Livetime;
    Tout1->SetBranchAddress("Livetime",&Livetime);
    
    vector<double> *MSDDe1Point = 0;
    TBranch *bMSDDe1Point = 0;
    Tout1->SetBranchAddress("MSDDe1Point", &MSDDe1Point, &bMSDDe1Point);
    
    vector<int> *TWChargePoint = 0;
    TBranch *bTWChargePoint = 0;
    Tout1->SetBranchAddress("TWChargePoint", &TWChargePoint, &bTWChargePoint);
    
    vector<double> *MSDXPoint = 0;
    TBranch *bMSDXPoint = 0;
    Tout1->SetBranchAddress("MSDXPoint", &MSDXPoint, &bMSDXPoint);
    
    vector<double> *MSDYPoint = 0;
    TBranch *bMSDYPoint = 0;
    Tout1->SetBranchAddress("MSDYPoint", &MSDYPoint, &bMSDYPoint);
    
    vector<double> *TWDe1Point = 0;
    TBranch *bTWDe1Point = 0;
    Tout1->SetBranchAddress("TWDe1Point", &TWDe1Point, &bTWDe1Point);
    

    TH1F *hTWPoints = new TH1F("hTWPoints","TW points", 10,0,10);
    TH1F *hMSDPoints = new TH1F("hMSDPoints","MSD points", 20,0,20);
    
    TH1F *hMSDP3L[3];
    hMSDP3L[0] = new TH1F("hMSDP3L","MSD points in layer 0",20,0,20);
    hMSDP3L[1] = new TH1F("hMSDP3L","MSD points in layer 1",20,0,20);
    hMSDP3L[2] = new TH1F("hMSDP3L","MSD points in layer 2",20,0,20);
    
    TH1F *hMSDDe1 = new TH1F("hMSDDe1","MSDDe", 500, 90,6000);
    
    TH1F *hlivetime = new TH1F("hlivetime","livetime", 600, 0, 2400);
    TH2F *hMSDDe_ltime = new TH2F("hMSDDe_ltime", "MSDDe vs Livetime",400, 90,6000,100,0,2400);
    
    TF1 *retta = new TF1("retta", "x", 0,2);
    
    
    for(Int_t i = 0; i < entries1; i++){
        Long64_t tentry = Tout1->LoadTree(i);
        Tout1->GetEntry(i);
        Tout2->GetEntry(i);
        // quì si fa il fill di un istogramma cumulativo dei punti visti dall'MSD con il vincolo che il TW abbia visto un solo punto////////////////////
        /*if(TWPoints==1){
            hTWPoints->Fill(TWPoints);
            int sum = 0;
            for(Int_t j = 0; j < MSDPoints->size(); j++){
                sum = sum + MSDPoints->at(j);
            }
            hMSDPoints->Fill(sum);
        }*/
        
        // quì si fa il fill di 3 istogrammi, uno per ogni layer dell'MSD, con il vincolo che il TW abbia visto un solo punto////////////////////
        /*if(TWPoints==1){
            hTWPoints->Fill(TWPoints);
            for(Int_t j = 0; j < MSDPoints->size(); j++){
                hMSDP3L[j]->Fill(MSDPoints->at(j));
            }
        }*/
        
        // quì si fa il fill di 3 istogrammi, uno per ogni layer, prenendo tutti gli eventi////////////////////
        /*for(Int_t j = 0; j < MSDPoints->size(); j++){
            hMSDP3L[j]->Fill(MSDPoints->at(j));
        }*/
        // //////////////////// Se c'è stata frammentazione quanti punti hanno visto le MSD?
        /*if(Frag==1){
            int sum = 0;
            for(Int_t j = 0; j < MSDPoints->size(); j++){
                sum = sum + MSDPoints->at(j);
            }
            hMSDPoints->Fill(sum);
        }*/
        // //////////////////// Cosa hanno visto le MSD quando il TW ha visto 0 nel pileup ?
        /*if((SCPileup==1) && (TWPoints == 0)){
            int sum = 0;
            for(Int_t j = 0; j < MSDPoints->size(); j++){
                sum = sum + MSDPoints->at(j);
            }
            hMSDPoints->Fill(sum);
            bMSDDe1Point->GetEntry(tentry);
            for(Int_t j = 0; j < MSDDe1Point->size(); j++){
                hMSDDe1->Fill(MSDDe1Point->at(j));
            }
        }*/
        
        // //////////////energia persa nelle MSD con specifica della carica////////////////
        /*bTWChargePoint->GetEntry(tentry);
        bMSDDe1Point->GetEntry(tentry);
        int counter=0;
        for(Int_t j = 0; j < MSDDe1Point->size(); j++){
            
            for(Int_t k = 0; k < TWChargePoint->size(); k++){
        
                if((TWChargePoint->at(k))==7){
                    counter++;
                    hMSDDe1->Fill(MSDDe1Point->at(j));
                }
            
            }
        }
        std::cout<<counter<<" ";*/
        
        // /////////////////geometria//////////////////////
        /*bTWChargePoint->GetEntry(tentry);
        bMSDPoints->GetEntry(tentry);
        bMSDXPoint->GetEntry(tentry);
        bMSDYPoint->GetEntry(tentry);
        if(TWPoints==1){
            for(Int_t j = 0; j < TWChargePoint->size(); j++){
                if(TWChargePoint->at(j)==8){
                        if((MSDPoints->at(0)==1)&&(MSDPoints->at(1)==1)&&(MSDPoints->at(2)==1)){
                            for(Int_t k = 0; k< MSDXPoint->size(); k++){
                                std::cout<< k <<": " << MSDXPoint->at(k) <<"\n";
                            }
                        }
                }
                
            }
        }
        std::cout << "\n";*/
        
        // ///////////////energia persa per carica nelle MSD con selezione di 1pt nel TW e 3 nelle MSD//////////////
        /*bTWChargePoint->GetEntry(tentry);
        bMSDDe1Point->GetEntry(tentry);
        bMSDPoints->GetEntry(tentry);
        if((TWPoints==1)&&(MSDPoints->at(0)==1)&&(MSDPoints->at(1)==1)&&(MSDPoints->at(2)==1)){
            for(Int_t j = 0; j < MSDDe1Point->size(); j++){
                for(Int_t k = 0; k < TWChargePoint->size(); k++){
                    if((TWChargePoint->at(k))==8){
                    hMSDDe1->Fill(MSDDe1Point->at(j));
                    }
                }
            }
        }*/
        
        // /////////////energia persa nelle MSD vs Livetime///////////////////////
        /*bTWChargePoint->GetEntry(tentry);
        bMSDDe1Point->GetEntry(tentry);
        int counter=0;
        for(Int_t j = 0; j < MSDDe1Point->size(); j++){
            for(Int_t k = 0; k < TWChargePoint->size(); k++){
                if((TWChargePoint->at(k))==1){
                    hMSDDe1->Fill(MSDDe1Point->at(j));
                    hMSDDe_ltime->Fill(MSDDe1Point->at(j),Livetime);
                }
            
            }
        }*/
        // ////////////////////////
        bMSDDe1Point->GetEntry(tentry);
        bTWDe1Point->GetEntry(tentry);
        for(Int_t j = 0; j < TWDe1Point->size(); j++){
            if((Frag==1)&&(TWChargePoint->at(j)<52)){
                for(Int_t k = 0; k < MSDDe1Point->size(); k++){
                    hMSDDe1->Fill(MSDDe1Point->at(k));
                }
            }
        }
    }
    
    
    
    // /////creazione delle Canvas//////////////////////////////
    /*TCanvas *cNumPoints = new TCanvas("cNumPoints", "",10,20,500,600);
    cNumPoints->Divide(2,1);
    cNumPoints->cd(1);
    gPad->SetLogy();
    hTWPoints->Draw();
    cNumPoints->cd(2);
    gPad->SetLogy();
    hMSDPoints->Draw();*/
    
    /*TCanvas *cNP3L = new TCanvas("cNP3L", "",10,20,500,600); //istogrammi dei punti per ogni layer dell'MSD
    cNP3L->Divide(2,2);
    cNP3L->cd(1);
    hTWPoints->Draw();
    cNP3L->cd(2);
    hMSDP3L[0]->Draw();
    cNP3L->cd(3);
    hMSDP3L[1]->Draw();
    cNP3L->cd(4);
    hMSDP3L[2]->Draw();*/
    
    /*TCanvas *c1 = new TCanvas("c1", "",10,20,500,600); //confronto dei punti nei tre layer
    hMSDP3L[0]->SetLineColor(2);
    hMSDP3L[1]->SetLineColor(3);
    hMSDP3L[2]->SetLineColor(4);
    hMSDP3L[0]->Draw();
    hMSDP3L[1]->Draw("same");
    hMSDP3L[2]->Draw("same");*/
    
    /*TCanvas *cPile1TW0 = new TCanvas("cPile1TW0", "",10,20,500,600);
    cPile1TW0->Divide(2,1);
    cPile1TW0->cd(1);
    gPad->SetLogy();
    hMSDDe1->Draw();
    cPile1TW0->cd(2);
    hMSDPoints->Draw();*/
    
    hMSDDe1->Draw();
    //hlivetime->Draw();
    //hMSDDe_ltime->Draw("COLZ");
    
}
