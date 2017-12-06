#include "amschain.h"
#include "root.h"
#include "SpFold.h"
#include "TFile.h"
#include "TKey.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include "TMath.h"
#include <iostream>
#include <map>
#include "mdst.h"
#include "commons.h"
#include "TH2.h"
#include "TF1.h"
#include "TROOT.h"
#include <cstdlib>
#include "TrCharge.h"

#include "TrExtAlignDB.h"
#include "Tofrec02_ihep.h"
#include "TrdKCluster.h"
#include "GammaFit.h"
#include "TrRecon.h"
#include "EcalH.h"
#include "bcorr.h"
#include "tkdcards.h"

#include<algorithm>

//using namespace std;
Int_t fError = 0;

bool debug = 1;


class DSTFill:public DST{
    public:
        TChain *fCh;
        TFile *fFile;
        TTree *fTree;
        Int_t fEntry;
        Int_t fMode;
        Int_t IsMC;

        Int_t    fMCpl;
        Double_t fMCW;
        Double_t fRcut;

        AMSChain *fAch;
        Int_t     fBrun;
        Int_t     fEofs;

        //ecal
        bool fIsMIP;
        bool fIsMIPL[18];
        Int_t Nhitl;
        Double_t fEdep;
        Int_t icell[18];
        Int_t NIsMIPL;

        AMSEventR * fpev;
        Double_t    fRig;
        TrTrackR *fTrk;
        Double_t fEkpn;

        bool fbtof,fbtrd,fbrich;

        DSTFill(TChain *ch, TString ofn, Int_t mode):fCh(ch),fMCW(1),fMCpl(1),fIsMIP(1),fRig(0.),
                                                     fAch(0),fBrun(0),fEofs(0),fEdep(0){

            fMode = mode%10;
            IsMC  = mode/10;

            fFile = new TFile(ofn, "recreate");
            if(!fFile || !fFile->IsOpen()){
                cout<<"I/O Error"<<endl;
                exit(0);
            }

            SetAddress(ch);
            //ch->SetBranchStatus("part", 0);
            TH1F *he = 0;
            TObjArray *ar = ch->GetListOfFiles();
            for(int i=0;i<ar->GetEntries();i++){
                TString hfn = ar->At(i)->GetTitle();
                TFile *hf = TFile::Open(hfn);
                if(!hf)continue;
                TList *keys = hf->GetListOfKeys();
                for(int j=0;j<keys->GetSize();j++){
                    TKey *key = (TKey*)keys->At(j);
                    if(!key)continue;
                    TString shn = key->GetName();
                    if(shn=="hist00"){
                        fFile->cd();
                        TH1F *htmp = (TH1F*)hf->Get("hist00");
                        if(!htmp)continue;
                        cout<<"Add hist00 from "<<hfn<<endl;
                        if(!he)he =  (TH1F*)htmp->Clone("hexp");
                        else he->Add(htmp);
                    }
                }
                hf->Close();
                delete hf;
            }
            Init();
            EventLoop();
        }


        void LoadTree(Int_t entry){fEntry = fCh->LoadTree(entry); Clear();}
        void GetEntry(const TString bname){
            TBranch *br = fCh->FindBranch(bname);
            if(br)br->GetEntry(fEntry);
            else{ if(fError++<20)cout<<"Error Branch: "<<bname<<endl;}
        }
        void Init();
        void EventLoop();
        void Fill(TString,Float_t,Float_t,Float_t);

        TrTrackR * GetTrack(void);
        int ecalmap_MIP(AMSEventR *pev, int ic_ecal);
        int read_ecalmap(AMSEventR *pev, int ic_ecal, int k);
        void NewRun(Int_t run);
        void Failed(const char *errinfo);
        int EventSelection();
        void BinLogX(TString,TString,TString,int);
        void SetTitles(TString hname,TString xtt,TString ytt);

        int BasicSelection();
        int TRDTOFRICHSelection(int);
        int TOFstudy(int);
        int RICHstudy(int);
        int TRDstudy(int);
        int TRACKERstudy(int);
        int ECALstudy(int);
};

void DSTFill::Fill(TString hname, Float_t x, Float_t y=1, Float_t w = 1){
    TObject *obj = fFile->Get(hname);
    if(!obj){
        cout<<Form("Hist %s not found",hname.Data())<<endl;
        return ;
    }
    TString htype = obj->ClassName();
    if(htype.Contains("TH1")){
        TH1 *h = (TH1 *)obj;
        h->Fill(x,y);
    }
    else{
        TH2 *h = (TH2 *)obj;
        h->Fill(x,y,w);
    }
}

void DSTFill::BinLogX(TString hname,TString xtt,TString ytt, int k) 
{
    TH2 *h = (TH2 *)fFile->Get(hname);
    if(!h){
        cout<<Form("Hist %s not found",hname.Data())<<endl;
        return ;
    }
    h->SetXTitle(xtt);
    h->SetYTitle(ytt);
    if(!k) return;

   TAxis *axis = h->GetXaxis(); 
   int bins = axis->GetNbins();

   Axis_t from = axis->GetXmin();
   Axis_t to = axis->GetXmax();
   Axis_t width = (to - from) / bins;
   Axis_t *new_bins = new Axis_t[bins + 1];

   for (int i = 0; i <= bins; i++) {
     new_bins[i] = TMath::Power(10, from + i * width);
   } 
   axis->Set(bins, new_bins); 
   delete new_bins; 
}

void DSTFill::SetTitles(TString hname,TString xtt,TString ytt) 
{
    TObject *obj = fFile->Get(hname);
    if(!obj){
        cout<<Form("Hist %s not found",hname.Data())<<endl;
        return ;
    }
    TString htype = obj->ClassName();
    if(htype.Contains("TH1")){
        TH1 *h = (TH1 *)obj;
        h->SetXTitle(xtt);
        h->SetYTitle(ytt);
    }
    else{
        TH2 *h = (TH2 *)obj;
        h->SetXTitle(xtt);
        h->SetYTitle(ytt);
    }
}


void DSTFill::Init(){
    Int_t    nbin = 0; 
    Double_t *bin = 0;

    TString sbn0 = "/afs/cern.ch/ams/Offline/AMSDataDir/v5.00/phe_bin2.root";
    TFile *fb = TFile::Open(sbn0); if(!fb) exit(-1);

    TH1F *hbin = (TH1F *)fb->Get("hist1");
    nbin = hbin->GetNbinsX();
    bin  = hbin->GetXaxis()->GetXbins()->fArray;
    cout << "Use new binning" << endl;
    fFile->cd();

    new TH1F("hist11", "Exposure (RTI V3)",          nbin, bin);

    TString shn;
    TString sthn,stht;
    for(int i=1;i<=2;i++){
        //for TOF
        shn = Form("hist%d1", i);
        new TH2F(shn+"1", "TOF Beta VS Rigidity",  300,0,10   ,300,0.3,1.4);
        new TH2F(shn+"10","TOF Beta VS Rigidity",  300,-1,1   ,300,0.3,1.4);
        SetTitles(shn+"1","Rigidity[GV]","TOF Beta");
        BinLogX(shn+"10","Rigidity[GV]","TOF Beta",1);
//estim
/*        new TH2F(shn+"10_estim", "Estimated Beta VS Rigidity Layer 0",  300,0,10   ,300,0.3,1.4);
        new TH2F(shn+"11_estim", "Estimated Beta VS Rigidity Layer 1",  300,0,10   ,300,0.3,1.4);
        new TH2F(shn+"12_estim", "Estimated Beta VS Rigidity Layer 2",  300,0,10   ,300,0.3,1.4);
        new TH2F(shn+"13_estim", "Estimated Beta VS Rigidity Layer 3",  300,0,10   ,300,0.3,1.4);
        SetTitles(shn+"10_estim","Rigidity[GV]","Estimated Beta");
        SetTitles(shn+"11_estim","Rigidity[GV]","Estimated Beta");
        SetTitles(shn+"12_estim","Rigidity[GV]","Estimated Beta");
        SetTitles(shn+"13_estim","Rigidity[GV]","Estimated Beta");*/
        new TH2F(shn+"10_estim", "Estimated Beta VS Rigidity Layer 0",  300,-1,1   ,300,0.3,1.4);
        new TH2F(shn+"11_estim", "Estimated Beta VS Rigidity Layer 1",  300,-1,1   ,300,0.3,1.4);
        new TH2F(shn+"12_estim", "Estimated Beta VS Rigidity Layer 2",  300,-1,1   ,300,0.3,1.4);
        new TH2F(shn+"13_estim", "Estimated Beta VS Rigidity Layer 3",  300,-1,1   ,300,0.3,1.4);
        BinLogX(shn+"10_estim","Rigidity[GV]","Estimated Beta",1);
        BinLogX(shn+"11_estim","Rigidity[GV]","Estimated Beta",1);
        BinLogX(shn+"12_estim","Rigidity[GV]","Estimated Beta",1);
        BinLogX(shn+"13_estim","Rigidity[GV]","Estimated Beta",1);


        new TH2F(shn+"2", "TOF 1/Beta VS Rigidity",  300,0,10   ,300,0,3);
        new TH2F(shn+"20","TOF 1/Beta VS Rigidity",  300,-1,1   ,300,0,3);
        SetTitles(shn+"2","Rigidity[GV]","TOF 1/Beta");
        BinLogX(shn+"20","Rigidity[GV]","TOF 1/Beta",1);

        new TH2F(shn+"3", "TOF 1/Beta^2 VS Rigidity",  300,0,10   ,400,0,6.5);
        new TH2F(shn+"30","TOF 1/Beta^2 VS Rigidity",  300,-1,1   ,400,0,6.5);
        SetTitles(shn+"3","Rigidity[GV]","TOF 1/Beta^2");
        BinLogX(shn+"30","Rigidity[GV]","TOF 1/Beta^2",1);

        new TH2F(shn+"4", "TOF Beta*Gamma VS Rigidity",  300,0,10   ,400,0,10);
        new TH2F(shn+"40","TOF Beta*Gamma VS Rigidity",  300,-1,1   ,400,0,10);
        SetTitles(shn+"4","Rigidity[GV]","TOF Beta*Gamma");
        BinLogX(shn+"40","Rigidity[GV]","TOF Beta*Gamma",1);

        new TH2F(shn+"5",  "TOF dE/dx VS Rigidity",         300,0,10,  500,0,30);
        new TH2F(shn+"500","TOF dE/dx VS Rigidity",         300,-1,1,  500,0,30);
        new TH2F(shn+"501","TOF dE/dx VS Rigidity(Truncated Mean Method)",300,-1,1,  500,0,30);
        new TH2F(shn+"50", "TOF dE/dx VS Rigidity Layer 0", 300,-1,1,  500,0,30);
        new TH2F(shn+"51", "TOF dE/dx VS Rigidity Layer 1", 300,-1,1,  500,0,30);
        new TH2F(shn+"52", "TOF dE/dx VS Rigidity Layer 2", 300,-1,1,  500,0,30);
        new TH2F(shn+"53", "TOF dE/dx VS Rigidity Layer 3", 300,-1,1,  500,0,30);
        SetTitles(shn+"5","Rigidity [GV]","dE/dx [MeV]");
        BinLogX(shn+"500","Rigidity [GV]","dE/dx [MeV]",1);
        BinLogX(shn+"501","Rigidity [GV]","dE/dx [MeV]",1);
        BinLogX(shn+"50","Rigidity [GV]","dE/dx [MeV]",1);
        BinLogX(shn+"51","Rigidity [GV]","dE/dx [MeV]",1);
        BinLogX(shn+"52","Rigidity [GV]","dE/dx [MeV]",1);
        BinLogX(shn+"53","Rigidity [GV]","dE/dx [MeV]",1);
//estim
        new TH2F(shn+"5_estim",  "TOF dE/dx VS Rigidity",         300,0,10,  500,0,30);
        new TH2F(shn+"500_estim","TOF dE/dx VS Rigidity",         300,-1,1,  500,0,30);
        new TH2F(shn+"501_estim","TOF dE/dx VS Rigidity(Truncated Mean Method)",300,-1,1,  500,0,30);
        SetTitles(shn+"5_estim","Rigidity [GV]","dE/dx [MeV]");
        BinLogX(shn+"500_estim","Rigidity [GV]","dE/dx [MeV]",1);
        BinLogX(shn+"501_estim","Rigidity [GV]","dE/dx [MeV]",1);


        new TH2F(shn+"6",  "TOF dE/dx VS Beta*Gamma",         300,0,10,  400,0,30);
        new TH2F(shn+"600","TOF dE/dx VS Beta*Gamma",         300,-1,1,  400,0,30);
        new TH2F(shn+"601","TOF dE/dx VS Beta*Gamma(Truncated Mean Method)",300,-1,1,  400,0,30);
        new TH2F(shn+"60", "TOF dE/dx VS Beta*Gamma Layer 0", 300,0,10,  400,0,30);
        new TH2F(shn+"61", "TOF dE/dx VS Beta*Gamma Layer 1", 300,0,10,  400,0,30);
        new TH2F(shn+"62", "TOF dE/dx VS Beta*Gamma Layer 2", 300,0,10,  400,0,30);
        new TH2F(shn+"63", "TOF dE/dx VS Beta*Gamma Layer 3", 300,0,10,  400,0,30);
        SetTitles(shn+"6","Beta*Gamma","dE/dx [MeV]");
        BinLogX(shn+"600","Beta*Gamma","dE/dx [MeV]",1);
        BinLogX(shn+"601","Beta*Gamma","dE/dx [MeV]",1);
        SetTitles(shn+"60","Beta*Gamma","dE/dx [MeV]");
        SetTitles(shn+"61","Beta*Gamma","dE/dx [MeV]");
        SetTitles(shn+"62","Beta*Gamma","dE/dx [MeV]");
        SetTitles(shn+"63","Beta*Gamma","dE/dx [MeV]");

        new TH2F(shn+"7",  "TOF dE/dx VS Beta",         300,0.3,1.4,  400,0,30);
        new TH2F(shn+"701","TOF dE/dx VS Beta(Truncated Mean Method)",300,0.3,1.4,  400,0,30);
        new TH2F(shn+"70", "TOF dE/dx VS Beta Layer 0", 300,0.3,1.4,  400,0,30);
        new TH2F(shn+"71", "TOF dE/dx VS Beta Layer 1", 300,0.3,1.4,  400,0,30);
        new TH2F(shn+"72", "TOF dE/dx VS Beta Layer 2", 300,0.3,1.4,  400,0,30);
        new TH2F(shn+"73", "TOF dE/dx VS Beta Layer 3", 300,0.3,1.4,  400,0,30);
        SetTitles(shn+"7","Beta","dE/dx [MeV]");
        BinLogX(shn+"701","Beta","dE/dx [MeV]",0);
        SetTitles(shn+"70","Beta","dE/dx [MeV]");
        SetTitles(shn+"71","Beta","dE/dx [MeV]");
        SetTitles(shn+"72","Beta","dE/dx [MeV]");
        SetTitles(shn+"73","Beta","dE/dx [MeV]");
//estim
        new TH2F(shn+"7_estim",  "TOF dE/dx VS Estimated Beta",         300,0.3,1.4,  400,0,30);
        new TH2F(shn+"701_estim","TOF dE/dx VS Estimated Beta(Truncated Mean Method)",300,0.3,1.4,  400,0,30);
        new TH2F(shn+"70_estim", "TOF dE/dx VS Estimated Beta Layer 0", 300,0.3,1.4,  400,0,30);
        new TH2F(shn+"71_estim", "TOF dE/dx VS Estimated Beta Layer 1", 300,0.3,1.4,  400,0,30);
        new TH2F(shn+"72_estim", "TOF dE/dx VS Estimated Beta Layer 2", 300,0.3,1.4,  400,0,30);
        new TH2F(shn+"73_estim", "TOF dE/dx VS Estimated Beta Layer 3", 300,0.3,1.4,  400,0,30);
        SetTitles(shn+"7_estim","Estimated Beta","dE/dx [MeV]");
        BinLogX(shn+"701_estim","Estimated Beta","dE/dx [MeV]",0);
        SetTitles(shn+"70_estim","Estimated Beta","dE/dx [MeV]");
        SetTitles(shn+"71_estim","Estimated Beta","dE/dx [MeV]");
        SetTitles(shn+"72_estim","Estimated Beta","dE/dx [MeV]");
        SetTitles(shn+"73_estim","Estimated Beta","dE/dx [MeV]");

/*
        new TH2F(shn+"8",  "Mass VS Rigidity",     300,0,10,  400,0,8);
        new TH2F(shn+"81", "1/Mass VS Rigidity",   300,0,10,  400,0,3);
        new TH2F(shn+"82", "1/Mass^2 VS Rigidity", 300,0,10,  400,0,3);
        SetTitles(shn+"8","Rigidity[GV]","Reconstructed Mass");
        SetTitles(shn+"81","Rigidity[GV]","Reconstructed Mass");
        SetTitles(shn+"82","Rigidity[GV]","Reconstructed Mass");
*/
        new TH2F(shn+"8",  "Mass VS Rigidity",     300,-1,1,  400,0,8);
        new TH2F(shn+"81", "1/Mass VS Rigidity",   300,-1,1,  400,0,3);
        new TH2F(shn+"82", "1/Mass^2 VS Rigidity", 300,-1,1,  400,0,3);
        BinLogX(shn+"8","Rigidity[GV]","Reconstructed Mass",1);
        BinLogX(shn+"81","Rigidity[GV]","Reconstructed Mass",1);
        BinLogX(shn+"82","Rigidity[GV]","Reconstructed Mass",1);


        new TH2F(shn+"9",  "1/Mass VS Kinetic Energy",     600,0,50,  300,0,3);
        SetTitles(shn+"9","Kinetic Energy/nuc [GeV/n]","Reconstructed Mass");


        //for RICH  j 0:Agl  1:NaF
        for(int j=0;j<=1;j++){
            shn = Form("hist%d2%d", i,j);
            if(j==0){
                new TH2F(shn+"1", "RICH Beta VS Rigidity",  300,0,10   ,300,0.93,1.03);
                new TH2F(shn+"10","RICH Beta VS Rigidity",  300,-1,1   ,300,0.93,1.03);
                new TH2F(shn+"2", "RICH 1/Beta VS Rigidity",  300,0,10   ,300,0.96,1.08);
                new TH2F(shn+"20","RICH 1/Beta VS Rigidity",  300,-1,1   ,300,0.96,1.08);
                new TH2F(shn+"3", "RICH 1/Beta^2 VS Rigidity",  300,0,10   ,300,0.9,1.2);
                new TH2F(shn+"30","RICH 1/Beta^2 VS Rigidity",  300,-1,1   ,300,0.9,1.2);
                new TH2F(shn+"4", "RICH Beta*Gamma VS Rigidity",  300,0,10   ,300,3,10);
                new TH2F(shn+"40","RICH Beta*Gamma VS Rigidity",  300,-1,1   ,300,3,10);
            }
            else if(j==1){
                new TH2F(shn+"1", "RICH Beta VS Rigidity",  300,0,10   ,300,0.75,1.05);
                new TH2F(shn+"10","RICH Beta VS Rigidity",  300,-1,1   ,300,0.75,1.05);
                new TH2F(shn+"2", "RICH 1/Beta VS Rigidity",  300,0,10   ,300,0.9,1.35);
                new TH2F(shn+"20","RICH 1/Beta VS Rigidity",  300,-1,1   ,300,0.9,1.35);
                new TH2F(shn+"3", "RICH 1/Beta^2 VS Rigidity",  300,0,10   ,300,0.8,1.8);
                new TH2F(shn+"30","RICH 1/Beta^2 VS Rigidity",  300,-1,1   ,300,0.8,1.8);
                new TH2F(shn+"4", "RICH Beta*Gamma VS Rigidity",  300,0,10   ,300,1,10);
                new TH2F(shn+"40","RICH Beta*Gamma VS Rigidity",  300,-1,1   ,300,1,10);
            }
            SetTitles(shn+"1","Rigidity[GV]","RICH Beta");
            BinLogX(shn+"10","Rigidity[GV]","RICH Beta",1);
            SetTitles(shn+"2","Rigidity[GV]","RICH 1/Beta");
            BinLogX(shn+"20","Rigidity[GV]","RICH 1/Beta",1);
            SetTitles(shn+"3","Rigidity[GV]","RICH 1/Beta^2");
            BinLogX(shn+"30","Rigidity[GV]","RICH 1/Beta^2",1);
            SetTitles(shn+"4","Rigidity[GV]","RICH Beta*Gamma");
            BinLogX(shn+"40","Rigidity[GV]","RICH Beta*Gamma",1);

            new TH2F(shn+"5", "RICH/TOF Beta Consistency Verification",  300,0.3,1.4,  300,-0.4,0.4 );
            new TH2F(shn+"6","RICH/TOF Beta Consistency Verification VS Rigidity",  300,-1,1, 300,-1,1);
            SetTitles(shn+"5","TOF Beta","(RICHBeta - TOFBeta)/ TOFBeta");
            BinLogX(shn+"6","Rigidity[GV]","RICHBeta/TOFBeta - 1)",1);

            new TH2F(shn+"7",  "Mass VS Rigidity",     300,0,10,  300,0,8);
            new TH2F(shn+"71", "1/Mass VS Rigidity",   300,0,10,  400,0,3);
            new TH2F(shn+"72", "1/Mass^2 VS Rigidity", 300,0,10,  400,0,3);
            new TH2F(shn+"73", "Mass^2 VS Rigidity",   300,0,10,  300,0,30);
            SetTitles(shn+"7","Rigidity[GV]","Reconstructed Mass");
            SetTitles(shn+"71","Rigidity[GV]","Reconstructed Mass^-1");
            SetTitles(shn+"72","Rigidity[GV]","Reconstructed Mass^-2");
            SetTitles(shn+"72","Rigidity[GV]","Reconstructed Mass^2");

            new TH2F(shn+"8",  "1/Mass VS Kinetic Energy",     500,0,40,  300,0,3);
            SetTitles(shn+"8","Kinetic Energy/nuc [GeV/n]","Reconstructed Mass");
        }



        //trd
        shn = Form("hist%d3", i);
        new TH2F(shn+"10","TRD dE/dx VS Rigidity",         300,0,10,  3000,0,3000);
        new TH2F(shn+"11","TRD dE/dx VS Rigidity",         300,-1,1,  3000,0,3000);
        new TH2F(shn+"12","TRD dE/dx VS Rigidity(Truncated Mean Method)",300,0,10,  3000,0,3000);
        new TH2F(shn+"13","TRD dE/dx VS Rigidity(Truncated Mean Method)",300,-1,1,  3000,0,3000);
        SetTitles(shn+"10","Rigidity [GV]","dE/dx ");
        BinLogX(shn+"11","Rigidity [GV]","dE/dx",1);
        SetTitles(shn+"12","Rigidity [GV]","dE/dx ");
        BinLogX(shn+"13","Rigidity [GV]","dE/dx",1);

        new TH2F(shn+"20","TRD dE/dx VS TOF Beta",         300,0.3,1.4,  3000,0,3000);
        new TH2F(shn+"21","TRD dE/dx VS Rigidity(Truncated Mean Method)",300,0.3,1.4,  3000,0,1500);
        BinLogX(shn+"20","Rigidity [GV]","dE/dx",0);
        BinLogX(shn+"21","Rigidity [GV]","dE/dx",0);

        TString stn,stt;
        for(int j=0;j<20;j++){
            stn = Form("%s1%02d",shn.Data(),j);
            stt = Form("TRD dE/dx VS Rigidity Layer %d",j);//lin->log
            new TH2F(stn, stt, 300,-1,1,  3000,0,3000);
            BinLogX(stn,"Rigidity [GV]","dE/dx",1);//lin->log
        }
        for(int j=0;j<20;j++){
            stn = Form("%s2%02d",shn.Data(),j);
            stt = Form("TRD dE/dx VS TOF Beta Layer %d",j);
            new TH2F(stn, stt, 300,0.3,1.4,  3000,0,3000);
            BinLogX(stn,"Rigidity [GV]","dE/dx",0);
        }



        //ecal
        shn = Form("hist%d4", i);
        if(i==1){
            new TH2F(shn+"10","Ecal Mean dE/dx vs Rigidity",                        300,0,10, 300,0,20);
            new TH2F(shn+"11","Ecal Mean dE/dx vs Rigidity(Truncated Mean Method)", 300,0,10, 300,0,10);
            new TH2F(shn+"12","Ecal Mean dE/dx vs Rigidity(Truncated Mean Method)", 300,0,10, 200,0,10);
            new TH2F(shn+"13","Ecal Mean dE/dx vs Rigidity(Truncated Mean Method)", 300,0,10, 200,0,10);
            new TH2F(shn+"20","Ecal Mean dE/dx vs Beta",                                     300,0.4,1.4, 200,0,10);
            new TH2F(shn+"21","Ecal Mean dE/dx vs Beta(Truncated Mean Method)",              300,0.4,1.4, 200,0,10);
            new TH2F(shn+"22","Ecal Mean dE/dx vs Beta(Truncated Mean Method)",              300,0.4,1.4, 200,0,10);
            new TH2F(shn+"23","Ecal Mean dE/dx vs Beta(Truncated Mean Method)",              300,0.4,1.4, 200,0,10);
            new TH2F(shn+"201","Ecal Mean dE/dx vs Beta(Truncated Mean Method)(1.0~1.1 GV)", 300,0.4,1.4, 200,0,10);
            new TH2F(shn+"202","Ecal Mean dE/dx vs Beta(Truncated Mean Method)(1.5~1.6 GV)", 300,0.4,1.4, 200,0,10);
            new TH2F(shn+"203","Ecal Mean dE/dx vs Beta(Truncated Mean Method)(2.0~2.1 GV)", 300,0.4,1.4, 200,0,10);
            new TH2F(shn+"202e","Ecal Mean dE/dx vs Estimated Beta(1.5~1.6 GV)", 300,0.4,1.4, 200,0,10);
            new TH2F(shn+"203e","Ecal Mean dE/dx vs Estimated Beta(2.0~2.1 GV)", 300,0.4,1.4, 200,0,10);
            new TH2F(shn+"231","Ecal Mean dE/dx vs Beta(Truncated Mean Method)(0.90~0.95 GV)",300,0.4,1.4, 200,0,10);
            new TH2F(shn+"232","Ecal Mean dE/dx vs Beta(Truncated Mean Method)(1.80~1.89 GV)",300,0.4,1.4, 200,0,10);
            new TH2F(shn+"30","Ecal Mean dE/dx vs Beta",                                     300,0.7,1.1, 200,0,10);
            new TH2F(shn+"31","Ecal Mean dE/dx vs Beta(Truncated Mean Method)",              300,0.7,1.1, 200,0,10);
            new TH2F(shn+"32","Ecal Mean dE/dx vs Beta(Truncated Mean Method)",              300,0.7,1.1, 200,0,10);
            new TH2F(shn+"33","Ecal Mean dE/dx vs Beta(Truncated Mean Method)",              300,0.7,1.1, 200,0,10);
            new TH2F(shn+"301","Ecal Mean dE/dx vs Beta(Truncated Mean Method)(1.0~1.1 GV)", 300,0.7,1.1, 200,0,10);
            new TH2F(shn+"302","Ecal Mean dE/dx vs Beta(Truncated Mean Method)(1.5~1.6 GV)", 300,0.7,1.1, 200,0,10);
            new TH2F(shn+"303","Ecal Mean dE/dx vs Beta(Truncated Mean Method)(2.0~2.1 GV)", 300,0.7,1.1, 200,0,10);
            new TH2F(shn+"40", "Mass VS Rigidity",      300,0,10,  400,0,8);
            new TH2F(shn+"41", "1/Mass VS Rigidity",    300,0,10,  400,0,3);
            new TH2F(shn+"42", "Mass^2 VS Rigidity",    300,0,10,  400,0,3);
            new TH2F(shn+"401", "Mass VS TOF Beta(1.0~1.1 GV)",      300,0,10,  400,0,8);
            new TH2F(shn+"411", "1/Mass VS TOF Beta(1.0~1.1 GV)",    300,0,10,  400,0,3);
            new TH2F(shn+"421", "Mass^2 VS TOF Beta(1.0~1.1 GV)",    300,0,10,  400,0,12);
            new TH2F(shn+"402", "Mass VS TOF Beta(1.5~1.6 GV)",      300,0,10,  400,0,8);
            new TH2F(shn+"412", "1/Mass VS TOF Beta(1.5~1.6 GV)",    300,0,10,  400,0,3);
            new TH2F(shn+"422", "Mass^2 VS TOF Beta(1.5~1.6 GV)",    300,0,10,  400,0,12);
            new TH2F(shn+"403", "Mass VS TOF Beta(2.0~2.1 GV)",      300,0,10,  400,0,8);
            new TH2F(shn+"413", "1/Mass VS TOF Beta(2.0~2.1 GV)",    300,0,10,  400,0,3);
            new TH2F(shn+"423", "Mass^2 VS TOF Beta(2.0~2.1 GV)",    300,0,10,  400,0,12);
        }
        if(i==2){
            new TH2F(shn+"10","Ecal Mean dE/dx vs Rigidity",                        200,0,10, 300,0,30);
            new TH2F(shn+"11","Ecal Mean dE/dx vs Rigidity(Truncated Mean Method)", 200,0,10, 300,0,30);
            new TH2F(shn+"12","Ecal Mean dE/dx vs Rigidity(Truncated Mean Method)", 200,0,10, 200,0,30);
            new TH2F(shn+"13","Ecal Mean dE/dx vs Rigidity(Truncated Mean Method)", 200,0,10, 200,0,30);
            new TH2F(shn+"20","Ecal Mean dE/dx vs Beta",                                     300,0.4,1.4, 200,0,30);
            new TH2F(shn+"21","Ecal Mean dE/dx vs Beta(Truncated Mean Method)",              300,0.4,1.4, 200,0,30);
            new TH2F(shn+"22","Ecal Mean dE/dx vs Beta(Truncated Mean Method)",              300,0.4,1.4, 200,0,30);
            new TH2F(shn+"23","Ecal Mean dE/dx vs Beta(Truncated Mean Method)",              300,0.4,1.4, 200,0,30);
            new TH2F(shn+"201","Ecal Mean dE/dx vs Beta(Truncated Mean Method)(1.5~1.6 GV)", 300,0.4,1.4, 200,0,30);
            new TH2F(shn+"202","Ecal Mean dE/dx vs Beta(Truncated Mean Method)(2.0~2.1 GV)", 300,0.4,1.4, 200,0,30);
            new TH2F(shn+"203","Ecal Mean dE/dx vs Beta(Truncated Mean Method)(1.8~1.9 GV)", 300,0.4,1.4, 200,0,30);
            new TH2F(shn+"202e","Ecal Mean dE/dx vs Estimated Beta(1.5~1.6 GV)", 300,0.4,1.4, 200,0,30);
            new TH2F(shn+"203e","Ecal Mean dE/dx vs Estimated Beta(2.0~2.1 GV)", 300,0.4,1.4, 200,0,30);
            new TH2F(shn+"231","Ecal Mean dE/dx vs Beta(Truncated Mean Method)(1.51~1.62 GV)",300,0.4,1.4, 200,0,30);
            new TH2F(shn+"232","Ecal Mean dE/dx vs Beta(Truncated Mean Method)(2.01~2.15 GV)",300,0.4,1.4, 200,0,30);
            new TH2F(shn+"30","Ecal Mean dE/dx vs Beta",                                     300,0.7,1.1, 200,0,30);
            new TH2F(shn+"31","Ecal Mean dE/dx vs Beta(Truncated Mean Method)",              300,0.7,1.1, 200,0,30);
            new TH2F(shn+"32","Ecal Mean dE/dx vs Beta(Truncated Mean Method)",              300,0.7,1.1, 200,0,30);
            new TH2F(shn+"33","Ecal Mean dE/dx vs Beta(Truncated Mean Method)",              300,0.7,1.1, 200,0,30);
            new TH2F(shn+"301","Ecal Mean dE/dx vs Beta(Truncated Mean Method)(1.5~1.6 GV)", 300,0.7,1.1, 200,0,10);
            new TH2F(shn+"302","Ecal Mean dE/dx vs Beta(Truncated Mean Method)(2.0~2.1 GV)", 300,0.7,1.1, 200,0,10);
            new TH2F(shn+"303","Ecal Mean dE/dx vs Beta(Truncated Mean Method)(2.5~2.6 GV)", 300,0.7,1.1, 200,0,10);
            new TH2F(shn+"40", "Mass VS Rigidity",      300,0,10,  400,0,8);
            new TH2F(shn+"41", "1/Mass VS Rigidity",    300,0,10,  400,0,3);
            new TH2F(shn+"42", "Mass^2 VS Rigidity",    300,0,10,  400,0,3);
            new TH2F(shn+"401", "Mass VS TOF Beta(1.5~1.6 GV)",      300,0,10,  400,0,8);
            new TH2F(shn+"411", "1/Mass VS TOF Beta(1.5~1.6 GV)",    300,0,10,  400,0,3);
            new TH2F(shn+"421", "Mass^2 VS TOF Beta(1.5~1.6 GV)",    300,0,10,  400,0,12);
            new TH2F(shn+"402", "Mass VS TOF Beta(2.0~2.1 GV)",      300,0,10,  400,0,8);
            new TH2F(shn+"412", "1/Mass VS TOF Beta(2.0~2.1 GV)",    300,0,10,  400,0,3);
            new TH2F(shn+"422", "Mass^2 VS TOF Beta(2.0~2.1 GV)",    300,0,10,  400,0,12);
            new TH2F(shn+"403", "Mass VS TOF Beta(2.5~2.6 GV)",      300,0,10,  400,0,8);
            new TH2F(shn+"413", "1/Mass VS TOF Beta(2.5~2.6 GV)",    300,0,10,  400,0,3);
            new TH2F(shn+"423", "Mass^2 VS TOF Beta(2.5~2.6 GV)",    300,0,10,  400,0,12);
        }
        SetTitles(shn+"10","Rigidity[GV]","Mean dE/dx [MeV]");
        SetTitles(shn+"11","Rigidity[GV]","Mean dE/dx [MeV]");
        SetTitles(shn+"12","Rigidity[GV]","Mean dE/dx [MeV]");
        SetTitles(shn+"13","Rigidity[GV]","Mean dE/dx [MeV]");
        SetTitles(shn+"20","TOF Beta","Mean dE/dx [MeV]");
        SetTitles(shn+"21","TOF Beta","Mean dE/dx [MeV]");
        SetTitles(shn+"22","TOF Beta","Mean dE/dx [MeV]");
        SetTitles(shn+"23","TOF Beta","Mean dE/dx [MeV]");
        SetTitles(shn+"201","TOF Beta","Mean dE/dx [MeV]");
        SetTitles(shn+"202","TOF Beta","Mean dE/dx [MeV]");
        SetTitles(shn+"203","TOF Beta","Mean dE/dx [MeV]");
        SetTitles(shn+"202e","Estimated Beta","Mean dE/dx [MeV]");
        SetTitles(shn+"203e","Estimated Beta","Mean dE/dx [MeV]");
        SetTitles(shn+"231","TOF Beta","Mean dE/dx [MeV]");
        SetTitles(shn+"232","TOF Beta","Mean dE/dx [MeV]");
        SetTitles(shn+"30","RICH Beta","Mean dE/dx [MeV]");
        SetTitles(shn+"31","RICH Beta","Mean dE/dx [MeV]");
        SetTitles(shn+"32","RICH Beta","Mean dE/dx [MeV]");
        SetTitles(shn+"33","RICH Beta","Mean dE/dx [MeV]");
        SetTitles(shn+"301","RICH Beta","Mean dE/dx [MeV]");
        SetTitles(shn+"302","RICH Beta","Mean dE/dx [MeV]");
        SetTitles(shn+"303","RICH Beta","Mean dE/dx [MeV]");
        SetTitles(shn+"40","Rigidity[GV]","Mass[GeV]");
        SetTitles(shn+"41","Rigidity[GV]","1/Mass[GeV^-1]");
        SetTitles(shn+"42","Rigidity[GV]","Mass^2[GeV^2]");
        SetTitles(shn+"401","Rigidity[GV]","Mass[GeV]");
        SetTitles(shn+"411","Rigidity[GV]","1/Mass[GeV^-1]");
        SetTitles(shn+"421","Rigidity[GV]","Mass^2[GeV^2]");
        SetTitles(shn+"402","Rigidity[GV]","Mass[GeV]");
        SetTitles(shn+"412","Rigidity[GV]","1/Mass[GeV^-1]");
        SetTitles(shn+"422","Rigidity[GV]","Mass^2[GeV^2]");
        SetTitles(shn+"403","Rigidity[GV]","Mass[GeV]");
        SetTitles(shn+"413","Rigidity[GV]","1/Mass[GeV^-1]");
        SetTitles(shn+"423","Rigidity[GV]","Mass^2[GeV^2]");


        TString sthn;
        for(int k=0;k<18;k++){
            sthn = Form("%s10%02d",shn.Data(),k);
            stht = Form("Ecal dE/dx vs Rigidity (Layer %d)",k);
            if(i==1)new TH2F(sthn, stht, 300,0,10, 300,0,20);
            if(i==2)new TH2F(sthn, stht, 300,0,10, 300,0,25);
            SetTitles(sthn,"Rigidity [GV]","dE/dx [MeV]");
        }
        for(int k=0;k<18;k++){
            sthn = Form("%s20%02d",shn.Data(),k);
            stht = Form("Ecal dE/dx vs TOF Beta (Layer %d)",k);
            if(i==1)new TH2F(sthn, stht, 200,0.4,1.4, 200,0,20);
            if(i==2)new TH2F(sthn, stht, 200,0.4,1.4, 200,0,25);
            SetTitles(sthn,"TOF Beta","dE/dx [MeV]");
        }
        for(int k=0;k<18;k++){
            sthn = Form("%s30%02d",shn.Data(),k);
            stht = Form("Ecal dE/dx vs RICH Beta (Layer %d)",k);
            if(i==1)new TH2F(sthn, stht, 200,0.7,1.1, 200,0,20);
            if(i==2)new TH2F(sthn, stht, 200,0.7,1.1, 200,0,25);
            SetTitles(sthn,"RICH Beta","dE/dx [MeV]");
        }



        shn = Form("hh%d1", i);
        new TH2F(shn+"1", "Mass VS Rigidity",     300,0,10,  300,0,8);
        new TH2F(shn+"10","Mass VS Rigidity",     300,-1,1,  300,0,8);
        new TH2F(shn+"2", "Mass VS Rigidity",     300,0,10,  300,0,8);
        new TH2F(shn+"20","Mass VS Rigidity",     300,-1,1,  300,0,8);
        SetTitles(shn+"1","Rigidity[GV]","TOF Reconstructed Mass");
        BinLogX(shn+"10","Rigidity[GV]","TOF Reconstructed Mass",1);
        SetTitles(shn+"2","Rigidity[GV]","RICH Reconstructed Mass");
        BinLogX(shn+"20","Rigidity[GV]","RICH Reconstructed Mass",1);

        shn = Form("hh%d2", i);
        new TH2F(shn+"10", "TRD Likelihoodratio VS Rigidity",     300,0,10,  300,0,3);
        new TH2F(shn+"20", "TRD Likelihoodratio VS Rigidity",     300,0,10,  300,0,3);
        new TH2F(shn+"30", "TRD Likelihoodratio VS Rigidity",     300,0,10,  300,0,3);
        SetTitles(shn+"10","Rigidity[GV]","Likelihood Ratio  llr");
        SetTitles(shn+"20","Rigidity[GV]","Likelihood Ratio  llrt");
        SetTitles(shn+"30","Rigidity[GV]","Likelihood Ratio  llre");
        new TH2F(shn+"11", "TRD Likelihoodratio VS Rigidity",     300,0,10,  300,0,3);
        new TH2F(shn+"21", "TRD Likelihoodratio VS Rigidity",     300,0,10,  300,0,3);
        new TH2F(shn+"31", "TRD Likelihoodratio VS Rigidity",     300,0,10,  300,0,3);
        SetTitles(shn+"11","Rigidity[GV]","Likelihood Ratio  llr");
        SetTitles(shn+"21","Rigidity[GV]","Likelihood Ratio  llrt");
        SetTitles(shn+"31","Rigidity[GV]","Likelihood Ratio  llre");


    }


    

}





void DSTFill::EventLoop(){
    Int_t Nent = fCh->GetEntries();
    cout<<"Ntr, Nent= "<<fCh->GetNtrees()<<' '<<fCh->GetEntries()<<endl;
    cout<<"IsMC= "<<IsMC <<endl;

    fFile->cd();
    fTree = fCh->CloneTree(0);

    Int_t Load = Nent/10;
    Int_t nfil = 0;
    Int_t fBtime =0;

    Int_t nrun = 0;
    Int_t nShower = 1;

    TStopwatch timer;
    timer.Start();

    AMSSetupR::sethead(new AMSSetupR);
    // Latest RTI database(pass6) 
    AMSSetupR::RTI::UseLatest(6);
    // Latest alignment  
    TkDBc::UseFinal();
    // Disabling reading settings from file  
    TRMCFFKEY.ReadFromFile = 0; 
    TRFITFFKEY.ReadFromFile = 0; 

    for(Int_t i=0;i<Nent;i++){
        if( (i+1)%Load==0 || (i+1)==Nent){
            Double_t rtm = timer.RealTime();
            Double_t ctm = timer.CpuTime();
            timer.Continue();
            cout << Form("%8d %8d (%5.1f%%) %4.0f sec (%5.1f kHz) CPU: %5.1f%%",
                    nfil, i+1, 100.*(i+1)/Nent, rtm, (i+1)/rtm*1e-3, ctm/rtm*100)<<endl;
        }
        LoadTree(i);


        GetEntry("header");
        GetEntry("status");
        //RTI and Status Cut
        if(!IsMC){
            Int_t stcut = BadTrig | InSAA | BadRTI;
            if( fStatus.ustat & stcut  ) continue;
            GetEntry("rti");
            if(fRTI.dl1>35 || fRTI.dl9>45) continue;
            fRcut = TMath::Abs(fRTI.cfi);
            if(fStatus.ustat & FirstRTI){
                TH1F *hist1 = (TH1F *)fFile->Get("hist11");
                if (hist1) {
                    TAxis *ax = hist1->GetXaxis();
                    Int_t  ib = ax->FindBin(fRcut*1.2);
                    fRcut = ax->GetBinLowEdge(ib+1);
                    for (Int_t i = ib; i <= hist1->GetNbinsX(); i++)
                        hist1->Fill(ax->GetBinCenter(i), fRTI.lf);
                }
            }
        }

        if (!IsMC && fRcut == 0) fRcut = fRTI.cfi*1.2;
        fRig=0.0;
        //ForMC pl1 selection
        fMCpl = 0;
        if(IsMC){
            GetEntry("mcinfo");
            fMCpl = 1;
            if(fMCinfo.coo[2]<194.999||fMCinfo.dir[2]>-0.7)fMCpl=2;
        }

        //fMCW = !IsMC?1:GetMCW(fMCinfo.rgt);
        fMCW = 1;//1114      
        //FROMHERE

        fpev= 0;

        GetEntry("trd");
        GetEntry("trdk");
        GetEntry("ecal");
        GetEntry("track");
        GetEntry("betah");
        GetEntry("beta");
        GetEntry("rich");

        int qsel = BasicSelection();

//continue;

        //if(qsel) nfil++;
        if( !(qsel==1 || qsel==2) ) continue;

        if (fHeader.run != fBrun) {
            nrun++;
            NewRun(fHeader.run);
            Double_t rtm = timer.RealTime();
            Double_t crm = timer.CpuTime();
            timer.Continue();
            cout << Form("%9d %9d (%5.1f%%) %4.0f sec (%4.f Hz) CPU: %5.1f%%",
                    nfil+1, i+1, 100.*(i+1)/Nent, rtm, 
                    (i+1)/rtm, crm/rtm*100) <<"   "<<nrun<<endl;
            fBrun = fHeader.run;
            fEofs = 0;
        }

        AMSEventR *pev = fAch->GetEvent(fHeader.ient+fEofs); //Get AMSEventR Object From AMSROOT Files
        Int_t ntry = 0;
        while (pev->Event() != fHeader.event) { //Check if we get the correct event
          if (ntry++ > 10) Failed("Cannot get the correct event!");
          fEofs += fHeader.event-pev->Event();
          pev = fAch->GetEvent(fHeader.ient+fEofs);
        }

        if(!pev || pev->NTrTrack()!=1) continue;  //only 1 trtrack
        fpev = pev;
        ParticleR * part = pev->pParticle(0);
if (!part) continue;
//TrTrackR * ptrk = part->pTrTrack();
//if (!ptrk) { fpev = pev;fTrk = GetTrack(); }
fTrk = part->pTrTrack();
if (!fTrk) fTrk = GetTrack(); 
/*
        fbtof = TOFstudy(qsel);
        //if( fbtof ) nfil++;
//if(!fbtof) continue;
        //if( RICHstudy(qsel) ) printf("RICH pass! \n");
        fbtrd = TRDstudy(qsel);
//if(!fbtrd) continue;
        fbrich = RICHstudy(qsel);
        //if( TRDstudy(qsel) ) nfil++;
//if(!fbrich) continue;
*/


ECALstudy(qsel);


if(! TRDTOFRICHSelection(qsel) ) continue;
fbtof = TOFstudy(qsel);
fbtrd = TRDstudy(qsel);
fbrich = RICHstudy(qsel);

        
        //TRACKERstudy(qsel);



    }//loopend
    fFile->cd();
    fFile->Write();
    fFile->Close();
}

//bool debug = 1;
/*
if(debug) printf("");
*/
//     1234
// histXXXX
// [1]: 1-proton    2-helium
// [2]: 1-TOF       2-RICH      3-TRACKER   4-ECAL
int DSTFill::TOFstudy(int qsel){
    //if(fBetaH.beta<=0.) return 0;
    TString shn = Form("hist%d1",qsel);
    TString sthn;

//for estimation 
//Double_t paras[2][12]={1.115872,0.949851,1.020059,1.184924,1.280588,1.099210,3.713922,3.063059,3.286361,4.241415,4.343575,3.707412,
//                       0.869103,0.960841,0.884911,0.708448,0.597497,0.871934,3.800462,4.490082,4.172249,2.808940,2.580988,3.731034};
Double_t paras[2][12]={
    1.121106, 0.958049, 1.027169, 1.189485, 1.286724, 1.102656, 3.786875, 3.138636, 3.389632, 4.314282, 4.472087, 3.788039, 
    0.867178, 0.954456, 0.876853, 0.701494, 0.587990, 0.870846, 3.788834, 4.501558, 4.096568, 2.793567, 2.448003, 3.700383};

Double_t beta_estim;

    Double_t gamma = 1.0 / TMath::Sqrt(1-fBetaH.beta*fBetaH.beta);
    Double_t betagamma = fBetaH.beta * gamma;
//if(debug) printf("TOF:  fRig:%f  \n",fRig);

    Fill(shn+"1",fRig,  fBetaH.beta);
    Fill(shn+"10",fRig,  fBetaH.beta);
    Fill(shn+"2",fRig,1./fBetaH.beta);
    Fill(shn+"20",fRig,1./fBetaH.beta);
    Fill(shn+"3",fRig,1./fBetaH.beta/fBetaH.beta);
    Fill(shn+"30",fRig,1./fBetaH.beta/fBetaH.beta);
    Fill(shn+"4",fRig,betagamma);
    Fill(shn+"40",fRig,betagamma);

    ParticleR * part = fpev->pParticle(0);
    if (!part) return 0;
        BetaHR *bth = (IsMC) ? part->pBetaH() : 0;
        for (Int_t i = 0; !bth && i < fpev->NBetaH(); i++){
//        if(fpev->pBetaH(i)->GetBeta()==fBetaH.beta) {bth=fpev->pBetaH(i);fTrk=fpev->pBetaH(i)->pTrTrack();}
        if(fpev->pBetaH(i)->GetBeta()==fBetaH.beta) {bth=fpev->pBetaH(i);if(!fTrk)fTrk=fpev->pBetaH(i)->pTrTrack();}

        //if(debug) printf("fBetaH.beta:%f   i:[%d]:%f \n",fBetaH.beta,i,fpev->pBetaH(i)->GetBeta());
        //if(fBetaH.beta - fpev->pBetaH(i)->GetBeta() > 1e-6) printf("fBetaH.beta:%f   i:[%d]:%f \n",fBetaH.beta,i,fpev->pBetaH(i)->GetBeta());
    }


    Double_t edep_tof[4] = {0.0};
    for(int k=0;k<4;k++) {
        edep_tof[k] = bth->GetEdepL(k);
        if(edep_tof[k] == 0.0) continue;
        sthn = Form("%s5%d",shn.Data(),k);
        Fill(sthn,fRig,edep_tof[k]);
        sthn = Form("%s6%d",shn.Data(),k);
        Fill(sthn,betagamma,edep_tof[k]);
        sthn = Form("%s7%d",shn.Data(),k);
        Fill(sthn,fBetaH.beta,edep_tof[k]);
//esitmation
        beta_estim = TMath::Sqrt(  paras[0][6*qsel+k-5]/(edep_tof[k]-paras[1][6*qsel+k-5])  );//proton:k+1   helium:k+1+6   -> k+1+6(qsel-1)=6*qsel+k-5
        sthn = Form("%s1%d_estim",shn.Data(),k);
        Fill(sthn,fRig,beta_estim);
        sthn = Form("%s7%d_estim",shn.Data(),k);
        Fill(sthn,beta_estim,edep_tof[k]);
    }

    Double_t edep_tof_aver = (edep_tof[0]+edep_tof[1]+edep_tof[2]+edep_tof[3])/4.0 ;
    Fill(shn+"5",fRig,edep_tof_aver);
    Fill(shn+"6",betagamma,edep_tof_aver);
    Fill(shn+"7",fBetaH.beta,edep_tof_aver);
    Fill(shn+"500",fRig,edep_tof_aver);
    Fill(shn+"600",betagamma,edep_tof_aver);
//estim
    beta_estim = TMath::Sqrt(  paras[0][qsel*6-6]/(edep_tof_aver-paras[1][qsel*6-6])  );
    Fill(shn+"5_estim",fRig,edep_tof_aver);
    Fill(shn+"500_estim",fRig,edep_tof_aver);
    Fill(shn+"7_estim",beta_estim,edep_tof_aver);

    sort(edep_tof,edep_tof+4);
    edep_tof_aver = (edep_tof[1]+edep_tof[2])/2.0 ;
    Fill(shn+"501",fRig,edep_tof_aver);
    Fill(shn+"601",betagamma,edep_tof_aver);
    Fill(shn+"701",fBetaH.beta,edep_tof_aver);
//estim
    beta_estim = TMath::Sqrt(  paras[0][qsel*6-1]/(edep_tof_aver-paras[1][qsel*6-1])  );
    Fill(shn+"501_estim",fRig,edep_tof_aver);
    Fill(shn+"701_estim",beta_estim,edep_tof_aver);

    Double_t charge = qsel;
    Double_t mass_tof = fRig * charge * TMath::Sqrt(1./fBetaH.beta/fBetaH.beta-1);
    Fill(shn+"8",fRig,mass_tof);
    Fill(shn+"81",fRig,1./mass_tof);
    Fill(shn+"82",fRig,1./mass_tof/mass_tof);


    //Double_t cl = 2.99792458e8; 
//    if(qsel==1)      fEkpn = (gamma-1)*0.938*(cl*cl);
//    else if(qsel==2) fEkpn = (gamma-1)*0.936*(cl*cl);
    if(qsel==1)      fEkpn = (gamma-1.)*0.938;
    else if(qsel==2) fEkpn = (gamma-1.)*0.936;

    Fill(shn+"9",fEkpn,1./mass_tof);
//if(debug) printf("q:%d    fEkpn:%f    1./mass:%d   ",qsel,fEkpn,1./mass_tof);
//if(debug) printf("1/mass:%d   \n",1/mass_tof);

    //if(debug) printf("==========================TOF==========================\n");

    return 1;

}


int DSTFill::RICHstudy(int qsel){ 
//TString shn = Form("hist%d2",qsel);

  ParticleR *part = fpev->pParticle(0);
  if(!part) return 0;
  RichRingR *rich = part->pRichRing();
  if(!rich) return 0;

/*//        GetEntry("betah");
if(debug) printf("===============TOF beta:%f==========RICH:prob:%f==============================\n",fBetaH.beta,fRich.prob);
if(debug) printf("===run:%d  event:%d  ient:%d   =========================\n",fHeader.run,fHeader.event,fHeader.ient);
if(debug) {printf("RICH Q:%f Qconsistency:%f  nPMTs:%d   ",rich->getCharge2Estimate(),rich->getPMTChargeConsistency(),rich->getPMTs());
           printf("IsNaf:%d  Beta:%f   IsClean:%d        nhits:%d  Prob:%f  \n",rich->IsNaF(),rich->getBeta(),rich->IsClean(),rich->getHits(),rich->getProb() );}
if(debug) {printf("mdst Q:%f tile:%f          nPMTs:%d",fRich.q,fRich.tile,fRich.pmts);
           printf("           Beta:%f  TOFBeta:%f  nhits:%d  Prob:%f  \n",fRich.beta,fBetaH.beta,fRich.nhits[1],fRich.prob);}
//if(debug) printf("======================RICH==========================\n");
//return 1;
*/
//if(debug) printf("======================RICH nPMTs:%d      fRich.pmts:%d==========================\n",rich->getPMTs(),fRich.pmts);


/*    //trdtofrichselection
    //Double_t q_rich = rich->getCharge2Estimate();//fRich.q
    Double_t q_rich = fRich.q;
    if(qsel==1){  if( q_rich<0.7||q_rich>1.4 ) return 0;}
    if(qsel==2){  if( q_rich<1.7||q_rich>2.4 ) return 0;}

    if(! rich->IsClean()) return 0;
    if(fRich.pmts<=3 ) return 0;//if(rich->getPMTs()<=3 ) return 0;
    //if(rich->getPMTChargeConsistency()<=10 )return 0;
*/


    //Double_t beta_rich = rich->getBeta();//fRich.beta
TString shn = Form("hist%d2%d",qsel,rich->IsNaF());

    if(fRich.beta<=0.) return 0;

    Double_t gamma = 1.0 / TMath::Sqrt(1-fRich.beta*fRich.beta);
    Double_t betagamma = fRich.beta * gamma;

    Fill(shn+"1",fRig,  fRich.beta);
    Fill(shn+"10",fRig,  fRich.beta);
    Fill(shn+"2",fRig,1./fRich.beta);
    Fill(shn+"20",fRig,1./fRich.beta);
    Fill(shn+"3",fRig,1./fRich.beta/fRich.beta);
    Fill(shn+"30",fRig,1./fRich.beta/fRich.beta);
    Fill(shn+"4",fRig,betagamma);
    Fill(shn+"40",fRig,betagamma);

    Fill(shn+"5",fBetaH.beta,(fRich.beta-fBetaH.beta)/fBetaH.beta);
    Fill(shn+"6",fRig,fRich.beta/fBetaH.beta-1);

    Double_t charge = qsel;
    Double_t mass_rich = fRig * charge * TMath::Sqrt(1./fRich.beta/fRich.beta-1);
    Fill(shn+"7",fRig,mass_rich);
    Fill(shn+"71",fRig,1./mass_rich);
    Fill(shn+"72",fRig,1./mass_rich/mass_rich);
    Fill(shn+"73",fRig,mass_rich*mass_rich);

    //Double_t cl = 2.99792458e8; 
    Double_t Ekpn_rich;
//    if(qsel==1)      Ekpn_rich = (gamma-1)*0.938*(cl*cl);
//    else if(qsel==2) Ekpn_rich = (gamma-1)*0.936*(cl*cl);
    if(qsel==1)      Ekpn_rich = (gamma-1)*0.938;
    else if(qsel==2) Ekpn_rich = (gamma-1)*0.936;
    Fill(shn+"8",Ekpn_rich,1.0/mass_rich);
//if(debug) printf("RICH q:%d    Ekpn_rich:%f    1./mass:%d    ",qsel,Ekpn_rich,1.0/mass_rich);
//if(debug) printf("  1/mass:%d   \n",1/mass_rich);
    return 1;
}

int DSTFill::TRDstudy(int qsel){
    TString shn = Form("hist%d3",qsel);
/*   //trdtofrichselection
    if (fTrdK.nhits <= 8) return 0;
    Double_t q_trd = fTrdK.q;
    if(qsel==1){  if( q_trd<0.7||q_trd>1.4 ) return 0;}
    if(qsel==2){  if( q_trd<1.7||q_trd>2.4 ) return 0;}

    int ntrdcluster = fpev->NTrdCluster();
//    if(debug) printf("ntrdcluster:%d  \n",ntrdcluster);
//if(debug) printf("ntrdtrack:%d  \n",fpev->NTrdTrack());
    //int ntrdtrack = fpev->NTrdTrack();
    if(fpev->NTrdTrack() != 1) return 0 ;

//    GetEntry("trhit");
//if(debug) for(int k=0;k<25;k++)printf("[%d]:layer:%f   plen:%f   amp:%f   plmc:%f  trfr:%f  \n",fTrdHit.layer[k],fTrdHit.plen[k],fTrdHit.amp[k],fTrdHit.plmc[k],fTrdHit.trfr[k]);
//printf("---------------------------------------------------------------------");

    if(!fTrk  ) return 0;
//if(!fTrk || !fpev->pParticle(0)->pTrdTrack() ) return 0;
//int fid = fTrk->iTrTrackPar(1,0,20);
*/

    Double_t ampsk,pathsk;
    Int_t namphitk=0;
    Double_t ampl[20]={0.};
    Double_t lenl[20]={0.};
    Double_t dEdxl[20]={0.};
    ampsk = pathsk = 0.0;
    TrdKCluster tkcl = TrdKCluster(fpev,fTrk,fTrk->Gettrdefaultfit());
    Int_t nhit = 0;
    Double_t llr[3] = { -1, -1, -1 }, ll[3];
    if (tkcl.IsReadAlignmentOK == 2 && tkcl.IsReadCalibOK == 1) {
        if(tkcl.GetLikelihoodRatio_TrTrack(15, &llr[0], nhit) == 1){
            AMSPoint pnt;
            AMSDir   dir;
            if( tkcl.GetTrTrackExtrapolation(pnt,dir)>=0){
                for(int k=0;k<tkcl.NHits();k++){
                    TrdKHit *hit = tkcl.GetHit(k);
                    if (!hit) continue;
	                Double_t len = hit->Tube_Track_3DLength(&pnt, &dir);
                    if(len<=0.02) continue;
	                Double_t amp = hit->TRDHit_Amp;
                    if(amp<=0) continue;
                    if(amp>32767) amp =32767;
                    ampsk += (amp/len);
                    pathsk+= len;
                    namphitk++;
                    int ilay = hit->TRDHit_Layer;
                    if(amp <= ampl[ilay]) continue;
                    ampl[ilay] = amp;
                    lenl[ilay] = len;
                    dEdxl[ilay]= amp/len;
                }
            }
        }
    }
    TString sthn;
    Double_t dEdx_aver;
    for(int k=0;k<20;k++) {
        dEdx_aver+=dEdxl[k];
        sthn = Form("%s1%02d",shn.Data(),k);
        Fill(sthn,fRig,dEdxl[k]);
        sthn = Form("%s2%02d",shn.Data(),k);
        Fill(sthn,fBetaH.beta,dEdxl[k]);
        //if(debug) {printf("[%d]:%f  ",k,dEdxl[k]);if(k==9 || k==19)printf("\n");}
    }

    dEdx_aver = dEdx_aver/namphitk;
    Fill(shn+"10",fRig,dEdx_aver);
    Fill(shn+"11",fRig,dEdx_aver);
    Fill(shn+"20",fBetaH.beta,dEdx_aver);
//if(debug) printf("dEdx_aver:[1]%f  ",dEdx_aver);
    sort(dEdxl,dEdxl+20);
    for(int k=2;k<namphitk-2;k++) dEdx_aver+=dEdxl[k];
    dEdx_aver = dEdx_aver/(namphitk-4);
    Fill(shn+"12",fRig,dEdx_aver);
    Fill(shn+"13",fRig,dEdx_aver);
    Fill(shn+"21",fBetaH.beta,dEdx_aver);
//if(debug) printf("[2]truncated met: %f  \n",dEdx_aver);
    return 1;

}



/*
if(debug) printf("");
*/





int DSTFill::TRACKERstudy(int qsel){
    TString shn = Form("hist%d5",qsel);
    
  ParticleR *part = fpev->pParticle(0);
  if(!part) return 0;
    

    float edeplay_offtrack[9]={0.};
    float edeplay_neighbor[9]={0.};
    int nhitlay_offtrack[9]={0};
    int nhitlay_neighbor[9]={0};


    for(int ihit =0;ihit<fTrk->GetNhits();ihit++ ){

        TrRecHitR *phit = fTrk->GetHit(ihit);
        int layJ = phit->GetLayerJ();
        int icly = phit->GetYClusterIndex();
        int idhit = fTrk->iTrRecHit(ihit);  //index in hit vector

        for(int jhit =0;jhit< fpev->nTrRecHit();jhit++){
            TrRecHitR * phit_offtrack = fpev->pTrRecHit(jhit);
            if(phit_offtrack->GetLayerJ() == layJ  
                && phit_offtrack->iTrCluster('y') != phit->iTrCluster('y')
                && idhit != jhit ){
                float distx = phit_offtrack->HitDist(*phit,0);
                float disty = phit_offtrack->HitDist(*phit,1);
                nhitlay_offtrack[layJ-1]++;
                edeplay_offtrack[layJ-1] += phit_offtrack->GetEdep(1);  //y side edep
                if( fabs(disty)<5. && fabs(distx)<5. ){  
                    nhitlay_neighbor[layJ-1]++;
                    edeplay_neighbor[layJ-1] += phit_offtrack->GetEdep(1);
                }
            }
        }

    }
printf("===================================================================================\n");
printf("run:%d  event:%d  ient:%d   \n",fHeader.run,fHeader.event,fHeader.ient);
printf("fRig:%f   fBetaH.beta:%f   fRich.beta:%f   fTrack.bith:%X   \n",fRig,fBetaH.beta,fRich.beta,fTrack.bith);
printf("fTrk->GetNhits()=%d     fpev->nTrRecHit()=%d   \n",fTrk->GetNhits(),fpev->nTrRecHit());
printf("ilay  nhitlay_offtrack  edeplay_offtrack   |  nhitlay_neighbor   edeplay_neighbor  \n");
for(int k=0;k<9;k++){
    printf("[%d]       %2d                %f                     %2d               %f  \n",
              k+1,nhitlay_offtrack[k],edeplay_offtrack[k],nhitlay_neighbor[k],edeplay_neighbor[k]);
if(nhitlay_neighbor[k]!=0)  printf("******************************************************************\n");
}
/*
printf("");






Double_t edepxyl[2][9]={0.};

//float eltr[2]={0.};
for(int il=0;il<9;il++){
    TrRecHitR * tkhit = fTrk->GetHitLJ(il+1);
    TrClusterR * tkcl[2] = {0};
    for(int ixy=0;ixy<2;ixy++){
        if(!tkhit)continue;
        TrClusterR * cl = (ixy==0)? tkhit->GetXCluster():tkhit->GetYCluster();
        if(cl){tkcl[ixy]=cl;edepxyl[ixy][il] =cl->GetEdep(); }
        if(debug) printf("il=%d  eltr[%d]:%f  \n",il,ixy,eltr[ixy]);
    }
}

if(debug) printf("================================================\n");
*/


    return 1;
}




int DSTFill::ECALstudy(int qsel){
    TString shn = Form("hist%d4",qsel);

    Double_t ecal_entry_z = -143.2;
    Double_t ecal_exit_z  = -158.92;
    AMSPoint _pnt;
    AMSDir   _dir;

    fTrk->Interpolate(ecal_entry_z, _pnt, _dir);
    Double_t L9theta = _dir.gettheta()*TMath::RadToDeg() ;
    if(  L9theta > 6  ) return 0;

    const Double_t Zecal[18]  = {
	    -143.215-0.005, -144.135+0.005, -145.065-0.005, -145.985+0.005, -146.915-0.005, -147.835+0.005,
	    -148.765-0.005, -149.685+0.005, -150.615-0.005, -151.535+0.005, -152.465-0.005, -153.385+0.005,
	    -154.315-0.005, -155.235+0.005, -156.165-0.005, -157.085+0.005, -158.015-0.005, -158.935+0.005
    };

    Double_t Pecalhit[3];
    //interpolate to each ecal layer   
    //X side: (xi-0.13)/36      Y side: (yi+0.07)/36      [0,71]  
    for(int k=0;k<18;k++){
        fTrk->Interpolate(Zecal[k], _pnt, _dir);
        for(int j=0;j<3;j++) Pecalhit[j]=_pnt[j];
        int iilay =  int(k/2.0)  % 2 ;
        if(  iilay == 0  )  icell[k] = int( (Pecalhit[1]+0.07)/0.9  +36 ); //Y side
        if(  iilay == 1  )  icell[k] = int( (Pecalhit[0]-0.13)/0.9  +36 ); //X side
        if(icell[k]>71 || icell[k]<0 ) icell[k]=-1;
    }


    fIsMIP= true;
    int iecals = fpev->nEcalShower();
    if( iecals != 1) return 0;

    for(int j=0;j<18;j++) fIsMIPL[j]=true;

    ecalmap_MIP(fpev,0);
    if(fRig<1.0)  { if(NIsMIPL<=3) return 0;}
    if(fRig>=1.0) { if(NIsMIPL<=3) return 0;}
    if(fRig>=2.0) { if(NIsMIPL<=4) return 0;}
    if(fRig>=2.0) { if(NIsMIPL<=5) return 0;}

    EcalShowerR * eshw = fpev->pEcalShower (0);
    Nhitl = eshw->NbLayerX + eshw->NbLayerY;

    Double_t dx = 0.9 / _dir.gettheta(); 


    if(qsel ==1){

        //if(fEcal.bdt > -0.8) continue;//1117

        Double_t edep[18],edepl[18];
        Int_t nedep0 = 0;
        Double_t edepthresh = 0.0;

//printf("====================================Proton==========================================\n");
        for(int k=0;k<18;k++){
            edepl[k] = edep[k] = eshw->DepositedEnergyInMeVPerLayer(k);
            if(edepl[k] == 0 ) nedep0++;
//printf("[%d]:%f  ",k,edepl[k]);
//if(k==8)printf("\n");
        }
//printf("\n---#0:[%d]-----------ABOVE: before resort     BELOW:after resort----------- \n",nedep0);
        sort(edep,edep+18);
//for(int k=0;k<18;k++) {printf("[%d]:%f  ",k,edep[k]);if(k==8)printf("\n");}
//        Int_t ithresh = nedep0 + TMath::Floor((18-nedep0)/2 + 0.5);
        Int_t ithresh = nedep0 + TMath::Floor((18-nedep0)*2/3.0 + 0.5);
        edepthresh = edep[ithresh-1];

        if(fIsMIPL[0]&&fIsMIPL[1]&&fIsMIPL[2]&&fIsMIPL[3]&&fIsMIPL[4]&&fIsMIPL[5]&&fIsMIPL[6]
                &&fIsMIPL[7]&&fIsMIPL[8]&&fIsMIPL[9]&&fIsMIPL[10]&&fIsMIPL[11]&&fIsMIPL[12]
                &&fIsMIPL[13]&&fIsMIPL[14]&&fIsMIPL[15]&&fIsMIPL[16]&&fIsMIPL[17]){

            Double_t edep_aver=0.0;
            Int_t nlayuse;

            for(int k=nedep0;k<18;k++)  edep_aver += edep[k];
            edep_aver = edep_aver / (18-nedep0);
            Fill(shn+"10", fRig, edep_aver/dx);
            Fill(shn+"20", fBetaH.beta, edep_aver/dx);
            if(fRich.beta!=0)Fill(shn+"30", fRich.beta, edep_aver/dx);
//printf("\n-------------No Truncation: Mean:%f    N:%d ------------------------\n",edep_aver,18-nedep0);

            nlayuse = TMath::Floor((18-nedep0)*2/3.0 + 0.5);
            for(int k=nedep0;k<nedep0+nlayuse;k++)  edep_aver += edep[k];
            edep_aver = edep_aver / nlayuse;
            Fill(shn+"11", fRig, edep_aver/dx);
            Fill(shn+"21", fBetaH.beta, edep_aver/dx);
            if(fRich.beta!=0)Fill(shn+"31", fRich.beta, edep_aver/dx);
//printf("-------------Truncation Mean (2/3):%f    N:%d ------------------------\n",edep_aver,nlayuse);

            nlayuse = TMath::Floor((18-nedep0)/2.0 + 0.5);
            for(int k=nedep0;k<nedep0+nlayuse;k++)  edep_aver += edep[k];
            edep_aver = edep_aver / nlayuse;
            Fill(shn+"12", fRig, edep_aver/dx);
            Fill(shn+"22", fBetaH.beta, edep_aver/dx);
            if(fRich.beta!=0)Fill(shn+"32", fRich.beta, edep_aver/dx);
//printf("-------------Truncation Mean (1/2):%f    N:%d ------------------------\n",edep_aver,nlayuse);

            //truncate both beginning and ending
            nlayuse = TMath::Floor((18-nedep0)/2.0 + 0.5) - 2;
            for(int k=nedep0+2;k<nedep0+2+nlayuse;k++)  edep_aver += edep[k];
            edep_aver = edep_aver / nlayuse;
            Fill(shn+"13", fRig, edep_aver/dx);
            Fill(shn+"23", fBetaH.beta, edep_aver/dx);
            if(fRich.beta!=0)Fill(shn+"33", fRich.beta, edep_aver/dx);


                                  //Fill(shn+"20",  fBetaH.beta,edep_aver/dx ,fMCW);
            if(1.0<fRig && fRig<1.1){
                Fill(shn+"201",  fBetaH.beta,edep_aver/dx );
                if(fRich.beta!=0)Fill(shn+"301",  fRich.beta,edep_aver/dx );
            }
            if(1.5<fRig && fRig<1.6){
                Fill(shn+"202",  fBetaH.beta,edep_aver/dx );
                if(fRich.beta!=0)Fill(shn+"302",  fRich.beta,edep_aver/dx );
            }
            if(2.0<fRig && fRig<2.1){
                Fill(shn+"203",  fBetaH.beta,edep_aver/dx ,fMCW);
                if(fRich.beta!=0)Fill(shn+"303",  fRich.beta,edep_aver/dx );
            }

            Double_t beta_estm;
            if(1.5<fRig && fRig<1.6){
                beta_estm = TMath::Power( 0.18 / (edep_aver/dx - 0.4832),1/4.0 );
                Fill(shn+"202e",  beta_estm,fBetaH.beta);
//printf("a==== fRig=%f   dE/dx=%f    [1]beta_estm=%f   ",fRig,edep_aver/dx, beta_estm);
            }
            if(2.0<fRig && fRig<2.1){
                beta_estm = TMath::Power( 0.0932 / (edep_aver/dx - 0.594), 1/5.0 );    // x^4:0.459   0.185
                Fill(shn+"203e",  beta_estm,fBetaH.beta);
//printf("b==== fRig=%f   dE/dx=%f    [1]beta_estm=%f   ",fRig,edep_aver/dx, beta_estm);
            }
//beta_estm = TMath::Power( 0.18 / (edep_aver/dx - 0.4832),1/4.0 );
//printf("==== fRig=%f   dE/dx=%f    [1]beta_estm=%f   ",fRig,edep_aver/dx, beta_estm);
//beta_estm = pow( 0.18 / (edep_aver/dx - 0.4832),1/4.0 );
//printf("[2] %f   \n",beta_estm);

            if(0.90<fRig && fRig<0.95)Fill(shn+"231",fBetaH.beta,edep_aver/dx);
            if(1.80<fRig && fRig<1.89)Fill(shn+"232",fBetaH.beta,edep_aver/dx);       

            Double_t charge = qsel;
            Double_t mass_tof = fRig * charge * TMath::Sqrt(1./fBetaH.beta/fBetaH.beta-1);

                                   Fill(shn+"40",fRig,mass_tof);
                                   Fill(shn+"41",fRig,1./mass_tof);
                                   Fill(shn+"42",fRig,mass_tof*mass_tof);
            if(1.0<fRig && fRig<1.1){
                Fill(shn+"401",fBetaH.beta,mass_tof);
                Fill(shn+"411",fBetaH.beta,1./mass_tof);
                Fill(shn+"421",fBetaH.beta,mass_tof*mass_tof);
            }
            if(1.5<fRig && fRig<1.6){
                Fill(shn+"402",fBetaH.beta,mass_tof);
                Fill(shn+"412",fBetaH.beta,1./mass_tof);
                Fill(shn+"422",fBetaH.beta,mass_tof*mass_tof);
            }
            if(2.0<fRig && fRig<2.1){
                Fill(shn+"403",fBetaH.beta,mass_tof);
                Fill(shn+"413",fBetaH.beta,1./mass_tof);
                Fill(shn+"423",fBetaH.beta,mass_tof*mass_tof);
            }
        }

        TString sthn;
//printf("\nedepthresh:%f  \n",edepthresh);
//printf("===================================================================================\n");
        bool IsMIPL=1;
        for(int k=0;k<18;k++){
            if(edepl[k] != 0 ){
                //truncation mean method
                if(edepl[k] > edepthresh) continue;
                if(k==0){ 
                    if(!fIsMIPL[1]) continue;}
                else{
                    for(int m=0;m<k;m++) IsMIPL = (IsMIPL && fIsMIPL[m]);
                    if(k!=17)IsMIPL = ( IsMIPL && fIsMIPL[k+1] );
                    if(!IsMIPL) continue;
                }
                sthn = Form("%s10%02d",shn.Data(),k);//
                Fill(sthn,   fRig, edepl[k]/dx);
                sthn = Form("%s20%02d",shn.Data(),k);//
                Fill(sthn,  fBetaH.beta,edepl[k]/dx );
                if(fRich.beta!=0){
                    sthn = Form("%s30%02d",shn.Data(),k);//
                    Fill(sthn,  fRich.beta,edepl[k]/dx );
                }
            }
        }
    }
    else if(qsel == 2){
//printf("===================================Helium==========================================\n");
        Double_t edep[18],edepl[18];
        Int_t nedep0 = 0;
        Double_t edepthresh = 0.0;

        for(int k=0;k<18;k++){
            edepl[k] = edep[k] = eshw->DepositedEnergyInMeVPerLayer(k);
            if(edepl[k] == 0 ) nedep0++;
//printf("[%d]:%f  ",k,edepl[k]);
//if(k==8)printf("\n");
        }

//printf("\n---#0:[%d]-----------ABOVE: before resort     BELOW:after resort----------- \n",nedep0);
        sort(edep,edep+18);
//for(int k=0;k<18;k++) {printf("[%d]:%f  ",k,edep[k]);if(k==8)printf("\n");}
//        Int_t ithresh = nedep0 + TMath::Floor((18-nedep0)/2 + 0.5);
        Int_t ithresh = nedep0 + TMath::Floor((18-nedep0)*2/3.0 + 0.5);
        edepthresh = edep[ithresh-1];


        if(fIsMIPL[0]&&fIsMIPL[1]&&fIsMIPL[2]&&fIsMIPL[3]&&fIsMIPL[4]&&fIsMIPL[5]&&fIsMIPL[6]
                &&fIsMIPL[7]&&fIsMIPL[8]&&fIsMIPL[9]&&fIsMIPL[10]&&fIsMIPL[11]&&fIsMIPL[12]
                &&fIsMIPL[13]&&fIsMIPL[14]&&fIsMIPL[15]&&fIsMIPL[16]&&fIsMIPL[17]){

            Double_t edep_aver=0.0;
            Int_t nlayuse;

            for(int k=nedep0;k<18;k++)  edep_aver += edep[k];
            edep_aver = edep_aver / (18-nedep0);
            Fill(shn+"10", fRig, edep_aver/dx);
            Fill(shn+"20", fBetaH.beta, edep_aver/dx);
            if(fRich.beta!=0)Fill(shn+"30", fRich.beta, edep_aver/dx);
//printf("\n-------------No Truncation: Mean:%f    N:%d ------------------------\n",edep_aver,18-nedep0);

            nlayuse = TMath::Floor((18-nedep0)*2/3.0 + 0.5);
            for(int k=nedep0;k<nedep0+nlayuse;k++)  edep_aver += edep[k];
            edep_aver = edep_aver / nlayuse;
            Fill(shn+"11", fRig, edep_aver/dx);
            Fill(shn+"21", fBetaH.beta, edep_aver/dx);
            if(fRich.beta!=0)Fill(shn+"31", fRich.beta, edep_aver/dx);
//printf("-------------Truncation Mean (2/3):%f    N:%d ------------------------\n",edep_aver,nlayuse);

            nlayuse = TMath::Floor((18-nedep0)/2.0 + 0.5);
            for(int k=nedep0;k<nedep0+nlayuse;k++)  edep_aver += edep[k];
            edep_aver = edep_aver / nlayuse;
            Fill(shn+"12", fRig, edep_aver/dx);
            Fill(shn+"22", fBetaH.beta, edep_aver/dx);
            if(fRich.beta!=0)Fill(shn+"32", fRich.beta, edep_aver/dx);
//printf("-------------Truncation Mean (1/2):%f    N:%d ------------------------\n",edep_aver,nlayuse);

            //truncate both beginning and ending
            nlayuse = TMath::Floor((18-nedep0)/2.0 + 0.5) - 2;
            for(int k=nedep0+2;k<nedep0+2+nlayuse;k++)  edep_aver += edep[k];
            edep_aver = edep_aver / nlayuse;
            Fill(shn+"13", fRig, edep_aver/dx);
            Fill(shn+"23", fBetaH.beta, edep_aver/dx);
            if(fRich.beta!=0)Fill(shn+"33", fRich.beta, edep_aver/dx);


                                  //Fill(shn+"20",  fBetaH.beta,edep_aver/dx ,fMCW);
            if(1.5<fRig && fRig<1.6){
                Fill(shn+"201",  fBetaH.beta,edep_aver/dx );
                if(fRich.beta!=0)Fill(shn+"301",  fRich.beta,edep_aver/dx );
            }
            if(2.0<fRig && fRig<2.1){
                Fill(shn+"202",  fBetaH.beta,edep_aver/dx );
                if(fRich.beta!=0)Fill(shn+"302",  fRich.beta,edep_aver/dx );
            }
            if(2.5<fRig && fRig<2.6){
                Fill(shn+"203",  fBetaH.beta,edep_aver/dx ,fMCW);
                if(fRich.beta!=0)Fill(shn+"303",  fRich.beta,edep_aver/dx );
            }

            Double_t beta_estm;
            if(1.5<fRig && fRig<1.6){
                beta_estm = TMath::Power( 0.18 / (edep_aver/dx - 0.4832),1/4.0 );
                Fill(shn+"202e",  beta_estm,fBetaH.beta);
//printf("a==== fRig=%f   dE/dx=%f    [1]beta_estm=%f   ",fRig,edep_aver/dx, beta_estm);
            }
            if(2.0<fRig && fRig<2.1){
                beta_estm = TMath::Power( 0.0932 / (edep_aver/dx - 0.594), 1/5.0 );    // x^4:0.459   0.185
                Fill(shn+"203e",  beta_estm,fBetaH.beta);
//printf("b==== fRig=%f   dE/dx=%f    [1]beta_estm=%f   ",fRig,edep_aver/dx, beta_estm);
            }
//beta_estm = TMath::Power( 0.18 / (edep_aver/dx - 0.4832),1/4.0 );
//printf("==== fRig=%f   dE/dx=%f    [1]beta_estm=%f   ",fRig,edep_aver/dx, beta_estm);
//beta_estm = pow( 0.18 / (edep_aver/dx - 0.4832),1/4.0 );
//printf("[2] %f   \n",beta_estm);

            if(1.51<fRig && fRig<1.62)Fill(shn+"231",fBetaH.beta,edep_aver/dx);
            if(2.01<fRig && fRig<2.15)Fill(shn+"232",fBetaH.beta,edep_aver/dx);       

            Double_t charge = qsel;
            Double_t mass_tof = fRig * charge * TMath::Sqrt(1./fBetaH.beta/fBetaH.beta-1);

                                   Fill(shn+"40",fRig,mass_tof);
                                   Fill(shn+"41",fRig,1./mass_tof);
                                   Fill(shn+"42",fRig,mass_tof*mass_tof);
            if(1.5<fRig && fRig<1.6){
                Fill(shn+"401",fBetaH.beta,mass_tof);
                Fill(shn+"411",fBetaH.beta,1./mass_tof);
                Fill(shn+"421",fBetaH.beta,mass_tof*mass_tof);
            }
            if(2.0<fRig && fRig<2.1){
                Fill(shn+"402",fBetaH.beta,mass_tof);
                Fill(shn+"412",fBetaH.beta,1./mass_tof);
                Fill(shn+"422",fBetaH.beta,mass_tof*mass_tof);
            }
            if(2.5<fRig && fRig<2.6){
                Fill(shn+"403",fBetaH.beta,mass_tof);
                Fill(shn+"413",fBetaH.beta,1./mass_tof);
                Fill(shn+"423",fBetaH.beta,mass_tof*mass_tof);
            }
        }

//printf("\nedepthresh:%f  \n",edepthresh);
//printf("===================================================================================\n");
        TString sthn;
        bool IsMIPL=1;
        for(int k=0;k<18;k++){
            //truncation mean method
            if(edepl[k] > edepthresh) continue;
            if(edepl[k] != 0 ){
                if(k==0){ 
                    if(!fIsMIPL[1]) continue;}
                else{
                    //if(eshw->DepositedEnergyInMeVPerLayer(k-1)< 30 ) continue;
                    for(int m=0;m<k;m++) IsMIPL = (IsMIPL && fIsMIPL[m]);
                    if(k!=17)IsMIPL = ( IsMIPL && fIsMIPL[k+1] );
                    if(!IsMIPL) continue;
                }
                sthn = Form("%s10%02d",shn.Data(),k);//
                Fill(sthn,   fRig, edepl[k]/dx);
                sthn = Form("%s20%02d",shn.Data(),k);//
                Fill(sthn,  fBetaH.beta,edepl[k]/dx );
                if(fRich.beta!=0){
                    sthn = Form("%s30%02d",shn.Data(),k);//
                    Fill(sthn,  fRich.beta,edepl[k]/dx );
                }
            }
        }
    }//qsel==2



    return 1;
}









int DSTFill::BasicSelection(){

    Int_t stsel = IsInL1 | HasL2 | TrdInL1 | TrdInTr ;
    if ((fStatus.ustat&stsel) != stsel) return 0;
    if (!(fHeader.phpat&0x3e)) return 0;

    Double_t qlm ;
    Double_t qin = fTrack.qin;
    Int_t qsel=0;
    if ( 0.7 < qin && qin < 1.4) {qsel=1;qlm=1;}
    if ( 1.7 < qin && qin < 2.4) {qsel=2;qlm=2;}
    if( !(qsel==1 || qsel==2) ) return 0; 

    if (fBetaH.beta <= 0.4 || fBeta.pattern > 4) return 0;
    if (fBetaH.clsn[0]+fBetaH.clsn[2]>4)   return 0;

    //At least one Y hit in all the inner planes
    if (    ! (            (fTrack.bith&(1<<1)) 
    && (fTrack.bith&(1<<2)||fTrack.bith&(1<<3))
    && (fTrack.bith&(1<<4)||fTrack.bith&(1<<5))
    && (fTrack.bith&(1<<6)||fTrack.bith&(1<<7))     )     )
        return 0;


    Double_t xl1 = fTrack.coox[0], yl1 = fTrack.cooy[0],
             xl9 = fTrack.coox[2], yl9 = fTrack.cooy[2];
    // LD's L1 fiducial
    if (xl1*xl1+yl1*yl1 > 62.14*62.14 ||
            TMath::Abs(yl1) > 47.4) return 0;
    //if(qsel==1 && fabs(xl9)<=33

    Double_t cin = fTrack.csqy[0];
    if (cin > 10 || cin < 0 ) return 0;

    Double_t qtu = (fBetaH.ql[0]+fBetaH.ql[1])/2.0;
    Double_t qtl = (fBetaH.ql[2]+fBetaH.ql[3])/2.0;
    if (fBetaH.ql[1] <= 0 && fBetaH.ql[0] > 0) qtu = fBetaH.ql[0];
    if (fBetaH.ql[0] <= 0 && fBetaH.ql[1] > 0) qtu = fBetaH.ql[1];
    if (fBetaH.ql[3] <= 0 && fBetaH.ql[2] > 0) qtl = fBetaH.ql[2];
    if (fBetaH.ql[2] <= 0 && fBetaH.ql[3] > 0) qtl = fBetaH.ql[3];
    if( qsel==1 && (qtu<0.5||qtu>2.5) ) return 0;//proton
    if( qsel==2 && (qtu<=1.25) )        return 0;//helium

    Double_t ql1 = fTrack.ql1;
    Double_t ql9 = fTrack.ql9;
    if(! ((qlm-0.4)<ql1 && ql1<(qlm+0.9) && (qlm-0.4)<ql9 && ql9<(qlm+0.9))  ) return 0;

    fRig = fTrack.rgt[1];
    if(!IsMC) fRig = AMSEventR::GetCorrectedRigidity(fRig,fHeader.utime ,0,1);
    if(fRig<=0 || fRig > 10) return 0;

//if(debug) printf("fRig:%f  \n",fRig);
    Double_t ra = TMath::Abs(fRig);
    if(1/fBetaH.beta < TMath::Sqrt(1+1.5*1.5/(ra+0.5)/(ra+0.5))-0.25) return 0;


    Double_t charge = qsel;
    Double_t massh;
    massh = fRig * charge * TMath::Sqrt(1./fBetaH.beta/fBetaH.beta-1);
    TString shhn;
    shhn = Form("hh%d11",qsel);
    Fill(shhn,fRig,massh);
    shhn = Form("hh%d110",qsel);
    Fill(shhn,fRig,massh);

    massh = fRig * charge * TMath::Sqrt(1./fRich.beta/fRich.beta-1);
    shhn = Form("hh%d12",qsel);
    Fill(shhn,fRig,massh);
    shhn = Form("hh%d120",qsel);
    Fill(shhn,fRig,massh);

    if(fTrdK.llr[0] + fTrdK.llr[1] + fTrdK.llr[2] > 0.){
        shhn = Form("hh%d210",qsel);
        if(qsel == 1)Fill(shhn,fRig,fTrdK.llr[0]);
        if(qsel == 2)Fill(shhn,fRig,fTrdK.llr[1]);
        shhn = Form("hh%d211",qsel);
        Fill(shhn,fRig,fTrdK.llr[2]);
    }
    if( fTrdK.llrt[0] + fTrdK.llrt[1] + fTrdK.llrt[2] > 0. ){
        shhn = Form("hh%d220",qsel);
        if(qsel == 1)  Fill(shhn,fRig,fTrdK.llrt[0]);
        if(qsel == 2)  Fill(shhn,fRig,fTrdK.llrt[1]);
        shhn = Form("hh%d221",qsel);
        Fill(shhn,fRig,fTrdK.llrt[2]);
    }
    if(fTrdK.llre[0] + fTrdK.llre[1] || fTrdK.llre[2] > 0.){
        shhn = Form("hh%d230",qsel);
        if(qsel == 1)Fill(shhn,fRig,fTrdK.llre[0]);
        if(qsel == 2)Fill(shhn,fRig,fTrdK.llre[1]);
        shhn = Form("hh%d231",qsel);
        Fill(shhn,fRig,fTrdK.llre[2]);
    }


if(qsel==1 && fTrdK.llr[0]<0.5 && fTrdK.llr[0]>0 ) return 0;
if(qsel==2 && fTrdK.llr[1]<0.7 && fTrdK.llr[1]>0 ) return 0;
if(qsel==1 && fEcal.bdt >-0.8) return 0;


    return qsel;
}



int DSTFill::TRDTOFRICHSelection(int qsel){
//TRD
    if (fTrdK.nhits <= 8) return 0;
    Double_t q_trd = fTrdK.q;
    if(qsel==1){  if( q_trd<0.7||q_trd>1.4 ) return 0;}
    if(qsel==2){  if( q_trd<1.7||q_trd>2.4 ) return 0;}

    int ntrdcluster = fpev->NTrdCluster();
    if(fpev->NTrdTrack() != 1) return 0 ;

    if(!fTrk  ) return 0;

//tof

//rich
  ParticleR *part = fpev->pParticle(0);
  if(!part) return 0;
  RichRingR *rich = part->pRichRing();
  if(!rich) return 0;

    Double_t q_rich = fRich.q;
    if(qsel==1){  if( q_rich<0.7||q_rich>1.4 ) return 0;}
    if(qsel==2){  if( q_rich<1.7||q_rich>2.4 ) return 0;}

    if(! rich->IsClean()) return 0;
    if(fRich.pmts<=3 ) return 0;//if(rich->getPMTs()<=3 ) return 0;
    //if(rich->getPMTChargeConsistency()<=10 )return 0;



    return 1;


}







int DSTFill::EventSelection(){
/*printf("fHeader.phpat: before:%X  ",fHeader.phpat);
GetEntry("header");
printf("fHeader.phpat: after:%X  ",fHeader.phpat);
printf("fBetaH.beta: before:%f  ",fBetaH.beta);
GetEntry("betah");
printf("fBetaH.beta: after:%f  \n",fBetaH.beta);
return 0;*/

    //GetEntry("header");
    if (!(fHeader.phpat&0x3e)) return 0;
//    Int_t stsel = HasL2 | IsInL1 | IsInL9;
    Int_t stsel = HasL2 | IsInL9;
    if ((fStatus.ustat&stsel) != stsel) return 0;

    Double_t rl9 = fTrack.rgt[2];
    if(!IsMC) rl9 = AMSEventR::GetCorrectedRigidity(rl9,fHeader.utime ,0,1);
    if(rl9<0 || rl9 > 10) return 0;

    GetEntry("betah");
    GetEntry("beta");
//mass cut reject electron and pion in low energy
    Double_t ra = TMath::Abs(rl9);
//printf("rl9=%f    Beta cut:%f     beta:%f  \n",rl9,1/( TMath::Sqrt(1+1.5*1.5/(ra+0.5)/(ra+0.5))-0.25  ), fBetaH.beta  );
    if(1/fBetaH.beta < TMath::Sqrt(1+1.5*1.5/(ra+0.5)/(ra+0.5))-0.25) return 0;

    Double_t qlm ;
    Double_t qin = fTrack.qin;
    Int_t qsel=0;
    if (qin < 0.7 || 1.4 < qin) {qsel=1;qlm=1;}
    if (qin < 1.7 || 2.4 < qin) {qsel=2;qlm=2;}
//use loose cuts here
//qsel = TMath::Floor(qin+0.5);
//qlm = qsel;

    if( !(qsel==1 || qsel==2) ) return 0; 
    if(qsel==2 && rl9>5) return 0;


    bool bcut1=1,bcut2=1,bcut3=1,bcut4=1,bcut5=1,bcut6=1,bcut7=1,bcut8=1;

    //GetEntry("betah");
    //GetEntry("beta");
//    if (fBetaH.chi2t > 10 || fBetaH.chi2c > 10)  bcut1=0;
    if (fBetaH.beta <= 0.3 || fBeta.pattern > 4) bcut1=0;
//    if (fBetaH.beta <= 0.1 || fBeta.pattern > 4) bcut1=0;
    if (fBetaH.type != 1) bcut1=0;//0502 added BetaHR->GetBuildType==1 
    if (fBetaH.clsn[0]+fBetaH.clsn[2]>4)   bcut1=0;

/*
    GetEntry("rti");
    Double_t rcut = (!IsMC) ? fRTI.cfi*1.2 : fMCinfo.rgt;
//    Double_t rfs = fTrack.rgt[3];
//    if(!IsMC) rfs = AMSEventR::GetCorrectedRigidity(rfs,fHeader.utime ,0,1);
//    if(rfs < rcut) bcut2=0;
    if(rl9 < rcut) bcut2=0;//11.12
*/
    if ( fTrack.qrms >= 0.4) bcut3=0;//11.12

    //At least one Y hit in all the inner planes
    if (    ! (            (fTrack.bith&(1<<1)) 
    && (fTrack.bith&(1<<2)||fTrack.bith&(1<<3))
    && (fTrack.bith&(1<<4)||fTrack.bith&(1<<5))
    && (fTrack.bith&(1<<6)||fTrack.bith&(1<<7))     )     )
        bcut4=0;

    Double_t xl1 = fTrack.coox[0], yl1 = fTrack.cooy[0],
             xl9 = fTrack.coox[2], yl9 = fTrack.cooy[2];
    // LD's L1 fiducial
    if (xl1*xl1+yl1*yl1 > 62.14*62.14 ||
            TMath::Abs(yl1) > 47.4) bcut5=0;
    if(qsel==1 && fabs(xl9)>=33 ) bcut5=0;   //IsInEcal

 if(!(bcut1 && bcut2 && bcut3 && bcut4 && bcut5) ) return 0;

    Double_t cin = fTrack.csqy[0];
    if (cin > 10 || cin < 0 ) bcut6=0;

    Double_t qtu = (fBetaH.ql[0]+fBetaH.ql[1])/2.0;
    Double_t qtl = (fBetaH.ql[2]+fBetaH.ql[3])/2.0;
    if (fBetaH.ql[1] <= 0 && fBetaH.ql[0] > 0) qtu = fBetaH.ql[0];
    if (fBetaH.ql[0] <= 0 && fBetaH.ql[1] > 0) qtu = fBetaH.ql[1];
    if (fBetaH.ql[3] <= 0 && fBetaH.ql[2] > 0) qtl = fBetaH.ql[2];
    if (fBetaH.ql[2] <= 0 && fBetaH.ql[3] > 0) qtl = fBetaH.ql[3];
    if( qsel==1 && (qtu<0.5||qtu>2.5) ) bcut7=0;//proton
    if( qsel==2 && (qtu<=1.25) )        bcut7=0;//helium

    Double_t ql1 = fTrack.ql1;
    Double_t ql9 = fTrack.ql9;
    bcut8 = ((qlm-0.4)<ql1 && ql1<(qlm+0.9) && (qlm-0.4)<ql9 && ql9<(qlm+0.9));

    if(!(bcut1 && bcut2 && bcut3 && bcut4 && bcut5 && bcut6 && bcut7 && bcut8) ) return 0;

    return qsel;
}










int DSTFill::ecalmap_MIP(AMSEventR *pev, int ic_ecal){

    Int_t ecal_nhit[18] = {0};
    float ecal_elay[18] = {0.0};
    int idxlay,idxcell;
    float edepcell;

    float ecalmap[18][72] = {0.0};

    EcalShowerR *pecal = pev->pEcalShower( ic_ecal );
    //if(pecal->Nhits)   //add total hits cuts

    for(int i2dcl=0; i2dcl<pecal->NEcal2DCluster(); i2dcl++)
    {
        Ecal2DClusterR *p2dcl = pecal->pEcal2DCluster(i2dcl);
        for(int i1dcl=0; i1dcl<p2dcl->NEcalCluster(); i1dcl++)
        {
            EcalClusterR *p1dcl = p2dcl->pEcalCluster(i1dcl);
            for(int ihit=0; ihit<p1dcl->NEcalHit(); ihit++)
            {
                EcalHitR *phit = p1dcl->pEcalHit(ihit);
                idxlay = phit->Plane;
                idxcell= phit->Cell;
                edepcell = phit->Edep;
                //ecalmap[ phit->Plane ][ phit->Cell ] = phit->Edep;
                ecalmap[ idxlay ][ idxcell ] = edepcell;
                ecal_nhit[idxlay]++;
                ecal_elay[idxlay] += edepcell;
            }   
        }   
    }


float er[18]={0.0};

    int NhitlIs2=0;
    int NhitlIs3=0;
    float eratio = 0.0;
    NIsMIPL = 0;
    for(int k=0;k<18;k++){
        if( ecal_nhit[k]>=4 ) {for(int j=k;j<18;j++) fIsMIPL[j]=false;  break;}
        if( ecal_nhit[k]==3 ) {
            NhitlIs3++;
            if( NhitlIs3>=2 ) {for(int j=k;j<18;j++) fIsMIPL[j]=false;  break;}
        }
        if( ecal_nhit[k]>=2 ) {
            NhitlIs2++;
            if( NhitlIs2>=4 ) {
                for(int j=k;j<18;j++) fIsMIPL[j]=false;
                if(pecal->Nhits < 25 && Nhitl>=15)   for(int j=k;j<18;j++) fIsMIPL[j]=true;
                else break;
            }
        }
        if(icell[k]<0) {for(int j=k;j<18;j++) fIsMIPL[j]=false;  break;}
        if(ecal_elay[k] >1e-4 ){
                idxcell = icell[k];
                if     (idxcell==0)  eratio = (ecalmap[k][idxcell]+ecalmap[k][idxcell+1] ) / ecal_elay[k];
                else if(idxcell==71) eratio = (ecalmap[k][idxcell]+ecalmap[k][idxcell-1] ) / ecal_elay[k];
                else                 eratio = (ecalmap[k][idxcell]+ecalmap[k][idxcell+1]+ecalmap[k][idxcell-1]) / ecal_elay[k];
                if(eratio < 0.9) fIsMIPL[k]=false;
        er[k] = eratio;
        }
    }

    for(int k=0;k<18;k++){ if(fIsMIPL[k]==true) NIsMIPL++; }


    return 1;
}



TrTrackR * DSTFill::GetTrack(void)
{
  ParticleR *p = fpev->pParticle(0); if (!p) return 0;
  if (p->pTrTrack()) return p->pTrTrack();

  BetaR *b = p->pBeta();
  if (b && b->Beta > 0.3 && b->pTrTrack()) return b->pTrTrack();

  BetaHR *bh = p->pBetaH();
  if (bh && bh->GetBeta() > 0.3 && bh->pTrTrack()) return bh->pTrTrack();

  VertexR *v = p->pVertex(); if (!v) return 0;
  TrTrackR *tr = 0;
  for (Int_t i = 0; i < v->NTrTrack(); i++) {
    TrTrackR *t = v->pTrTrack(i);
    if (!tr || TMath::Abs(tr->GetRigidity()) < TMath::Abs(t->GetRigidity()))
      tr = t;
  }
  return tr;
}




Int_t mdfil(TChain &ch, const char *oname, Int_t mode)
{
    DSTFill dst(&ch,oname,mode);
    return 0;
}

/* fname: Filelist
 * oname: Output Name
 *  mode: MC and Particle
 *   rid: Select Runs in Filelist
 */
Int_t mdfil(const char *fname, const char *oname, Int_t mode, Int_t rid){
//    AMSChain ca;
    ifstream fin(fname);
    if(!fin){
        cout<<"File list not found: "<<fname<<endl;
        return -1;
    }
    TString srm = "root://eosams.cern.ch//eos/ams/user/t/tracker/PG/mdst/p6/";
    TString stk = Form("mdst_%d", rid);
    TString stn = "tree";
    TString sfn = fname;

    if(sfn.Contains("pr822")){
        srm = "root://eosams.cern.ch//eos/ams/user/t/tracker/PG/mdst/pr822/";
    }
//helium
    if(sfn.Contains("He1081")){
        srm = "root://eosams.cern.ch//eos/ams/user/h/hehuang/HeliumFlux/Helium_MC_mdst/He.B1081_24000/";  
        stk = Form("mdst_mc_%d", rid);
    }
    if(sfn.Contains("He1036")){
        srm = "root://eosams.cern.ch//eos/ams/user/h/hehuang/HeliumFlux/Helium_MC_mdst/He.B1036_216000/";
        stk = Form("mdst_mc_%d", rid);
    }
    cout<<fname<<" "<<srm<<endl;
    TChain ch(stn);
    TString str;
    while(str.ReadLine(fin)){
        if(str.Contains(stk)){
            ch.Add(srm+str);
            cout<<"Added: "<<srm+str<<endl;
        }
    }

    return mdfil(ch, oname, mode);
}


#include <cstdlib>
int main(int argc, char *argv[]){
    if (argc < 4) {
        cout << "mdfil [fname] [oname] [mode]"       << endl;
        cout << "mdfil [fname] [oname] [mode] [rid]" << endl;
        cout << "mode: 3:Li +10: ismc" << endl;
        return 1;
    }

    Int_t ret = -1;
//    if (argc == 4) ret = mdfil(argv[1], argv[2], atoi(argv[3]));
    if (argc >= 5) ret = mdfil(argv[1], argv[2], atoi(argv[3]),atoi(argv[4]));
    if(ret==0)cout<<"Done"<<endl;
    return ret;
}


void DSTFill::Failed(const char *errinfo="0")
{
  cout << "Failed:" << errinfo << endl;
  exit(-1);
}




void DSTFill::NewRun(Int_t run)
{
  TString sfn = "list/list.ISS.B950.pass6";
if     (IsMC==1) sfn = "list/list.raw.MC.Pr.B1082";
else if(IsMC==2) sfn = "list/list.raw.MC.He.B1081";
//printf("Here 000000 =======================================\n");
  ifstream fin(sfn);
  if (!fin) {
    cout << "File not found: " << sfn.Data() << endl;
    if(!IsMC)Failed("list/list.ISS.B950.pass6 not found.");
    else if(IsMC==1)     Failed("list/list.raw.MC.Pr.B1082 not found.");
    else if(IsMC==2)     Failed("list/list.raw.MC.He.B1081 not found.");
  }

//printf("Here 000000 ===============IsMC:%d=====sfn:%s========================\n",IsMC,sfn.Data());
  if (fAch) delete fAch;
  fAch = new AMSChain;

  TString str;
  TString srn = Form("/%d", run);

  while (str.ReadLine(fin))
    if (str.Contains(srn)) {
      fAch->Add(str);
      cout << "Add: " << str.Data() << endl;
    }

  cout << "New RUN " << run << " "
       << "Ntr,Nent= " << fAch->GetNtrees() << " "
                       << fAch->GetEntries() << endl;
}










int DSTFill::read_ecalmap(AMSEventR *pev, int ic_ecal, int k)
{
    float ecalmap[18][72] = {0};
fEdep = 0;
Int_t NhitlIs2=0;
Int_t NhitlIs3=0;
Int_t Lmax;
//printf("---------------------------------------------------------------------k=%d \n",k);
    EcalShowerR *pecal = pev->pEcalShower( ic_ecal );
//printf("pecal->NEcal2DCluster()=%d   \n",pecal->NEcal2DCluster());    //1010 pecal->NEcal2DCluster()==2
if(pecal->NEcal2DCluster() > 2) return 0;
    //for(int i2dcl=0; i2dcl<pecal->NEcal2DCluster(); i2dcl++)
    for(int i2dcl=pecal->NEcal2DCluster()-1; i2dcl>=0; i2dcl--)
    {
        Lmax=0;
        Ecal2DClusterR *p2dcl = pecal->pEcal2DCluster(i2dcl);
//printf("NEcalCluster: %d:%d \n",i2dcl,p2dcl->NEcalCluster());
if(i2dcl==1) { if(p2dcl->NEcalCluster()>11) {fIsMIP=false; return 0;} }
if(i2dcl==0) { if(p2dcl->NEcalCluster()>9)  {fIsMIP=false; return 0;} }

//        printf("p2dcl->NEcalCluster()=%d  \n",p2dcl->NEcalCluster());
        for(int i1dcl=0; i1dcl<p2dcl->NEcalCluster(); i1dcl++)
        {
            EcalClusterR *p1dcl = p2dcl->pEcalCluster(i1dcl);
            //if(p2dcl->NEcalCluster()>9) printf("i1dcl:%d  fNhitL[%d]:%d  \n",i1dcl,i1dcl,fNhitL[i1dcl]);
//            printf("p1dcl->NEcalHit()=%d ",p1dcl->NEcalHit());
            
            if(p1dcl->NEcalHit() > 3 )  {fIsMIP = false;/* return 0;*/ continue;}
            if(p1dcl->NEcalHit() == 2 ) NhitlIs2 +=1;
            if(p1dcl->NEcalHit() == 3 ) NhitlIs3 +=1;
            if(NhitlIs3>1) return 0;
            /*    if(p1dcl->NEcalHit() == 2 ) {
                    NhitlIs2 +=1;
                    if(NhitlIs2 > 2) {fIsMIP = false; return 0;}
                }
            */


//1015_2
//============================================================================
    if( p1dcl->NEcalHit()>1 ) fIsMIPL[p1dcl->pEcalHit(0)->Plane] = false;
//    printf("      [%s]  phit->Plane:%d  \n",i2dcl==1?"Y":"X",p1dcl->pEcalHit(0)->Plane);

    if(p1dcl->pEcalHit(0)->Plane >= Lmax) Lmax=p1dcl->pEcalHit(0)->Plane;
    else fIsMIPL[p1dcl->pEcalHit(0)->Plane] = false;

            //printf("p1dcl->NEcalHit()=%d  fIsMIP:%d \n",p1dcl->NEcalHit(),fIsMIP);
            for(int ihit=0; ihit<p1dcl->NEcalHit(); ihit++)
            {
                //printf("[1] \n");
                EcalHitR *phit = p1dcl->pEcalHit(ihit);
                ecalmap[ phit->Plane ][ phit->Cell ] = phit->Edep;
                fEdep += phit->Edep;
            }   
        }   
    }

/*
printf("fIsMIPL 0~17:\n");
for(int j=0;    j<18;j++) printf("%d  ",fIsMIPL[j]);
printf("\n");
*/

    if(      Nhitl>=0  && Nhitl<=1 ) 
        {if( NhitlIs2 > 0 ) fIsMIP=false;}
    else if( Nhitl>=2  && Nhitl<=7 )
        {if( NhitlIs2 > 1 ) fIsMIP=false;}
    else if( Nhitl>=8  && Nhitl<=13)
        {if( NhitlIs2 > 2 ) fIsMIP=false;}
    else if( Nhitl>=14 && Nhitl<=17)
        {if( NhitlIs2 > 3 ) fIsMIP=false;}



    if(k!=-1 && k < 1000){
    //if(k!=-1 ){
    //if(k >1000 ){
    //k=k-1000;
        TString hn = Form("h_ecalmap_%d",k);
        TH2D* h_ecalmap = (TH2D *)fFile->Get(hn);
        for(int m = 0; m<18 ;m++){
            for(int n=0; n<72; n++){
                //h_ecalmap->Fill(m,n,ecalmap[m][n]);
                h_ecalmap->Fill(n,17-m,ecalmap[m][n]);
            }
        } 
    return 2;
//}

}
    //printf("==================================  X:%f  Y:%f     =========================\n",h_ecalmap->GetMean(1),h_ecalmap->GetMean(2));
    return 1;
}






