#define Ks_info_cxx

//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Oct  9 12:22:21 2015 by ROOT version 5.34/28
// from TTree Ks_info/track N-Tuple example
// found on file: KsKl_20000.root
//////////////////////////////////////////////////////////

#ifndef Ks_info_h
#define Ks_info_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <vector>
using std::vector;
class TH1D;

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class Ks_info {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Double_t        rec_truth_Mks;
   Double_t        rec_truth_Pks;
   Double_t        rec_truth_Eks;
   Int_t           indexmc;
   Int_t           pdgid[100];   //[indexmc]
   Int_t           motheridx[100];   //[indexmc]
   Double_t        m_chisq0;
   Double_t        m_chisqvtxnd;
   Double_t        m_chisq1c;
   Double_t        m_decayL;
   Double_t        m_decayLerr;
   Double_t        m_ctau;
   Double_t        Mpippim;
   Double_t        Ppippim;
   Double_t        Epippim;
   Double_t        Thepippim;
   Double_t        Phipippim;
   Double_t        Ppip;
   Double_t        Ppim;
   Double_t        Thepip;
   Double_t        Thepim;
   Double_t        Phipip;
   Double_t        Phipim;
   Double_t        pippx_mc;
   Double_t        pippy_mc;
   Double_t        pippz_mc;
   Double_t        pipe_mc;
   Double_t        pipp_mc;
   Double_t        pimpx_mc;
   Double_t        pimpy_mc;
   Double_t        pimpz_mc;
   Double_t        pime_mc;
   Double_t        pimp_mc;
   Double_t        pippx;
   Double_t        pippy;
   Double_t        pippz;
   Double_t        pipe;
   Double_t        pipp;
   Double_t        pimpx;
   Double_t        pimpy;
   Double_t        pimpz;
   Double_t        pime;
   Double_t        pimp;
   Int_t           ntof1;
   Int_t           toflayer1[5];   //[ntof1]
   Double_t        tof1[5];   //[ntof1]
   Double_t        t01[5];   //[ntof1]
   Int_t           ntof2;
   Int_t           toflayer2[5];   //[ntof2]
   Double_t        tof2[5];   //[ntof2]
   Double_t        t02[5];   //[ntof2]
   Int_t           mustatus;
   Int_t           emcstatus;
   Double_t        emctrk1;
   Double_t        emctrk2;
   Double_t        epratio1;
   Double_t        epratio2;

   // List of branches
   TBranch        *b_rec_truth_Mks;   //!
   TBranch        *b_rec_truth_Pks;   //!
   TBranch        *b_rec_truth_Eks;   //!
   TBranch        *b_indexmc;   //!
   TBranch        *b_pdgid;   //!
   TBranch        *b_motheridx;   //!
   TBranch        *b_m_chisq0;   //!
   TBranch        *b_m_chisqvtxnd;   //!
   TBranch        *b_m_chisq1c;   //!
   TBranch        *b_m_decayL;   //!
   TBranch        *b_m_decayLerr;   //!
   TBranch        *b_m_ctau;   //!
   TBranch        *b_Mpippim;   //!
   TBranch        *b_Ppippim;   //!
   TBranch        *b_Epippim;   //!
   TBranch        *b_Thepippim;   //!
   TBranch        *b_Phipippim;   //!
   TBranch        *b_Ppip;   //!
   TBranch        *b_Ppim;   //!
   TBranch        *b_Thepip;   //!
   TBranch        *b_Thepim;   //!
   TBranch        *b_Phipip;   //!
   TBranch        *b_Phipim;   //!
   TBranch        *b_pippx_mc;   //!
   TBranch        *b_pippy_mc;   //!
   TBranch        *b_pippz_mc;   //!
   TBranch        *b_pipe_mc;   //!
   TBranch        *b_pipp_mc;   //!
   TBranch        *b_pimpx_mc;   //!
   TBranch        *b_pimpy_mc;   //!
   TBranch        *b_pimpz_mc;   //!
   TBranch        *b_pime_mc;   //!
   TBranch        *b_pimp_mc;   //!
   TBranch        *b_pippx;   //!
   TBranch        *b_pippy;   //!
   TBranch        *b_pippz;   //!
   TBranch        *b_pipe;   //!
   TBranch        *b_pipp;   //!
   TBranch        *b_pimpx;   //!
   TBranch        *b_pimpy;   //!
   TBranch        *b_pimpz;   //!
   TBranch        *b_pime;   //!
   TBranch        *b_pimp;   //!
   TBranch        *b_ntof1;   //!
   TBranch        *b_toflayer1;   //!
   TBranch        *b_tof1;   //!
   TBranch        *b_t01;   //!
   TBranch        *b_ntof2;   //!
   TBranch        *b_toflayer2;   //!
   TBranch        *b_tof2;   //!
   TBranch        *b_t02;   //!
   TBranch        *b_mustatus;   //!
   TBranch        *b_emcstatus;   //!
   TBranch        *b_emctrk1;   //!
   TBranch        *b_emctrk2;   //!
   TBranch        *b_epratio1;   //!
   TBranch        *b_epratio2;   //!

   Ks_info(TTree *tree=0);
   virtual ~Ks_info();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(TFile* outfile, vector<TH1D*> h1col, const char* inname, vector<int> evtlist);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Ks_info_cxx
Ks_info::Ks_info(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("KsKl_20000.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("KsKl_20000.root");
      }
      f->GetObject("Ks_info",tree);

   }
   Init(tree);
}

Ks_info::~Ks_info()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Ks_info::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Ks_info::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Ks_info::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("rec_truth_Mks", &rec_truth_Mks, &b_rec_truth_Mks);
   fChain->SetBranchAddress("rec_truth_Pks", &rec_truth_Pks, &b_rec_truth_Pks);
   fChain->SetBranchAddress("rec_truth_Eks", &rec_truth_Eks, &b_rec_truth_Eks);
   fChain->SetBranchAddress("indexmc", &indexmc, &b_indexmc);
   fChain->SetBranchAddress("pdgid", pdgid, &b_pdgid);
   fChain->SetBranchAddress("motheridx", motheridx, &b_motheridx);
   fChain->SetBranchAddress("m_chisq0", &m_chisq0, &b_m_chisq0);
   fChain->SetBranchAddress("m_chisqvtxnd", &m_chisqvtxnd, &b_m_chisqvtxnd);
   fChain->SetBranchAddress("m_chisq1c", &m_chisq1c, &b_m_chisq1c);
   fChain->SetBranchAddress("m_decayL", &m_decayL, &b_m_decayL);
   fChain->SetBranchAddress("m_decayLerr", &m_decayLerr, &b_m_decayLerr);
   fChain->SetBranchAddress("m_ctau", &m_ctau, &b_m_ctau);
   fChain->SetBranchAddress("Mpippim", &Mpippim, &b_Mpippim);
   fChain->SetBranchAddress("Ppippim", &Ppippim, &b_Ppippim);
   fChain->SetBranchAddress("Epippim", &Epippim, &b_Epippim);
   fChain->SetBranchAddress("Thepippim", &Thepippim, &b_Thepippim);
   fChain->SetBranchAddress("Phipippim", &Phipippim, &b_Phipippim);
   fChain->SetBranchAddress("Ppip", &Ppip, &b_Ppip);
   fChain->SetBranchAddress("Ppim", &Ppim, &b_Ppim);
   fChain->SetBranchAddress("Thepip", &Thepip, &b_Thepip);
   fChain->SetBranchAddress("Thepim", &Thepim, &b_Thepim);
   fChain->SetBranchAddress("Phipip", &Phipip, &b_Phipip);
   fChain->SetBranchAddress("Phipim", &Phipim, &b_Phipim);
   fChain->SetBranchAddress("pippx_mc", &pippx_mc, &b_pippx_mc);
   fChain->SetBranchAddress("pippy_mc", &pippy_mc, &b_pippy_mc);
   fChain->SetBranchAddress("pippz_mc", &pippz_mc, &b_pippz_mc);
   fChain->SetBranchAddress("pipe_mc", &pipe_mc, &b_pipe_mc);
   fChain->SetBranchAddress("pipp_mc", &pipp_mc, &b_pipp_mc);
   fChain->SetBranchAddress("pimpx_mc", &pimpx_mc, &b_pimpx_mc);
   fChain->SetBranchAddress("pimpy_mc", &pimpy_mc, &b_pimpy_mc);
   fChain->SetBranchAddress("pimpz_mc", &pimpz_mc, &b_pimpz_mc);
   fChain->SetBranchAddress("pime_mc", &pime_mc, &b_pime_mc);
   fChain->SetBranchAddress("pimp_mc", &pimp_mc, &b_pimp_mc);
   fChain->SetBranchAddress("pippx", &pippx, &b_pippx);
   fChain->SetBranchAddress("pippy", &pippy, &b_pippy);
   fChain->SetBranchAddress("pippz", &pippz, &b_pippz);
   fChain->SetBranchAddress("pipe", &pipe, &b_pipe);
   fChain->SetBranchAddress("pipp", &pipp, &b_pipp);
   fChain->SetBranchAddress("pimpx", &pimpx, &b_pimpx);
   fChain->SetBranchAddress("pimpy", &pimpy, &b_pimpy);
   fChain->SetBranchAddress("pimpz", &pimpz, &b_pimpz);
   fChain->SetBranchAddress("pime", &pime, &b_pime);
   fChain->SetBranchAddress("pimp", &pimp, &b_pimp);
   fChain->SetBranchAddress("ntof1", &ntof1, &b_ntof1);
   fChain->SetBranchAddress("toflayer1", toflayer1, &b_toflayer1);
   fChain->SetBranchAddress("tof1", tof1, &b_tof1);
   fChain->SetBranchAddress("t01", t01, &b_t01);
   fChain->SetBranchAddress("ntof2", &ntof2, &b_ntof2);
   fChain->SetBranchAddress("toflayer2", toflayer2, &b_toflayer2);
   fChain->SetBranchAddress("tof2", tof2, &b_tof2);
   fChain->SetBranchAddress("t02", t02, &b_t02);
   fChain->SetBranchAddress("mustatus", &mustatus, &b_mustatus);
   fChain->SetBranchAddress("emcstatus", &emcstatus, &b_emcstatus);
   fChain->SetBranchAddress("emctrk1", &emctrk1, &b_emctrk1);
   fChain->SetBranchAddress("emctrk2", &emctrk2, &b_emctrk2);
   fChain->SetBranchAddress("epratio1", &epratio1, &b_epratio1);
   fChain->SetBranchAddress("epratio2", &epratio2, &b_epratio2);
   Notify();
}

Bool_t Ks_info::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Ks_info::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Ks_info::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Ks_info_cxx


#include <TH2.h>
#include <TH1D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <iostream>
#include <fstream>
using namespace std;

#include "CommonFunc.h"

void Ks_info::Loop(TFile* outfile, vector<TH1D*> h1col, const char* inname, vector<int> evtlist)
{
//   In a ROOT session, you can do:
//      Root > .L Ks_info.C
//      Root > Ks_info t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   
   double Ebeam = getEne(inname);
   cout<<"Beam energy is "<< Ebeam << "GeV"<< endl;

   char name[1000];
   sprintf(name,"hAngpip_%s",getPureName2(inname));
   TH1D *hAngpip = new TH1D(name,name,360,0,180);
   sprintf(name,"hAngpim_%s",getPureName2(inname));
   TH1D *hAngpim = new TH1D(name,name,360,0,180);
   
   sprintf(name,"hPpippim_%s",getPureName2(inname));
   TH1D *hPpippim = new TH1D(name,name,200,0,1.5*Ebeam/2.0);
   sprintf(name,"hMpippim_%s",getPureName2(inname));
   TH1D *hMpippim = new TH1D(name,name,200,0,1);
   sprintf(name,"hEpippim_%s",getPureName2(inname));
   TH1D *hEpippim = new TH1D(name,name,200,0,1.5*Ebeam/2.0);

   sprintf(name,"hAngInee_%s",getPureName2(inname));
   TH1D *hAngInee = new TH1D(name,name,360,0,180);
   sprintf(name,"hPmiss_%s",getPureName2(inname));
   TH1D *hPmiss = new TH1D(name,name,200,0,1.5*Ebeam/2.0);

   sprintf(name,"hAngInKs_%s",getPureName2(inname));
   TH1D *hAngInKs = new TH1D(name,name,360,0,180);
   sprintf(name,"hPpip_%s",getPureName2(inname));
   TH1D *hPpip = new TH1D(name,name,200,0,1.5*Ebeam/2.0);
   sprintf(name,"hPpim_%s",getPureName2(inname));
   TH1D *hPpim = new TH1D(name,name,200,0,1.5*Ebeam/2.0);

   sprintf(name,"hMKs_%s",getPureName2(inname));
   TH1D *hMKs = new TH1D(name,name,100,0.45,0.55);

   sprintf(name,"hchi2_%s",getPureName2(inname));
   TH1D *hchi2 = new TH1D(name,name,200, 0,200);

   sprintf(name,"hpKs_%s",getPureName2(inname));
   TH1D *hpKs = new TH1D(name,name,200,0,1.5*Ebeam/2.0);
   sprintf(name,"hpKsL_%s",getPureName2(inname));
   TH1D *hpKsL = new TH1D(name,name,200,0,1.5*Ebeam/2.0);
   sprintf(name,"hpKsU_%s",getPureName2(inname));
   TH1D *hpKsU = new TH1D(name,name,200,0,1.5*Ebeam/2.0);
   
   TLorentzVector pip, pim;
   double mpi = 0.13957;
   double mk0 = 0.497614;
   double pexp = sqrt(pow(Ebeam/2,2)-pow(mk0,2));
   cout<<"expected p Ks "<< pexp << endl;
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      hchi2->Fill(m_chisq1c);

      pip.SetVectMag(TVector3(pippx,pippy,pippz),mpi);
      pim.SetVectMag(TVector3(pimpx,pimpy,pimpz),mpi);
      double angpip = pip.Theta()/TMath::Pi()*180;
      double angpim = pim.Theta()/TMath::Pi()*180;

      // boost to ee cms
      pip.Boost(-0.011,0,0);
      pim.Boost(-0.011,0,0);
      double ang = pip.Vect().Angle(pim.Vect())/TMath::Pi()*180;
      double pmiss = (pip+pim).Rho();

      // boost to Ks cms
      TLorentzVector Ks = pip+pim;
      Ks.SetRho(pexp);
      Ks.SetVectMag(Ks.Vect(),mk0);
      TVector3 Vboost = Ks.BoostVector();
      pip.Boost(-Vboost);
      pim.Boost(-Vboost);
      double angInKs = pip.Vect().Angle(pim.Vect())/TMath::Pi()*180;
      double ppip = pip.Rho();
      double ppim = pim.Rho();

      hAngpip->Fill(angpip);
      hAngpim->Fill(angpim);

      hMpippim->Fill(Mpippim);
      hPpippim->Fill(Ppippim);
      hEpippim->Fill(Epippim);

      hAngInee->Fill(ang);
      hPmiss->Fill(pmiss);
      hAngInKs->Fill(angInKs);
      hPpip->Fill(ppip);
      hPpim->Fill(ppim);

      if (Ppippim<pexp+5*0.008 && Ppippim>pexp-5*0.008){
        if (m_decayL/m_decayLerr>5){
	  if (ang>40 && ang<120 && angInKs>160) {
	    hMKs->Fill(Mpippim);
	  }
	}
      }
      if (m_decayL/m_decayLerr<5) continue;
      if (ang<40 || ang>120 || angInKs<160) continue;
      if (Mpippim>0.490 && Mpippim<0.505) hpKs->Fill(Ppippim);
      if (Mpippim>0.475 && Mpippim<0.490) hpKsL->Fill(Ppippim);
      if (Mpippim>0.505 && Mpippim<0.520) hpKsU->Fill(Ppippim);
   }

   outfile->cd();
   double par[10];
   double parerr[10];
   double &Nsig = par[0];
   double &Nsigerr = parerr[0];
   int Npar = FitHist(hMKs,mk0,0,par,parerr,getPureName2(inname));
   ofstream pars("pars.txt",ios::app);
   pars << getPureName2(inname) << "\t" << Nsig << "\t" << Nsigerr << endl;

   outfile->WriteTObject(hAngpip);
   outfile->WriteTObject(hAngpim);
   outfile->WriteTObject(hPpippim);
   outfile->WriteTObject(hMpippim);
   outfile->WriteTObject(hEpippim);
   outfile->WriteTObject(hAngInee);
   outfile->WriteTObject(hPmiss);
   outfile->WriteTObject(hAngInKs);
   outfile->WriteTObject(hPpip);
   outfile->WriteTObject(hPpim);
    
   outfile->WriteTObject(hchi2);
   outfile->WriteTObject(hpKs);
   outfile->WriteTObject(hpKsL);
   outfile->WriteTObject(hpKsU);
   
   h1col.push_back(hAngpip);
   h1col.push_back(hAngpim);
   h1col.push_back(hPpippim);
   h1col.push_back(hMpippim);
   h1col.push_back(hEpippim);
   h1col.push_back(hAngInee);
   h1col.push_back(hPmiss);
   h1col.push_back(hAngInKs);
   h1col.push_back(hPpip);
   h1col.push_back(hPpim);
   
   outfile->WriteTObject(hMKs);

}


int main(int argc, char** argv)
{
    if (argc<2) {
      cout<<"No input file!"<<endl;
      return -1;
    }

    TFile *outfile = new TFile("output.root","recreate");
    for (int datai=1; datai<argc; datai++){
      char *filename = argv[datai];
      TFile *infile = new TFile(filename);
      if (infile==0) continue;
      TTree *tree = (TTree*)infile->Get("Ks_info");
      if (tree==0) continue;
      Ks_info data(tree);
      
      vector<TH1D*> histcol;
      vector<int> evtlist;
      data.Loop(outfile,histcol,filename,evtlist);
    }
    return 0;
}

