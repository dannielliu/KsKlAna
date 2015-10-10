#ifndef CommonFunc_h
#define CommonFunc_h
const int MAXSTRSIZE=10000;

double GetEnergy(int run)
{
  if (run>=39335 && run<=39618) return 3.08;
  if (run>=39711 && run<=39738) return 3.02;
  if (run>=39680 && run<=39710) return 3.00;
  if (run>=39651 && run<=39679) return 2.981;
  if (run>=39619 && run<=39650) return 2.95;
  
  if (run>=39775 && run<=40069) return 2.90;
  if (run>=40128 && run<=40296) return 2.6444;
  if (run>=40300 && run<=40435) return 2.6464;
  if (run>=40436 && run<=40439) return 2.70;
  if (run>=40440 && run<=40443) return 2.80;
  if (run>=40459 && run<=40769) return 2.396;
  if (run>=40771 && run<=40776) return 2.5;
  if (run>=40777 && run<=40804) return 2.6444;//separated beam
  if (run>=40806 && run<=40951) return 2.3864;
  if (run>=40989 && run<=41121) return 2.2;
  if (run>=41122 && run<=41239) return 2.2324;
  if (run>=41240 && run<=41411) return 2.3094;
  if (run>=41416 && run<=41532) return 2.175;
  if (run>=41533 && run<=41570) return 2.15;
  if (run>=41588 && run<=41727) return 2.1;

  if (run>=41729 && run<=41909) return 2.0;
  if (run>=41911 && run<=41958) return 2.05;
  if (run>=41959 && run<=41999) return 2.2324; // separated beam
  return -1;
}

int GetStringLength(const char* str, const char* flag=0)
{
  if (str!=0){
    int i=0;
    for (; i<MAXSTRSIZE; i++){
      if (flag!=0 && str[i]==flag[0]) return i;
      if (str[i]=='\0') return i;
      else continue;
    }
    if (i==MAXSTRSIZE) return -2;
  }
  return -1;
}

bool CompareStr(const char* str1, const char* str2, const char* flag=0)
{
  int size1 = GetStringLength(str1,flag);
  int size2 = GetStringLength(str2,flag);
  if (size1!=size2) {
    std::cout<<"Waring: string sizes are not equal"<<std::endl;
    return 0;
  }
  if (size1<0) return false;
  if (strncmp(str1,str2,size1)==0) return true;
  else return false;
}


const char* getPureName(const char* name)
{
  int pos=0;
  for (int i=0;;i++){
    if (name[i]=='\0') return &name[pos];
    if (name[i]=='/') pos = i+1;
    if (i>10000) return NULL;
  } 
}

const char* getPureName2(const char* name)
{
  char *name1 = new char[1000];
  int pos=0;
  int ppos=0;

  for (int i=0;;i++){
    if (name[i]=='\0') break; //return &name[pos];
    if (name[i]=='/') pos = i+1;
    if (name[i]=='.') ppos = i;
    if (i>1000) return NULL;
  }
  //std::cout<<"/ pos is "<< pos <<", . pos is "<< ppos <<std::endl;
  for (int i=0;i<ppos-pos;i++){
    name1[i] = name[pos+i];
  }
  name1[ppos-pos]='\0';

  return name1;
}


double getEne(const char* name)
{
  name = getPureName(name);
  std::string str = name;
  int _pos=0;
  int ppos=0;
  for (int i=0;;i++){
    if (name[i]=='\0') break;
    if (name[i]=='_') _pos=i;
    if (name[i]=='.') ppos=i;
  }
  int len = ppos-_pos-1;
  string vstr=str.substr(_pos+1,len);
  if (atof(vstr.c_str())/10000>1) return atof(vstr.c_str())/10000;

  else {// find index number in _*_
    int _pos2[2]={0,0};
    int j=0;
    for (int i=0;;i++){
      if (name[i]=='_') {_pos2[j]=i; j++;}
      if (j>=2) break;
    }
    if (j!=2) return -1;
    string vstr = str.substr(_pos2[0]+1,3);
    int idx = atoi(vstr.c_str());
    std::cout<<vstr<<'\t'<<idx<<std::endl;
    switch (idx) {
      case 1:    return  2.0;   
      case 2:    return  2.05;  
      case 3:    return  2.1;   
      case 4:    return  2.15;  
      case 5:    return  2.175; 
                                
      case 6:    return  2.2;   
      case 7:    return  2.2324;
      case 8:    return  2.3094;
      case 9:    return  2.3864;
      case 10:   return  2.396; 
      case 11:   return  2.5;   
      case 12:   return  2.6444;
      case 13:   return  2.6464;
      case 14:   return  2.7;   
      case 15:   return  2.8;   
      case 16:   return  2.9;   
      case 17:   return  2.95;  
      case 18:   return  2.981; 
      case 19:   return  3.0;   
      case 20:   return  3.02;  
                                
      case 21:   return  3.08;  
    }
  }
  return -2;
}

#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TF1.h>
#include "TPaveText.h"
#include "TMath.h"
#include "RooFit.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooPolynomial.h"
#include "RooCBShape.h"
#include "RooGenericPdf.h"
#include "RooChebychev.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooPlot.h"
#include "RooMsgService.h"
using namespace RooFit;

int FitHist(TH1D *&hist, double mean=0, double sigma=0, double* par=0, double *pare=0, const char* sfx=0)
{
    if (sfx==0) sfx="nosfx";
    double minX = hist->GetXaxis()->GetXmin();
    double maxX = hist->GetXaxis()->GetXmax();
    if (mean==0) mean = (minX+maxX)/2.0;
    if (sigma==0) sigma = (maxX-minX)/10.0;
    RooRealVar x("x","mass",mean,minX,maxX,"GeV");
    RooRealVar R_mean("R_mean","mean of gaussian",mean,minX,maxX);
    RooRealVar R_sigma("R_sigma","width of gaussian",sigma,sigma/10,sigma*4);
    RooGaussian F_sig("F_sig","gauss(x,m,s)",x,R_mean,R_sigma);
    
    RooRealVar co1("co1","coefficient #1",   0 ,-100.,100.);
    RooRealVar co2("co2","coefficient #1",   0 ,-100.,100.);
    RooChebychev F_bck("F_bck","background",x,RooArgList(co1,co2));
    
    RooRealVar R_sig("R_sig"," ",1000,0,100000000);//event number
    RooRealVar R_bck("R_bck"," ",1000,0,100000000);

    RooAddPdf sum("sum","sum",RooArgList(F_sig,F_bck),RooArgList(R_sig,R_bck));
    RooPlot *xframe = x.frame();;
    int Npar = 6;
    RooDataHist datahist("datahist","data",x,hist);
    
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit
    TCanvas *c1 = new TCanvas("","",800,600);
    
    sum.fitTo(datahist);
    datahist.plotOn(xframe);
    sum.plotOn(xframe,Components(F_sig),LineStyle(2),LineColor(2));
    sum.plotOn(xframe,Components(F_bck),LineStyle(2),LineColor(3));
    sum.plotOn(xframe);
    xframe->Draw();

    TPaveText *pt = new TPaveText(0.60,0.5,0.90,0.90,"BRNDC");
  //pt->SetBorderSize(0);
  //pt->SetFillStyle(4000);
  //pt->SetTextAlign(12);
  //pt->SetTextFont(42);
  //pt->SetTextSize(0.035);
    char item[100];
    sprintf(item,"#mu = %1.6f #pm %1.6f",R_mean.getVal(),R_mean.getError());
    pt->AddText(item);
    sprintf(item,"#sigma = %1.6f #pm %1.6f",R_sigma.getVal(),R_sigma.getError());
    pt->AddText(item);
    sprintf(item,"N_{sig} = %.2f #pm %.2f",R_sig.getVal(),R_sig.getError());
    pt->AddText(item);
    sprintf(item,"#chi^{2}/%d = %5.3f",Npar,xframe->chiSquare(Npar));
    pt->AddText(item);
    pt->Draw();
    sprintf(item,"spectrum_%s",sfx);
    c1->SetName(item);
    c1->Write();
    
    delete xframe;
    return Npar;
}


#endif
