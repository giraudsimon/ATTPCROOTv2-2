/*******************************************************************
* Class for Monte Carlo Minimization                               *
* Log: Class started 29-07-2016                                    *
* Author: Y. Ayyad and W. Mittig (NSCL)                            *
********************************************************************/

#ifndef ATMCQMINIMIZATION_H
#define ATMCQMINIMIZATION_H

//ATTPCROOT
#include "ATMinimization.hh"
#include "ATHit.hh"
#include "ATDigiPar.hh"

#include "AtTpcMap.h"

// FairRoot classes
#include "FairRuntimeDb.h"
#include "FairRun.h"

//ROOT
#include "TRotation.h"
#include "TMatrixD.h"
#include "TArrayD.h"
#include "TVector3.h"
#include "TH2Poly.h"

//System
#include <iostream>  //wm
#include <fstream>   // file stream wm
#include <cstdlib>  //wm
#include <vector>  //wm
using std::vector;  //wm
using namespace std; //wm


#define cRED "\033[1;31m"
#define cYELLOW "\033[1;33m"
#define cNORMAL "\033[0m"
#define cGREEN "\033[1;32m"

class ATMCQMinimization : public ATMinimization{

      public:
	       ATMCQMinimization();
        ~ATMCQMinimization();

        void AddELossFunc(std::function<Double_t(Double_t,std::vector<Double_t>&)>& func);
        void AddRtoEFunc(std::function<Double_t(Double_t,std::vector<Double_t>&)>& func);
        void AddELossPar(std::vector<Double_t> (&par)[10]);
        void AddRtoEPar(std::vector<Double_t> (&par)[10]);
        void AddParticle(std::vector<std::pair<Int_t,Int_t>>& ptcl);
        void SetEntTB(Int_t value);


        Bool_t Minimize(Double_t* parameter,ATEvent *event);
        Bool_t MinimizeOpt(Double_t* parameter,ATEvent *event);
        Bool_t MinimizeOptMap(Double_t* parameter,ATEvent *event, TH2Poly* hPadPlane);
        Bool_t MinimizeOptMapAmp(Double_t* parameter,ATEvent *event, TH2Poly* hPadPlane,const multiarray& PadCoord);

        template <typename T, typename R>
        bool  MinimizeGen(Double_t* parameter,const T* event ,const std::function<std::vector<R>*()>& func,TH2Poly* hPadPlane,const multiarray& PadCoord);

        std::vector<Double_t> GetPosXMin();
        std::vector<Double_t> GetPosYMin();
        std::vector<Double_t> GetPosZMin();
        std::vector<Double_t> GetPosXExp();
        std::vector<Double_t> GetPosYExp();
        std::vector<Double_t> GetPosZExp();
        std::vector<Double_t> GetPosXInt();
        std::vector<Double_t> GetPosYInt();
        std::vector<Double_t> GetPosZInt();
        std::vector<Double_t> GetPosXBack();
        std::vector<Double_t> GetPosYBack();
        std::vector<Double_t> GetPosZBack();
        Int_t GetMinimization();
        std::vector<ATHit> GetTBHitArray(Int_t TB,std::vector<ATHit> *harray);
        std::vector<std::function<Double_t(Double_t,std::vector<Double_t>&)>> *GetELossFunctionArray();

        void ResetParameters();




      protected:
       void GetEnergy(Double_t M,Double_t IZ,Double_t BRO,Double_t &E);
       void GetBro(Double_t M,Double_t IZ,Double_t &BRO,Double_t E);
       Double_t GetSimThetaAngle(TVector3* pos, TVector3* posforw);
       void BackwardExtrapolation();
       void SetMap(AtTpcMap* map);

       Int_t GetTBHit(Int_t TB,std::vector<ATHit> *harray);
       TVector3 TransformIniPos(Double_t x,Double_t y, Double_t z); //Transforms initial position from Pad plane to Lab frame
       TVector3 InvTransIniPos(Double_t x,Double_t y, Double_t z); //Transforms lab frame to pad plane

       void PrintParameters(Int_t index);

       // New MC Amplitude functions
          void MCvar( double* parameter, int & modevar,int & iconvar,double & x0MC, double & y0MC, double & z0MC,
                      double & aMC, double & phiMC, double & Bmin, double & dens, double & romin,
                      double & x0MCv,  double & y0MCv, double & z0MCv, double & aMCv, double & phiMCv, double & Bminv,
                      double & densv, double & rominv);

          void QMCsim(double* parameter, double* Qsim,double *zsimq,double & QMCtotal,
                                    double x0MC,double y0MC, double z0MC, double phiMCv,double aMCv,
                                    double Bminv, double densv, double rominv,double & e0sm, multiarray PadCoord,TH2Poly *padplane);

          void Chi2MC(double  Qtrack[10000],double  ztrackq[10000],double & Qtracktotal,
                                    double  Qsim[10000],double  zsimq[10000],double & QMCtotal,
                                    double & Chi2fit, double & sigmaq, double & sigmaz);

              std::vector<ATHit> fHitTBArray;
             std::vector<ATHit> *fHitArray;

             Double_t fThetaMin;
             Double_t fEnerMin;
             TVector3 fPosMin;
             Double_t fBrhoMin;
             Double_t fBMin;
             Double_t fPhiMin;
             Double_t fDensMin;
             Double_t fVertexEner;
             TVector3 fVertexPos;
             Double_t fDens;
             Double_t fPressure;


             std::vector<Double_t> fPosXmin;
             std::vector<Double_t> fPosYmin;
             std::vector<Double_t> fPosZmin;
             std::vector<Int_t>    fPosTBmin;
             std::vector<Double_t> fPosXexp;
             std::vector<Double_t> fPosYexp;
             std::vector<Double_t> fPosZexp;
             std::vector<Double_t> fPosXinter;
             std::vector<Double_t> fPosYinter;
             std::vector<Double_t> fPosZinter;
             std::vector<Int_t>    fPosTBinter;

             std::vector<Double_t> fPosXBack;
             std::vector<Double_t> fPosYBack;
             std::vector<Double_t> fPosZBack;

             std::vector<Double_t> fPosTBexp;

             Double_t fDriftVelocity;
             Int_t fTBTime;
             Double_t fBField;
             Double_t fTiltAng;
             Double_t fThetaPad;
             Double_t fThetaLorentz;
             Double_t fThetaRot;
             Double_t fZk;
             Int_t fEntTB; //Beam entrance Time Bucket

             //TRotation* fPadtoDetRot;

             //!AtTpcMap *fAtMapPtr;
             //!TH2Poly* fPadPlane;

             std::vector<std::function<Double_t(Double_t,std::vector<Double_t>&)>> fEloss_func_array; //!
             std::vector<std::function<Double_t(Double_t,std::vector<Double_t>&)>> fRtoE_func_array; //!
             std::vector<Double_t> fELossPar_array[10];
             std::vector<Double_t> fRtoEPar_array[10];
             std::vector<std::pair<Int_t,Int_t>> fParticleAZ;

             Bool_t kDebug;
             Bool_t kVerbose;

             //Global variables
             Double_t sm1;
             Double_t m;
             Double_t dzstep;
             Int_t    integrationsteps;
             Double_t restmass;
             Double_t esm;
             Double_t iz1;
             Double_t z1;
             Double_t B0;
             Double_t B;




           ClassDef(ATMCQMinimization, 1);

};

template <typename T, typename R>
bool  ATMCQMinimization::MinimizeGen(Double_t* parameter,const T* event,const std::function<std::vector<R>*()>& func,TH2Poly* hPadPlane,const multiarray& PadCoord)
{

          TH2Poly* fPadPlane = new TH2Poly();
          fPadPlane = hPadPlane;

          std::cout<<cGREEN<<" ============================"<<std::endl;
          std::cout<<" Starting Monte Carlo event  "<<cNORMAL<<std::endl;

          //TODO:: Enable indexing to select the function and the parameters depending on the particle

           if(fParticleAZ.size()>0){
                 m                  = fParticleAZ.at(0).first;
                 iz1                = fParticleAZ.at(0).second;
                 sm1                = m;
                 z1                 = iz1;
                 restmass           = sm1*931.49432;
                 esm                = z1*1.75879e-3*0.510998918/restmass;// ![e/m electron cm**2/(Volt*nsec**2] this is not the energy/mass but charge/mass
                 std::cout<<cGREEN<<" Particle  A : "<<m<<"  -  Z : "<<iz1<<cNORMAL<<std::endl;
            }else std::cerr<<cRED<<" ATMCQMinimization::MinimizeOptMapAmp -  Warning ! Particle (A,Z) not found. Using A : "<<m<<"  -  Z : "<<iz1<<cNORMAL<<std::endl;

          if(fELossPar_array[0].size()>0 && fRtoEPar_array[0].size()>0)
            PrintParameters(0);
          else std::cerr<<cRED<<" ATMCQMinimization::MinimizeOptMapAmp - Function parameters not found ! "<<cNORMAL<<std::endl;




          double Qtrack[10240]={0.}; //simulated amplitude for track to analyse
          double ztrackq[10240]={0.}; //simulated amplitude for track to analyse *ztrack fo find center of gravity
          std::vector<R>* HitArray = func();



          for(Int_t i=0;i<HitArray->size();i++){
             ATHit hit = HitArray->at(i);
             TVector3 position = hit.GetPosition();
             Int_t hitTB = hit.GetTimeStamp();
             Int_t hitPadNum = hit.GetHitPadNum();
             Int_t hitcharge = hit.GetCharge();
             //std::cout<<cGREEN<<" Hit number : "<<i<<" X : "<<position.X()<<" Y : "<<position.Y()<<" Z : "<<position.Z()<<" TB : "<<hitTB<<" Pad Number : "<<hitPadNum<<" charge "<<hitcharge<<cNORMAL<<std::endl;
             double xhelp = position.X();
             double yhelp = position.Y();
             Int_t iplot= 10*sqrt(hitcharge/500.);
             int ipl;
             if (iplot>10)iplot=10;
             //for (ipl=0;ipl<(iplot);++ipl){
            //   fout1<<"  "<<yhelp;
            // }
            // fout1<< std::endl;
             Qtrack[hitPadNum] = hitcharge;
             ztrackq[hitPadNum] = position.Z() ; // values to be used in chi2

             //For visualization
             fPosXexp.push_back(position.X());
             fPosYexp.push_back(position.Y());
             fPosTBexp.push_back(hit.GetTimeStamp());
             fPosXinter.push_back(position.X());
             fPosYinter.push_back(position.Y());
             fPosZinter.push_back(position.Z());
             fPosTBinter.push_back(hit.GetTimeStamp());


           }


           Double_t ymax            = 0.0;
           Double_t e0sm            = 0.0;
           Double_t smprot          = 931.49432;
           Double_t Bmin            = B;
           Double_t chi2min         = 1E10;

           Double_t xmin;
           Double_t ymin;
           Double_t zmin;  //Micromegas Origin  at 1000 mm of the entrance
           Double_t TBmin; // Absolute TB to compare between exp and sim
           Double_t phi0;
           Double_t bro;// !Tm*/
           Double_t theta0;
           ////////////////////////////////////////////
           double x0MC  = xmin ; //mm
           double y0MC  = ymin;
           double z0MC  = zmin;
           double aMC   = theta0; //radians
           double phiMC = phi0;

           double romin;
           double Bminv;//B //after MC variation
           double densv  ; //gas density after MC variation
           double rominv; //magnetic radius after MC variation
           double xpad[10240]={0}; //these are the partial integral results in small steps
           double ypad[10240]={0};
           double zpad[10240]={0};
           double Qpad[10240]={0}; //wm
           double Qsim[10240]={0.}; //simulated amplitude
           double zsimq[10240]={0.};
           double Qtracktotal;
           double QMCtotal;
           double CHi2fit;
           int modevar=1;
           double rangeMC= 200.;  //simulation of track
           int iconvar;
           double x0MCv= x0MC ;
           double y0MCv =y0MC ;
           double z0MCv = z0MC;
           double aMCv = aMC;
           double phiMCv = phiMC;
           double rangeMCv = rangeMC ;  //simulation of track values for MC variation
           double sigmaq=0.2 ;  //defined as the fraction of sum Qsim+Qtrack
           double sigmaz=4.0  ;  //defined as deviation of the center of gravity in mm modified from 5.4 on june 10 wm
           int imc1=0;
           int imc1max=10;//10
           iconvar=imc1;
           int imc2=0;
           int imc2max=50;//100
           int icontrol=1;

           MCvar(parameter, icontrol,iconvar,x0MC, y0MC, z0MC,aMC,phiMC, Bmin, fDens,romin,x0MCv, y0MCv,z0MCv,aMCv, phiMCv, Bminv, densv, rominv); // for initialisation

               std::cout<<std::endl;
               std::cout<<cGREEN<<" X : "<<x0MC<<" cm  - Y : "<<y0MC<<" cm - Z : "<<z0MC<<" cm "<<std::endl;
               std::cout<<" Brho : "<<(Bmin*romin)<<" Tm "<<std::endl;
               std::cout<<" Magnetic field : "<<Bmin<<" T "<<std::endl;
               std::cout<<" Radius of curvature : "<<parameter[5]<<" mm "<<std::endl;
               std::cout<<" Scattering Angle : "<<aMC*180.0/TMath::Pi()<<" deg "<<std::endl;
               std::cout<<" Azimutal Angle : "<<phiMC*180.0/TMath::Pi()<<" deg "<<std::endl;
               std::cout<<" Lenght of the experimental data : "<<parameter[7]<<cNORMAL<<std::endl;

               double Chimin=10000000. ;




               for (imc1=0;imc1<imc1max;imc1++){
                 iconvar = imc1; //this controls the step of MC in MCvar
                   for (imc2=0;imc2<imc2max;imc2++){
                     icontrol=2;
                       //fPadPlane->Reset(0);
                       MCvar(parameter,icontrol,iconvar,x0MC, y0MC, z0MC, aMC,phiMC, Bmin, fDens,romin,x0MCv, y0MCv,z0MCv, aMCv, phiMCv, Bminv, densv, rominv); // for MC variation with same starting value as before
                       QMCsim(parameter, Qsim, zsimq, QMCtotal,x0MCv, y0MCv, z0MCv,  aMCv, phiMCv, Bminv, densv,rominv,e0sm,PadCoord,fPadPlane);
                      //std::cout<<cRED<<" After QMCsim x "<<x0MCv<<" y "<< y0MCv<< " z "<< z0MCv<<" theta "<< aMCv<<" phi "<< phiMCv<<" B " <<Bminv<<" dens "<< densv<< " e0sm "<<e0sm<<" ro "<<rominv<<cNORMAL;
                       Chi2MC(Qtrack,ztrackq,Qtracktotal,Qsim,zsimq,QMCtotal,CHi2fit,sigmaq,sigmaz);   //Chi2 to compare track and MC



                     if(CHi2fit<Chimin){

                         Chimin=CHi2fit;

                           icontrol=3;
                           MCvar(parameter,icontrol,iconvar,x0MC, y0MC, z0MC, aMC,phiMC, Bmin,fDens,romin,x0MCv, y0MCv,z0MCv, aMCv, phiMCv, Bminv, densv, rominv);

                     }


                   }//imc2 loop
               }//imc1 loop

               //fPadPlane->Draw("zcol");

               QMCsim(parameter, Qsim, zsimq, QMCtotal,x0MC, y0MC, z0MC,  aMC, phiMC, Bmin,fDens,romin,e0sm,PadCoord,fPadPlane); // simulation with Chimin parameters

               FitParameters.sThetaMin = aMC;
               FitParameters.sEnerMin=e0sm;
               FitParameters.sPosMin.SetXYZ(x0MC,y0MC,z0MC);
               FitParameters.sBrhoMin=Bmin*romin;
               FitParameters.sBMin=Bmin;
               FitParameters.sPhiMin=phiMC;


               //std::cout<<cRED<<" final fit x "<<x0MC <<" y "<< y0MC << " z "<< z0MC <<" theta "<< aMC <<" phi "<< phiMC <<" B " <<Bmin <<" dens "<< dens <<" e0sm "<<e0sm<<" ro "<<romin <<cNORMAL<<std::endl;
               Chi2MC(Qtrack,ztrackq,Qtracktotal,Qsim,zsimq,QMCtotal,CHi2fit,sigmaq,sigmaz);   //Chi2 to compare track and MC for Chimin parameters

               FitParameters.sChi2Min=CHi2fit;
               //FitParameters.sNumMCPoint=num_MC_Point;
               //FitParameters.sNormChi2=chi2min/num_MC_Point;

               BackwardExtrapolation();

               std::cout<<cYELLOW<<" Minimization result : "<<std::endl;
               std::cout<<" Scattering Angle : "<<aMC*180.0/TMath::Pi()<<std::endl;
               std::cout<<" Azimutal angle : "<<phiMC*180.0/TMath::Pi()<<std::endl;
               std::cout<<" B : "<<Bmin<<std::endl;
               std::cout<<" Brho : "<<FitParameters.sBrhoMin<<std::endl;
               std::cout<<" Energy : "<<e0sm<<std::endl;
               std::cout<<" Vertex Position - X : "<<x0MC<<" - Y : "<<y0MC<<" - Z : "<<z0MC<<std::endl;
               std::cout<<" Vertex Position (Backward extrapolation) - X : "<<fVertexPos.X()<<" - Y : "<<fVertexPos.Y()<<" - Z : "<<fVertexPos.Z()<<std::endl;
               std::cout<<" Vertex Energy : "<<fVertexEner<<" MeV "<<std::endl;
               //std::cout<<" Reduced chi2 : "<<chi2min/FitParameters.sNumMCPoint<<std::endl;
               std::cout<<" Minimum chi2 : "<<CHi2fit<<cNORMAL<<std::endl;

               fPadPlane = NULL;
               delete fPadPlane;

               return kTRUE;


}

#endif
