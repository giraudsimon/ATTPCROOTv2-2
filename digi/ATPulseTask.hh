/*****************************************************************/
/*    ATPulseTask: Simulates the ionized electrons which are */
/*    drifted and amplified by avalanche when reaching the mesh. */
/*    Log: Created in 24-10-2016     				 */
/*    Author: Nathan Watwood (NSCL) 				 */
/*    ayyadlim@nscl.msu.edu                                      */
/*****************************************************************/
#ifndef ATPulseTask_H
#define ATPulseTask_H

#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"

#include "FairTask.h"
#include "FairMCPoint.h"

#include "TClonesArray.h"
#include "ATDigiPar.hh"
#include "TH2Poly.h"
#include "AtTpcMap.h"
#include "ATRawEvent.hh"
#include "ATGas.hh"



class ATPulseTask : public FairTask
{
  public:
     ATPulseTask();
     ~ATPulseTask();

    void SetPersistence(Bool_t val) { fIsPersistent = val; }
    virtual InitStatus Init();        //!< Initiliazation of task at the beginning of a run.
    virtual void Exec(Option_t* opt); //!< Executed for each event.
    virtual void SetParContainers();  //!< Load the parameter container from the runtime database.

   private:
    ATGas*     fGas;                     //!< Gas parameter container.
    ATDigiPar* fPar; //!< Base parameter container.
    Int_t fEventID;                      //!< EventID
    Double_t fGain;                      //!< Gain.
    Bool_t fIsPersistent;                //!< If true, save container
    TClonesArray* fDriftedElectronArray; //!< drifted electron array (input)
    TClonesArray* fRawEventArray;        //!< Raw Event array(only one)
    ATRawEvent* fRawEvent;               //!< Raw Event Object
    TH2Poly *fPadPlane;                  //!< pad plane
    AtTpcMap *fMap;                      //!<ATTPC map
    Int_t fInternalID;                   //!<Internal ID


     ClassDef(ATPulseTask,1);

};

#endif
