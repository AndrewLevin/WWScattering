//==============================================================================
//
//                        Logon file
//
//==============================================================================
{
//#include <iomanip.h>
//#include <time.h>
 
  // Look for include files
  //gSystem->SetIncludePath("-I$BS_ANALYSIS/include -I$CDFSOFT2_DIR/include");
 
  // Load in ROOT physics vectors and event generator libraries
  gSystem->Load("$ROOTSYS/lib/libPhysics.so");
  gSystem->Load("$ROOTSYS/lib/libEG.so");

  gInterpreter->ExecuteMacro("/home/ceballos/root/macros/MITStyle.C");
  gROOT->Macro("/home/ceballos/releases/CMSSW_4_2_8_patch4/src/auxiliar/LeptonScaleLookup.cc+");

  // This line reports the process ID which simplifies debugging
  gInterpreter->ProcessLine(".! ps | grep root.exe");
   
 END:;
}
