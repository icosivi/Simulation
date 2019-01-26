#include <Riostream.h>
#include <TObject.h>
#include <TSystem.h>
#include "TStopwatch.h"
#include "TMath.h"
#include "TAxis.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <TTree.h>
#include <TBranch.h>
#include "TClonesArray.h"
#include "TRandom3.h"
#include "Punto.h"
#include "Retta.h"
#include "Neutron.h"
#include "Propagatore.h"
#include "Generatore.h"
#include "Rivelatore.h"


void Compilatore(){
  gSystem->CompileMacro("Punto.cxx","kgf");
  gSystem->CompileMacro("Retta.cxx","kgf");
  gSystem->CompileMacro("Neutron.cxx","kgf");
  gSystem->CompileMacro("Generatore.cxx","kgf");
  gSystem->CompileMacro("Propagatore.cxx","kgf");
  gSystem->CompileMacro("Rivelatore.cxx","kgf");
  gSystem->CompileMacro("Simulatore.C","kgf");

}
