#include <Riostream.h>
#include "TObject.h"
#include "TMath.h"
#include "TAxis.h"
#include "TFile.h"
#include "TH1.h"
#include "TRandom3.h"
#include "Punto.h"
#include "Retta.h"
#include "Neutron.h"
#include "Generatore.h"



ClassImp(Generatore)



Generatore::Generatore(int N_particles,double E_iniziale,double dy_beam, double dz_beam, double y0, double z0):TObject(),
  Np(N_particles),
  Ep(E_iniziale),
  dy_p(dy_beam),
  dz_p(dz_beam),
  y0p(y0),
  z0p(z0){
}


Generatore::~Generatore(){
  // distruttore
 
}



//-------------------------------------------


Neutron* Generatore::Genera_neutrone(){


  double y0=this->y0p+(this->dy_p/2-this->dy_p*gRandom->Rndm());
  double z0=this->z0p+(this->dz_p/2-this->dz_p*gRandom->Rndm()); 
  double x0=0; //per ora partiamo contro lo slab in x=0

  Punto *p=new Punto(x0,y0,z0);
  Retta *r=new Retta(p,TMath::Pi()*0.5,0);
  Neutron *n=new Neutron(p,r,this->Ep);

  //delete p;////
  //delete r;////
  
  return n;


}
