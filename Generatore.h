#ifndef GENERATORE_H
#define GENERATORE_H

#include "TObject.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TH1.h"
#include "Retta.h"
#include "Punto.h"
#include "Neutron.h"





class Generatore : public TObject
{

 public:

  //ho definito un generatore con x0 già settato a zero
  Generatore(int N_particles,double E_iniziale,double dy_beam, double dz_beam, double y0, double z0);
  

  virtual ~Generatore();

  
  Neutron* Genera_neutrone();
  int GetParticles() const {return Np;}
  double GetEn() const {return Ep;}
  double Getdy() const {return dy_p;}
  double Getdz() const {return dz_p;}
  double Gety0() const {return y0p;}
  double Getz0() const {return z0p;}
  
  
 private:

  int Np;
  double Ep;
  double dy_p;
  double dz_p;
  double y0p;
  double z0p;
 

  
  ClassDef(Generatore,1)
  
};

#endif
