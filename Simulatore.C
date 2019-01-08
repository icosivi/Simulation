//#if !defined(__CINT__) || defined(__MAKECINT__))
#include <Riostream.h>
#include "TObject.h"
#include "TSystem.h"
#include "TMath.h"
#include "TAxis.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "TRandom3.h"
#include "Punto.h"
#include "Retta.h"
#include "Neutron.h"
#include "Propagatore.h"
#include "Generatore.h"
#include "Rivelatore.h"


void Simulatore(){

  gSystem->CompileMacro("Punto.cxx","kg");
  gSystem->CompileMacro("Retta.cxx","kg");
  gSystem->CompileMacro("Neutron.cxx","kg");
  gSystem->CompileMacro("Propagatore.cxx","kg");
  gSystem->CompileMacro("Generatore.cxx","kg");
  gSystem->CompileMacro("Rivelatore.cxx","kg");

  
  TFile *f0=new TFile("Cabs.root","READ");
  TFile *f1=new TFile("Habs.root","READ");
  TFile *f2=new TFile("Cscat.root","READ");
  TFile *f3=new TFile("Hscat.root","READ");

   TH1F *h_Cabs=dynamic_cast<TH1F*>(f0->Get("C_absorption"));
   TH1F *h_Habs=dynamic_cast<TH1F*>(f1->Get("H_absorption"));
   TH1F *h_Cscatt=dynamic_cast<TH1F*>(f2->Get("C_scattering"));
   TH1F *h_Hscatt=dynamic_cast<TH1F*>(f3->Get("H_scattering"));
   
  

   const double polietilene_density=1.2e23;
   double Estart;
   int Nstart;
   double lateral_size;
   double thick;
   double beam_size;
   double y_start=0;// y_start e z_start sono le posizioni del centro del fascio, che per ora lascio a zero
   double z_start=0;
   double radius;
   double shield_sphere_dist;
   double y_centro=0;//sono le coordinate del centro della sfera, per ora sono a zero
   double z_centro=0;

   double seed=0.1;
   gRandom->SetSeed(seed);
   

   cout << "Initial Beam Energy (eV)" << endl;
   cin >> Estart;
   cout << "Number of neutrons " << endl;
   cin >> Nstart;
   cout << "Shield lateral size (mm) " << endl;
   cin >> lateral_size;
   cout << "Shield thickness (mm) " << endl;
   cin >> thick;
   cout << "Beam lateral size (mm) " << endl;
   cin >> beam_size;
   cout << "Radius of the spherical detector (mm) " << endl;
   cin >> radius;
   cout << "Distance shield-sphere centre (mm) " << endl;
   cin >> shield_sphere_dist;

   
   TFile file("neutroni.root","RECREATE");
   TTree *tree=new TTree("tree","tree di neutroni uscenti");
   TClonesArray *ptr=new TClonesArray("Neutron",100);
   TClonesArray &nptr=*ptr;

   tree->Branch("neutron out",&nptr);

   
  
  Propagatore *prop=new Propagatore(h_Cscatt,h_Hscatt,h_Cabs,h_Habs,lateral_size,thick,polietilene_density);
  Generatore *gen=new Generatore(Nstart,Estart,beam_size,beam_size,y_start,z_start);
  Rivelatore *riv=new Rivelatore(radius,shield_sphere_dist+prop->GetTargetThick(),0,0);
  //Punto *p=new Punto(0,0,0);
  //Retta *r=new Retta(p,0,0);
  Neutron *n;
  Neutron *n_out;
  
  for(int i=0;i<gen->GetParticles();i++){
   
    n=gen->Genera_neutrone();
    n_out=prop->Propagation(n);

    //new(nptr[i]) Neutron(n_out->GetPunto(),n_out->GetRetta(),n_out->GetEnergy());
    
    riv->Intersezione(n_out);
    
 		  
  }

  //tree->Fill();
  file.Write();
  file.Close();

  cout<<" "<<endl;
  cout<<"Fluence per Starting Particle: "<<riv->GetFluence()/Nstart<<endl;
  cout<<" "<<endl;
  cout<<" "<<endl;


}




//DA AGGIUNGERE: fluence spectrum, distribuzione in energia neutroni
