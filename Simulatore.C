//#if !defined(__CINT__) || defined(__MAKECINT__))
#include <Riostream.h>
#include "TObject.h"
#include "TSystem.h"
#include "TMath.h"
#include "TAxis.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TCanvas.h"
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
   
  

   const double hydrogen_density=8e22;
   const double carbon_density=4e22;
   double Estart;
   double Nstart;
   double lateral_size;
   double thick;
   double beam_size;
   double y_start=0;// y_start e x_start sono le posizioni del centro del fascio, che per ora lascio a zero
   double x_start=0;
   double radius;
   double shield_sphere_dist;
   double y_centro=0;//sono le coordinate del centro della sfera, per ora sono a zero
   double x_centro=0;

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

   
   //TFile file("neutroni","RECREATE");
   //TTree *tree=new TTree("tree","tree di neutroni uscenti");

  
  Propagatore *prop=new Propagatore(h_Cscatt,h_Hscatt,h_Cabs,h_Habs,lateral_size,thick,hydrogen_density,carbon_density);
  Generatore *gen=new Generatore(Nstart,Estart,beam_size,beam_size,x_start,y_start);
  Rivelatore *riv=new Rivelatore(radius,0,0,shield_sphere_dist+prop->GetTargetThick());
  Neutron *n;
  Neutron *n_out;

  //Neutron &nn_in=*n;
  //Neutron &nn_out=*n_out;

  //tree->Branch("nn_in",&nn_in);
  //tree->Branch("nn_out",&nn_out);
  
  TH1D *spectrum=new TH1D("spectrum","spectrum",1e6,0,1e6);
  
  for(int i=0;i<gen->GetParticles();i++){
   
    n=gen->Genera_neutrone();
    n_out=prop->Propagation(n);
    
     double length=riv->Intersezione(n_out);
     if((n_out->GetEnergy()>0) && !(n_out->GetAbsorption()) && length!=0 ) {spectrum->Fill(n_out->GetEnergy(),length/riv->GetVolume());}
    
    
    //tree->Fill();
  
 
  }
 
  spectrum->GetXaxis()->SetTitle("Energy (eV)");
  spectrum->GetYaxis()->SetTitle("Fluence Spectrum #Delta#phi/#Delta E");
  spectrum->Draw();
  
  //file.Write();
  //file.Close();

  cout<<" "<<endl;
  cout<<"Fluence per Starting Particle: "<<riv->GetFluence()/Nstart<<endl;
  cout<<" "<<endl;
  //cout<<"Fluence Spectrum per Starting Particle: "<<spectrum->Integral()/Nstart<<endl;

}




//DA AGGIUNGERE: fluence spectrum, distribuzione in energia neutroni
