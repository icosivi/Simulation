//#if !defined(__CINT__) || defined(__MAKECINT__))
#include <Riostream.h>
#include "TObject.h"
#include "TSystem.h"
#include "TStopwatch.h"
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

  gSystem->CompileMacro("Punto.cxx","kgf");
  gSystem->CompileMacro("Retta.cxx","kgf");
  gSystem->CompileMacro("Neutron.cxx","kgf");
  gSystem->CompileMacro("Propagatore.cxx","kgf");
  gSystem->CompileMacro("Generatore.cxx","kgf");
  gSystem->CompileMacro("Rivelatore.cxx","kgf");

  
  TFile *f0=new TFile("CrossSections/Cabs.root","READ");
  TFile *f1=new TFile("CrossSections/Habs.root","READ");
  TFile *f2=new TFile("CrossSections/Cscat.root","READ");
  TFile *f3=new TFile("CrossSections/Hscat.root","READ");

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
  
  Propagatore *prop=new Propagatore(h_Cscatt,h_Hscatt,h_Cabs,h_Habs,lateral_size,thick,hydrogen_density,carbon_density);
  Generatore *gen=new Generatore(Nstart,Estart,beam_size,beam_size,x_start,y_start);
  Rivelatore *riv=new Rivelatore(radius,0,0,shield_sphere_dist+prop->GetTargetThick());

  double x_vec[461];


  for(int i=0;i<=100;i++){
    x_vec[i]=i;
  }

  for(int j=1;j<=90;j++){
	x_vec[j+100]=100+j*10;
      }

  for(int j=1;j<=90;j++){
	x_vec[j+190]=1000+j*100;
      }

   for(int j=1;j<=90;j++){
	x_vec[j+280]=10000+j*1000;
      }

    for(int j=1;j<=90;j++){
	x_vec[j+370]=100000+j*10000;
	}

  
    Neutron *n_in=new Neutron();
    Neutron *n_out=new Neutron();
    
    //Neutron &nn_in=*n;
    // Neutron &nn_out=*n_out;

    TFile *outFile=new TFile("neutronTree.root","RECREATE");
    TTree *tree=new TTree("tree","tree di neutroni uscenti");
    tree->Branch("n_in",&n_in,32000,2);
    tree->Branch("n_out",&n_out,32000,2);

    TCanvas *c1=new TCanvas("c1","c1",1200,800);
    TCanvas *c2=new TCanvas("c2","c2",1200,800);
  
    TH1D *spectrum=new TH1D("spectrum","spectrum",460,x_vec);  //x_vec deve avere dimensione nbins+1 !!!!
    TH1D *en_distrib=new TH1D("en distrib","en distrib",460,x_vec); 

    TStopwatch *watch=new TStopwatch();
    double time=0;
   
  
  for(int i=0;i<gen->GetParticles();i++){
    
    watch->Start();
    
    gen->Genera_neutrone(n_in);
    *n_out=Neutron(n_in);
    prop->Propagation(n_out);

    watch->Stop();
    time+=watch->RealTime();
    
    
     double length=riv->Intersezione(n_out);
     
     if((n_out->GetEnergy()>0) && !(n_out->GetAbsorption()) && length!=0 ) {

       int bin=spectrum->FindBin(n_out->GetEnergy());
       double deltaE=spectrum->GetBinWidth(bin);
       spectrum->Fill(n_out->GetEnergy(),((length/riv->GetVolume())/deltaE)/Nstart); //fluence spectrum per starting particle
       en_distrib->Fill(n_out->GetEnergy(),1./deltaE);

     }
    
    
    tree->Fill();
    n_in->Reset();
  
 
  }

  tree->Write();
  delete tree;
  delete n_in;
  delete n_out;
  outFile->Close();

  // c1->cd();
  // spectrum->GetXaxis()->SetTitle("Energy (eV)");
  // spectrum->GetYaxis()->SetTitle("Fluence Spectrum #Delta#phi/#Delta E");
  // spectrum->Draw("HIST");

  // c2->cd();
  // en_distrib->SetTitle("Distribuzione Energia neutroni uscenti");
  // en_distrib->GetXaxis()->SetTitle("eV");
  // en_distrib->GetYaxis()->SetTitle("Counts/#Delta E");
  // en_distrib->Draw("HIST");

  // c1->Close();
  // c2->Close();
  
  //file.Write();
  //file.Close();

  cout<<" "<<endl;
  cout<<"Fluence per Starting Particle: "<<riv->GetFluence()/Nstart<<endl;
  cout<<" "<<endl;
  cout<<"Mean number of collisions: "<<prop->GetNcoll()/Nstart<<endl;
  cout<<" "<<endl;
  cout<<"Mean CPU time per neutron: "<<time/Nstart<<endl;
  cout<<" "<<endl;
 

}



//DA AGGIUNGERE: dose, aggiungere le uncertainties, confronto con MCNP con S(alfa,beta) turned on and off (?)

//DA SISTEMARE: sezioni d'urto con bin vuoti

//Opzione aggiuntiva: riempire un tree coi neutroni uscenti e utilizzarlo per l'analisi

//Deadline: 10 Febbraio (report con nostri risultati)
