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
#include "Propagatore.h"



ClassImp(Propagatore)



Propagatore::Propagatore(TH1 *Cscatt, TH1 *Hscatt, TH1 *Cabs, TH1 *Habs, double target_side, double target_thickness, double density):TObject(),
  Cscattp(Cscatt),
  Hscattp(Hscatt),
  Cabsp(Cabs),
  Habsp(Habs),
  ptarget_side(target_side),
  ptarget_thickness(target_thickness),
  pdensity(density){
}


Propagatore::~Propagatore(){
  // distruttore
  delete Cscattp;
  delete Hscattp;;
  delete Cabsp;
  delete Habsp;
}



//_____________________________________________________________________________
//_____________________________________________________________________________




double Propagatore::x_interaction(double En){
  
  const TAxis *xaxis0=this->Cscattp->GetXaxis();
  int binx0=xaxis0->FindBin(En);
  double sigma_csc=this->Cscattp->GetBinContent(binx0)*this->pdensity*(1e-24);  

  const TAxis *xaxis1=this->Hscattp->GetXaxis();
  int binx1=xaxis1->FindBin(En);
  double sigma_hsc=this->Hscattp->GetBinContent(binx1)*this->pdensity*(1e-24);

  const TAxis *xaxis2=this->Cabsp->GetXaxis();
  int binx2=xaxis2->FindBin(En);
  double sigma_cab=this->Cabsp->GetBinContent(binx2)*this->pdensity*(1e-24);

  const TAxis *xaxis3=this->Habsp->GetXaxis();
  int binx3=xaxis3->FindBin(En);
  double sigma_hab=this->Habsp->GetBinContent(binx3)*this->pdensity*(1e-24);
 

  double sigma_tot=sigma_csc+sigma_hsc+sigma_cab+sigma_hab;

  //double xint=-(1./sigma_tot)*log(1-gRandom->Rndm())*10;
  // cout<<xint<<endl;

  return -(1./sigma_tot)*log(1-gRandom->Rndm())*10;  //moltiplico per 10 per passare da cm a mm (la sigma macroscopica è in 1/cm )

 
}



//_____________________________________________________________________________



Neutron* Propagatore::scattering(Neutron *n){ //Ricevo un warning perchè potrei non ritornare nulla, ma non dovrebbe essere un problema 
  
  const TAxis *xaxis0=this->Cscattp->GetXaxis();
  int binx0=xaxis0->FindBin(n->GetEnergy());
  double sigma_csc=this->Cscattp->GetBinContent(binx0)*this->pdensity*(1e-24);  

  const TAxis *xaxis1=this->Hscattp->GetXaxis();
  int binx1=xaxis1->FindBin(n->GetEnergy());
  double sigma_hsc=this->Hscattp->GetBinContent(binx1)*this->pdensity*(1e-24);

  const TAxis *xaxis2=this->Cabsp->GetXaxis();
  int binx2=xaxis2->FindBin(n->GetEnergy());
  double sigma_cab=this->Cabsp->GetBinContent(binx2)*this->pdensity*(1e-24);

  const TAxis *xaxis3=this->Habsp->GetXaxis();
  int binx3=xaxis3->FindBin(n->GetEnergy());
  double sigma_hab=this->Habsp->GetBinContent(binx3)*this->pdensity*(1e-24);

  double sigma_tot=sigma_csc+sigma_hsc+sigma_cab+sigma_hab;
  double sigma_abs=sigma_cab+sigma_hab;

  double absorption_limit=sigma_abs/sigma_tot;
  double carbon_scatt=sigma_csc/sigma_tot;
  double hidrogen_scatt=sigma_hsc/sigma_tot;


  double alfa_cm=1-2*gRandom->Rndm();  //angolo theta nel CM, uniforme tra -1 ed 1
  double beta_cm=1-2*gRandom->Rndm();  //angolo phi nel CM, uniforme tra -1 ed 1

  double interaction_type=gRandom->Rndm();

  double new_energy;
  double delta_th;
  double new_theta;
  double delta_ph;
  double new_phi;

  if(interaction_type<=absorption_limit){

    n->SetAbsorption();
    return n;

  }

  if(interaction_type>absorption_limit && interaction_type<=(absorption_limit+hidrogen_scatt)){  //il nucleo scelto è l'idrogeno 

    new_energy=n->GetEnergy()*(0.5)*(1.+alfa_cm);
    
    delta_th=TMath::ACos((alfa_cm+1.)/(TMath::Sqrt(2.+(2.)*alfa_cm)));//definisco l'angolo delta teta nel SR-LAB
    new_theta=delta_th+n->GetTheta();

    delta_ph=TMath::ACos((beta_cm+1.)/(TMath::Sqrt(2.+(2.)*beta_cm)));//definisco l'angolo delta phi nel SR-LAB
    new_phi=delta_ph+n->GetPhi();
 
    n->SetEnergia(new_energy);
   
    n->GetRetta()->SetTheta(new_theta);
    n->GetRetta()->SetPhi(new_phi);
    
    return n;
    
  }

  if(interaction_type>(absorption_limit+hidrogen_scatt)){  //il nucleo è il carbonio
  
    new_energy=n->GetEnergy()*(1./169.)*(145.+24*alfa_cm);
 
    delta_th=TMath::ACos(((12.)*alfa_cm+1.)/(TMath::Sqrt(145.+(24.)*alfa_cm)));
    new_theta=delta_th+n->GetTheta();
   

    delta_ph=TMath::ACos(((12.)*beta_cm+1.)/(TMath::Sqrt(145.+(24.)*beta_cm)));
    new_phi=delta_ph+n->GetPhi();
    

    n->SetEnergia(new_energy);
    n->GetRetta()->SetTheta(new_theta);
    n->GetRetta()->SetPhi(new_phi);
    
    return n;
     
  }


}




Neutron* Propagatore::Propagation(Neutron *n){

  bool contr=kTRUE;
  
  while(!(n->GetAbsorption()) && contr){

    /*cout<<"x "<<n->GetX()<<endl;
    cout<<"y "<<n->GetY()<<endl;
    cout<<"z "<<n->GetZ()<<endl;
    cout<<"th "<<n->GetTheta()<<endl;
    cout<<"ph "<<n->GetPhi()<<endl;
    cout<<" "<<endl;*/

      n->SetNuovoPunto(this->x_interaction(n->GetEnergy()));

      
      if( n->GetX()>=0 && n->GetX()<=this->ptarget_thickness && n->GetY()>=-(0.5)*this->ptarget_side &&  n->GetY()<=(0.5)*this->ptarget_side && n->GetZ()>=-(0.5)*this->ptarget_side && n->GetZ()<=(0.5)*this->ptarget_side ){
	
	this->scattering(n);

      }else{contr=kFALSE; }


 }

 return n;


}


//_____________________________________________________________________________


//matrice per passare da CM a LAB
/*void Rivelatore::Rotate(double theta, double phi, double theta0, double phi0, double *cd){

 
  double mr[3][3];

  mr[0][0]=-TMath::Sin(phi);
  mr[1][0]=TMath::Cos(phi);
  mr[2][0]=0.;
  mr[0][1]=-TMath::Cos(phi)*TMath::Cos(theta);
  mr[1][1]=-TMath::Cos(theta)*TMath::Sin(phi);
  mr[2][1]=TMath::Sin(theta);
  mr[0][2]=TMath::Sin(theta)*TMath::Cos(phi);
  mr[1][2]=TMath::Sin(theta)*TMath::Sin(phi);
  mr[2][2]=TMath::Cos(theta);

  double cdp[3];
  cdp[0]=TMath::Sin(theta0)*TMath::Cos(phi0);
  cdp[1]=TMath::Sin(theta0)*TMath::Sin(phi0);
  cdp[2]=TMath::Cos(theta0);

  for(int i=0;i<3;i++){
    cd[i]=0.;
    for(int j=0;j<3;j++){
      cd[i]+=mr[i][j]*cdp[j];
    }
    }*/



