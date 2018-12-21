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
#include "Rivelatore.h"



ClassImp(Rivelatore)



Rivelatore::Rivelatore(double raggio, double x_centro, double y_centro, double z_centro):TObject(),
  rp(raggio),
  xp(x_centro),
  yp(y_centro),
  zp(z_centro),
  tp(0){
  
}


Rivelatore::~Rivelatore(){
  // distruttore
 
}



//-------------------------------------------




void Rivelatore::Intersezione(Neutron *n){

  
      if( !(n->GetAbsorption()) && n->GetX()>=0 ){

    
          double vx=TMath::Sin(n->GetTheta())*TMath::Cos(n->GetPhi());
	  double vy=TMath::Sin(n->GetTheta())*TMath::Sin(n->GetPhi());
	  double vz=TMath::Cos(n->GetTheta());

	  double beq=2*(vx*(n->GetX()-this->xp)+vy*(n->GetY()-this->yp)+vz*(n->GetZ()-this->zp));
	  double ceq=TMath::Power(n->GetX(),2)+TMath::Power(n->GetY(),2)+TMath::Power(n->GetZ(),2)+TMath::Power(this->xp,2)-2*n->GetX()*this->xp+TMath::Power(this->yp,2)-2*n->GetY()*this->yp+TMath::Power(this->zp,2)-2*n->GetZ()*this->zp-(this->rp*this->rp);

	  double delta=TMath::Power(beq,2)-4*ceq;
	  double t_1=(0.5)*(-beq+TMath::Sqrt(delta));
	  double t_2=(0.5)*(-beq-TMath::Sqrt(delta));
	  cout<<delta<<endl;
	  
	   if(delta>0){ this->tp=this->tp+TMath::Abs(t_1-t_2); }
	   
		      

  }
}
