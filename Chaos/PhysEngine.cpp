#include <iostream>
#include <cmath> 
#include <fstream>
#include <iomanip>

#include "ConfigFile.hpp"
#include "ConfigFile.tpp"
#include "PhysEngine.hpp"

using namespace std;

/** Constructor
 */
double pi(acos(-1));
PhysEngine::PhysEngine(const ConfigFile configFile){
  outputFile=NULL;
  try{    
    Ig=configFile.get<double>("Ig");
    mu=configFile.get<double>("mu");
    B0=configFile.get<double>("B0");
    B1=configFile.get<double>("B1");
    nu=configFile.get<double>("nu");   
    omega= configFile.get<double>("omega");
    tFin=configFile.get<double>("tFin");   
    dt=configFile.get<double>("dt");
    t=0;  
    Q=configFile.get<double>("theta0");  
    P=configFile.get<double>("vtheta0");  
   
    
    double w0=sqrt(mu*B0/Ig);
    double w1=sqrt(mu*B1/Ig);
    
    if(dt<0){
      double nDtParT= configFile.get<double>("nDtParT");   
      //Frequency always omega (Poincare map)
      dt=2*pi/nDtParT/omega;
      tFin=2*pi*tFin/omega;
    }    
    
    energyMec=0.5*Ig*pow(P,2)-mu*B0*cos(Q);//manque l'énerge potentielle
    powerNC=mu*B1*sin(omega*t)*sin(Q)*P+nu*pow(P,2); 
    
    std::string path=configFile.get<std::string>("outputPath");    
    outputFile = new std::ofstream();
    outputFile->open( path.c_str());
    if (!outputFile->is_open())
    {
      delete outputFile;
      outputFile=NULL;
    }
    
    cout << "wRef=" << omega << ", w0=" << w0 << ", w1=" << w1 ;
    cout << ", omega=" << omega << ", dt=" << dt << ", tFin=" << tFin << endl;
    cout << "nu/I=" << nu/Ig << ", w0^2=" << w0*w0 ;
    cout << ", A=w1^2=" << w1*w1 << ", omega/w1=" << omega/w1 << endl;
    
  }catch(std::string e){
    cerr << "Exception: " << e.c_str() << endl;
    if(outputFile!=NULL){outputFile->close();delete outputFile;}
  }
}

/** Destructor
 */
PhysEngine::~PhysEngine(){
  outputFile->close();
  delete outputFile;
}
  
int PhysEngine::run(){
  
  //Si il y a eu un probleme a la construction, on sort.
  if(outputFile==NULL){
    return -1;
  }
  cout << "Start the physics engine" << endl;
  
  printout(); 
  while(t<tFin-dt*0.5){
	  
      step();
      t=t+dt;//Current
        
   
    energyMec=0.5*Ig*pow(P,2)-mu*B0*cos(Q);
    powerNC = -(mu*B1)*sin(omega*t)*sin(Q)*P-nu*pow(P,2);
    Energylost+=powerNC*dt;
    printout();
  }
  cout << "Stop the physics engine" << endl;
  
return 0;
} 

/** Print the data in the output file
 * */
void PhysEngine::printout(){
  *outputFile << setprecision(14) << t << " " << Q<< " " << P ;
  *outputFile << " " << energyMec << " " << powerNC << " "<< (powerNC*dt) << " " <<Energylost <<endl;
}

/** Return the angle acceleration
 */
double PhysEngine::acceleration(const double q, const double p, const double t){
  return (-mu*(B0+B1*sin(omega*t))*sin(q)-nu*p)/Ig;
}

///////////////////////////////////////////////
//
// Implementation of the virtual step() method
//
//////////////////////////////////////////////

void PhysEngineEuler::step(){
    double Pprev=P;
    P+=dt*acceleration(Q,P,t);
    Q+=dt*Pprev;
} 

void PhysEngineEulerCromerA::step(){        
double Pbis;
    Pbis=P+acceleration(Q,P,t)*0.5*dt;
    Q=Q+Pbis*0.5*dt;
} 

void PhysEngineEulerStormerVerlet::step(){
    double Qprev(Q);
    double Pbis;
    Q+=P*dt+0.5*acceleration(Q,P,t)*pow(dt,2);
    Pbis=P+0.5*acceleration(Qprev,P,t)*dt,
    P+=(acceleration(Qprev,Pbis,t)+acceleration(Q,Pbis,t+dt))*0.5*dt;
}
	 
