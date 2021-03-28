#ifndef TEST_PHYSENGINE_H_INCLUDED
#define TEST_PHYSENGINE_H_INCLUDED

#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <valarray>


#include "ConfigFile.hpp"
//#include "ConfigFile.tpp"




typedef std::valarray<double> Vecteur;


class PhysEngine
{
public:
  
  PhysEngine(const ConfigFile configFile);
  
  
    void set_dt (const double& d);
    double get_distance (const unsigned int& astre /* 1 = terre, 2 = asteroide, 3 = lune*/) const;
    double get_dt() const;
    

    
    virtual int run() =0;
    
    virtual ~PhysEngine();
protected:
    double dt,tfin,G,Mt,Ma,Ml,Rt,Rl,alpha,epsilon;
    Vecteur ra;
    Vecteur rt;
    Vecteur rl;
    Vecteur va;
    Vecteur vt;
    Vecteur vl;
    double t;
    bool Jeanbob;
    double R,V,Emec;
  
  void printout(double h);  
  void acceleration(Vecteur& A, Vecteur r0a,Vecteur r0t,Vecteur r0l,Vecteur v0a, Vecteur v0t, Vecteur v0l);
 
  /** step() is a virtual method which has to be implemented in each 
   *  children classes
   */
    
  virtual void step() =0 ;
  virtual void adaptateur() =0;

    std::ofstream *outputFile;

};

class PhysengineRungeKutta4 : public PhysEngine {

    public :
    
    PhysengineRungeKutta4(const ConfigFile configfile) : PhysEngine(configfile){}
    virtual int run() override;


    
        virtual void adaptateur() override;

    virtual void step() override;
    
    virtual ~PhysengineRungeKutta4() {std::cout << "destructeur Physengine4" << std::endl;}

};
#endif //TEST_PHYSENGINEs_H_INCLUDED
