#ifndef TEST_PHYSENGINE_H_INCLUDED
#define TEST_PHYSENGINE_H_INCLUDED

using namespace std;
class PhysEngine
{
public:
  
  PhysEngine(const ConfigFile configFile);  
  ~PhysEngine();
  
  int run();
  
protected:
  double dt,tFin,omega,B0,B1,mu,nu,Ig;
  double t,Q,P,powerNC,energyMec,Energylost;
  
  void printout();  
  double acceleration(const double q, const double p, const double t);
 
  /** step() is a virtual method which has to be implemented in each 
   *  children classes
   */
  virtual void step() =0;
  
private:
  std::ofstream *outputFile;
  
};

class PhysEngineEuler : public PhysEngine{
  public:
  
  PhysEngineEuler(const ConfigFile configFile):PhysEngine(configFile){};  
  
  protected:
    void step();
};

class PhysEngineEulerCromerA : public PhysEngine{
  public:
  
  PhysEngineEulerCromerA(const ConfigFile configFile):PhysEngine(configFile){};  
  
  protected:
    void step();
};

class PhysEngineEulerStormerVerlet : public PhysEngine{
  public:
  
  PhysEngineEulerStormerVerlet(const ConfigFile configFile):PhysEngine(configFile){};  
  
  protected:
    void step();
    double Qprev,Pprev;
};

#endif //TEST_PHYSENGINEs_H_INCLUDED
