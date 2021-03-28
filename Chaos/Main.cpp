#include <iostream>

#include "ConfigFile.hpp"
#include "PhysEngine.hpp"

int main(int argc, char* argv[]) 
{  
  std::string inputPath = "configuration.in";
  if(argc>1){
    inputPath=argv[1];
  } 
  
  ConfigFile configFile(inputPath);
  
  for(int i=2;i<argc;i++){
    configFile.process(argv[i]);
  }
  
  std::string solver = configFile.get<std::string>("solver");
  
  PhysEngine* engine;
  if (solver.compare("Euler")==0){
      engine = new PhysEngineEuler(configFile);
  }else if(solver.compare("EulerCromerA")==0){
      engine = new PhysEngineEulerCromerA(configFile);
  }else if(solver.compare("StormerVerlet")==0){
      engine = new PhysEngineEulerStormerVerlet(configFile);
  }else{
      engine=NULL;
      cerr << "Unkown solver "<< solver <<" !\n";
  }
  
  configFile.printOut(configFile.get<std::string>("outputPath").append(".in"));
    
  if(engine!=NULL){
    engine->run();
    delete engine;    
  }
}