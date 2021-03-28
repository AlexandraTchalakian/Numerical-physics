#include <iostream>

#include "PhysEngine.hpp"
#include "ConfigFile.hpp"

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
    PhysEngine* engine;

    engine = new PhysengineRungeKutta4(configFile);
  
    
  configFile.printOut(configFile.get<std::string>("outputPath").append(".in"));
    if(engine!=NULL){
        engine->run();
        delete engine;
    }

}
