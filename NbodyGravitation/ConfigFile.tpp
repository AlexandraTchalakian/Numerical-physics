#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <fstream>

#include "ConfigFile.hpp"

ConfigFile::ConfigFile(const std::string& fileName){
  std::ifstream file;
  file.open(fileName.c_str());
  if (!file){
    throw "[ConfigFile] Error in openening file "+fileName;
  }
  try{
    std::string lineread;
    while(std::getline(file, lineread)){
      process(lineread);
    }
    file.close();   
  }catch(std::string e){
    file.close();   
    throw e;
  }
}

ConfigFile::~ConfigFile(){}

void ConfigFile::printOut(std::string path){  
/*  std::ofstream outputFile();
  outputFile.open(path.c_str());*/ 
  std::ofstream outputFile(path.c_str());
  if (outputFile.is_open())
  {
    outputFile << toString() << endl;
  }
  outputFile.close();
}

std::string ConfigFile::toString() {
        std::map<std::string, std::string>::iterator iter;
        std::string strToReturn; //This is no longer on the heap

        for (iter = configMap.begin(); iter != configMap.end(); ++iter) {
           strToReturn.append(iter->first); //Not a method call
           strToReturn.append("=");
           strToReturn.append(iter->second);
           strToReturn.append("\n");
           //....
           // Make sure you don't modify table here or the iterators will not work as you expect
        }
        //...
        return strToReturn;
}

void ConfigFile::process(const std::string& lineread) {
  //cout << "Process: "<< lineread <<endl;
  size_t commentPosition=lineread.find('#',0);
      if(commentPosition!=0){
	size_t spacePosition=lineread.find(' ',0);
	if(spacePosition!=string::npos){
	  cerr << "ConfigFile object do not support character ' ' at position " << spacePosition << " in line '" << lineread << endl; 
	}
	size_t equalPosition = lineread.find('=',1);
	if(equalPosition==string::npos){
	  cerr << "Line without '=' " << lineread << endl; 
	}else{
	  string key=lineread.substr(0,equalPosition);
	  string value=lineread.substr(equalPosition+1,lineread.length());  
	  std::map<std::string, std::string>::const_iterator val = configMap.find(key);
	  if ( val != configMap.end() ) {
	    configMap.erase(key);
	  }
	  configMap.insert( std::pair<std::string, std::string>(key,value) );      
	}
      }
  
}


template<typename T> T ConfigFile::get(const std::string& key) const{
   std::map<std::string, std::string>::const_iterator val = configMap.find(key);
   if ( val != configMap.end() ) {
     T out;
     std::istringstream iss(val->second);
     iss >> out;
     cout << "\t" << key << "=" << out << endl;
     return out;
   }else{
     std::string msg = "[Phys input] Key missing from file, key=";
     throw msg+ key;    
   }
 }

 template<> bool ConfigFile::get<bool>(const std::string& key) const{
  std::istringstream iss(configMap.find(key)->second);
  bool result(false);
  iss >> result;
  if (iss.fail()){
    iss.clear();
    iss >> std::boolalpha >> result;
  }
    return result;
}
