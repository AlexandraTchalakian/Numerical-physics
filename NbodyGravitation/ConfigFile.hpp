#ifndef TEST_CONFIGFILE_H_INCLUDED
#define TEST_CONFIGFILE_H_INCLUDED

#include <map>

using namespace std;

class ConfigFile{
  
  public:
      ConfigFile(const std::string& filename);
      ~ConfigFile();

      template<typename T> T get(const std::string& key) const;
      
      void process(const std::string& lineread);
      
      std::string toString();
      
      void printOut(std::string path);
      
  private:
      std::map<std::string, std::string> configMap;
};

#endif //TEST_CONFIGFILE_H_INCLUDED