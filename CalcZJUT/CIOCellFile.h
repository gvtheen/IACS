#ifndef CIOCELLFILE_H
#define CIOCELLFILE_H
#include <string>
#include <vector>
#include "CIOBase.h"

namespace CATAZJUT{
   class CConfigurationBase;;
}
namespace CALCZJUT{

class CIOCellFile:public CIOBase
{
    public:
        CIOCellFile(CATAZJUT::CConfigurationBase*);
        virtual ~CIOCellFile();

        void output(const std::string& file);
        void  input(std::string file="");
        Bitset input(std::string file,CALCZJUT::CParameter::SIMULATION_MODE);

    protected:

    private:
        std::map<std::string,std::string> m_pseudoPotentialFiles;
};


}
#endif // CIOCELLFILE_H
