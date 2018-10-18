#ifndef CIOCELLFILE_H
#define CIOCELLFILE_H
#include <string>
#include "CIOBase.h"

namespace CATAZJUT{
   class CPeriodicFramework;;
}
namespace CALCZJUT{

class CIOCellFile:public CIOBase
{
    public:
        CIOCellFile(CATAZJUT::CPeriodicFramework*);
        virtual ~CIOCellFile();

        void output(const std::string& file);
        void  input(std::string file="");
        Bitset input(std::string file,CALCZJUT::CParameter::SIMULATION_MODE);

    protected:

    private:
};


}
#endif // CIOCELLFILE_H
