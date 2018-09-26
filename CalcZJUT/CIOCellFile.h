#ifndef CIOCELLFILE_H
#define CIOCELLFILE_H
#include <string>
#include "Cios.h"

namespace CATAZJUT{
   class CPeriodicFramework;;
}
namespace CALCZJUT{

class CIOCellFile:public Cios
{
    public:
        CIOCellFile(CATAZJUT::CPeriodicFramework*);
        virtual ~CIOCellFile();

        void output(std::string& file);
        void  input(std::string file="");
        Bitset input(std::string file,CALCZJUT::CParameter::SIMULATION_MODE);

    protected:

    private:
};


}
#endif // CIOCELLFILE_H
