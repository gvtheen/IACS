#ifndef CIOMOL_H
#define CIOMOL_H
#include <string>
#include "CIOBase.h"
namespace CATAZJUT{
   class CPeriodicFramework;
}

namespace CALCZJUT{
class CIOBase;

class CIOMol:public CIOBase
{
    public:
        CIOMol(CATAZJUT::CPeriodicFramework*);
        virtual ~CIOMol();

        void output(const std::string& file);
        void  input(std::string  file="");
        Bitset input(std::string file,CALCZJUT::CParameter::SIMULATION_MODE);
    protected:

    private:

};


}
#endif // CIOMOL_H
