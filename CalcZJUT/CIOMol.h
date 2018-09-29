#ifndef CIOMOL_H
#define CIOMOL_H
#include <string>
#include "Cios.h"
namespace CATAZJUT{
   class CPeriodicFramework;
}

namespace CALCZJUT{
class Cios;

class CIOMol:public Cios
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
