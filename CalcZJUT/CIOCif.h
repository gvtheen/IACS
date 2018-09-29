#ifndef CIOCIF_H
#define CIOCIF_H

#include <string>
#include "Cios.h"
#include "../Util/Bitset.h"

using util::Bitset;

namespace CATAZJUT{
   class CPeriodicFramework;;
}
namespace CALCZJUT{

class CIOCif:public Cios
{
    public:
        CIOCif(CATAZJUT::CPeriodicFramework*);
        virtual ~CIOCif();

        void output(const std::string& file);
        void  input(std::string& file="");
        Bitset input(std::string file,CALCZJUT::CParameter::SIMULATION_MODE mode);
    protected:

};


}
#endif // CIOCIF_H
