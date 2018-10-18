#ifndef CIOCAR_H
#define CIOCAR_H

#include <string>
#include "CIOBase.h"
#include "../Util/Bitset.h"

using util::Bitset;

namespace CATAZJUT{
   class CPeriodicFramework;
   class CAtom;
}
namespace CALCZJUT{

class CIOCar:public CIOBase
{
    public:
        CIOCar(CATAZJUT::CPeriodicFramework*);
        virtual ~CIOCar();

        void output(const std::string& file);
         void input(std::string file="");
         Bitset input(std::string file,CALCZJUT::CParameter::SIMULATION_MODE);

    protected:

    private:
};

}
#endif // CIOCAR_H
