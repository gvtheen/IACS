#ifndef CIOCAR_H
#define CIOCAR_H

#include <string>
#include "Cios.h"
#include "Bitset.h"
namespace CATAZJUT{
   class CPeriodicFramework;
   class CAtom;
}
namespace CALCZJUT{


class CIOCar:public Cios
{
    public:
        CIOCar(CATAZJUT::CPeriodicFramework*);
        virtual ~CIOCar();

        void output(std::string& file);
         void input(std::string file="");
         Bitset input(std::string file,CALCZJUT::CParameter::SIMULATION_MODE);

    protected:

    private:
};

}
#endif // CIOCAR_H
