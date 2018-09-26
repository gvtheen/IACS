#ifndef CIOS_H
#define CIOS_H
#include<string>
#include "Bitset.h"
#include "CParameter.h"
namespace CATAZJUT{
   class CPeriodicFramework;;
}
namespace CALCZJUT{

class Cios
{
    public:
        Cios(CATAZJUT::CPeriodicFramework*);
        virtual ~Cios();

        virtual void output(std::string& file)=0;
        virtual void input(std::string  file="")=0;
        virtual Bitset input(std::string file,CALCZJUT::CParameter::SIMULATION_MODE);
        void setConfiguration(CATAZJUT::CPeriodicFramework*);

    public:
        CATAZJUT::CPeriodicFramework*  m_pPeriodicFramework;

};


}
#endif // COUTPUTSTRUCTURE_H
