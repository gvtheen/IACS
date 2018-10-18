#ifndef CIOS_H
#define CIOS_H
#include<string>
#include "../Util/Bitset.h"
#include "CParameter.h"
using util::Bitset;

namespace CATAZJUT{
   class CPeriodicFramework;;
}
namespace CALCZJUT{

class CIOBase
{
    public:
        CIOBase(CATAZJUT::CPeriodicFramework*);
        virtual ~CIOBase();

        virtual void output(const std::string& file)=0;
        virtual void input(std::string  file="")=0;
        virtual Bitset input(std::string file,CALCZJUT::CParameter::SIMULATION_MODE);
        void setConfiguration(CATAZJUT::CPeriodicFramework*);

    public:
        CATAZJUT::CPeriodicFramework*  m_pPeriodicFramework;

};


}
#endif // COUTPUTSTRUCTURE_H
