#ifndef CIOS_H
#define CIOS_H
#include<string>
#include "../Util/Bitset.h"
#include "CParameter.h"
using util::Bitset;

namespace CATAZJUT{
   class CConfigurationBase;;
}
namespace CALCZJUT{

class CIOBase
{
    public:
        CIOBase(CATAZJUT::CConfigurationBase*);
        virtual ~CIOBase();

        virtual void output(const std::string& file)=0;
        virtual void input(std::string  file="")=0;
        virtual Bitset input(std::string file,CALCZJUT::CParameter::SIMULATION_MODE);
        void setConfiguration(CATAZJUT::CConfigurationBase*);

    public:
        CATAZJUT::CConfigurationBase*  m_pPeriodicFramework;

};


}
#endif // COUTPUTSTRUCTURE_H
