#ifndef CIOXYZ_H
#define CIOXYZ_H

#include <string>
#include "CIOBase.h"
#include "../Util/Bitset.h"

using util::Bitset;

namespace CATAZJUT{
   class CConfigurationBase;;
}

namespace CALCZJUT{

class CIOXyz:public CIOBase
{
    public:
        CIOXyz(CATAZJUT::CConfigurationBase*);
        virtual ~CIOXyz();

        void output(const std::string& file);
        void  input(std::string file="");
        Bitset input(std::string file,CALCZJUT::CParameter::SIMULATION_MODE mode);
    private:
};


}
#endif // CIOXYZ_H
