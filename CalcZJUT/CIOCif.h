#ifndef CIOCIF_H
#define CIOCIF_H

#include <string>
#include "CIOBase.h"
#include "../Util/Bitset.h"

using util::Bitset;

namespace CATAZJUT{
   class CConfigurationBase;;
}
namespace CALCZJUT{

class CIOCif:public CIOBase
{
    public:
        CIOCif(CATAZJUT::CConfigurationBase*);
        virtual ~CIOCif();

        void output(const std::string& file);
        void  input(std::string file="");
        Bitset input(std::string file,CALCZJUT::CParameter::SIMULATION_MODE mode);
    protected:

    private:


};


}
#endif // CIOCIF_H
