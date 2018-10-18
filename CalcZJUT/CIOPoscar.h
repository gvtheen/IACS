#ifndef CIOPOSCAR_H
#define CIOPOSCAR_H
//#include "CatalystUniverseDefine.h"
#include "CIOBase.h"
#include "Bitset.h"

namespace CATAZJUT{
   class CPeriodicFramework;
}

namespace CALCZJUT{
class CIOBase;

class CIOPoscar:public CIOBase
{
    public:
        CIOPoscar(CATAZJUT::CPeriodicFramework*);
        virtual ~CIOPoscar();

        void output(const std::string& file);
        void  input(std::string file="");
        Bitset input(std::string file,CALCZJUT::CParameter::SIMULATION_MODE);

    protected:

    private:
        size_t atomicIndex(std::vector<size_t>& mht,size_t& atomNum);

    private:



};




}
#endif // CIOPOSCAR_H
