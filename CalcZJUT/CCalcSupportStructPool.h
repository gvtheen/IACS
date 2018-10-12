#ifndef CCALCSUPPORTSTRUCTPOOL_H
#define CCALCSUPPORTSTRUCTPOOL_H
#include <string>
#include "CCalcStructureBasePool.h"

namespace CATAZJUT{
  class CPeriodicFramework;
  class CElement;
}

namespace CALCZJUT{

class CCalcSupportStructPool:public CCalcStructureBasePool
{
    public:
        CCalcSupportStructPool(CParameter*);
        virtual ~CCalcSupportStructPool();

        void init();
        void GeneVARRange(std::vector<GeneVAR>&);


    protected:
        void getIO(std::string &file_name,CATAZJUT::CPeriodicFramework* currentPeriodicFramework)
    private:
};



}//namespace
#endif // CCALCSUPPORTSTRUCTPOOL_H
