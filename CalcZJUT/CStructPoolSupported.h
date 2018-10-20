#ifndef CCALCSUPPORTSTRUCTPOOL_H
#define CCALCSUPPORTSTRUCTPOOL_H
#include <string>
#include "CStructPoolBase.h"

namespace CATAZJUT{
  class CPeriodicFramework;
  class CElement;
}

namespace CALCZJUT{

class CStructPoolSupported:public CStructPoolBase
{
    public:
        CStructPoolSupported(CParameter*);
        virtual ~CStructPoolSupported();

        void init();
        void GeneVARRange(std::vector<GeneVAR>&);


    protected:
        void getIO(std::string &file_name,CATAZJUT::CPeriodicFramework* currentPeriodicFramework);
    private:
        CATAZJUT::CPeriodicFramework*  copy_pPeriodicFramework;
};



}//namespace
#endif // CCALCSUPPORTSTRUCTPOOL_H
