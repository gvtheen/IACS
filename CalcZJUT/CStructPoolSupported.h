#ifndef CCALCSUPPORTSTRUCTPOOL_H
#define CCALCSUPPORTSTRUCTPOOL_H
#include <string>
#include "CStructPoolBase.h"

namespace CATAZJUT{
  class CConfigurationBase;
  class CElement;
}

namespace CALCZJUT{

class CStructPoolSupported:public CStructPoolBase
{
    public:
        CStructPoolSupported(CParameter*);
        virtual ~CStructPoolSupported();

        void init();
        void VarRangeStructRange(std::vector<VarRangeStruct>&);


    protected:
        void getIO(std::string &file_name,CATAZJUT::CConfigurationBase* currentPeriodicFramework);
    private:
        CATAZJUT::CConfigurationBase*  copy_pPeriodicFramework;
};



}//namespace
#endif // CCALCSUPPORTSTRUCTPOOL_H
