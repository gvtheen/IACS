#ifndef CIOGJF_H
#define CIOGJF_H

#include<vector>
#include<string>
#include "Cios.h"
#include "../Util/Bitset.h"

using util::Bitset;

namespace CATAZJUT{
   class CPeriodicFramework;
}

namespace CALCZJUT{

class CIOGjf:public Cios
{
    public:
        CIOGjf(CATAZJUT::CPeriodicFramework*);
        virtual ~CIOGjf();

        void output(const std::string& file);
         void input(std::string file="");
         Bitset input(std::string file,CALCZJUT::CParameter::SIMULATION_MODE);
    protected:

    private:
        std::vector<std::string> m_Commandline;
        std::vector<int> m_Constranit;  // 0,-1;

};


}
#endif // CIOGJF_H
