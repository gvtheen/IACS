#ifndef CIOGJF_H
#define CIOGJF_H

#include<vector>
#include<string>
#include "CIOBase.h"
#include "../Util/Bitset.h"

using util::Bitset;

namespace CATAZJUT{
   class CConfigurationBase;
}

namespace CALCZJUT{

class CIOGjf:public CIOBase
{
    public:
        CIOGjf(CATAZJUT::CConfigurationBase*);
        virtual ~CIOGjf();

        void output(const std::string& file);
         void input(std::string file="");
         Bitset input(std::string file,CALCZJUT::CParameter::SIMULATION_MODE);
         void setCommandlines(const std::vector<std::string>&);
         std::vector<std::string> Commandlines() const;

    protected:

    private:
        std::vector<std::string> m_Commandline;
        std::vector<int> m_Constranit;  // 0,-1;

};


}
#endif // CIOGJF_H
