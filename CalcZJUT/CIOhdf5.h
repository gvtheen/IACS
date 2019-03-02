#ifndef CIOHDF5_H
#define CIOHDF5_H

#include<string>

namespace CATAZJUT{
   class CConfigurationBase;
   class CAtom;
}
namespace CALCZJUT{

class CIOhdf5
{
    public:
        CIOhdf5(CATAZJUT::CConfigurationBase*);
        virtual ~CIOhdf5();

        void input(std::string file);
        void output(std::string file);

    protected:

    private:
        CATAZJUT::CConfigurationBase* currentPeriodicFramework;
};


}
#endif // CIOHDF5_H
