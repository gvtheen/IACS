#ifndef CIOCIF_H
#define CIOCIF_H

namespace CATAZJUT{
   class CPeriodicFramework;;
}
namespace CALCZJUT{

class CIOCif:public Cios
{
    public:
        CIOCif(CATAZJUT::CPeriodicFramework*);
        virtual ~CIOCif();

        void output(std::string& file);
        void  input(std::string& file="");
        Bitset input(std::string file,CALCZJUT::CParameter::SIMULATION_MODE mode);
    protected:

};


}
#endif // CIOCIF_H
