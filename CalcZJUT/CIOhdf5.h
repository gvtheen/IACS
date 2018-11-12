#ifndef CIOHDF5_H
#define CIOHDF5_H

namespace CATAZJUT{
   class CPeriodicFramework;
   class CAtom;
}
namespace CALCZJUT{

class CIOhdf5
{
    public:
        CIOhdf5(CATAZJUT::CPeriodicFramework*);
        virtual ~CIOhdf5();

    protected:

    private:
};


}
#endif // CIOHDF5_H
