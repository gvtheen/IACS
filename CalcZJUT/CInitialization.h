#ifndef CINITIALIZATION_H
#define CINITIALIZATION_H


namespace CATAZJUT{
  class CPeriodicFramework;
}

namespace CALCZJUT{

class CCalcModeStruct;
class CParameter;
class Cios;

class CInitialization
{
    public:
        CInitialization(CParameter* mpara);
        virtual ~CInitialization();

        CCalcModeStruct* constructCalcModeStruct();

        void Init2DSupport();
        void InitCluster();
        void InitClusterSupport();

    protected:
        void getIO(std::string &file_name,CATAZJUT::CPeriodicFramework* currentPeriodicFramework);

        void InitClusterFromChemFormula();
    private:
        std::vector<CCalcModeStruct*>  m_PoolCalcModeStruct;
                     CParameter       *m_Parameter;
                     Cios             *m_IO;
};


}
#endif // CINITIALIZATION_H
