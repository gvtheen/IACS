#ifndef CCATACALCADAPTER_H
#define CCATACALCADAPTER_H

namespace CALCZJUT{

class CExeFitnessInterface;
class CParameter;
class CStructPoolBase;

class CCataCalcAdapter:public CAdapterBase
{
    public:
        CCataCalcAdapter(CALCZJUT::CParameter*);
        virtual ~CCataCalcAdapter();

    public:
        void engineToComputation(std::vector<double>& _input,
                                 size_t _index,
                                 bool& _state,
                                 double& _output);
        void computationToEngine(std::vector<IACSZJUT::VarRangeStruct>& _mvalue);
        void computationToEngine(std::vector<double>& _mvalue);

    private:
        double convertOrigValueToRawvalue(const double&);

    private:
        std::vector<CALCZJUT::CExeFitnessInterface*>  m_FitnessCalculatorPool;
                    CALCZJUT::CExeFitnessInterface   *m_currentEvaluator;
        CALCZJUT::CStructPoolBase                    *m_pStructurePool;
	    CALCZJUT::CParameter                         *m_pParameter;
	    std::vector<IACSZJUT::VarRangeStruct>         m_VarRangeStruct;
	    double m_
};


}
#endif // CCATACALCADAPTER_H
