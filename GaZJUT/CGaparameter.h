#ifndef CGAPARAMETER_H
#define CGAPARAMETER_H
#include <string>
#include <map>
#include <vector>
#include "GaDeclaration.h"
namespace CALCZJUT{
   class CParameter;
}
namespace GAZJUT{

class CGaparameter
{
    public:

        CGaparameter();
        virtual ~CGaparameter();
        CGaparameter(std::vector <GAZJUT::GeneVAR>&  myVar);
        CGaparameter( CGaparameter& other);
        //operator
        CGaparameter& operator=( CGaparameter& other);

        void   defaultInit();

        size_t    GenerationNum();
        size_t    PopNum();
        size_t    CrossNum();
        void   add_Curr_Generation();
        //void   SetPopNum(int);
        double CrossProb();
        //void   SetCrossProb(double);
        double MutaProb();
        //void   SetMutaProb(double);

        std::vector <GeneVAR>& GeneVAR();
        void                   setGeneVAR(std::vector<GAZJUT::GeneVAR>&);
        void                   checkGeneVAR();
        std::string            GeneFile();
        std::string            InputFile();
        //void   SetGeneFile(string);
        GAZJUT::E_GA_TYPE      SearchType();
        //void   SetSearchType(E_GA_TYPE);
        E_SELECT_OPERATOR      SelectMode();

        E_CROSSOVER_OPERATOR   CrossMode();
        //void   SetCrossMode(E_CROSSOVER_OPERATOR);
        E_MUTATE_OPERATOR      MutateMode();
        //void   SetMutateMode(E_MUTATE_OPERATOR);
        E_GENEFORMATION_TYPE   InitGenMode();
        //void   SetIntilGenMode(E_GENEFORMATION_TYPE);
        E_EVALUATOR_EXE        EvaluateEXE();
        //void   SetEvaluateEXE(E_EVALUATOR_EXE);
        E_SCALING_TYPE         ScalingMode();
        //void   SetScalingMode(E_SCALING_TYPE);
        E_CODE_TYPE            CodeMode();
        //void   SetCodeMode(E_CODE_TYPE);
        //overload operator
        std::string& operator[](std::string key_name);
        //const parameter
        //scalling
        size_t       Curr_Generation;
        // scale parameters
        const double ScaleLinearMultiplier     = 1.2;
        const double ScaleSigmaTruncMultiplier = 2.0;
        const double ScalePowerLawFactor       = 1.0005;
        const double ScaleBoltzMinTemp         = 1.0;
        const double ScaleBoltzFactor          = 0.05;
        const double ScaleBoltzStart           = 40.0;

        const double UnifArithmCrossConstant   = 0.25;
        const double NoUnifArithmCrossConstant = 0.10;

        std::map<std::string, std::string>      *m_mapCmdString;
    private:
        friend class CALCZJUT::CParameter;

        std::vector<GAZJUT::GeneVAR>     m_GeneVARofPopulation;
//        size_t                   m_PopNum;
//        size_t                   m_GenerationNum;
//        size_t                   m_CrossNum;
//        double                   m_CrossProb;
//	    double                   m_MutaProb;
//	    std::string              m_GeneFile;
//        std::string              m_SearchType;
//	    std::string              m_CrossMode;
//	    std::string              m_MutateMode;
//	    std::string              m_SelectMode;
//	    std::string              m_IntilGenMode;
//	    std::string              m_EvaluateEXE;
//	    std::string              m_ScalingMode;
//	    std::string              m_CodeMode;
};

}
#endif // CGAPARAMETER_H
