#ifndef CGAPARAMETER_H
#define CGAPARAMETER_H
#include<string>
#include<map>
#include<vector>
#include "GaDeclaration.h"
namespace GAZJUT{

class CGaparameter
{
    public:

        CGaparameter();
        virtual ~CGaparameter();
        CGaparameter(std::vector <GeneVAR>* myVar);
        CGaparameter(const CGaparameter& other);
        CGaparameter& operator=(const CGaparameter& other);

        void   defaultInit();
        void   OutputTofile();
        int    GenerationNum();
        int    PopNum();
        int    CrossNum();
        void   add_Curr_Generation();
        //void   SetPopNum(int);
        double CrossProb();
        //void   SetCrossProb(double);
        double MutaProb();
        //void   SetMutaProb(double);

        std::vector <GeneVAR>* GeneVAR();
        void                   setGeneVAR(std::vector <GeneVAR>*);
        void                   checkGeneVAR();
        std::string            GeneFile();
        std::string            InputFile();
        //void   SetGeneFile(string);
        E_GA_TYPE              SearchType();
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
        int          Curr_Generation;
        const double ScaleLinearMultiplier     = 1.2;
        const double ScaleSigmaTruncMultiplier = 2.0;
        const double ScalePowerLawFactor       = 1.0005;
        const double ScaleBoltzMinTemp         = 1.0;
        const double ScaleBoltzFactor          = 0.05;
        const double ScaleBoltzStart           = 40.0;

        const double UnifArithmCrossConstant   = 0.25;
        const double NoUnifArithmCrossConstant = 0.10;
    protected:
        std::map<std::string, std::string>      *m_mapCmdString;
    private:
        std::vector <GeneVAR>   *m_pGeneVARofPopulation;
        int                      m_PopNum;
        int                      m_GenerationNum;
        int                      m_CrossNum;
        double                   m_CrossProb;
	    double                   m_MutaProb;
	    std::string              m_GeneFile;
        std::string              m_SearchType;
	    std::string              m_CrossMode;
	    std::string              m_MutateMode;
	    std::string              m_SelectMode;
	    std::string              m_IntilGenMode;
	    std::string              m_EvaluateEXE;
	    std::string              m_ScalingMode;
	    std::string              m_CodeMode;
};

}
#endif // CGAPARAMETER_H
