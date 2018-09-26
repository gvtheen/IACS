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
        typedef enum EM01{
                         MIN=0x1000,
                         MAX=0x1001 }            E_GA_TYPE;
        typedef enum EM02{SINGLE=0x1010,
                          MULTIPLE=0x1011,
                          UNIFORM_C=0x1100,
                          ARITHMETIC=0x1101,
                          UNARITHMETIC=0x1110 }  E_CROSSOVER_OPERATOR;
        typedef enum EM03{UNIFORM_M=0x1111,
                          BOUNDARY=0x10000,
                          NOUNIFORM=0x10001,
                          GAUSSIAN_M=0x10010  }  E_MUTATE_OPERATOR;
        typedef enum EM04{ROULETTE_WHEEL=0x10011,
                          RANDOM=0x10100,
                          TOURNAMENT=0x10101,
                          MIXED=0x10110  }       E_SELECT_OPERATOR;
        typedef enum EM05{RANDOM_FORMATION=0x10111,
                          FILE_INPUT=0x11000 }   E_GENEFORMATION_TYPE;
        typedef enum EM06{VASP=0x11001,
                          GAUSSIAN=0x11010,
                          DMOL=0x11011,
                          LAMMPS=0x11100,
                          CASTEP=0x100011 }      E_EVALUATOR_EXE;
        typedef enum EM07{LINEAR=0x11101,
                          SIGMA=0x11110,
                          POWER=0x11111}         E_SCALING_TYPE;
        typedef enum EM08{BINARY=0x100000,
                          GRAY=0x100001,
                          REAL=0x100010}         E_CODE_TYPE;

        CGaparameter();
        virtual ~CGaparameter();
        CGaparameter(std::vector <GENEVAR>* myVar);
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

        std::vector <GENEVAR>* GeneVar();
        void                   setGeneVar(std::vector <GENEVAR>*);
        void                   checkGeneVar();
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
        std::vector <GENEVAR>    *m_pGeneVarofPopulation;
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
