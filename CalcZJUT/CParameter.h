#ifndef CPARAMETER_H
#define CPARAMETER_H
/*   format of input.dat
#:   comment label;

*/
#include <string>
#include <vector>
#include <map>
#include "CGaparameter.h"
namespace CALCZJUT{

    namespace GAZJUT{
       class CGaparameter;
    }

class CParameter
{
    public:
        typedef void (CParameter::*cmdFun)(std::string);

        typedef enum EM_P_01{ENERGY=0x101,
                              FORCE=0x102} EVALUATOR_CRIT;
        typedef enum EM_P_02{CLUSTER=0x103,
                             MOL_CLUSTER=0x104,
                             MOL_2DMATERIAL=0x105} SIMULATION_MODE;

        CParameter(char* filename);
        virtual ~CParameter();

        void input();
        GAZJUT::CGaparameter* GaParameter();

        //Parameter of general setting
        std::string sysName;
        std::vector<std::string>   runCmd;

        std::vector<std::pair<std::string,std::string>> excludeBond;
        std::pair<double,double> bondToleranceFactor;
        SIMULATION_MODE simulationMode;
        EVALUATOR_CRIT  evaluatorCriterion;

        std::string supportStructFile;
        std::string adsorbentStructFile;
        std::string adso_supp_Struct;

        std::string output_struct_format;
    protected:
        void setSysName(std::string);
        void setPopSize(std::string);
        void setPm(std::string);
        void setPc(std::string);
        void setGenNum(std::string);
        void setRunCmd(std::string);
        void setBond_Tolerance_Factor(std::string);
        void setExclude_Bond(std::string);
        void setEvaluator_Code(std::string);
        void setEvaluator_Criterion(std::string);
        void setSimulationMode(std::string);
        void setScaling_Mode(std::string);
        void setMutation_Mode(std::string);
        void setGene_Code(std::string);
        void setCross_Mode(std::string);
        void setSelect_Mode(std::string);
        void setGene_Formation_Mode(std::string);
        void setOutput_struct_format(std::string);

        bool checkIsValidParameters();

    private:
        std::string *m_pfile;
        GAZJUT::CGaparameter* m_pGAParameter;
        std::map<std::string,cmdFun> m_mapCmdFunc;
};


}
#endif // CPARAMETER_H
