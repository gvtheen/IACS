#ifndef CPARAMETER_H
#define CPARAMETER_H
/*   format of input.dat
#:   comment label;

*/
#include <string>
#include <vector>
#include <map>
#include "../GaZJUT/CGaparameter.h"
namespace GAZJUT{
       class CGaparameter;
}
namespace CALCZJUT{

class CParameter
{
    public:
        //pointer of function
        typedef void (CParameter::*cmdFun)(std::string);

        typedef enum EM_P_01{ ENERGY  =0x101,
                              FORCE   =0x102,
                              BAND_GAP=0x103} EVALUATOR_CRITERION;

        typedef enum EM_P_02{CLUSTER       =0x201,
                             PERIODIC      =0x202,
                             MOL_CLUSTER   =0x203,
                             MOL_2DMATERIAL=0x204,
                      MOL_CLUSTER2DMATERIAL=0x205} SIMULATION_MODE;

        CParameter(char* filename);
        virtual ~CParameter();

        void input();
        void output();
        GAZJUT::CGaparameter* GaParameter();
        size_t currentGenerationNum();
        size_t popNum();

    public:
        //Parameter of general setting
        std::string sysName;
        std::vector<std::string>   runCmd;

        std::vector<std::pair<std::string*,std::string*>> excludeBond;
        std::pair<double,double> bondToleranceFactor;
        SIMULATION_MODE simulationMode;
        EVALUATOR_CRITERION  evaluatorCriterion;
        // for supported catalyst
        std::string supportStructFile;
        std::string adsorbentStructFile;
        std::vector<std::string*> adso_supp_Input_File;
        // for pure cluster
        std::string cluster_Formula;
        std::vector<std::string*> cluster_Input_File;


        std::string output_struct_format;

        double optimal_gap_value;

//        size_t PopSize()const;
//        double Pm()const;
//        double Pc()const;
//        size_t GeneNum();

    protected:
        // system command
        void setSysName(std::string);
        void setRunCmd(std::string);
        void setSimulationMode(std::string);
        void setEvaluator_Code(std::string);
        void setEvaluator_Criterion(std::string);
        // structure command
        void setBond_Tolerance_Factor(std::string);
        void setExclude_Bond(std::string);
        // GA parameters
        void setPopSize(std::string);
        void setPm(std::string);
        void setPc(std::string);
        void setGenNum(std::string);
        void setScaling_Mode(std::string);
        void setMutation_Mode(std::string);
        void setGene_Code(std::string);
        void setCross_Mode(std::string);
        void setSelect_Mode(std::string);
        void setGene_Formation_Mode(std::string);
        //
        void setSupport_Structure(std::string);
        void setAdsorbent_Structure(std::string);
        void setAdsorbent_Support_Structure(std::string);
        void setCluster_Formula(std::string);
        void setCluster_Input_File(std::string);

        void setOutput_struct_format(std::string);

    protected:

        bool checkIsValidParameters();

    private:
        std::string *m_pfile;
        GAZJUT::CGaparameter* m_pGAParameter;
        std::map<std::string,cmdFun> m_mapCmdFunc;
};


}
#endif // CPARAMETER_H
