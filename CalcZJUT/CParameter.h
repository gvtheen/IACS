#ifndef CPARAMETER_H
#define CPARAMETER_H
/*   format of input.dat
#:   comment label;

*/
#include <string>
#include <vector>
#include <map>
#include <boost/filesystem.hpp>
#include "../CATAZJUT/CBondPrivate.h"
#include "../GaZJUT/CGaparameter.h"
namespace GAZJUT{
   class CGaparameter;
}
namespace CATAZJUT{
   class CBondPrivate;
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
                             MOL_CLUSTER   =0x202,
                             MOL_2DMATERIAL=0x203,
                      MOL_CLUSTER2DMATERIAL=0x204,
                                   PERIODIC=0x205} SIMULATION_MODE;
        typedef enum EM_P_03{WORK      = 0x301,
                             PARAMETER = 0x302,
                             SCRATCH   = 0x303,
                             MAINPATH  = 0x304}PATH_TYPE;
        //minmum or maxmum
        typedef enum EM_P_04{MIN = 0x401,
                             MAX = 0x402}SEARCH_MODE;

        typedef enum EM_P_05{GA = 0x501,
                             PSO = 0x502,
                             ABC = 0x503}OPT_ENGINE;

        //Parameter of general setting
        std::string sysName;
        std::vector<std::string>   runCmd;

        std::vector<std::pair<std::string*,std::string*>> excludeBond;
        std::pair<double,double> bondToleranceFactor;
        std::vector<CATAZJUT::CBondPrivate*> bondTolerance;
        SIMULATION_MODE simulationMode;
        EVALUATOR_CRITERION  evaluatorCriterion;
        SEARCH_MODE searchMode;

        // for supported catalyst
        std::string supportStructFile;
        std::string adsorbentStructFile;
        std::vector<std::string*> adso_supp_Input_File;
        std::pair<double,double> DistRange_Adsorbent_Support;
        // for pure cluster
        std::string cluster_Formula;
        std::vector<std::string*> cluster_Input_File;
        // for periodic structure
        std::string periodic_structure_Formula;
        std::vector<std::string*> periodic_structure_Input_File;

        std::vector<std::vector<char*>> evaluatorParameterFile;
        std::string output_struct_format;
        const char* Backup_Parameter_Folder;
        double optimal_gap_value;

    public:
        CParameter(std::string filename);
        virtual ~CParameter();


        void input();
        void output();

        void GetEvaluateEXE(std::vector<std::string>& res);

        //
        GAZJUT::CGaparameter* GaParameter();
        size_t currentGenerationNum();
        size_t popNum();
        // set working environment
        void initWorkEnvironment();
        void setCurrentWorkPathAt(size_t,size_t);
        void setCurrentWorkPathAt(PATH_TYPE);
        std::string currentWorkPath();
        std::string rootPath();
        std::string scratchPath();
        std::string parameterPath();
        void checkExeNecessaryFiles(std::vector<std::string*>& files);
        void moveFileToPath(const std::string& file, const std::string& dir);
        void copyFileToPath(const std::string& file, const std::string& dir);

    private:
        // system command
        void setSysName(std::string);
        void setRunCmd(std::string);
        void setSimulationMode(std::string);
        void setEvaluator_Code(std::string);
        void setEvaluator_Criterion(std::string);
        void setSearch_Mode(std::string);
        void setOptEngineMethod(std::string);
        // structure command
        void setBond_Tolerance_Factor(std::string);   // Cu  O
        void setExclude_Bond(std::string);
        void setBond_Tolerance(std::string);  // Cu  O  1.21   2.1
        // GA parameters
        void setPopSize(std::string);
        void setPm(std::string);
        void setPc(std::string);
        void setGenNum(std::string);
        void setScaling_Mode(std::string);
        void setMutation_Mode(std::string);
        void setGene_Code(std::string);
        void setCross_Mode(std::string);
        void setCross_Num(std::string);
        void setSelect_Mode(std::string);
        void setGene_Formation_Mode(std::string);
        //
        void setSupport_Structure(std::string);
        void setAdsorbent_Structure(std::string);
        void setAdsorbent_Support_Structure(std::string);
        void setCluster_Formula(std::string);
        void setCluster_Input_File(std::string);
        void setDistanceRangeAdsorbentOnSupport(std::string);
        void setOutputStructureFormat(std::string);

        void setOutput_struct_format(std::string);

        bool checkIsValidParameters();


    private:
        std::string *m_pfile;
        GAZJUT::CGaparameter* m_pGAParameter;
        std::map<std::string,cmdFun> m_mapCmdFunc;

        std::vector<std::string> m_StrEvaluateEXE;

        boost::filesystem::path m_root_WorkingPath;
};


}
#endif // CPARAMETER_H
