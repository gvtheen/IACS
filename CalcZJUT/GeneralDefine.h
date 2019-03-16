#ifndef GENERALIZEDEFINE_H
#define GENERALIZEDEFINE_H


namespace CALCZJUT{

typedef enum EMMMM01{VASP_CALC=0x6011,
           DMOL_CALC==0x6012,
           GAUSSIAN_CALC==0x6013}FIT_CALCULATOR;

typedef enum EMMMM02{CLUSTER=0x6021,
             ADSORBENT_CLUSTER=0x6022,
             ADSORBENT_2D=0x6023}  WORK_MODEL;

typedef enum EM06{VASP=0x6001,
                  GAUSSIAN=0x6002,
                  DMOL=0x6003,
                  LAMMPS=0x6004,
                  CASTEP=0x6005,
                  DFTB=0x6006} E_EVALUATOR_EXE;


}




#endif

