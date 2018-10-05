#ifndef UTILITY_H
#define UTILITY_H
//#define NDEBUG
//#include "assert.h"

namespace GAZJUT{


typedef enum EM01{
                  MIN=0x1000,
                  MAX=0x1001 }           E_GA_TYPE;

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

/*
Varible of gene
*/
typedef struct str_001{
	double min;
	double max;
	double accuracy;
}GeneVAR;

}


#endif
