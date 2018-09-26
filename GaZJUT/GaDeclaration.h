#ifndef UTILITY_H
#define UTILITY_H
//#define NDEBUG
#include "assert.h"

namespace GAZJUT{

//typedef enum EM01{MIN, MAX}                               E_GA_TYPE;
//typedef enum EM02{SINGLE,MULTIPLE,UNIFORM_C,ARITHMETIC,UNARITHMETIC}  E_CROSSOVER_OPERATOR;
//typedef enum EM03{UNIFORM_M,BOUNDARY,NOUNIFORM,GAUSSIAN_M}  E_MUTATE_OPERATOR;
//typedef enum EM04{ROULETTE_WHEEL,RANDOM,TOURNAMENT,MIXED} E_SELECT_OPERATOR;
//typedef enum EM05{RANDOM_FORMATION,FILE_INPUT}            E_GENEFORMATION_TYPE;
//typedef enum EM06{VASP,GAUSSIAN,DMOL,LAMMPS}              E_EVALUATOR_EXE;
//typedef enum EM07{LINEAR,SIGMA,POWER}                     E_SCALING_TYPE;
//typedef enum EM08{BINARY,GRAY,REAL}                       E_CODE_TYPE;

//macro definition

// struct definition

/*
Varible of gene
*/
typedef struct str_001{
	double min;
	double max;
	double accuracy;
}GENEVAR;

}


#endif
