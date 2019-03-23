#ifndef UTILITY_H
#define UTILITY_H

namespace GAZJUT{


typedef enum EM01{
                  MIN=0x1000,
                  MAX=0x1001 }           E_GA_TYPE;

typedef enum EM02{SINGLE=0x2002,
                  MULTIPLE=0x2003,
                  UNIFORM_C=0x2004,
                  ARITHMETIC=0x2005,
                  UNARITHMETIC=0x2001 }  E_CROSSOVER_OPERATOR;

typedef enum EM03{UNIFORM_M=0x3001,
                  BOUNDARY=0x3002,
                  NOUNIFORM=0x3003,
                  GAUSSIAN_M=0x3004  }  E_MUTATE_OPERATOR;

typedef enum EM04{ROULETTE_WHEEL=0x4001,
                  RANDOM=0x4002,
                  TOURNAMENT=0x4003,
                  MIXED=0x4004  }       E_SELECT_OPERATOR;

typedef enum EM05{RANDOM_FORMATION=0x5001,
                  FILE_INPUT=0x5002 }   E_GENEFORMATION_TYPE;

typedef enum EM07{LINEAR=0x7001,
                  SIGMA=0x7002,
                  POWER=0x7003}         E_SCALING_TYPE;
typedef enum EM08{BINARY=0x8001,
                  GRAY=0x8002,
                  REAL=0x8003}          E_CODE_TYPE;


}


#endif
