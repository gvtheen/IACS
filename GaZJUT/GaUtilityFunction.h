#ifndef GAUTILITYFUNCTION_H_INCLUDED
#define GAUTILITYFUNCTION_H_INCLUDED
#include<string>
#include<vector>
#include "GaDeclaration.h"
#include "../Util/Bitset.h"

using util::Bitset;

namespace GAZJUT{

double binaryDecode(Bitset&,GeneVAR);
int    calcBitNum(GeneVAR);
void   grayTobit(Bitset&);
void   bitTogray(Bitset&);

//void   ERROR_OUTPUT(std::string );
//void   ERROR_OUTPUT(std::string,std::string );
//void   ERROR_OUTPUT(std::string,std::string,std::string );
//void   ERROR_OUTPUT(std::string,std::string,std::string,std::string);

}


#endif // GAUTILITYFUNCTION_H_INCLUDED
