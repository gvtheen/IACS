#ifndef GAUTILITYFUNCTION_H_INCLUDED
#define GAUTILITYFUNCTION_H_INCLUDED
#include<string>
#include<vector>
#include "GaDeclaration.h"
#include "../Util/Bitset.h"
#include "../IACS.h"

using util::Bitset;
using IACSZJUT::VarRangeStruct;
namespace GAZJUT{


double binaryDecode(Bitset&,VarRangeStruct);
int    calcBitNum(VarRangeStruct);
void   grayTobit(Bitset&);
void   bitTogray(Bitset&);


}


#endif // GAUTILITYFUNCTION_H_INCLUDED
