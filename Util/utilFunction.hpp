#ifndef UTILFUNCTION_H
#define UTILFUNCTION_H

#include <cstddef>
#include "Point-Vector.h"
#include "../CataZJUT/Geometry.h"
#include "../GaZJUT/GaDeclaration.h"
#include "Bitset.h"

using GAZJUT::GeneVAR;

namespace util{
//
Vector4 SphereEquationFromPoints(const std::vector<Point3>& coordinate);
double binaryDecode(const Bitset & myCode,GeneVAR myGenVar);
int calcBitNum(GeneVAR myGeneVAR);
void grayTobit(Bitset& data);
void bitTogray(Bitset& data);



}



#endif
