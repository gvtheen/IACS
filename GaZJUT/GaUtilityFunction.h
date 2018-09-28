#ifndef GAUTILITYFUNCTION_H_INCLUDED
#define GAUTILITYFUNCTION_H_INCLUDED
#include<string>
#include<vector>
#include "GaDeclaration.h"

namespace GAZJUT{

double binaryDecode(std::vector<unsigned int>*,GENEVAR);
int    calcBitNum(GENEVAR);
void   grayTobit(std::vector<unsigned int>*);
void   bitTogray(std::vector<unsigned int>* );

//void   ERROR_OUTPUT(std::string );
//void   ERROR_OUTPUT(std::string,std::string );
//void   ERROR_OUTPUT(std::string,std::string,std::string );
//void   ERROR_OUTPUT(std::string,std::string,std::string,std::string);

}


#endif // GAUTILITYFUNCTION_H_INCLUDED
