#ifndef CINITIALIZE_H
#define CINITIALIZE_H
#include "GaDeclaration.h"
#include "GaUtilityFunction.h"
namespace GAZJUT{


class CInitialize
{
	public:
	    typedef void (*cmdFun)(std::string);
		CInitialize();
		~CInitialize();
	protected:
	    std::map<std::string,cmdFun> m_mapCmdFunc;
};



}
#endif
