#ifndef CCOMPUTATEFITNESSINTERFACE_H
#define CCOMPUTATEFITNESSINTERFACE_H
#include<vector>
#include<iostream>
#include<string>
#include "../Util/Bitset.h"

namespace CATAZJUT{
  class CPeriodicFramework;
}

using util::Bitset;
namespace CALCZJUT{

class CCalcModeStruct;
class Cios;
class CParameter;

class CCalcFitnessInterface
{
	public:
		 CCalcFitnessInterface(CParameter*);
		 virtual ~CCalcFitnessInterface();
		//virtual function:  interface function
		 virtual CCalcFitnessInterface* clone()=0;

		 virtual void init()=0;
		 virtual double CalcuRawFit(std::vector<double>* RealValueOfGenome,size_t& pop_index, bool& isNormalexist)=0;
         virtual void   ConvOrigToRawScore(std::vector<double>*)=0;

         void setCalcModeStruct(CCalcModeStruct* Temp_calcModeStruct);
         CCalcModeStruct* calcModeStruct();

         void setIO(Cios* m_IO);
         Cios* IO()const;

         Cios* getIO(std::string &file_name,CATAZJUT::CPeriodicFramework* currentPeriodicFramework);

//		 virtual void   outputResult()=0;

	protected:
	     // cluster, molecule/cluster, molecule/2D_support
	               CCalcModeStruct  *m_pCalcModeStruct;
	                          Cios  *m_pIO;
	                    CParameter  *m_Parameter;
	                    Bitset       pop_run_state;
};


}
#endif
