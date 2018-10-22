#ifndef CCOMPUTATEFITNESSINTERFACE_H
#define CCOMPUTATEFITNESSINTERFACE_H
#include<vector>
#include<iostream>
#include<string>


namespace CATAZJUT{
  class CPeriodicFramework;
}


namespace CALCZJUT{

class CModelBase;
class CIOBase;
class CParameter;

class CExeFitnessInterface
{
	public:
		 CExeFitnessInterface(CParameter*);
		 virtual ~CExeFitnessInterface();
		//virtual function:  interface function
		 virtual CExeFitnessInterface* clone();

		 virtual void init();
		 virtual double CalcuRawFit(std::vector<double>& RealValueOfGenome,size_t& pop_index, bool& isNormalexist);
         virtual void   ConvOrigToRawScore(std::vector<double>&);

         virtual char* ExeName()=0;

         void setCalcModeStruct(CModelBase* Temp_calcModeStruct);
         CModelBase* calcModeStruct();

         void setIO(CIOBase* m_IO);
         CIOBase* IO()const;

         CIOBase* getIO(std::string &file_name,CATAZJUT::CPeriodicFramework* currentPeriodicFramework);

//		 virtual void   outputResult()=0;

	protected:
	     // cluster, molecule/cluster, molecule/2D_support
	               CModelBase  *m_pCalcModeStruct;
                      CIOBase  *m_pIO;
                   CParameter  *m_Parameter;

};


}
#endif
