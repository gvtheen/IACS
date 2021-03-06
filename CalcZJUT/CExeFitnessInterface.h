#ifndef CCOMPUTATEFITNESSINTERFACE_H
#define CCOMPUTATEFITNESSINTERFACE_H
#include<vector>
#include<iostream>
#include<string>


namespace CATAZJUT{
  class CConfigurationBase;
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
		 virtual CExeFitnessInterface* clone()=0;
		 virtual void init();
		 virtual double CalcuRawFit(std::vector<double>& RealValueOfGenome,size_t& pop_index, bool& isNormalexist);
         virtual void ConvOrigToRawScore(std::vector<double>&);
         virtual std::string ExeName()=0;

         void setCalcModeStruct(CModelBase* Temp_calcModeStruct);
         CModelBase* calcModeStruct();

         void setIO(CIOBase* m_IO);
         CIOBase* IO()const;

         void getIO(std::string &file_name,CATAZJUT::CConfigurationBase* currentPeriodicFramework,CIOBase* ios);

//		 virtual void   outputResult()=0;

	protected:
	     // cluster, molecule/cluster, molecule/2D_support
	               CModelBase  *m_pCalcModeStruct;
                      CIOBase  *m_pIO;
                   CParameter  *m_Parameter;

};


}
#endif
