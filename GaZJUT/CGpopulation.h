#ifndef CGPOPULATION_H
#define CGPOPULATION_H
#include <string>
#include <iostream>
#include <vector>
#include "CGaparameter.h"
#include "CGenome.h"
#include "GaDeclaration.h"

namespace GAZJUT{

/*
   The class have two points:
   1.   construct m_Gpopulation array, which is CGenome object.
*/
class CGpopulation
{
	public:
		CGpopulation();
		CGpopulation(CGaparameter*);
		CGpopulation(CGpopulation&);
		~CGpopulation();

		//member function
	//	void init();
		CGpopulation* clone();
		void copyBestGenome();
		void updatePopulation();

		void raw_statistic();
		void fitness_statistic();
		void asscendSort();
		void descendSort();
		//input output
		size_t  popNum();
        void setPopNum(size_t);

        std::vector <GeneVAR>& GeneVARArray();
        void setGeneVARArray(std::vector <GeneVAR>&);

        void modifyPopulation(std::vector <CGenome*>,size_t pos);
		//operator function
		CGenome* operator[](const size_t);             //return   pointer of genome of x th
        double   operator[](std::string );          //return   minraw, maxraw, minfit,maxfit of population
    public:
        CGaparameter            *m_pObjGaparameter;
             CGenome            *m_pMinRawGenome;      // pointer of address of m_pGpopulation
		     CGenome            *m_pMaxRawGenome;
             CGenome            *m_pMinFitGenome;
		     CGenome            *m_pMaxFitGenome;
		     CGenome            *m_pMaxOriGenome;
		     CGenome            *m_pMinOriGenome;

		std::vector <CGenome*>   m_Gpopulation;
    protected:
		// analyse the score of current Population
		     double              m_MinRawScore;
		     double              m_MaxRawScore;
		     double              m_AvgRawScore;
		     double              m_DevrawScore;
		// analyse the fitness of current Population
		     double              m_MinFitness;
		     double              m_MaxFitness;
		     double              m_AvgFitness;

		     double              m_MinOriValue;
		     double              m_MaxOriValue;
};

}
#endif
