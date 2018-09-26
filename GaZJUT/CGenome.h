#ifndef CGENOME_H
#define CGENOME_H

#include <vector>
#include <string>
#include "CGenebase.h"
#include "GaDeclaration.h"


namespace GAZJUT{

class CGenome
{
	public:
		//constructor
		CGenome();
		CGenome(std::vector <GENEVAR>*, E_CODE_TYPE );
		CGenome(CGenome&);
		~CGenome();
		//operations function
		CGenome* clone();
		CGenebase* extractGene(int);
		void   insertGeneToGenome(CGenebase*);
		void   init(std::vector <GENEVAR>*);
		void   init();
		//input output of varible
		double fitness();
		void   setFitness(const double);

		double rawscore();
		void   setRawScore(const double);

		double origValue();
		void   setOrigValue(const double);

        void   setCumufitness(const double);
		double cumufitness();

		double relativefitness();
		void   setRelativefitness(const double);

		std::vector <GENEVAR>*     geneVariable();       //return varible
		std::vector <double>       *getDecValue();

		void updateDecValueGene(std::vector <double>*);              // after Genetor operatoration, it will be updated.
        void   updateTotalGeneToIndividualGene();

		std::vector <unsigned int> *totalbitGene();
		void setTotalbitGene(std::vector <unsigned int>*);

		std::vector <double>       *totalrealGene();
		void setTotalrealGene(std::vector <double>*);

		int                    geneNum();
		int                    totalbitNum();

		bool   isNormalFinish()const;
		void   setFinishState(bool);
		//operator~
		bool   operator == (CGenome&);         // equal label
		bool   operator ^= (CGenome&);         // Approximately equal
		double operator [] (std::string);      // obtained value;
		// varible
	protected:
		std::vector <CGenebase*>   *m_pGenome;
		std::vector <unsigned int> *m_ptotalgeneofGenome;   // for bit  gray gene;
		std::vector <double>       *m_ptotalRealofGenome;   // for real
		std::vector <GENEVAR>      *m_pgeneVarofGenome;
		          E_CODE_TYPE       m_codeType;
                          int       m_totalbitNum;
		                  int       m_geneNum;
		                  int       m_varNumofGenome;
		// three value for evaluation of evaluator
		                double      m_origValue;      // it is obtained from energy of DFT
		                                         // Importantly, only this value can be compared with those from parents, grandparents
		                double      m_rawscore;       // it is obtained from conversion of m_origValue by using some formula.  Usually
		                                         // this conversion formula differed among the generations.
		                double      m_fitness;        // m_fitness= the value after scalling treatment of m_rawscore
		                                         // statistic treatment of fitness, and it is used for selector, crosser, mutator
		                double      m_relativefitness;
                        double      m_cumufitness;
		                  bool      FinishState;

};

}
#endif
