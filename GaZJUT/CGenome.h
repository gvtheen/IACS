#ifndef CGENOME_H
#define CGENOME_H

#include <vector>
#include <string>
#include "CGenebase.h"
#include "GaDeclaration.h"
#include "../Util/Bitset.h"

using util::Bitset;

namespace GAZJUT{

class CGenome
{
	public:
		//constructor
		CGenome();
		CGenome(std::vector <GeneVAR>*, E_CODE_TYPE,size_t);
		CGenome(CGenome&);
		~CGenome();
		//operations function
		CGenome* clone();
		CGenebase* extractGene(int);
		void   insertGeneToGenome(CGenebase*);
		void   init(std::vector <GeneVAR>&);
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

		std::vector <GeneVAR>&     GeneVARiable();       //return varible
        void getDecValue(std::vector <double>&);

		void updateDecValueGene(std::vector <double>&);              // after Genetor operatoration, it will be updated.
        void updateTotalGeneToIndividualGene();

		Bitset& totalbitGene();
		void setTotalbitGene(const Bitset);

		std::vector <double>& totalrealGene();
		void setTotalrealGene(std::vector <double>&);

		size_t geneNum();
		size_t totalbitNum();

		size_t index();

		bool   isNormalFinish()const;
		void   setFinishState(bool);
		//operator~
		bool   operator == (CGenome&);         // equal label
		bool   operator ^= (CGenome&);         // Approximately equal
		double operator [] (std::string);      // obtained value;
		// varible
    private:
        void concatenateBitsets(Bitset& first, const Bitset& second);
	protected:
	                      size_t    m_index;
		std::vector <CGenebase*>    m_Genome;
		Bitset                      m_totalgeneofGenome;   // for bit  gray gene;
		std::vector <double>        m_totalRealofGenome;   // for real
		std::vector <GeneVAR>       m_GeneVARofGenome;
		          E_CODE_TYPE       m_codeType;
                       size_t       m_totalbitNum;
                       size_t       m_geneNum;
                       size_t       m_varNumofGenome;
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
