#include "CSelector.h"
#include "CGaparameter.h"
#include "CGenome.h"
#include "GaDeclaration.h"
#include "../Util/log.hpp"
#include "../GACatalyst.h"
#include "../Util/CRandomgenerator.h"
using util::Log;
using namespace GAZJUT;

namespace GAZJUT{

//constructor
CSelector::CSelector()
:CGaOperatorBase()
{
    this->mix_index=0;
}
CSelector::CSelector(const CSelector& myselector)
{
    //this->CSelector();
}
//deconstructor
CSelector::~CSelector()
{
}
//property function
CSelector* CSelector::clone()
{
	return new CSelector(*this);
}

void CSelector::run(CGpopulation* PCurrentPopulation)
{
	assert(PCurrentPopulation);

    if(PCurrentPopulation==nullptr)
    {
          Log::Error<<"Error: CurrentPopulation is null! CSelector_run!\n";
          boost::throw_exception(std::runtime_error("Error: CurrentPopulation is null! CSelector_run!\n"));
    }
    // run statistic calculation
	PCurrentPopulation->fitness_statistic();
	//dispatch selector
	switch ((int)(PCurrentPopulation->m_pObjGaparameter->SearchType()))
	{
		case ROULETTE_WHEEL:
		    this->roulette_Wheel_Select(PCurrentPopulation);
		    break;
		case RANDOM:
			this->rondom_Select(PCurrentPopulation);
			break;
		case TOURNAMENT:
			this->tournament_Select(PCurrentPopulation);
			break;
		case MIXED:
		    this->mix_index +=1;
		    if(mix_index%3==0)
               this->roulette_Wheel_Select(PCurrentPopulation);
            else if(mix_index%3==1)
               this->rondom_Select(PCurrentPopulation);
            else if(mix_index%3==2){
               this->tournament_Select(PCurrentPopulation);
               this->mix_index=2;
            }
			break;
		default:
			this->roulette_Wheel_Select(PCurrentPopulation);
		    break;
	}
}

void CSelector::roulette_Wheel_Select(CGpopulation* P_CurrentPopulation)
{
	// copy original population to temp population for further treatment.
	// After selector process, current population (m_PCurrentPopulation) modified is new population

	CGpopulation* P_copy_of_CurrentPopulation = P_CurrentPopulation->clone();
	std::vector <CGenome*>& R_copy_currentPopulation = P_copy_of_CurrentPopulation->m_Gpopulation;

	std::vector <CGenome*>::iterator it;
	double rndNum;
    util::CRandomgenerator *Rndgenerator=new util::CRandomgenerator();
    unsigned int seedNum = 0;

    std::vector<double> tempDoubleVect;

	for (it=P_CurrentPopulation->m_Gpopulation.begin();it<P_CurrentPopulation->m_Gpopulation.end();it++)
	{
		rndNum=Rndgenerator->uniformRandom01(seedNum++);
		if( rndNum < R_copy_currentPopulation[0]->cumufitness()){
            // Get Dec value from copy population
            R_copy_currentPopulation[0]->getDecValue(tempDoubleVect);
            // Set Dec value into current Population
           (*it)->updateDecValueGene(tempDoubleVect);
        }else{
			for(size_t i=0; i < R_copy_currentPopulation.size() - 1; i++){
                    if(rndNum > R_copy_currentPopulation[i]->cumufitness() &&  \
                       rndNum < R_copy_currentPopulation[i+1]->cumufitness()){
                        // Get Dec value from copy population
                        R_copy_currentPopulation[i]->getDecValue(tempDoubleVect);
                        // Set Dec value into current Population
                        (*it)->updateDecValueGene(tempDoubleVect);
				    }
			}
		}
	}
	delete Rndgenerator;
	delete P_copy_of_CurrentPopulation;
}
void CSelector::rondom_Select(CGpopulation* P_CurrentPopulation)
{
	// copy original population to temp population for further treatment.
	// After selector process, current population (m_PCurrentPopulation) modified is new population
	CGpopulation* P_copy_of_currentPopulation = P_CurrentPopulation->clone();
	std::vector <CGenome*>& R_p_copy_currentPopulation = P_copy_of_currentPopulation->m_Gpopulation;

	double pointerDistance;
	int pointer;
    util::CRandomgenerator *Rndgenerator=new util::CRandomgenerator();

	pointerDistance=1.0/P_copy_of_currentPopulation->popNum();
	pointer=pointerDistance*Rndgenerator->uniformRandom01(100);

	std::vector<double> tempDoubleVect;
	for (size_t i=0;i<P_CurrentPopulation->popNum();i++)
	{
		pointer=pointer+i*pointerDistance;
		if(pointer < R_p_copy_currentPopulation[0]->cumufitness()){
             R_p_copy_currentPopulation[0]->getDecValue(tempDoubleVect);
            (*P_CurrentPopulation)[i]->updateDecValueGene(tempDoubleVect);
		}else{
			for(size_t j=0;j<R_p_copy_currentPopulation.size() - 1;j++){
				if(pointer > R_p_copy_currentPopulation[j]->cumufitness() && \
				   pointer < R_p_copy_currentPopulation[j+1]->cumufitness() ){
				      R_p_copy_currentPopulation[j]->getDecValue(tempDoubleVect);
				      (*P_CurrentPopulation)[i]->updateDecValueGene(tempDoubleVect);
				   }
			}
		}

	}
	delete P_copy_of_currentPopulation;
	delete Rndgenerator;

}

void CSelector::tournament_Select(CGpopulation* P_CurrentPopulation)
{
	CGpopulation* P_copy_of_currentPopulation = P_CurrentPopulation->clone();
	std::vector <CGenome*>& R_p_copy_currentPopulation = P_copy_of_currentPopulation->m_Gpopulation;

	int p_n1,p_n2;
    util::CRandomgenerator *Rndgenerator=new util::CRandomgenerator();
    unsigned int seedNum=1;
    size_t popNum=P_CurrentPopulation->popNum();

    std::vector<double> tempDoubleVect;
	for(size_t i=0;i<popNum;i++)
	{
		p_n1=(int)(popNum*(Rndgenerator->uniformRandom01(seedNum++)));
		p_n2=(int)(popNum*(Rndgenerator->uniformRandom01(seedNum++)));

		// comparing value is fitness, but not cumufitness;
		if(R_p_copy_currentPopulation[p_n1]->fitness() > R_p_copy_currentPopulation[p_n2]->fitness()){
            R_p_copy_currentPopulation[p_n1]->getDecValue(tempDoubleVect);
            (*P_CurrentPopulation)[i]->updateDecValueGene(tempDoubleVect);

		}else{
		    R_p_copy_currentPopulation[p_n2]->getDecValue(tempDoubleVect);
            (*P_CurrentPopulation)[i]->updateDecValueGene(tempDoubleVect);
		}
	}
	delete P_copy_of_currentPopulation;
	delete Rndgenerator;

}



}///namespace
