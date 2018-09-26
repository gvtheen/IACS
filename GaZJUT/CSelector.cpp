#include "CSelector.h"
namespace GAZJUT{

//constructor
CSelector::CSelector()
{
}
CSelector::CSelector(const CSelector& myselector)
{
    return new CSelector();
}
//deconstructor
CSelector::~CSelector()
{
}
//property function
void CSelector::clone()
{
	return new CSelector(*this);
}

void CSelector::run(CGpopulation* PCurrentPopulation)
{
	assert(PCurrentPopulation);

    if(PCurrentPopulation==nullptr)
    {
          ERROR_OUTPUT("Error: CurrentPopulation is null","CSelector","run");
          abort();
    }
    // run statistic calculation
	PCurrentPopulation->fitness_statistic();
	//dispatch selector
	switch (PCurrentPopulation->m_pObjGaparameter->SearchT())
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
			break;
		default:
			this->roulette_Wheel_Select(PCurrentPopulation);
		    break;
	}
}

void CSelector::roulette_Wheel_Select(CGpopulation* m_PCurrentPopulation)
{
	// copy original population to temp population for further treatment.
	// After selector process, current population (m_PCurrentPopulation) modified is new population

	CGpopulation* copy_of_currentPopulation = m_PCurrentPopulation->clone();
	std::vector <CGenome*> p_copy_currentPopulation=copy_of_currentPopulation->m_pGpopulation;

	std::vector <CGenome*>::iterator it;
	double rndNum;
    CRandomgenerator *Rndgenerator=new CRandomgenerator();
    unsigned int seedNum = 0;

	for (it=m_PCurrentPopulation->m_pGpopulation->begin();it<m_PCurrentPopulation->m_pGpopulation->end();it++)
	{
		rndNum=Rndgenerator->uniformRandom01(seedNum++);
		if( rndNum < p_copy_currentPopulation[0]->cumufitness()
		    (*it)->updateDecValueGene(p_copy_currentPopulation[0]->getDecValue());
		else
		{
			for(int i=0; i < p_copy_currentPopulation.size() - 1; i++)
            {
                    if(rndNum > p_copy_currentPopulation[i]->cumufitness() &&  \
                       rndNum < p_copy_currentPopulation[i+1]->cumufitness())
                    {
                        (*it)->updateDecValueGene(p_copy_currentPopulation[i]->getDecValue());
				    }
			}
		}
	}
	delete Rndgenerator;
	delete copy_of_currentPopulation;
}
void CSelector::rondom_Select(CGpopulation* m_PCurrentPopulation)
{
	// copy original population to temp population for further treatment.
	// After selector process, current population (m_PCurrentPopulation) modified is new population
	CGpopulation* copy_of_currentPopulation = m_PCurrentPopulation->clone();
	std::vector <CGenome*> p_copy_currentPopulation=copy_of_currentPopulation->m_pGpopulation;

    std::vector <CGenome*>::iterator it;
	double pointerdistance;
	int pointer;
    CRandomgenerator *Rndgenerator=new CRandomgenerator();

	pointerdistance=1.0/m_PCurrentPopulation->popNum();
	pointer=pointerdistance*Rndgenerator->uniformRandom01(100);

	for (it=m_PCurrentPopulation->m_pGpopulation->begin();it<m_PCurrentPopulation->m_pGpopulation->end();it++)
	{
		pointer=pointer+i*pointerDistance;
		if(pointer < p_copy_currentPopulation[0]->cumufitness())
		   (*it)->updateDecValueGene(p_copy_currentPopulation[0]->getDecValue());
		else{
			for(int j=0;j<p_copy_currentPopulation.size() - 1;j++)
			{
				if(pointer > p_copy_currentPopulation[j]->cumufitness() && \
				   pointer < p_copy_currentPopulation[j+1]->cumufitness() )
				   (*it)->updateDecValueGene(p_copy_currentPopulation[i]->getDecValue());
			}
		}

	}
	delete copy_of_currentPopulation;
	delete Rndgenerator;
	p_copy_currentPopulation.clear();
}

void CSelector::tournament_Select(CGpopulation* m_PCurrentPopulation)
{
	CGpopulation* copy_of_currentPopulation = m_PCurrentPopulation->clone();
	std::vector <CGenome*> p_copy_currentPopulation = copy_of_currentPopulation->m_pGpopulation;
    std::vector <CGenome*>::iterator it;
    std::vector <CGenome*> p_currpopu = m_PCurrentPopulation->m_pGpopulation;
	int p_n1,p_n2;
    CRandomgenerator *Rndgenerator=new CRandomgenerator();
    unsigned int seedNum=100;
	for(it=p_currpopu.begin();it<p_currpopu.end();it++)
	{
		p_n1=(int)(m_PCurrentPopulation->popNum()*(Rndgenerator->uniformRandom01(seedNum++));
		p_n2=(int)(m_PCurrentPopulation->popNum()*(Rndgenerator->uniformRandom01(seedNum++));

		// comparing value is fitness, but not cumufitness;
		if(p_copy_currentPopulation[p_n1]->fitness() > p_copy_currentPopulation[p_n2]->fitness())
		   (*it)->updateDecValueGene(p_copy_currentPopulation[p_n1]->getDecValue());
		else
		   (*it)->updateDecValueGene(p_copy_currentPopulation[p_n2]->getDecValue());
	}
	delete copy_of_currentPopulation;
	delete Rndgenerator;
	p_copy_currentPopulation.clear();
}



}///namespace
