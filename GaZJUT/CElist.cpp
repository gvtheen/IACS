#include "CElist.h"
#include "CGaparameter.h"
namespace GAZJUT{

CElist::CElist()
{
}

void CElist::run(CGpopulation* pCurrentPopulation)
{
    assert(pCurrentPopulation);

    int current_generation = pCurrentPopulation->m_pObjGaparameter->Curr_Generation;

    if(current_generation==0){
        this->m_OldPopulation = new CGpopulation(*pCurrentPopulation);
      this->m_MinOrignalValue = (*pCurrentPopulation)["minOri"];
          this->m_pBestGenome = new CGenome(*(pCurrentPopulation->m_pMinOriGenome));
    }else if( this->m_MinOrignalValue > (*pCurrentPopulation)["minOri"]){
       this->m_MinOrignalValue = (*pCurrentPopulation)["minOri"];
       delete this->m_OldPopulation;
       this->m_OldPopulation = new CGpopulation(*pCurrentPopulation);
    }else if((*m_OldPopulation)["minOri"]<(*pCurrentPopulation)["minOri"]){
       this->m_OldPopulation->descendSort();   //max   min  Fitness
          pCurrentPopulation->asscendSort();   //min   max  Fitness

       int popnum = pCurrentPopulation->popNum();
       for(int i=0;i<popnum;i++)
          if( (*m_OldPopulation)[i]->origValue() < (*pCurrentPopulation)[i]->origValue() )
            (*pCurrentPopulation)[i]->updateDecValueGene((*m_OldPopulation)[i]->getDecValue());
    }
}
CElist::~CElist()
{
    delete m_OldPopulation;
    delete m_pBestGenome;
}

}
