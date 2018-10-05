#include "CElist.h"
#include "CGaparameter.h"
#include <vector>
namespace GAZJUT{

CElist::CElist()
:CGaOperatorBase()
{
}

void CElist::run(CGpopulation* pCurrentPopulation)
{
    assert(pCurrentPopulation);

    size_t current_generation = pCurrentPopulation->m_pObjGaparameter->Curr_Generation;

    if(current_generation==0){
        this->m_OldPopulation = new CGpopulation(*pCurrentPopulation);
      this->m_MinOrignalValue = (*pCurrentPopulation)["minOri"];
          this->m_pBestGenome = new CGenome(*(pCurrentPopulation->m_pMinOriGenome));
    }else if( this->m_MinOrignalValue > (*pCurrentPopulation)["minOri"]){
       this->m_MinOrignalValue = (*pCurrentPopulation)["minOri"];
       delete this->m_OldPopulation;
       this->m_OldPopulation = new CGpopulation(*pCurrentPopulation);
    }else if((*m_OldPopulation)["minOri"]<(*pCurrentPopulation)["minOri"]){
                                            //  index= 0 1 2   ... N-1
       this->m_OldPopulation->descendSort();   //      max         min          Fitness
          pCurrentPopulation->ascendSort();    //      min         max          Fitness

       size_t popnum = pCurrentPopulation->popNum();
       std::vector<double> tempValue;
       for(size_t i=0;i<popnum;i++)
          if( (*m_OldPopulation)[i]->origValue() < (*pCurrentPopulation)[i]->origValue() ){
             (*m_OldPopulation)[i]->getDecValue(tempValue);
             (*pCurrentPopulation)[i]->updateDecValueGene(tempValue);
          }
    }
}
CElist::~CElist()
{
    delete m_OldPopulation;
    delete m_pBestGenome;
}

}
