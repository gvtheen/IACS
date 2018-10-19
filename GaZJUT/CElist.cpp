/******************************************************************************
**
** Copyright (C) 2019-2031 Dr.Gui-lin Zhuang <glzhuang@zjut.edu.cn>
** All rights reserved.
**
**
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
** OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**
******************************************************************************/
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
