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
#include <algorithm>
#include<vector>
#include<iostream>
#include<math.h>
#include "CCross.h"
#include "../Util/CRandomgenerator.h"
#include "GaDeclaration.h"
#include "../Util/Bitset.h"
#include "../Util/log.hpp"

using util::Bitset;
using util::Log;

namespace GAZJUT{

CCross::CCross()
:CGaOperatorBase()
{
}

CCross::~CCross()
{
}

void CCross::run(CGpopulation* CurrentPopulation)
{
   assert(CurrentPopulation);

   size_t first=0,one;
   std::vector <CGenome*>::iterator it;
   util::CRandomgenerator *Rndgenerator=new util::CRandomgenerator();

   for(size_t i=0;i<CurrentPopulation->popNum();i++)
        if(Rndgenerator->uniformRandom01(i) < this->m_Crossprob)
        {
            first = first + 1;
            if(first%2==0)
                X_Crossover(CurrentPopulation,one,i);
            else
                one=i;
        }
   delete Rndgenerator;
}
void CCross::X_Crossover(CGpopulation* CurrentPopulation,int one,int two)
{
    size_t pointer,totalbitNum;
    util::CRandomgenerator *Rndgenerator  = new util::CRandomgenerator();
    E_CODE_TYPE codeMode            = CurrentPopulation->m_pObjGaparameter->CodeMode();
    E_CROSSOVER_OPERATOR cross_Mode = CurrentPopulation->m_pObjGaparameter->CrossMode();
    if(codeMode != REAL )
    {
        switch(cross_Mode)
        {
             case  SINGLE:
              {  //binary code gray code
                totalbitNum=((*CurrentPopulation)[one])->totalbitNum();
                pointer = (int)((Rndgenerator->uniformRandom01(one*two))*(totalbitNum));
                Bitset& temp_one = ((*CurrentPopulation)[one])->totalbitGene();
                Bitset& temp_two = ((*CurrentPopulation)[two])->totalbitGene();
                size_t temp;
                for(size_t i=0;i<pointer;i++)
                {
                    temp = temp_one[i];
                    temp_one[i] = temp_two[i];
                    temp_two[i] = temp;
                }
              }
              break;
             case  MULTIPLE:
              {
                std::vector<int> *pointer = new (std::vector<int>);
                // using operator [index] return the pointer of index th Genome
                totalbitNum=((*CurrentPopulation)[one])->totalbitNum();
                size_t crossNum = CurrentPopulation->m_pObjGaparameter->CrossNum();
                for(size_t i=0;i<crossNum;i++)
                    pointer->push_back((int)((Rndgenerator->uniformRandom01(one*two*i*i+1))*(totalbitNum)));
                // sort the crossing points.
                std::sort(pointer->begin(),pointer->end(),[](int a,int b){return a<b;});

                size_t start,stop,temp;
                for(size_t i=0;i<crossNum;i++)
                {
                    if((i+1)%2==1)
                    {
                        start = pointer->at(i);
                        if((i+1)<crossNum)
                            stop = pointer->at(i+1);
                        else
                            stop = totalbitNum;
                        Bitset& temp_one = ((*CurrentPopulation)[one])->totalbitGene();
                        Bitset& temp_two = ((*CurrentPopulation)[two])->totalbitGene();
                        for(size_t j=start; j<stop; j++)
                        {
                            temp = temp_one[i];
                            temp_one[i] = temp_two[i];
                            temp_two[i] = temp;
                        }
                    }
                 }
               }
                break;
             case  UNIFORM_C:
             {
                 Bitset& temp_one = ((*CurrentPopulation)[one])->totalbitGene();
                 Bitset& temp_two = ((*CurrentPopulation)[two])->totalbitGene();
                 totalbitNum = ((*CurrentPopulation)[one])->totalbitNum();
                 size_t temp;
                 for(size_t i=0;i<totalbitNum;i++)
                     if(Rndgenerator->uniformRandom01(10*i) <= 0.5 )
                     {
                        temp = temp_one[i];
                        temp_one[i] = temp_two[i];
                        temp_two[i] = temp;
                     }
             }
              break;
             default:
                break;
         }
     }else{
         switch(cross_Mode)
         {
            case ARITHMETIC:
             {
                 //std::vector <GeneVAR> *varofGenome= ((*CurrentPopulation)[one])->GeneVARiable();
                 std::vector <double>& temp_one = ((*CurrentPopulation)[one])->totalrealGene();
                 std::vector <double>& temp_two = ((*CurrentPopulation)[two])->totalrealGene();
                 size_t geneNum = ((*CurrentPopulation)[one])->geneNum();
                 double tmp_Vone,tmp_Vtwo, a;
                 for(size_t i=0;i<geneNum;i++)
                 {
                   //   maxV = (varofGenome->at(i)).max;
                   //   minV = (varofGenome->at(i)).min;
                      tmp_Vone = temp_one[i];
                      tmp_Vtwo = temp_two[i];
                      a = CurrentPopulation->m_pObjGaparameter->UnifArithmCrossConstant;
                      /*
                        a = uniform random value of [-d,1+d], d=0.25;
                      */
                      a = -1*a + (1+2*a)* Rndgenerator->uniformRandom01(one*two + i);
                      temp_one[i] = tmp_Vone + a*( tmp_Vtwo - tmp_Vone );
                      temp_two[i] = tmp_Vtwo + a*( tmp_Vone - tmp_Vtwo );
                 }
             }
             break;
            case UNARITHMETIC:
             {
                // std::vector <GeneVAR> *varofGenome= ((*CurrentPopulation)[one])->GeneVARiable();
                 std::vector <double>& temp_one = ((*CurrentPopulation)[one])->totalrealGene();
                 std::vector <double>& temp_two = ((*CurrentPopulation)[two])->totalrealGene();
                 size_t geneNum = ((*CurrentPopulation)[one])->geneNum();
                 size_t currentGenerationNum = CurrentPopulation->m_pObjGaparameter->Curr_Generation;
                 size_t totalGenerationNum   = CurrentPopulation->m_pObjGaparameter->GenerationNum();
                 double tmp_Vone,tmp_Vtwo, a;
                 for(size_t i=0;i<geneNum;i++)
                 {
//                      maxV = (varofGenome->at(i)).max;
//                      minV = (varofGenome->at(i)).min;
                      tmp_Vone = temp_one[i];
                      tmp_Vtwo = temp_two[i];
                      a = CurrentPopulation->m_pObjGaparameter->NoUnifArithmCrossConstant;
                      /*
                        a = uniform random value of [-d,1+d], d=0.25;
                      */
                      a = std::exp(-a*totalGenerationNum /currentGenerationNum);
                      temp_one[i] = tmp_Vone + a*( tmp_Vtwo - tmp_Vone );
                      temp_two[i] = tmp_Vtwo + a*( tmp_Vone - tmp_Vtwo );
                 }
             }
               break;
            default:
               break;
         }
     }
     delete Rndgenerator;
}

}  //namespace

