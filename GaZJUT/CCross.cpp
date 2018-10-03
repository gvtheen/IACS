#include <algorithm>
#include<vector>
#include<iostream>
#include<math.h>
#include "CCross.h"
#include "CRandomgenerator.h"
#include "GaDeclaration.h"
namespace GAZJUT{

CCross::CCross()
{
}

CCross::~CCross()
{
}

void CCross::run(CGpopulation* CurrentPopulation)
{
   assert(CurrentPopulation);

   int first=0,one;
   std::vector <CGenome*>::iterator it;
   CRandomgenerator *Rndgenerator=new CRandomgenerator();

   for(int i=0;i<CurrentPopulation->popNum();i++)
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
    int pointer,totalbitNum;
    CRandomgenerator *Rndgenerator  = new CRandomgenerator();
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
                std::vector <unsigned int>*temp_one = ((*CurrentPopulation)[one])->totalbitGene();
                std::vector <unsigned int>*temp_two = ((*CurrentPopulation)[two])->totalbitGene();
                unsigned int temp;
                for(int i=0;i<pointer;i++)
                {
                    temp = temp_one->at(i);
                    temp_one->at(i) = temp_two->at(i);
                    temp_two->at(i) = temp;
                }
              }
              break;
             case  MULTIPLE:
              {
                std::vector<int> *pointer = new (std::vector<int>);
                // using operator [index] return the pointer of index th Genome
                totalbitNum=((*CurrentPopulation)[one])->totalbitNum();
                int crossNum = CurrentPopulation->m_pObjGaparameter->CrossNum();
                for(int i=0;i<crossNum;i++)
                    pointer->push_back((int)((Rndgenerator->uniformRandom01(one*two*i*i+1))*(totalbitNum)));
                // sort the crossing points.
                std::sort(pointer->begin(),pointer->end(),[](int a,int b){return a<b;});

                unsigned int start,stop,temp;
                for(int i=0;i<crossNum;i++)
                {
                    if((i+1)%2==1)
                    {
                        start = pointer->at(i);
                        if((i+1)<crossNum)
                            stop = pointer->at(i+1);
                        else
                            stop = totalbitNum;
                        std::vector <unsigned int>*temp_one = ((*CurrentPopulation)[one])->totalbitGene();
                        std::vector <unsigned int>*temp_two = ((*CurrentPopulation)[two])->totalbitGene();
                        for(unsigned int j=start; j<stop; j++)
                        {
                            temp = temp_one->at(i);
                            temp_one->at(i) = temp_two->at(i);
                            temp_two->at(i) = temp;
                        }
                    }
                 }
               }
                break;
             case  UNIFORM_C:
             {
                 std::vector <unsigned int>*temp_one = ((*CurrentPopulation)[one])->totalbitGene();
                 std::vector <unsigned int>*temp_two = ((*CurrentPopulation)[two])->totalbitGene();
                 totalbitNum = ((*CurrentPopulation)[one])->totalbitNum();
                 unsigned int temp;
                 for(int i=0;i<totalbitNum;i++)
                     if(Rndgenerator->uniformRandom01(10*i) <= 0.5 )
                     {
                        temp = temp_one->at(i);
                        temp_one->at(i) = temp_two->at(i);
                        temp_two->at(i) = temp;
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
                 std::vector <double>*temp_one = ((*CurrentPopulation)[one])->totalrealGene();
                 std::vector <double>*temp_two = ((*CurrentPopulation)[two])->totalrealGene();
                 int geneNum = ((*CurrentPopulation)[one])->geneNum();
                 double tmp_Vone,tmp_Vtwo, a;
                 for(int i=0;i<geneNum;i++)
                 {
                   //   maxV = (varofGenome->at(i)).max;
                   //   minV = (varofGenome->at(i)).min;
                      tmp_Vone = temp_one->at(i);
                      tmp_Vtwo = temp_two->at(i);
                      a = CurrentPopulation->m_pObjGaparameter->UnifArithmCrossConstant;
                      /*
                        a = uniform random value of [-d,1+d], d=0.25;
                      */
                      a = -1*a + (1+2*a)* Rndgenerator->uniformRandom01(one*two + i);
                      temp_one->at(i) = tmp_Vone + a*( tmp_Vtwo - tmp_Vone );
                      temp_two->at(i) = tmp_Vtwo + a*( tmp_Vone - tmp_Vtwo );
                 }
             }
             break;
            case UNARITHMETIC:
             {
                // std::vector <GeneVAR> *varofGenome= ((*CurrentPopulation)[one])->GeneVARiable();
                 std::vector <double>*temp_one = ((*CurrentPopulation)[one])->totalrealGene();
                 std::vector <double>*temp_two = ((*CurrentPopulation)[two])->totalrealGene();
                 int geneNum = ((*CurrentPopulation)[one])->geneNum();
                 int currentGenerationNum = CurrentPopulation->m_pObjGaparameter->Curr_Generation;
                 int totalGenerationNum   = CurrentPopulation->m_pObjGaparameter->GenerationNum();
                 double tmp_Vone,tmp_Vtwo, a;
                 for(int i=0;i<geneNum;i++)
                 {
//                      maxV = (varofGenome->at(i)).max;
//                      minV = (varofGenome->at(i)).min;
                      tmp_Vone = temp_one->at(i);
                      tmp_Vtwo = temp_two->at(i);
                      a = CurrentPopulation->m_pObjGaparameter->NoUnifArithmCrossConstant;
                      /*
                        a = uniform random value of [-d,1+d], d=0.25;
                      */
                      a = std::exp(-a*totalGenerationNum /currentGenerationNum);
                      temp_one->at(i) = tmp_Vone + a*( tmp_Vtwo - tmp_Vone );
                      temp_two->at(i) = tmp_Vtwo + a*( tmp_Vone - tmp_Vtwo );
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

