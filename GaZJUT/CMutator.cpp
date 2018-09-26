#include<cmath>
#include<iostream>
#include "CMutator.h"
namespace GAZJUT{

CMutator::CMutator()
{
}

CMutator::~CMutator()
{
}
void CMutator::run(const CGpopulation* CurrentPopulation)
{
    CRandomgenerator *Rndgenerator=new CRandomgenerator();
    E_CODE_TYPE codeMode         = CurrentPopulation->m_pObjGaparameter->CodeMode();
    E_MUTATE_OPERATOR mutateMode = CurrentPopulation->m_pObjGaparameter->MutateMode();
                       double Pm = CurrentPopulation->m_pObjGaparameter->MutaProb();
    for(int i=0;i<CurrentPopulation->popNum();i++)
    {
       if(codeMode!= REAL ){   //  binary gray code have only 0 and 1 value;
          int bitnum = CurrentPopulation->m_pGpopulation->at(i)->totalbitNum();
          std::vector <unsigned int> *currentGeneofGenome = CurrentPopulation->m_pGpopulation->at(i)->totalbitGene();
          for(int j=0;j<bitnum;j++)
             if ( Rndgenerator->uniformRandom01(i*j) > Pm )
                 switch(mutateMode)
                 {
                    case UNIFORM_M:
                         if (currentGeneofGenome->at(j)== 0)
                             currentGeneofGenome->at(j)= 1;
                         else
                             currentGeneofGenome->at(j)= 0;
                         break;
                    case BOUNDARY:
                         if(Rndgenerator->uniformRandom01(i*(j+2)) > 0.5)
                            currentGeneofGenome->at(j)= 1;
                         else
                            currentGeneofGenome->at(j)= 0;
                         break;
                    default:
                        break;
                 }
       }else{                      //real gene
           int geneNum = CurrentPopulation->m_pGpopulation->at(i)->geneNum();
           std::vector <double> *currentGeneofGenome = CurrentPopulation->m_pGpopulation->at(i)->totalrealGene();
           std::vector <GENEVAR> *varofGenome= CurrentPopulation->m_pGpopulation->at(i)->geneVariable();
           for(int j=0;j<geneNum;j++)
             if ( Rndgenerator->uniformRandom01(i*j) > Pm )
               switch(mutateMode)
              {
                case NOUNIFORM:
                    minV=varofGenome->at(j)->min;
                    maxV=varofGenome->at(j)->max;
                    double b=2.0;
                    double expon=b*(1-(double)(CurrentPopulation->m_pObjGaparameter->Curr_Generation) \
                                    / CurrentPopulation->m_pObjGaparameter->GenerationNum());
                    double currentValue=currentGeneofGenome->at(j);
                    if(Rndgenerator->uniformRandom01(i*j+j+10) < 0.5)
                       currentGeneofGenome->at(j)= currentValue + (maxV-currentValue)* \
                                                  (1-std::pow(Rndgenerator->uniformRandom01(i*j+10*j),expon));
                    else
                       currentGeneofGenome->at(j)= currentValue - (currentValue - minV )* \
                                                  (1-std::pow(Rndgenerator->uniformRandom01(i*j+10*j),expon));
                    break;
                case GAUSSIAN_M:
                    minV=varofGenome->at(j)->min;
                    maxV=varofGenome->at(j)->max;
                    double *rndNumber = new double(12);
                    double sumRnd=0.0;
                    for(int k=0;k<12;k++)
                        sumRnd = sumRnd + uniformRandom01(i*j+k+10);
                    currentGeneofGenome->at(j) =(minV+maxV)/2.0+(maxV-minV)*(sumRnd - 6.0)/6.0;
              }
       }
       CurrentPopulation->m_pGpopulation->at(i)->updateTotalGeneToIndividualGene();
    }

}


}
