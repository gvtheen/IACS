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
#include <boost/algorithm/string.hpp>
#include "CExeFitnessInterface.h"
#include "../CataZJUT/CPeriodicFramework.h"
#include "CParameter.h"
#include "../CataZJUT/CFragment.h"
#include "CModelCluster.h"
#include "CModelClusterSupport.h"
#include "CModel2DSupport.h"
#include "CIOBase.h"
#include "CIOCar.h"
#include "CIOMol.h"
#include "CIOPoscar.h"
#include "CIOGjf.h"
#include "CIOCif.h"
#include "CIOCellFile.h"
#include "../Util/Point-Vector.h"
#include "../GaZJUT/CGaparameter.h"
#include "../Util/log.hpp"
#include "../CataZJUT/CBondTolerance.h"

using util::Log;
using util::Bitset;

namespace CALCZJUT{

CExeFitnessInterface::CExeFitnessInterface(CParameter* mpara)
:m_Parameter(mpara)
{

}

CExeFitnessInterface::~CExeFitnessInterface()
{
}

void CExeFitnessInterface::init()
{

}

double CExeFitnessInterface::CalcuRawFit(std::vector<double>& RealValueOfGenome,size_t& pop_index, bool& isNormalexist)
{
    return 0;
}

void CExeFitnessInterface::ConvOrigToRawScore(std::vector<double>& othr)
{

}
void CExeFitnessInterface::setCalcModeStruct(CModelBase* Temp_calcModeStruct)
{
    this->m_pCalcModeStruct = Temp_calcModeStruct;
}
CModelBase* CExeFitnessInterface::calcModeStruct()
{
    return this->m_pCalcModeStruct;
}
void CExeFitnessInterface::setIO(CIOBase* m_IO)
{
    this->m_pIO=m_IO;
}
CIOBase* CExeFitnessInterface::IO()const
{
    return this->m_pIO;
}
void CExeFitnessInterface::getIO(std::string &file_name,CATAZJUT::CPeriodicFramework* currentPeriodicFramework,CIOBase* resultIO)
{
    std::vector<std::string> vectstr;
    boost::algorithm::split(vectstr,file_name,boost::algorithm::is_any_of("."),boost::algorithm::token_compress_on);
    boost::algorithm::trim(vectstr[1]);
    if(vectstr[1]=="mol")
        resultIO = new CIOMol(currentPeriodicFramework);
    else if(vectstr[1]=="car")
        resultIO = new CIOCar(currentPeriodicFramework);
    else if(vectstr[1]=="poscar")
        resultIO = new CIOPoscar(currentPeriodicFramework);
    else if(vectstr[1]=="gjf")
        resultIO = new CIOGjf(currentPeriodicFramework);
    else if(vectstr[1]=="cif")
        resultIO = new CIOCif(currentPeriodicFramework);
    else if(vectstr[1]=="cell")
        resultIO = new CIOCellFile(currentPeriodicFramework);
    else
        resultIO = nullptr;
}

}
