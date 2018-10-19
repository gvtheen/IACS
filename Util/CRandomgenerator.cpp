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
#include <boost/random.hpp>
#include<time.h>
#include <limits>
#include "CRandomgenerator.h"

namespace util{

unsigned int CRandomgenerator::classcalled=0;

CRandomgenerator::CRandomgenerator()
{
    classcalled++;
    //cout<<"new Ranm: " <<classcalled<<endl;
    m_seedInt=classcalled;
    if (classcalled==std::numeric_limits<unsigned int>::max())
        classcalled=0;
}
CRandomgenerator::~CRandomgenerator()
{
    //dtor
}
double CRandomgenerator::uniformRandom01(int seedNum)
{
    boost::mt19937 rng(time(0)+seedNum*seedNum+m_seedInt);
    boost::uniform_01<boost::mt19937&> u01(rng);
    boost::uniform_int<> ui(0,MNNN);
    double *myRandom=new double[MNNN];
    for(int i=0;i<MNNN;i++)
        myRandom[i]=u01();
    double res= myRandom[ui(rng)];
    delete myRandom;
    return res;
}
int CRandomgenerator::uniformRandomRandge(int start, int stop,int seedNum)
{
     boost::mt19937 rng(time(0)+seedNum*seedNum*10+m_seedInt*m_seedInt*10);
     boost::uniform_int<> ui(start,stop);
     return ui(rng);
}
Vector3 CRandomgenerator::randomVector01(int seedNum)
{
     Vector3 res;
     res<<uniformRandom01(seedNum*5+1),uniformRandom01(seedNum*8+2),uniformRandom01(seedNum*16+3);
     return res;
}

}
