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
#include<iostream>
#include<fstream>
#include<string>
#include<cmath>
#include<ctime>
#include "GaUtilityFunction.h"


namespace GAZJUT{

//
double binaryDecode(const Bitset & myCode,GeneVAR myGenVar)
{
     int i,m;
     double low,high,sum;

     low= myGenVar.min;
     high=myGenVar.max;
     sum=0.0;
     m=myCode.size();
     for(i=0;i<m;i++)
          sum=sum+(double)(myCode[i])*std::pow(2.0,m-i-1);
     return  low+(high-low)*sum/(std::pow(2.0,m)-1);
}
int calcBitNum(GeneVAR myGeneVAR)
{
    double c,m;
    m=(myGeneVAR.max - myGeneVAR.min)/myGeneVAR.accuracy;
    c=std::log10(m)/std::log10(2.0);
    return (int)c+1;
}
void grayTobit(Bitset& data)
{
     int num=data.size();
     Bitset temp= data;
     data[0]=temp[0];
     for(int i=1;i<num;i++)
        data[i]=(data[i-1])^(temp[i]);
     temp.clear();
}
void bitTogray(Bitset& data)
{
     int num=data.size();
     Bitset temp= data;

     data[0]=temp[0];
     for(int i=1;i<num;i++)
        data[i]=(temp[i])^(temp[i-1]);
     temp.clear();
}

}  //namespace

