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
#include "CGenebase.h"
#include "GaUtilityFunction.h"

namespace GAZJUT{
CGenebase::CGenebase()
{
    this->m_GeneVAR=nullptr;
}
CGenebase::CGenebase(GeneVAR* var)
{
    this->init(var);
}
double CGenebase::decode()
{
    return this->m_value;
}
void CGenebase::init(GeneVAR* var)
{
    m_GeneVAR=var;
}
void CGenebase::updatecode(double m)
{

}
Bitset& CGenebase::bitGene()
{
    Bitset a1;        // nothing
    return a1;
}
double CGenebase::realGene()
{
    return this->m_value;
}
double CGenebase::value()
{
    return this->m_value;
}
size_t CGenebase::bitNum()
{
    return 0;
}
CGenebase::~CGenebase()
{
    //delete this->m_GeneVAR;
}


}
