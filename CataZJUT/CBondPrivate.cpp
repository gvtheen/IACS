#include "CBondPrivate.h"
namespace CATAZJUT{

CBondPrivate::CBondPrivate()
: m_1stAtom(std::string("Dummy")),m_2ndAtom(std::string("Dummy")),m_MinBondLength(-1.0),m_MaxBondLength(-1.0)
{
}
CBondPrivate::CBondPrivate(std::string a1,std::string a2,double d1,double d2)
: m_1stAtom(a1),m_2ndAtom(a2),m_MinBondLength(d1),m_MaxBondLength(d2)
{

}
CBondPrivate::~CBondPrivate()
{
}


}
