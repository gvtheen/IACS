#ifndef CBONDPRIVATE_H
#define CBONDPRIVATE_H
namespace CATAZJUT{

class CBondPrivate
{
     public:
        CBondPrivate();
        virtual ~CBondPrivate();
     public:
       std::string m_1stAtom;
       std::string m_2ndAtom;
       double m_MinBondLength;
       double m_MaxBondLength;
};


}
#endif
