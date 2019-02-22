#ifndef CBONDTOLERANCE_H
#define CBONDTOLERANCE_H
#include<string>
#include<vector>
#include<map>
#include "CBond.h"

namespace CATAZJUT{
//declare
class CAtom;
class CElement;
class CConfigurationBase;
class CBondPrivate;

class CBondTolerance
{
    public:
        CBondTolerance(CConfigurationBase*);
        CBondTolerance(CConfigurationBase*,std::vector<std::pair<std::string*,std::string*>>);
        CBondTolerance(CConfigurationBase*,double,double);

        CBondTolerance(CBondTolerance&other);
        virtual ~CBondTolerance();
        //element of bond
        void addElement(const CElement&);
        void removeElement(CElement&);
        void setExcludeBond(std::vector<std::pair<std::string*,std::string*>>&);
        //bond property
        int IsExistBond(const std::string&,const std::string&);

        bool IsBond(const std::string&,const std::string&,const double&);
        bool IsBond(const size_t&,const size_t&,const double&);
        bool IsBond(const CAtom*,const CAtom*);

        void AddBondType(const std::string& _1stAtom,const std::string& _2ndAtom,const double& mind,const double& maxd);
        void AddBondType(const std::string&,const std::string&);
        void AddBondType(const CElement& Elem_1,const CElement& Elem_2,const double& mind,const double& maxd);
        void AddBondType(const CElement& Elem_1,const CElement& Elem_2);
        void RemoveBondType(const std::string&,const std::string&);

        void setTolerancefactor(const double,const double);
        void setTolerancefactor(std::pair<double,double> &mht);
        void OutputTofile(std::string& outfile);

    protected:
        bool isExcludeBond(const std::string&,const std::string&);
        bool isExcludeBond(CElement& Elem_1,CElement& Elem_2);
    private:
        CConfigurationBase*        m_pConfiguration;
        std::vector<std::pair<std::string*,std::string*>>  m_pExcludeBond;
        double lower_tolerance_factor; //BondCriteria
        double upper_tolerance_factor;

        std::vector<CBondPrivate*> m_pBondType;
        std::vector<CElement>      m_pElement;
};


}
#endif // CBONDTOLERANCE_H
