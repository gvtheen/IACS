#ifndef CMOLECULARSYMMETRY_H
#define CMOLECULARSYMMETRY_H


namespace CATAZJUT{

class CConfigurationBase;

class CMolecularSymmetry
{
    public:
        CMolecularSymmetry(CConfigurationBase*);
        virtual ~CMolecularSymmetry();

        char* GetPointGroup();

    private:
        CConfigurationBase    *m_pConfiguration;
        char* point_group;
};


}
#endif // CMOLECULARSYMMETRY_H
