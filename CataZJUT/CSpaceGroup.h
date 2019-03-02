#ifndef CSPACEGROUP_H
#define CSPACEGROUP_H

namespace CATAZJUT{

class CConfigurationBase;

class CSpaceGroup
{
    public:
        CSpaceGroup(CConfigurationBase*);
        virtual ~CSpaceGroup();

        char* GetSpaceGroup();
        bool  Find_primitive();
    protected:

    private:
        CConfigurationBase* m_pCConfigurationBase;
        char* m_pSpaceGroupName;
};


}
#endif // CSPACEGROUP_H
