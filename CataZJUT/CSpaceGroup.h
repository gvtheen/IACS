#ifndef CSPACEGROUP_H
#define CSPACEGROUP_H

namespace CATAZJUT{

class CPeriodicFramework;

class CSpaceGroup
{
    public:
        CSpaceGroup(CPeriodicFramework*);
        virtual ~CSpaceGroup();

        char* GetSpaceGroup();
    protected:

    private:
        CPeriodicFramework* m_pCPeriodicFramework;
        char* m_pSpaceGroupName;
};


}
#endif // CSPACEGROUP_H
