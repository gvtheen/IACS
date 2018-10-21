#ifndef CIOMONITOR_H
#define CIOMONITOR_H

namespace CALCZJUT{

class CParameter;

class CIOMonitor
{
    public:
        CIOMonitor(CParameter*);
        virtual ~CIOMonitor();

        void createWorkingPath();
        std::string& absoluteFilePath(size_t,size_t);
        std::string& currentPath();

        void initWorkingEnvironment(size_t,size_t);
    protected:

    private:
        CParameter* m_pParameter;
};


}
#endif // CIOMONITOR_H
