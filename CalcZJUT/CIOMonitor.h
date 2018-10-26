#ifndef CIOMONITOR_H
#define CIOMONITOR_H
#include <string>
namespace CALCZJUT{

class CParameter;

class CIOMonitor
{
    public:
        CIOMonitor(CParameter*);
        virtual ~CIOMonitor();

        void createWorkingPath();
        void absoluteFilePath(size_t,size_t,std::string&);
        std::string& currentPath();

        void initWorkingEnvironment(size_t,size_t);

    protected:

    private:
        CParameter* m_pParameter;
        std::string m_currentPath;
};


}
#endif // CIOMONITOR_H
