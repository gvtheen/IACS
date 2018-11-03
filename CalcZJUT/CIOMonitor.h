#ifndef CIOMONITOR_H
#define CIOMONITOR_H
#include <string>
#include <vector>
namespace CALCZJUT{

class CParameter;

class CIOMonitor
{
    public:
        CIOMonitor(CParameter*);
        virtual ~CIOMonitor();

        void createWorkingPath();
        void init();
        void absoluteFilePath(size_t,size_t,std::string&);
        std::string& currentPath();
        bool checkExeNecessaryFiles(std::vector<std::string>& files, size_t checkmode);
        void initWorkingEnvironment(size_t,size_t);

    protected:

    private:
        CParameter* m_pParameter;
        std::string m_root_WorkingPath;
};


}
#endif // CIOMONITOR_H
