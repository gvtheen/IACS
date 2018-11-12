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

        void initWorkEnvironment();
        void setCurrentWorkPathAt(size_t,size_t);
        std::string currentWorkPath();
        std::string currentInitPath();
        void checkExeNecessaryFiles(std::vector<std::string>& files);

        void moveFileToPath(const std::string& file, const std::string& dir);
        void copyFileToPath(const std::string& file, const std::string& dir);
    protected:

    private:
        CParameter* m_pParameter;
        boost::filesystem::path m_root_WorkingPath;
};


}
#endif // CIOMONITOR_H
