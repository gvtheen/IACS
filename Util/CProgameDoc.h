#ifndef CPROGAMEDOC_H
#define CPROGAMEDOC_H
#include <string>
namespace util{

class CProgameDoc
{
    public:
        CProgameDoc();
        virtual ~CProgameDoc();

        void output();

    protected:

    private:
        std::string versionText;
        std::string authorText;
        std::string addressText;
};

}
#endif // CPROGAMEDOC_H
