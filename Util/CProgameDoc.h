#ifndef CPROGAMEDOC_H
#define CPROGAMEDOC_H

namespace util{

class CProgameDoc
{
    public:
        CProgameDoc();
        virtual ~CProgameDoc();

    protected:

    private:
        std::string versionText;
        std::string authorText;
        std::string logoText;
};

}
#endif // CPROGAMEDOC_H
