#ifndef CADAPTERBASE_H
#define CADAPTERBASE_H

/*
1. transfer the initial variable range of computation to engine
2. transfer the setting value from engine to computation
3. convert the search mode(max,specific value) to min mode
*/
namespace IACSZJUT{

class CAdapterBase
{
    public:
        CAdapterBase();
        virtual ~CAdapterBase();

        virtual void engineToComputation(std::vector<double>& _input,size_t _index; bool& _state, double &_output)=0;
        virtual void computationToEngine()=0;
    protected:

    private:
};

}
#endif // CADAPTERBASE_H
