#include<iostream>
#include "../CalcZJUT/CIOXyz.h"
#include "../CalcZJUT/CParameter.h"
#include "../CataZJUT/CPeriodicFramework.h"
using namespace std;
int main()
{
    CALCZJUT::CParameter* CurrentPara = new CALCZJUT::CParameter("iacs.input");
    CATAZJUT::CPeriodicFramework* mol= new CATAZJUT::CPeriodicFramework(CurrentPara);
    CALCZJUT::CIOXyz* myoutput = new CALCZJUT::CIOXyz(mol);
    myoutput->input("c6h6.xyz");
    myoutput->output("c6h6-1.xyz");
    cout<<"hello"<<endl;
    delete myoutput;
    delete mol;
    delete CurrentPara;
    return 1;
}
