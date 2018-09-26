#ifndef CCONFIGURATIONPRIVATEDATA_H
#define CCONFIGURATIONPRIVATEDATA_H
#include<vector>
#include<string>

namespace CATAZJUT{

class CBond;
class CAtom;
class CConfigurationWatcher;
class CCoordinateSet;
class CFragment;

class CConfigurationPrivateData
{
public:
    CConfigurationPrivateData();
    ~CConfigurationPrivateData();

    std::string name;
    std::vector<CBond *> bonds;
   // std::vector<CConfigurationWatcher *> watchers;
    bool fragmentsPerceived;
    std::vector<CFragment *> fragments;
    std::vector<std::string> atomTypes;
    std::vector<std::pair<CAtom*, CAtom*> > bondAtoms;
    std::vector<std::vector<CBond *> > atomBonds;
    std::vector<double> partialCharges;
    std::vector<boost::shared_ptr<CCoordinateSet> > coordinateSets;
};






}

#endif
