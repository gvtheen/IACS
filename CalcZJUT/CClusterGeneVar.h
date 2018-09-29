#ifndef CCLUSTERGENEVAR_H
#define CCLUSTERGENEVAR_H

#include<vector>
#include "../Util/Point-Vector.h"

using util::Point3i;

namespace CALCZJUT{

class CClusterGeneVar
{
    typedef enum EEUUR0001{
                 DISTANCE,
                 ANGLE,
                 TORSIN
              }VAR_TYPE;

    public:
        CClusterGeneVar();
        CClusterGeneVar(VAR_TYPE, int, Point3i);

        virtual ~CClusterGeneVar();

        VAR_TYPE type();
        void setType(const VAR_TYPE&);

        std::vector<size_t> index();
        void setIndex(std::vector<size_t>&);
        void setIndex(int i,int j,int k,int h);
        void setIndex(size_t,Point3i);

        bool operator == (CClusterGeneVar&);
        void operator = (CClusterGeneVar&);

    protected:

    private:
                   VAR_TYPE   m_type;
        std::vector<size_t>   m_index;



};



}
#endif // CCLUSTERGENEVAR_H
