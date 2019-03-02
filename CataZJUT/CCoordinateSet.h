#ifndef CCOORDINATESET_H
#define CCOORDINATESET_H
#include "../Util/Point-Vector.h"
#include "CatalystUniverseDefine.h"

using util::Point3;
using CATAZJUT::DEFINED::CoordinateType;


namespace CATAZJUT{

class CCartesianCoordinates;
class CInternalCoordinates;
class CFractionCoordinates;

class CCoordinateSet
{
public:
    // construction and destruction
    CCoordinateSet();
    explicit CCoordinateSet(CCartesianCoordinates *coordinates);
    explicit CCoordinateSet(CInternalCoordinates  *coordinates);
    explicit CCoordinateSet(CFractionCoordinates  *coordinates);
    CCoordinateSet(const CCoordinateSet &other);
    ~CCoordinateSet();

    // properties
    CATAZJUT::DEFINED::CoordinateType type() const;
    size_t size() const;
    bool isEmpty() const;
    void setCoordinates(CCartesianCoordinates *coordinates);
    void setCoordinates(CInternalCoordinates  *coordinates);
    void setCoordinates(CFractionCoordinates  *coordinates);
    CCartesianCoordinates* cartesianCoordinates() ;//const;
    CInternalCoordinates*  internalCoordinates() ;//const;
    CFractionCoordinates*  fractionCoordinates() ;//const;
    void clear();

    Point3 position(size_t index) const;

    CCoordinateSet& operator=( CCoordinateSet &other );

private:
    CATAZJUT::DEFINED::CoordinateType m_type;
    union {
        CCartesianCoordinates *m_cartesianCordinates;
        CInternalCoordinates  *m_internalCoordinates;
        CFractionCoordinates  *m_fractionCoordinates;
    };
};


}
#endif // CCOORDINATESET_H
