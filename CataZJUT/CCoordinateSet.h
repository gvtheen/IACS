#ifndef CCOORDINATESET_H
#define CCOORDINATESET_H
#include "../Util/Point-Vector.h"

using util::Point3;

namespace CATAZJUT{

class CCartesianCoordinates;
class CInternalCoordinates;
class CFractionCoordinates;

class CCoordinateSet
{
public:
    // enumerations
    enum Type {
        None=0,
        Cartesian=1,
        Internal=2,
        Fraction=3
    };

    // construction and destruction
    CCoordinateSet();
    explicit CCoordinateSet(CCartesianCoordinates *coordinates);
    explicit CCoordinateSet(CInternalCoordinates  *coordinates);
    explicit CCoordinateSet(CFractionCoordinates  *coordinates);
    CCoordinateSet(const CCoordinateSet &other);
    ~CCoordinateSet();

    // properties
    Type type() const;
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
    Type m_type;
    union {
        CCartesianCoordinates *m_cartesianCordinates;
        CInternalCoordinates  *m_internalCoordinates;
        CFractionCoordinates  *m_fractionCoordinates;
    };
};


}
#endif // CCOORDINATESET_H
