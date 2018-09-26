#ifndef CINTERNALCOORDINATES_H
#define CINTERNALCOORDINATES_H
#include "Point-Vector.h"

namespace CATAZJUT{

class CCartesianCoordinates;
class CInternalCoordinatesPrivate;
class CConfigurationBase;

class CInternalCoordinates
{
public:
    // construction and destruction
    CInternalCoordinates();
    CInternalCoordinates(CConfigurationBase *coordinates,size_t size_m);
    CInternalCoordinates(const CInternalCoordinates &coordinates);
    ~CInternalCoordinates();

    // properties
    size_t size() const;
    bool isEmpty() const;

    // coordinates
    void setCoordinates(size_t row, double r, double theta = 0, double phi = 0);
    void setCoordinatesRadians(size_t row, double r, double theta = 0, double phi = 0);
    Point3 coordinates(size_t row) const;
    Point3 coordinatesRadians(size_t row) const;
    void setConnections(size_t row, int a, int b= -1, int c = -1);
    void addInternalCoordinates(size_t a,double r, int b=-1, double theta=0,\
                                int c=-1, double phi=0);
    void addInternalCoordinates(Point3i,Point3);

    Point3i connections(size_t row) const;
    CConfigurationBase* configuration()const;

    // conversions
    CCartesianCoordinates* toCartesianCoordinates() const;

    // operators
    CInternalCoordinates& operator=(const CInternalCoordinates &coordinates);

private:
    CInternalCoordinatesPrivate* const m_data;
             CConfigurationBase* m_pConfiguration;
};

}
#endif // CINTERNALCOORDINATES_H
