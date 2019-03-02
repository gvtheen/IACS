/******************************************************************************
**
** Copyright (C) 2019-2031 Dr.Gui-lin Zhuang <glzhuang@zjut.edu.cn>
** All rights reserved.
**
**
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
** OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**
******************************************************************************/
#include "CInternalCoordinates.h"
#include "CCartesianCoordinates.h"
#include "CConfigurationPrivateData.h"
#include "Constant.h"
namespace CATAZJUT{
/// Creates a new, empty set of internal coordinates.

class CInternalCoordinatesPrivate
{
public:
    CInternalCoordinatesPrivate();
    ~CInternalCoordinatesPrivate();
    //std::size_t size;
    std::vector<Point3i> connections;
    std::vector<Point3>  coordinates;
};
// definition of CInternalCoordinatesPrivate
//
//
//
//
CInternalCoordinatesPrivate::CInternalCoordinatesPrivate()
{

}
CInternalCoordinatesPrivate::~CInternalCoordinatesPrivate()
{

}
/*
definition of CInternalCoordinates



*/
CInternalCoordinates::CInternalCoordinates()
    : m_data(new CInternalCoordinatesPrivate)
{
//    m_data->connections = nullptr;
//    m_data->coordinates = nullptr;
}
CInternalCoordinates::CInternalCoordinates(CConfigurationBase *structure)
: m_data(new CInternalCoordinatesPrivate)
{
    this->m_pConfiguration=structure;
}

/// Creates a new internal coordinate set with \p size rows.
CInternalCoordinates::CInternalCoordinates(CConfigurationBase *m_conf,size_t size_m)
    : m_data(new CInternalCoordinatesPrivate),m_pConfiguration(m_conf)
{

    m_data->connections.insert(m_data->connections.begin(),size_m, Point3i().setZero());
    m_data->coordinates.insert(m_data->coordinates.begin(),size_m, Point3().setZero());

}

/// Creates a new internal coordinates object as a copy of
/// \p coordinates.
CInternalCoordinates::CInternalCoordinates(const CInternalCoordinates &coordinates)
    : m_data(new CInternalCoordinatesPrivate)
{
    //size_t tempsize = coordinates.m_data->connections.size();

    m_data->connections.assign(coordinates.m_data->connections.begin(),coordinates.m_data->connections.end());
    m_data->coordinates.assign(coordinates.m_data->coordinates.begin(),coordinates.m_data->coordinates.end());
    m_pConfiguration = coordinates.configuration();
}

/// Destroys the internal coordinates object.
CInternalCoordinates::~CInternalCoordinates()
{
    delete m_data;
}

// --- Properties ---------------------------------------------------------- //
/// Returns the number of rows of coordinates.
size_t CInternalCoordinates::size() const
{
    return m_data->connections.size();
}

/// Returns \c true if the internal coordinates object contains
/// no coordinates (i.e. size() == \c 0).
bool CInternalCoordinates::isEmpty() const
{
    return size() == 0;
}

// --- Coordinates --------------------------------------------------------- //
/// Sets the distance, angle, and torsion at \p row to \p r,
/// \p theta and \p phi respectively. The angles are in degrees.
void CInternalCoordinates::setCoordinates(size_t row, double r, double theta, double phi)
{
    assert(row < size());
    m_data->coordinates[row]<<r,theta,phi;
}

/// Sets the distance, angle, and torsion at \p row to \p r,
/// \p theta and \p phi respectively. The angles are in radians.
void CInternalCoordinates::setCoordinatesRadians(size_t row, double r, double theta, double phi)
{
    assert(row < size());

    double r1 = r;
    double theta1 = theta * CATAZJUT::constants::RadiansToDegrees;
    double phi1 = phi * CATAZJUT::constants::RadiansToDegrees;
    m_data->coordinates[row]<<r1,theta1,phi1;
}
void CInternalCoordinates::addInternalCoordinates(size_t a,double r, int b, double theta,int c, double phi)
{
    size_t row= size();
    setCoordinates(row,r,theta,phi);
    setConnections(row,a,b,c);
}
void CInternalCoordinates::addInternalCoordinates(Point3i _connection,Point3 _coordinate)
{
    size_t row= size();
    setCoordinates(row,_coordinate[0],_coordinate[1],_coordinate[2]);
    setConnections(row,_connection[0],_connection[1],_connection[2]);
}
/// Returns the distance, angle, and torsion coordinates at \p row.
/// The returned angles are in degrees.
Point3 CInternalCoordinates::coordinates(size_t row) const
{
    assert(row < size());

    return m_data->coordinates[row];
}
CConfigurationBase* CInternalCoordinates::configuration()const
{
    return this->m_pConfiguration;
}
/// Returns the distance, angle, and torsion coordinates at \p row.
/// The returned angles are in radians.
Point3 CInternalCoordinates::coordinatesRadians(size_t row) const
{
    assert(row < size());
    Point3 tmp=m_data->coordinates.at(row);
    tmp[1]=tmp[1]*CATAZJUT::constants::DegreesToRadians;
    tmp[2]=tmp[2]*CATAZJUT::constants::DegreesToRadians;
    return tmp;
}

/// Sets the connections for the coordinates at \p row to \p a, \p b
/// and \p c.
void CInternalCoordinates::setConnections(size_t row, int a, int b, int c)
{
    assert(row < size());

    m_data->connections[row]<<a,b,c;
}

/// Returns the connections for the coordinates at \p row.
Point3i CInternalCoordinates::connections(size_t row) const
{
    assert(row < size());

    return m_data->connections[row];
}

// --- Conversions --------------------------------------------------------- //
/// Converts the internal coordinates into cartesian coordinates.
CCartesianCoordinates* CInternalCoordinates::toCartesianCoordinates() const
{
    CCartesianCoordinates *cartesianCoordinates = new CCartesianCoordinates();

    // set positions for the first three atoms
    if(this->size() > 0){
        cartesianCoordinates->setPosition(0, Point3(0, 0, 0));

        if(this->size() > 1){
            double r1 = coordinates(1)(0);
            cartesianCoordinates->setPosition(1, Point3(r1, 0, 0));

            if(this->size() > 2){
                   double r2 = coordinates(2)(0);
                double theta = coordinates(2)(1);

                double x = r2 * cos((180.0 - theta) * CATAZJUT::constants::DegreesToRadians);
                double y = r2 * sin((180.0 - theta) * CATAZJUT::constants::DegreesToRadians);

                cartesianCoordinates->setPosition(2, Point3(r1 + x, y, 0));
            }
        }
    }

    // set positions for the rest of the atoms
    for(size_t i = 3; i < this->size(); i++){
        Point3 coordinates = this->coordinates(i);
        double r = coordinates(0);
        double theta = coordinates(1);
        double phi = coordinates(2);

        double sinTheta = sin(theta * CATAZJUT::constants::DegreesToRadians);
        double cosTheta = cos(theta * CATAZJUT::constants::DegreesToRadians);
        double sinPhi = sin(phi * CATAZJUT::constants::DegreesToRadians);
        double cosPhi = cos(phi * CATAZJUT::constants::DegreesToRadians);

        double x = r * cosTheta;
        double y = r * cosPhi * sinTheta;
        double z = r * sinPhi * sinTheta;

        Point3i connections = this->connections(i);

        const Point3 &a = cartesianCoordinates->position(connections(2));
        const Point3 &b = cartesianCoordinates->position(connections(1));
        const Point3 &c = cartesianCoordinates->position(connections(0));

        Vector3 ab = (b - a);
        Vector3 bc = (c - b).normalized();
        Vector3 n = ab.cross(bc).normalized();
        Vector3 ncbc = n.cross(bc);

        Eigen::Matrix<double, 3, 3> M;
        M << bc.x(), ncbc.x(), n.x(),
             bc.y(), ncbc.y(), n.y(),
             bc.z(), ncbc.z(), n.z();

        Point3 d = (M * Point3(-x, y, z)) + c;
        cartesianCoordinates->setPosition(i, d);
    }

    return cartesianCoordinates;
}

// --- Operators ----------------------------------------------------------- //
CInternalCoordinates& CInternalCoordinates::operator=(const CInternalCoordinates &coordinates)
{
    if(&coordinates == this){
        return *this;
    }

    //std:size_t tempsize = coordinates.m_data->connections.size();
//    m_data->connections = new std::vector<Point3i>();
//    m_data->coordinates = new std::vector<Point3>();

    m_data->connections.assign(coordinates.m_data->connections.begin(),coordinates.m_data->connections.end());
    m_data->coordinates.assign(coordinates.m_data->coordinates.begin(),coordinates.m_data->coordinates.end());

    return *this;
}


}


