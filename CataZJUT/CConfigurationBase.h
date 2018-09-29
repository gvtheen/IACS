#ifndef CCONFIGURATIONBASE_H
#define CCONFIGURATIONBASE_H
#include <string>
#include <map>
#include <vector>
#include <boost/range/iterator_range.hpp>
#include <boost/shared_ptr.hpp>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include "../Util/Point-Vector.h"
#include "../Util/Bitset.h"
#include "CCoordinateSet.h"
#include "../CalcZJUT/CParameter.h"
#include "CatalystUniverseDefine.h"
/*
In this class, main functions have 3 aspects:
1)
*/

using util::Bitset;
using util::Point3;
using util::Point3i;

namespace CALCZJUT{
      class CParameter;
      class CIOPoscar;
      class CIOGjf;
      class CIOMol;
      class CIOCellFile;
      class CCalcSupportBase;
      class CCalcClusterSupport;
      class CCalcVASP;
      class CCalcDMol;
      class CCalcGaussian;
      class CCalcLammps;
      class CCalc2DSupport;
}
namespace CATAZJUT{

class CCartesianCoordinates;
class CElement;
class CAtom;
class CBond;
class CFragment;
class CBondTolerance;
class CConfigurationPrivateData;
class CFractionCoordinates;

// declarazation of class
class CConfigurationBase
{
    public:
        typedef boost::iterator_range<std::vector<CAtom *>::const_iterator> AtomRange;
        typedef boost::iterator_range<std::vector<CBond *>::const_iterator> BondRange;
        typedef boost::iterator_range<std::vector<CFragment *>::const_iterator> FragmentRange;
        typedef boost::iterator_range<std::vector<boost::shared_ptr<CCoordinateSet> >::const_iterator> CoordinateSetRange;

        CConfigurationBase(CALCZJUT::CParameter*);
        CConfigurationBase(CConfigurationBase&);
        virtual ~CConfigurationBase();

        //chemical name
        std::vector<std::pair<std::string,size_t>>* composition();
        std::string formula();

        //interface function  from other objects
        //Cartesian coordinate
        void addAtom(std::string& elementName, Point3& position);
        void addAtom(std::string& elementName, double& x, double& y, double& z);
        //Internal Coordinates
        void addAtom(std::string& elementName, Point3i& connection,Point3& coordinate);
        void addAtom(std::string& elementName, int atom1=-1, double dist=-1, int atom2=-1, double angle=0.0, int atom3=-1, double twist=0.0);

        inline size_t size() const;
        inline bool isEmpty() const;
        inline void checkElement(CElement& );
        //
        CAtom* addAtom(const CElement &element);
        CAtom* addAtomCopy(const CAtom *atom);
        void removeAtom(CAtom *atom);
        void removeAtoms(const std::vector<CAtom *> &atoms);
        template<typename Range> void removeAtoms(Range range);
        inline CAtom* atom(size_t index) const;
        inline AtomRange atoms() const;
        inline size_t atomCount() const;
        size_t atomCount(const CElement &element) const;
        bool contains(const CAtom *atom) const;
        bool contains(const CElement &element) const;

        //bond
        void perceiveBonds();
        CBond* addBond(CAtom *a, CAtom *b);
        CBond* addBond(size_t a, size_t b);
        void removeBond(CBond *bond);
        void removeBond(CAtom *a, CAtom *b);
        void removeBond(size_t a, size_t b);
        void removeBonds(const std::vector<CBond *> &bonds);
        template<typename Range> void removeBonds(Range range);
        CBond* bond(size_t index) const;
        CBond* bond(const CAtom *a, const CAtom *b) const;
        CBond* bond(size_t a, size_t b) const;
        BondRange bonds() const;
        size_t bondCount() const;
        bool contains(const CBond *bond) const;
        void clear();

        CCartesianCoordinates* coordinates();
         CInternalCoordinates* Internalcoordinates();
         virtual CFractionCoordinates* Fractioncoordinates()=0;

        void addCoordinateSet(const boost::shared_ptr<CCoordinateSet> &coordinates);
        void addCoordinateSet(CCartesianCoordinates *coordinates);
        void addCoordinateSet(CInternalCoordinates *coordinates);
        void addCoordinateSet(CFractionCoordinates *coordinates);
        bool removeCoordinateSet(const boost::shared_ptr<CCoordinateSet> &coordinates);
        boost::shared_ptr<CCoordinateSet> coordinateSet(size_t index) const;
        boost::shared_ptr<CCoordinateSet> coordinateSet(CCoordinateSet::Type type) const;
        CoordinateSetRange coordinateSets() const;
        size_t coordinateSetCount() const;

        //constraint

        Bitset constraintBit()const;
        void setConstraintBit(const Bitset& othr);
        //fragment
        bool isFragments(std::vector<size_t>*);   //identify whether the molecule have fragments
        void perceiveFragments();
        CFragment* fragmentForAtom(const CAtom *atom);
        FragmentRange fragments();
        CFragment* fragment(size_t index);
        void setFragmentsPerceived(bool perceived);
        bool fragmentsPerceived()const;

        double distance(const CAtom *a, const CAtom *b);
        double bondAngle(const CAtom *a, const CAtom *b, const CAtom *c);
        double torsionAngle(const CAtom *a, const CAtom *b, const CAtom *c, const CAtom *d);
        double wilsonAngle(const CAtom *a, const CAtom *b, const CAtom *c, const CAtom *d);
        void setCenter(const Point3 &position);
        void setCenter(double x, double y, double z);
        Point3 center();
        //
        void setCoordinateType(const CATAZJUT::DEFINED::CoordinateType);
        CATAZJUT::DEFINED::CoordinateType coordinateType()const;
        void setDimensionalType(const CATAZJUT::DEFINED::DimensionalType);
        CATAZJUT::DEFINED::DimensionalType dimensionalType()const;

    protected:
        CConfigurationPrivateData *m_pData;
             CALCZJUT::CParameter *m_pParameter;
    private:
        friend class CAtom;
        friend class CBond;
        friend class CBondTolerance;
        friend class CCartesianCoordinates;
        friend class CSphere;
        friend class CFragment;
        friend class CPeriodicFramework;

        friend class CALCZJUT::CCalcClusterSupport;
        //output

        friend class CALCZJUT::CIOPoscar;
        friend class CALCZJUT::CIOGjf;
        friend class CALCZJUT::CIOMol;
        friend class CALCZJUT::CIOCellFile;
        friend class CALCZJUT::CCalcSupportBase;
        friend class CALCZJUT::CCalc2DSupport;

        friend class CALCZJUT::CCalcVASP;
        friend class CALCZJUT::CCalcDMol;
        friend class CALCZJUT::CCalcGaussian;
        friend class CALCZJUT::CCalcLammps;

        mutable CCartesianCoordinates     *m_pCartesian;
        std::vector<CElement>      m_Element;         // contain the chemical property of all atoms
        std::vector<CAtom*>        m_Atom;            // contain the behavior of all atoms
        CBondTolerance            *m_pBondEvaluator;   //
        // type of coordinate
        CATAZJUT::DEFINED::CoordinateType    m_CoordinateType;
        // type of dimensional
        CATAZJUT::DEFINED::DimensionalType   m_DimensionalType;

        Bitset constraintBits;   //0:no constraint; 1: constraint

};


}
#include "CConfigurationBase_sub.h"

#endif // CCONFIGURATIONBASE_H
