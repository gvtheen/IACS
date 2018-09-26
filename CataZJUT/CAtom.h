#ifndef CATOM_H
#define CATOM_H

#include "gacata.h"
#include <vector>
#include <boost/function.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/iterator/filter_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include "CElement.h"
#include "Point-Vector.h"

namespace CATAZJUT{

class CConfigurationBase;
class CBond;
class CFragment;

class CAtom
{
    public:
        typedef boost::iterator_range<std::vector<CBond *>::const_iterator> BondRange;
        typedef boost::iterator_range<boost::transform_iterator<  \
               boost::function<CAtom* (CBond*)>,std::vector<CBond*>::const_iterator> > NeighborRange;


        CAtom(CConfigurationBase*,size_t);
        virtual ~CAtom();

        CAtom* Clone();

        void SetPosition(const Point3 &position);
        void SetPosition(double x, double y, double z);
        Point3 position() const;
        double x() const;
        double y() const;
        double z() const;

//        void setMassNumber(MassNumberType massNumber);
        void setPartialCharge(double charge);
        double PartialCharge()const;
        // property of atom
        std::string Symbol() const;
        std::string Name() const;
        double Mass() const;
        double CovalentRadius() const;
        double VanDerWaalsRadius() const;
        double Electronegativity() const;
        double AtomNumber() const;
        int    ExpectedValence() const;
        int    CoorditionNum ()const;
        double distance(const CAtom *atom) const;
        size_t maxCoordinationNum()const;
        std::string valentConfigurationStr()const;

        inline CElement element() const;
        inline bool IsElement(const CElement &element) const;
        inline CConfigurationBase* Configuration() const;
        inline size_t index() const;

        //geometry
        CBond* bond(size_t index) const;
        BondRange bonds() const;
        size_t bondCount() const;
        CBond* bondTo(const CAtom *atom) const;
        CAtom* neighbor(size_t index) const;
        NeighborRange neighbors() const;
        size_t neighborCount() const;
        size_t neighborCount(const CElement &element) const;
        bool isBondedTo(const CAtom *atom) const;
        bool isBondedTo(const CElement &element) const;
        bool isTerminal() const;
        bool isTerminalHydrogen() const;
        CFragment* fragment() const;

        enum AtomName{
            Hydrogen = 1,
            Helium = 2,
            Lithium = 3,
            Beryllium = 4,
            Boron = 5,
            Carbon = 6,
            Nitrogen = 7,
            Oxygen = 8,
            Fluorine = 9,
            Neon = 10,
            Sodium = 11,
            Magnesium = 12,
            Aluminum = 13,
            Silicon = 14,
            Phosphorus = 15,
            Sulfur = 16,
            Chlorine = 17,
            Argon = 18,
            Potassium = 19,
            Calcium = 20,
            Scandium = 21,
            Titanium = 22,
            Vanadium = 23,
            Chromium = 24,
            Manganese = 25,
            Iron = 26,
            Cobalt = 27,
            Nickel = 28,
            Copper = 29,
            Zinc = 30,
            Gallium = 31,
            Germanium = 32,
            Arsenic = 33,
            Selenium = 34,
            Bromine = 35,
            Krypton = 36,
            Rubidium = 37,
            Strontium = 38,
            Yttrium = 39,
            Zirconium = 40,
            Niobium = 41,
            Molybdenum = 42,
            Technetium = 43,
            Ruthenium = 44,
            Rhodium = 45,
            Palladium = 46,
            Silver = 47,
            Cadmium = 48,
            Indium = 49,
            Tin = 50,
            Antimony = 51,
            Tellurium = 52,
            Iodine = 53,
            Xenon = 54,
            Cesium = 55,
            Barium = 56,
            Lanthanum = 57,
            Cerium = 58,
            Praseodymium = 59,
            Neodymium = 60,
            Promethium = 61,
            Samarium = 62,
            Europium = 63,
            Gadolinium = 64,
            Terbium = 65,
            Dysprosium = 66,
            Holmium = 67,
            Erbium = 68,
            Thulium = 69,
            Ytterbium = 70,
            Lutetium = 71,
            Hafnium = 72,
            Tantalum = 73,
            Tungsten = 74,
            Rhenium = 75,
            Osmium = 76,
            Iridium = 77,
            Platinum = 78,
            Gold = 79,
            Mercury = 80,
            Thallium = 81,
            Lead = 82,
            Bismuth = 83,
            Polonium = 84,
            Astatine = 85,
            Radon = 86,
            Francium = 87,
            Radium = 88,
            Actinium = 89,
            Thorium = 90,
            Protactinium = 91,
            Uranium = 92,
            Neptunium = 93,
            Plutonium = 94,
            Americium = 95,
            Curium = 96,
            Berkelium = 97,
            Californium = 98,
            Einsteinium = 99,
            Fermium = 100,
            Mendelevium = 101,
            Nobelium = 102,
            Lawrencium = 103,
            Rutherfordium = 104,
            Dubnium = 105,
            Seaborgium = 106,
            Bohrium = 107,
            Hassium = 108,
            Meitnerium = 109
        };
    private:
        friend class CCartesianCoordinates;
        friend class CConfigurationBase;
    private:
        CConfigurationBase *m_pConfiguration;
                    size_t  m_index;

};

}
#include "CAtom_sub.h"

#endif // CATOM_H
