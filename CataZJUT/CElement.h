#ifndef CCELEMENT_H
#define CCELEMENT_H
#include <string>
#include <map>
namespace CATAZJUT{


class CElement
{
public:
    // typedefs
    typedef unsigned char AtomicNumberType;

    typedef struct e0001{
                size_t n;
                char   l;
                size_t ele;
            }Atomic_Block;

    // construction and destruction
    inline CElement(AtomicNumberType atomicNumber = 0);
    CElement(const char *symbol);
    CElement(const std::string &symbol);
    CElement(const CElement&);

    CElement* clone();
    // properties
    inline void setAtomicNumber(AtomicNumberType atomicNumber);
    inline AtomicNumberType atomicNumber() const;
    std::string symbol() const;
    std::string name() const;
    int period() const;
    double mass() const;
    double exactMass() const;
    double ionizationEnergy() const;
    double electronAffinity() const;
    double electronegativity() const;
    double covalentRadius() const;
    double vanDerWaalsRadius() const;
    double boilingPoint() const;
    double meltingPoint() const;
    int expectedValence() const;
    bool isValid() const;

    bool isMetal() const;
    bool isTransitionMetal()const;
    bool is1stTransitionMetal()const;
    bool is2ndTransitionMetal()const;
    bool isNobelMetal()const;
    bool isMainGroupElement()const;
    bool isLanthanideMetal()const;
    bool isActinidesMetal()const;
    bool isNonmetal() const;

    bool isSblockElement()const;
    bool isPblockElement()const;
    bool isDblockElement()const;
    bool isFblockElement()const;
    bool isSPblockElement()const;

    bool isAlkaliMetal()const;
    bool isAlkalineEarthMeta()const;
    bool isBGroupElement()const;
    bool isCGroupElement()const;
    bool isNGroupElement()const;
    bool isOGroupElement()const;
    bool isHalogen()const;
    bool isZeroGroupElement()const;

    std::vector<Atomic_Block> valentConfiguration();
    std::string valentConfigurationStr();
    size_t maxCoordinationNum();
    //

    // operators
    inline bool operator==(const CElement &CElement) const;
    inline bool operator!=(const CElement &CElement) const;

    // static methods
    static CElement fromName(const std::string &name);
    static CElement fromName(const char *name);
    static CElement fromSymbol(const std::string &symbol);
    static CElement fromSymbol(const char *symbol);
    static CElement fromSymbol(const char *symbol, size_t length);
    static CElement fromSymbol(char symbol);
    static bool isValidAtomicNumber(AtomicNumberType atomicNumber);
    static bool isValidSymbol(const std::string &symbol);

private:
    AtomicNumberType m_atomicNumber;
    inline size_t ElectronNumberofShell(char p);
};

inline CElement::CElement(AtomicNumberType atomicNumber)
{
    m_atomicNumber = atomicNumber;
}

// --- Properties ---------------------------------------------------------- //
/// Sets the atomic number for the CElement to \p atomicNumber.
inline void CElement::setAtomicNumber(AtomicNumberType atomicNumber)
{
    m_atomicNumber = atomicNumber;
}

/// Returns the atomic number of the CElement.
inline CElement::AtomicNumberType CElement::atomicNumber() const
{
    return m_atomicNumber;
}

// --- Operators ----------------------------------------------------------- //
inline bool CElement::operator==(const CElement &CElement) const
{
    return m_atomicNumber == CElement.m_atomicNumber;
}

inline bool CElement::operator!=(const CElement &CElement) const
{
    return m_atomicNumber != CElement.m_atomicNumber;
}

inline size_t CElement::ElectronNumberofShell(char p)
{
    size_t res;
    switch ( p)
    {
       case 's':
           res=2;
           break;
       case 'p':
           res=6;
           break;
       case 'd':
           res=10;
           break;
       case 'f':
           res=14;
           break;
       default:
           res=0;
           break;
    }
    return res;
}



}   //namespace



#endif


