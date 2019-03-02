#ifndef CATASTUNIVERSEDEFINE_H
#define CATASTUNIVERSEDEFINE_H

namespace CATAZJUT{
namespace DEFINED{

typedef enum Coordinate_Type{
     Cartesian=0x0601,
     Internal=0x0602,
     Fraction=0x0603,
     None=0x0604,
}CoordinateType;

typedef enum Dimensional_Type{
     Periodic=0x0605,
     Molecule=0x0606
}DimensionalType;

typedef enum StructureFile_type{
     POSCAR=0x0607,
     MOL=0x0608,
     CIF=0x0609,
     GJF=0x0610,
     CAR=0x0611,
     XYZ=0x0612
}StructureFiletype;

typedef enum UNITCELL_TYPE{
     Reciprocal=0x0701,
     Bravais=0x0702
}UnitcellType;

typedef enum crystal_001{
           Cubic=0x40A1,
           Hexagonal=0x40A2,
           Trigonal=0x40A3,
           Tetragonal=0x40A4,
           Orthorhombic=0x40A5,
           Monoclinic=0x40A6,
           Triclinic=0x40A17,
           ERROR=0x40A38
}CrystalSystemType;



} // namespace DEFINED

} //namespace CATAZJUT


#endif
