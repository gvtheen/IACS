#ifndef CATASTUNIVERSEDEFINE_H
#define CATASTUNIVERSEDEFINE_H

namespace CATAZJUT{
namespace DEFINED{

enum Coordinate_Type{
     Cartesian,
     Internal,
     Fraction
};
enum Dimensional_Type{
     Periodic,
     Molecule
};
enum StructureFile_type{
     POSCAR,
     MOL,
     CIF,
     GJF,
     CAR
};
typedef enum Coordinate_Type  CoordinateType;
typedef enum Dimensional_Type DimensionalType;
typedef enum StructureFile_type StructureFiletype;

} // namespace DEFINED

} //namespace CATAZJUT


#endif
