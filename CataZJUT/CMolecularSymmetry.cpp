#include<iostream>
#include <string>
#include <string.h>
#include "CConfigurationBase.h"
#include "CAtom.h"
#include "../Util/foreach.h"
#include "CMolecularSymmetry.h"
extern "C"{
#include "../Molsymm/msym.h"
}

namespace CATAZJUT{

CMolecularSymmetry::CMolecularSymmetry(CConfigurationBase* mbf)
:m_pConfiguration(mbf)
{
    point_group = new char[6];

    m_AtomicBits.resize(mbf.size());
    m_AtomicBits.set();  //set all elements to 1;
}
CMolecularSymmetry(CConfigurationBase*  mbf,Bitset dat)
::m_pConfiguration(mbf),m_AtomicBits(dat)
{
    point_group = new char[6];
}
CMolecularSymmetry::~CMolecularSymmetry()
{
    if(point_group==NULL)
        delete []point_group;
}
char* CMolecularSymmetry::GetPointGroup()
{
    msym_error_t ret = MSYM_SUCCESS;
    msym_element_t *elements = NULL;

    char *error = NULL;
    double cm[3], radius = 0.0;

    /* Do not free these variables */
    msym_element_t *melements = NULL;
   // msym_basis_function_t *mbfs = NULL;
    /* these are not mutable */
//    const msym_symmetry_operation_t *msops = NULL;
    const msym_subgroup_t *msg = NULL;
//    const msym_subrepresentation_space_t *msrs = NULL;
//    const msym_character_table_t *mct = NULL;
//    const msym_equivalence_set_t *mes = NULL;
//    int mesl = 0;
//    double *irrep = NULL;
//    msym_basis_function_t *bfs = NULL;

    int msgl = 0,  mlength = 0;
    //msopsl = 0,msrsl = 0, mbfsl = 0, bfsl = 0;
    /* This function reads xyz files.
     * It initializes an array of msym_element_t to 0,
     * then sets the coordinates and name of the elements */
    int length;// = read_xyz(in_file, &elements);
    //if(length <= 0) return -1;
    length = this->m_AtomicBits.count();
    melements = new msym_element_t[length];

    CAtom* atom=nullptr;
   // foreach(CAtom* atom, this->m_pConfiguration->atoms()){
    size_t pos=m_AtomicBits.find_first();

    if(pos==Bitset::npos){
        Log::Error<<" Error in the computation of molecular symmetry! CMolecularSymmetry::GetPointGroup()!\n";
        boost::throw_exception(std::runtime_error(" Error in the computation of molecular symmetry! CMolecularSymmetry::GetPointGroup()"));
    }
     size_t index=0;
     while(pos!=Bitset::npos){
        atom=this->m_pConfiguration->atom(pos);
        strcpy( melements[index].name, atom->Symbol().c_str() );
        melements[index].v[0]=atom->position()[0];
        melements[index].v[1]=atom->position()[1];
        melements[index].v[2]=atom->position()[2];
        index++;
        pos = m_AtomicBits.find_next(pos);
    }

    /* Create a context */
    msym_context ctx = msymCreateContext();
    msym_thresholds_t *thresholds=NULL;

    if(NULL != thresholds){
        if(MSYM_SUCCESS != (ret = msymSetThresholds(ctx, thresholds))) goto err;
    }

    /* Use default thresholds otherwise call:
     * msymSetThresholds(msym_context ctx, msym_thresholds_t *thresholds); */

    /* Set elements */
    if(MSYM_SUCCESS != (ret = msymSetElements(ctx, length, elements))) goto err;
    /* Get elements msym elements */
    if(MSYM_SUCCESS != (ret = msymGetElements(ctx, &mlength, &melements))) goto err;
    /* These are no longer needed, internal versions of these are kept in the context,
     * They are indexed in the same way that they have been allocated.
     * I.e. during orbital symmetrization or when getting the symmetrized LCAO,
     * the coefficients will correspond to the same indexing as "orbitals",
     * this is the main reason for the two levels of indirection */
    delete elements;
    elements = NULL;

    /* Some trivial information */
    if(MSYM_SUCCESS != (ret = msymGetCenterOfMass(ctx,cm))) goto err;
    if(MSYM_SUCCESS != (ret = msymGetRadius(ctx,&radius))) goto err;
//    printf("Molecule has center of mass [%lf; %lf; %lf] "
//           "and a radius of %lf\n",cm[0],cm[1],cm[2],radius);
    /* Find molecular symmetry */
    if(MSYM_SUCCESS != (ret = msymFindSymmetry(ctx))) goto err;
    /* Get the point group name */
    if(MSYM_SUCCESS != (ret = msymGetPointGroupName(ctx, sizeof(char[6]), point_group))) goto err;
    if(MSYM_SUCCESS != (ret = msymGetSubgroups(ctx, &msgl, &msg))) goto err;

    return point_group;

err:
    return error;
}


}
