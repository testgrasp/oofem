/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2008   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

//
// FILE: sparsematrixtype.h
//

#ifndef sparsematrixtype_h
#define sparsematrixtype_h

/**
 * Enumerative type used to identify the sparse matrix type
 */
enum SparseMtrxType {
    SMT_Skyline,       // symmetric skyline
    SMT_SkylineU,      // unsymmetric skyline
    SMT_CompCol,       // compressed column
    SMT_DynCompCol,    // dynamically growing compressed column
    SMT_SymCompCol,    // symmetric compressed column
    SMT_DynCompRow,    // dynamically growing compressed row
    SMT_SpoolesMtrx,   // spooles sparse mtrx representation
    SMT_PetscMtrx,     // PETSc library mtrx representation
    SMT_DSS_sym_LDL, // Richard Vondracek's sparse direct solver
    SMT_DSS_sym_LL, // Richard Vondracek's sparse direct solver
    SMT_DSS_unsym_LU // Richard Vondracek's sparse direct solver
};

#endif // sparsematrixtype_h
