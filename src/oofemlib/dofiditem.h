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
// FILE: dofiditem.h
//

#ifndef dofiditemh
#define dofiditemh

#include "enumitem.h"


/* mask definning the physical meaning of particular DOF in node.
 * mask array are also used in elements, where these arrays
 * are determining required DOFs needed by element and which are then
 * requsted on particular nodes. They are used in function
 * Node::giveLOcationArray which returns equations numbers
 * corresponding to selected dofs.
 */
typedef char DofID;
/**
 * Type representing particular dof type. Values of this type describe the physical meaning of
 * available DOFs.
 * Note: implementation of Node::computeGNTransformation rely on D_u, D_v and D_w (R_u, R_v, R_w) order.
 * Do not change their order and do not insert any values between these values.
 */
#define DofIDItem_DEF \
    ENUM_ITEM_WITH_VALUE(Undef, 0) /* Erorr value */ \
    ENUM_ITEM_WITH_VALUE(D_u, 1) /* u-displacement (in direction of x-axis) */   \
    ENUM_ITEM_WITH_VALUE(D_v, 2) /* v-displacement (in direction of y-axis) */   \
    ENUM_ITEM_WITH_VALUE(D_w, 3) /* w-displacement (in direction of z-axis) */   \
    ENUM_ITEM_WITH_VALUE(R_u, 4) /* rotation around x-axis (right hand rule assumed) */   \
    ENUM_ITEM_WITH_VALUE(R_v, 5) /* rotation around y-axis */   \
    ENUM_ITEM_WITH_VALUE(R_w, 6) /* rotation around z-axis */   \
  \
    ENUM_ITEM_WITH_VALUE(V_u, 7) /* u-velocity (in direction of x-axis) */   \
    ENUM_ITEM_WITH_VALUE(V_v, 8) /* v-velocity (in direction of y-axis) */   \
    ENUM_ITEM_WITH_VALUE(V_w, 9) /* w-velocity (in direction of z-axis) */   \
  \
    ENUM_ITEM_WITH_VALUE(T_f, 10) /* temperature field */  \
    ENUM_ITEM_WITH_VALUE(P_f, 11) /* pressure field */  \
    ENUM_ITEM_WITH_VALUE(G_0, 12) /* DOF for gradient formulation no. 0 */  \
    ENUM_ITEM_WITH_VALUE(G_1, 13) /* DOF for gradient formulation no. 1 */  \
    ENUM_ITEM_WITH_VALUE(C_1, 14) /* mass concentration of the first constituent */

enum DofIDItem {
    DofIDItem_DEF
};

#undef ENUM_ITEM
#undef ENUM_ITEM_WITH_VALUE
#undef enumitem_h


char *__DofIDItemToString(DofIDItem _value);

// max length of text string with DofIdName + 1
// see Dof::giveDofIDName function
#define DofIdNameMaxLength 5

#endif // dofiditem_h
