//==============================================================================
//
//  System Simulation Group
//  Jochen Haerdtlein
//
//  University Erlangen-Nuremberg
//  Department of Computational Science
//  System Simulation Group
//  Cauerstrasse 6
//  91058 Erlangen
//  Germany
//
//------------------------------------------------------------------------------
//
//  Copyright (C) 2004 University Erlangen-Nuremberg
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
//------------------------------------------------------------------------------
//
//  $Id: Colsamm.h,v 1.29 2006/07/25 08:11:31 jochen Exp $
//
//==============================================================================
      //  _______  _______  __       _______  _______  ___  ____  ___  ____   // 
     //  / _____/ / ___  / / /      / _____/ / ___  / / _ |/_  / / _ |/_  /  //
    //  / /      / /  / / / /      / /____  / /__/ / / /|__// / / /|__// /  //
   //  / /      / /  / / / /      /____  / / ___  / / /    / / / /    / /  //
  //  / /____  / /__/ / / /____  _____/ / / /  / / / /    / / / /    / /  // 
 //  /______/ /______/ /______/ /______/ /_/  /_/ /_/    /_/ /_/    /_/  //
//                                                                      //
//==============================================================================
#ifndef COLSAMM_H
#define COLSAMM_H

     #include <iostream>
     #include <cmath>
     #include <typeinfo>
     #include <cassert>
     #include <complex>
     #include <vector>
     #include <typeinfo>
#define INITIALIZE_ARRAYS 
     namespace _COLSAMM_ 
        {
          #include "Loops.h"
          #include "Loops.C"

          #include "Polynom.h"
          #include "Polynom.C"

          #include "Gauss_Points.h"
          #include "Gauss_Points.C"

          #include "VertexComputations.h"
          #include "VertexComputations.C"

          #include "Mapping.h"
          #include "Mapping.C"

          #include "BasisFunctions.h"
          #include "BasisFunctions.C"
          #include "VectorBasisFunctions.C"
          #include "DerivativeBasisFunctions.C"

          #include "Domain.h"
          #include "Domain.C"
 
#ifdef FET_COLSAMM
          #include "Stencil.h"
          #include "Stencil.C"

          #include "CT_Manip.h"
          #include "CT_Manip.C"
#else 
          #include "Stencil_ET.h"
          #include "Stencil_ET.C"

          #include "CT_Manip_ET.h"
          #include "CT_Manip_ET.C"

          #include "Integrand_Vector.h"
          #include "Integrand_Vector.C"
#endif
          #include "Elements.h"
        }
#undef INITIALIZE_ARRAYS 
//-----------------------------------------------------------------------------
#endif // COLSAMM_H
//==============================================================================
//
//  Status:
//  -------
//
//  $Date: 2006/07/25 08:11:31 $
//  $Name:  $
//  $Author: jochen $
//  $Revision: 1.29 $
//
//------------------------------------------------------------------------------
//
//  History:
//  --------
//
//  $Log: Colsamm.h,v $
//  Revision 1.29  2006/07/25 08:11:31  jochen
//  fixed Haralds problem concenring the nan for boundary integrals, caused
//  by dividing the normal by the determinant, not the length of the normal!
//
//  Revision 1.28  2006/07/20 14:03:13  jochen
//  *** empty log message ***
//
//  Revision 1.27  2006/07/19 14:34:54  jochen
//  *** empty log message ***
//
//  Revision 1.26  2006/07/19 09:05:18  jochen
//
//  introduced the absolut value for testing whether the determinant is zero
//  or not !
//
//  Revision 1.25  2006/06/30 10:35:19  jochen
//  fixed problems concerning vector-vector multiplication on integrands
//
//  Revision 1.24  2006/06/30 09:27:34  jochen
//  fixed problems in FET implementation
//
//  Revision 1.23  2006/06/28 06:18:48  jochen
//  finishing the vector basis functions ...
//
//  Revision 1.22  2006/06/23 14:42:06  jochen
//  iiiiiiiii
//
//  Revision 1.21  2006/06/23 12:26:38  jochen
//  ooo
//
//  Revision 1.20  2006/06/23 10:12:29  jochen
//  .....
//
//  Revision 1.19  2006/06/09 10:12:46  jochen
//  splitted Mapping files ...
//
//  Revision 1.18  2006/06/09 10:01:26  jochen
//  splitted BasisFuncitons into separate files
//
//  Revision 1.17  2006/06/09 09:50:53  jochen
//  compressing the extension arrays
//
//  Revision 1.16  2006/05/31 09:41:54  jochen
//  faster version
//
//  Revision 1.15  2006/04/28 12:52:44  jochen
//  Introduced the easier generation of elments ...
//
//  Revision 1.14  2006/04/28 08:24:20  jochen
//  added CT_Manip.C
//
//  Revision 1.13  2006/04/28 08:06:12  jochen
//  added Stencil.C
//
//  Revision 1.12  2006/04/28 07:10:13  jochen
//  added Loops.C
//
//  Revision 1.11  2006/04/28 06:41:12  jochen
//  adding Gauss_Points.C
//
//  Revision 1.10  2006/04/28 06:07:43  jochen
//  adding Mapping.C
//
//  Revision 1.9  2006/04/26 12:35:27  jochen
//  easying the initialization of elements!
//
//  Revision 1.8  2006/04/21 14:01:31  jochen
//  correctimg the computation for complex problems
//
//  Revision 1.7  2006/02/16 17:02:29  jochen
//  Additional implementations for the solutiion of the dipol problem!
//
//  Revision 1.6  2006/02/13 08:32:15  jochen
//  Fixes for gcc version 4.0.2
//  introducing compile time data structures
//
//  Revision 1.5  2006/02/13 06:58:23  jochen
//  Update of many details, compile time vector and matrix data structures,
//  solved some compiling errors occuring via gcc 4.0.2
//
//  Revision 1.4  2005/08/25 12:33:49  jochen
//  TRAFO eingefügt
//
//  Revision 1.3  2005/08/19 07:52:20  jochen
//  Update und Verbesserungen
//
//  Revision 1.2  2005/06/27 19:16:32  jochen
//  All Compile-Time-Transfer
//
//  Revision 1.1  2005/05/19 17:14:04  jochen
//  Initial version
//
//
//==============================================================================

