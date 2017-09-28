/**
 * @file boundary_condition.cpp
 * @author Olivier Delestre <olivierdelestre41@yahoo.fr> (2008)
 * @author Christian Laguerre <christian.laguerre@math.cnrs.fr> (2012-2015)
 * @version 1.07.01
 * @date 2017-05-09
 *
 * @brief Boundary condition
 * @details 
 * Common part for all the boundary conditions.
 *
 * @copyright License Cecill-V2 \n
 * <http://www.cecill.info/licences/Licence_CeCILL_V2-en.html>
 *
 * (c) CNRS - Universite d'Orleans - BRGM (France)
 */
/* 
 *
 * This file is part of FullSWOF_2D software. 
 * <https://sourcesup.renater.fr/projects/fullswof-2d/> 
 *
 * FullSWOF_2D = Full Shallow-Water equations for Overland Flow, 
 * in two dimensions of space.
 * This software is a computer program whose purpose is to compute
 * solutions for 2D Shallow-Water equations.
 *
 * LICENSE
 *
 * This software is governed by the CeCILL license under French law and
 * abiding by the rules of distribution of free software. You can use,
 * modify and/ or redistribute the software under the terms of the CeCILL
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * <http://www.cecill.info>.
 *
 * As a counterpart to the access to the source code and rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty and the software's author, the holder of the
 * economic rights, and the successive licensors have only limited
 * liability.
 *
 * In this respect, the user's attention is drawn to the risks associated
 * with loading, using, modifying and/or developing or reproducing the
 * software by the user in light of its specific status of free software,
 * that may mean that it is complicated to manipulate, and that also
 * therefore means that it is reserved for developers and experienced
 * professionals having in-depth computer knowledge. Users are therefore
 * encouraged to load and test the software's suitability as regards their
 * requirements in conditions enabling the security of their systems and/or
 * data to be ensured and, more generally, to use and operate it in the
 * same conditions as regards security.
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL license and that you accept its terms.
 *
 ******************************************************************************/


#include "boundary_condition.hpp"

Boundary_condition::Boundary_condition(Parameters & par):NXCELL(par.get_Nxcell()),NYCELL(par.get_Nycell()){
  
  /**
   * @details   
   * Defines the number of cells. 
   * @param[in] par parameter, contains all the values from the parameters file.
   */

  Choice_Lbound= par.get_Lbound();
  Choice_Rbound= par.get_Rbound();
  Choice_Bbound= par.get_Bbound();
  Choice_Tbound= par.get_Tbound();
  
}

SCALAR Boundary_condition::get_hbound() const{
   
  /**
   * @details
   * @return Boundary_condition#hbound water height on the fictive cell.
   */
  
  return hbound;
}

SCALAR Boundary_condition::get_unormbound() const{
  
  /**
   * @details
   * @return Boundary_condition#unormbound normal velocity on the fictive cell.
   */
  
  return unormbound;
}

SCALAR Boundary_condition::get_utanbound() const{
  
  /**
   * @details
   * @return Boundary_condition#utanbound tangential velocity on the fictive cell.
   */
  
  return utanbound;
}


Boundary_condition::~Boundary_condition() {
}
