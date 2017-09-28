/**
 * @file minmod.cpp
 * @author Olivier Delestre <olivierdelestre41@yahoo.fr> (2008)
 * @author Christian Laguerre <christian.laguerre@math.cnrs.fr> (2012-2015)
 * @version 1.06.00
 * @date 2015-02-19
 *
 * @brief %Minmod limiter
 * @details 
 * Slope limiter: minmod.
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

#include "minmod.hpp"

Minmod::Minmod(){
}

void Minmod::calcul(SCALAR a, SCALAR b){
  
  /**
   * @details 
   * %Minmod function:
   * \f[ \mbox{minmod}(x,y)=\left\{\begin{array}{ll}
   *    \min(x,y)&\mbox{if }x,y\geq{0},\\
   *    \max(x,y)&\mbox{if }x,y\leq{0},\\
   *    0&\mbox{else}.    \end{array}\right. \f]
   * @param[in] a slope on the left of the cell.
   * @param[in] b slope on the right of the cell.
   * @par Modifies
   * Limiter#rec recontructed value.
   */
  
  if (a>=0. && b>=0.){
    rec=min(a,b);
  }
  else if (a<=0. && b<=0.){
    rec=max(a,b);
  }
  else{
    rec=0.;
  }
}

Minmod::~Minmod(){
}

