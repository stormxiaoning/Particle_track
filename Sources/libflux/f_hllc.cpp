/**
 * @file f_hllc.cpp
 * @author Olivier Delestre <olivierdelestre41@yahoo.fr> (2008)
 * @author Christian Laguerre <christian.laguerre@math.cnrs.fr> (2012-2015)
 * @version 1.07.01
 * @date 2016-12-22
 *
 * @brief HLLC flux
 * @details 
 * Numerical flux: Harten, Lax, van Leer formulation with restoration
 * of the Contact Surface.
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

#include "f_hllc.hpp"

F_HLLC::F_HLLC(){
}

void F_HLLC::calcul(SCALAR h_L,SCALAR u_L,SCALAR v_L,SCALAR h_R,SCALAR u_R,SCALAR v_R){
    
  /**
   * @details
   * The HLLC approximate Riemann solver is a modification of the basic HLL scheme to account for the contact and shear waves (see \cite Toro01).\n
   * If the water heights on the two sides are small or \f$c_1 \approx c_2 \approx 0\f$, there is no water.
   * Else, HLL formulation is used (see \cite Bouchut04):
   * \f[{\cal F}(U_L,U_R)=\left\{\begin{array}{lll} F(U_L) & \mbox{if}& 0<c_{1} (\leq c_2), \\[0.2cm]
   * \displaystyle\frac{c_{2}F(U_L)-c_{1}F(U_R)}{c_{2}-c_{1}}+\frac{c_{1}c_{2}}{c_{2}-c_{1}}(U_R-U_L)&\mbox{if}& c_{1}<0<c_{2}, \\[0.4cm]
   * F(U_R)&\mbox{if}& (c_1 \leq ) c_{2}<0, \end{array}\right.\f]
   * with \f[c_{1}={\inf\limits_{U=U_L,U_R}}\left({\inf\limits_{j\in\{1,2\}}}|\lambda_{j}(U)|\right) \mbox{ and }
   * c_{2}={\sup\limits_{U=U_L,U_R}}\left({\sup\limits_{j\in\{1,2\}}}|\lambda_{j}(U)|\right),\f]
   * where \f$\lambda_1(U)=u-\sqrt{gh}\f$ and \f$\lambda_2(U)=u+\sqrt{gh}\f$ are the eigenvalues of the Shallow Water system,
   * \f$U = {}^t (h,hu,hv)\f$ and \f$ F(U) = {}^t (hu, hu^2+gh^2/2, hv^2) \f$.
   * \n If we consider the approximate flux HLL
   * \f[F_{i+\frac{1}{2}}=\left( \begin{array}{lll} F^{1}_{i+\frac{1}{2}}  \\[0.2cm]
   * F^{2}_{i+\frac{1}{2}}\\[0.2cm]
   * F^{3}_{i+\frac{1}{2}} \end{array}\right) \f]

   * \n then to obtain the HLLC solver just add the following expression for the third component
   * \f[F^{3}_{i+\frac{1}{2}}=\left\{\begin{array}{lll} F^{1}_{i+\frac{1}{2}}*V_L & \mbox{if}& 0 \leq u_{*}, \\[0.2cm]
   * F^{1}_{i+\frac{1}{2}}*V_R&\mbox{if}&  u_{*}<0, \end{array}\right.\f]
   * Where \f[\displaystyle u_{*}=\frac{c_{1}h_{R}(u_{R}-c_{2})-c_{2}h_{L}(u_{L}-c_{1})}{h_{R}(u_{R}-c_{2})-h_{L}(u_{L}-c_{1})}\f]

   * @param[in] h_L water height at the left of the interface where the flux is calculated.
   * @param[in] u_L velocity (in the x direction) at the left of the interface where the flux is calculated.
   * @param[in] v_L velocity (in the y direction) at the left of the interface where the flux is calculated.
   * @param[in] h_R water height at the right of the interface where the flux is calculated.
   * @param[in] u_R velocity (in the x direction) at the right of the interface where the flux is calculated.
   * @param[in] v_R velocity (in the y direction) at the right of the interface where the flux is calculated.
   * @par Modifies
   * Flux#f1 first component of the numerical flux.\n
   * Flux#f2 second component of the numerical flux.\n
   * Flux#f3 third component of the numerical flux.\n
   * Flux#cfl value of the CFL.
   * @note Long double are used locally in the computation to avoid numerical approximations.
   */
  
  if (h_L<HE_CA && h_R<HE_CA){
    f1 = 0.;
    f2 = 0.;
    f3 = 0.;
    cfl = 0.;
  }
  else{
    SCALAR grav_h_L = GRAV*h_L;
    SCALAR grav_h_R = GRAV*h_R;
    SCALAR q_R = u_R*h_R;
    SCALAR q_L = u_L*h_L;
    SCALAR c1, c2;
    if(h_L < HE_CA) {
      c1 = u_R - 2*sqrt(GRAV*h_R);
    }else{
      c1 = min(u_L-sqrt(grav_h_L),u_R-sqrt(grav_h_R)); // as u-sqrt(grav_h) <= u+sqrt(grav_h)    
    }
    if(h_R < HE_CA) {
       c2 = u_L + 2*sqrt(GRAV*h_L);
    }else{
       c2 = max(u_L+sqrt(grav_h_L),u_R+sqrt(grav_h_R)); // as u+sqrt(grav_h) >= u-sqrt(grav_h)
    }

    //cfl is the velocity to calculate the real cfl=max(fabs(c1),fabs(c2))*tx with tx=dt/dx
    if (fabs(c1)<EPSILON && fabs(c2)<EPSILON){              //dry state
      f1=0.;
      f2=0.;
      f3=0.;
      cfl=0.; //max(fabs(c1),fabs(c2))=0
    }
    else if (c1>=EPSILON){ //supercritical flow, from left to right : we have max(abs(c1),abs(c2))=c2>0
      f1=q_L;
      f2=q_L*u_L+GRAV*h_L*h_L*0.5;
      f3=q_L*v_L;
      cfl=c2; //max(fabs(c1),fabs(c2))=c2>0
    }
    else if (c2<=-EPSILON){ //supercritical flow, from right to left : we have max(abs(c1),abs(c2))=-c1>0
      f1=q_R;
      f2=q_R*u_R+GRAV*h_R*h_R*0.5;
      f3=q_R*v_R;
      cfl=fabs(c1); //max(fabs(c1),fabs(c2))=fabs(c1)
    }
    else{ //subcritical flow
      SCALAR c_star = (c1*h_R *(u_R - c2) - c2*h_L *(u_L - c1))/(h_R *(u_R - c2) - h_L *(u_L - c1));
      SCALAR tmp = 1./(c2-c1);
      f1=(c2*q_L-c1*q_R)*tmp+c1*c2*(h_R-h_L)*tmp;
      // long double are used locally to avoid numerical approximations
      f2=(SCALAR) ((long double)c2*((long double)q_L*(long double)u_L+(long double)GRAV*(long double)h_L*(long double)h_L*0.5)-c1*((long double)q_R*(long double)u_R+(long double)GRAV*(long double)h_R*(long double)h_R*0.5))*(long double)tmp+(long double)c1*(long double)c2*((long double)q_R-(long double)q_L)*(long double)tmp;
      if(c_star > EPSILON) {
	f3=f1*v_L;
      }else{
	f3=f1*v_R;
      }
      cfl=max(fabs(c1),fabs(c2));
    }
  }
}

F_HLLC::~F_HLLC(){
}

