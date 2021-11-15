/*!
 * \file transport_model.hpp
 * \brief Headers of the main transport properties subroutines of the SU2 solvers.
 * \author S. Vitale, M. Pini, G. Gori, A. Guardone, P. Colonna
 * \version 6.0.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef TRANSPORT_MODEL_HPP_
#define TRANSPORT_MODEL_HPP_
#endif /* TRANSPORT_MODEL_HPP_ */
#pragma once

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <string>
#include <cmath>

#define LEN_COMPONENTS 32

#include "stdio.h"
#include "math.h"

#include "../../Common/include/datatype_structure.hpp"

using namespace std;


/*!
 * \class CViscosityModel
 * \brief Main class for defining the Transport-Physical Model
 * a child class for each particular Model (Power law, Sutherland, Chung, etc.)
 * \author S.Vitale, M.Pini
 * \version 1.0
 */
class CViscosityModel {
protected:
su2double      Mu,      /*!< \brief Dynamic viscosity. */
       dmudrho_T,   /*!< \brief DmuDrho_T. */
       dmudT_rho;   /*!< \brief DmuDT_rho. */
public:

    /*!
     * \brief Constructor of the class.
     */
    CViscosityModel(void);

    /*!
     * \brief Destructor of the class.
     */
    virtual ~CViscosityModel(void);

    /*!
     * \brief return viscosity value.
     */
    su2double GetViscosity(void);

    /*!
     * \brief return viscosity partial derivative value.
     */
    su2double Getdmudrho_T(void);

    /*!
     * \brief return viscosity partial derivative value.
     */
    su2double GetdmudT_rho(void);

    /*!
     * \brief Set Viscosity.
     */
    virtual   void SetViscosity(su2double T, su2double rho);

    /*!
     * \brief Set Viscosity Derivatives.
     */
    virtual   void SetDerViscosity(su2double T, su2double rho);

};


/*!
 * \class CConstantViscosity
 * \brief this class defines a constant viscosity
 * a child class for each particular Model (Power law, Sutherland, Chung, etc.)
 * \author S.Vitale, M.Pini
 * \version 1.0
 */
class CConstantViscosity : public CViscosityModel {
  
public:
  
  /*!
   * \brief Constructor of the class.
   */
  CConstantViscosity(void);
  
  /*!
   * \brief Constructor of the class.
   */
  CConstantViscosity(su2double mu_const);
  
  /*!
   * \brief Destructor of the class.
   */
  virtual ~CConstantViscosity(void);
  
  
};


/*!
 * \class CSutherland
 * \brief this class defines a constant viscosity
 * a child class for each particular Model (Power law, Sutherland, Chung, etc.)
 * \author S.Vitale, M.Pini
 * \version 1.0
 */
class CSutherland : public CViscosityModel {
protected:
  su2double      Mu_ref,    /*!< \brief Internal Energy. */
  T_ref,     /*!< \brief DpDd_e. */
  S;       /*!< \brief DpDe_d. */
  
public:
  
  /*!
   * \brief Constructor of the class.
   */
  CSutherland(void);
  
  /*!
   * \brief Constructor of the class.
   */
  CSutherland(su2double mu_ref, su2double t_ref, su2double s);
  
  /*!
   * \brief Destructor of the class.
   */
  virtual ~CSutherland(void);
  
  /*!
   * \brief Set Viscosity.
   */
  void SetViscosity(su2double T, su2double rho);
  
  /*!
   * \brief Set Viscosity Derivatives.
   */
  void SetDerViscosity(su2double T, su2double rho);
  
};


/*!
 * \class CThermalConductivityModel
 * \brief Main class for defining the Transport-Physical Model
 * a child class for each particular Model (Prandtl-based, etc.)
 * \author S. Vitale, M. Pini
 * \version 1.0
 */
class CConductivityModel {
protected:
su2double      Kt,      /*!< \brief Thermal conductivity. */
       dktdrho_T,   /*!< \brief DktDrho_T. */
       dktdT_rho;   /*!< \brief DktDT_rho. */
public:

    /*!
     * \brief Constructor of the class.
     */
    CConductivityModel(void);

    /*!
     * \brief Destructor of the class.
     */
    virtual ~CConductivityModel(void);

    /*!
     * \brief return viscosity value.
     */
    su2double GetConductivity(void);

    /*!
     * \brief return viscosity partial derivative value.
     */
    su2double Getdktdrho_T(void);

    /*!
     * \brief return viscosity partial derivative value.
     */
    su2double GetdktdT_rho(void);

    /*!
     * \brief Set Thermal conductivity.
     */
    virtual   void SetConductivity(su2double T, su2double rho, su2double mu, su2double cp);

    /*!
     * \brief Set Thermal conductivity derivatives.
     */
    virtual   void SetDerConductivity(su2double T, su2double rho, su2double dmudrho_T, su2double dmudT_rho, su2double cp);

};


/*!
 * \class CConstantPrandtl
 * \brief this class defines a constant thermal conductivity using a constant Prandtl's number
 * \author S.Vitale, M.Pini
 * \version 1.0
 */
class CConstantConductivity : public CConductivityModel {

public:

    /*!
     * \brief Constructor of the class.
     */
      CConstantConductivity(void);

    /*!
     * \brief Constructor of the class.
     */
      CConstantConductivity(su2double kt_const);

    /*!
     * \brief Destructor of the class.
     */
    virtual ~CConstantConductivity(void);

};


/*!
 * \class CConstantPrandtl
 * \brief this class defines a non-constant thermal conductivity using a constant Prandtl's number
 * \author S.Vitale, M.Pini
 * \version 1.0
 */
class CConstantPrandtl : public CConductivityModel {
protected:
  su2double      Pr_const;    /*!< \brief Prandtl's number. */

public:

    /*!
     * \brief Constructor of the class.
     */
      CConstantPrandtl(void);

    /*!
     * \brief Destructor of the class.
     */
    virtual ~CConstantPrandtl(void);

    /*!
     * \brief Constructor of the class.
     */
      CConstantPrandtl(su2double pr_const);

    /*!
     * \brief Set Thermal conductivity.
     * \brief par1 -> Cp.
     * \brief par2 -> Mu.
     */
    void SetConductivity(su2double T, su2double rho, su2double mu, su2double cp);

    /*!
     * \brief Set Thermal conductivity derivatives.
     */
    void SetDerConductivity(su2double T, su2double rho, su2double dmudrho_T, su2double dmudT_rho, su2double cp);

};


#include "transport_model.inl"
