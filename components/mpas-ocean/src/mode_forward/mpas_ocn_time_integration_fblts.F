! Copyright (c) 2021,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_integration_lts
!
!> \brief  MPAS barotropic ocean LTS Time integration scheme
!> \author Jeremy Lilly
!> \date   November 2022
!> \details
!>  This module contains the LTS init routine and the LTS
!>  barotropic ocean time integration scheme with splitting
!>  on the fast and slow tendency terms.
!
!-----------------------------------------------------------------------

module ocn_time_integration_fblts
       
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_dmpar
   use mpas_threading
   use mpas_vector_reconstruction
   use mpas_spline_interpolation
   use mpas_timer

   use ocn_constants
   use ocn_tendency
   use ocn_diagnostics
   use ocn_mesh
   use ocn_vmix
   use ocn_config
   use ocn_diagnostics_variables
   use ocn_equation_of_state
   use ocn_time_average_coupled
   use ocn_time_varying_forcing

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------
   public :: ocn_time_integrator_fblts, &
             ocn_time_integration_fblts_init

   contains


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_integrator_fblts
!
!> \brief  MPAS barotropic ocean FB_LTS time integration scheme
!> \author Jeremy Lilly
!> \date   October 2023
!> \details
!>  This routine integrates one timestep (dt) using an FB_LTS time
!>  integrator with a splitting of the fast and slow tendency terms
!
!-----------------------------------------------------------------------

   subroutine ocn_time_integrator_fblts(domain, dt)!{{{
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Advance model state forward in time by the specified time step using
   ! a local time stepping scheme with spltting of the fast and slow tendencies
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      !-----------------------------------------------------------------
      ! Input variables
      !-----------------------------------------------------------------

      real (kind=RKIND), intent(in) :: &
         dt  !< [in] time step (sec) to move forward

      !-----------------------------------------------------------------
      ! Input/output variables
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: &
         domain  !< [inout] model state to advance forward

   end subroutine ocn_time_integrator_fblts!}}}


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_integration_fblts_init
!
!> \brief  Initialization for MPAS ocean FB_LTS time integration scheme
!> \author Giacomo Capodaglio
!> \date   October 2023
!> \details
!>  This routine computes the LTS arrays
!
!-----------------------------------------------------------------------

   subroutine ocn_time_integration_fblts_init(domain)!{{{
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Sets the LTS arrays
   !
   ! Output: LTS arrays are written    
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      implicit none
      
      !-----------------------------------------------------------------
      ! Input/output variables
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: &
         domain  !< [inout] model state

   end subroutine ocn_time_integration_fblts_init!}}}


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_lts_tends
!
!> \brief  MPAS ocean LTS3 Time integration scheme
!> \author Giacomo Capodaglio
!> \date   October 2023 
!> \details
!>  This routine calculates the tendencies on different LTS regions
!
!-----------------------------------------------------------------------

   subroutine ocn_lts_tends(statePool, LTSPool, tendPool, timeLevelIn, computeOnFineBig, &
                            computeOnFineSmall, computeOnCoarse, computeOnInterface)!{{{
      
      !-----------------------------------------------------------------
      ! Input variables
      !-----------------------------------------------------------------
      
      integer, intent(in) :: &
         timeLevelIn, &
         computeOnFineBig, &
         computeOnFineSmall, &
         computeOnCoarse, &
         computeOnInterface

      type (mpas_pool_type), intent(in) :: &
         statePool          !< Input: state variables

      type (mpas_pool_type), intent(in) :: &
         LTSPool          !< Input: LTS data

      !-----------------------------------------------------------------
      ! Input/output variables
      !-----------------------------------------------------------------
      
      type (mpas_pool_type), intent(inout) :: &
         tendPool          !< Input: tendency variables
     
   end subroutine ocn_lts_tends!}}}


end module ocn_time_integration_fblts

