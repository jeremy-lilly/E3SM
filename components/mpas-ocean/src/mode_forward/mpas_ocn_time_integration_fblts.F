! Copyright (c) 2021,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_integration_fblts
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
!  ocn_lts_thick_tend
!
!> \brief  Calculate thickness tendencies for FB_LTS
!> \author Giacomo Capodaglio
!> \date   October 2023 
!> \details
!>  This routine calculates the tendencies on different LTS regions
!
!-----------------------------------------------------------------------

   subroutine ocn_lts_thick_tend(statePool, LTSPool, tendPool, timeLevelIn, computeOnFineBig, &
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
         statePool  !< Input: state variables

      type (mpas_pool_type), intent(in) :: &
         LTSPool  !< Input: LTS data

      !-----------------------------------------------------------------
      ! Input/output variables
      !-----------------------------------------------------------------
      
      type (mpas_pool_type), intent(inout) :: &
         tendPool  !< Input/output: tendency variables
    
      !-----------------------------------------------------------------
      ! Local variables
      !-----------------------------------------------------------------
      
      integer, dimension(:,:), pointer :: &
         nCellsInLTSRegion
      
      integer, dimension(:,:,:), pointer :: &
         cellsInLTSRegion
      
      real (kind=RKIND), dimension(:,:), pointer :: &
         layerThicknessTend, &
         normalVelocity, & 
         layerThickness

      integer :: &
         iEdge, cell1, cell2, k, ie, iRegion, nRegions, &
         ic, i, iCell, kmin, kmax
      real (kind=RKIND) :: &
         invdcEdge, flux 
      
      !-----------------------------------------------------------------
      ! Begin routine
      !-----------------------------------------------------------------
      
      nRegions = 2

      call mpas_pool_get_array(LTSPool, 'cellsInLTSRegion', cellsInLTSRegion)
      call mpas_pool_get_array(LTSPool, 'nCellsInLTSRegion', nCellsInLTSRegion)

      call mpas_pool_get_array(tendPool, 'layerThickness', layerThicknessTend)

      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, timeLevelIn)
      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevelIn)

      ! calculate layerThickEdgeFlux
      do iEdge = 1, nEdgesAll
         kmin = minLevelEdgeBot(iEdge)
         kmax = maxLevelEdgeTop(iEdge)
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         do k = 1,nVertLevels
            ! initialize layerThicknessEdgeMean to avoid divide by
            ! zero and NaN problems.
            layerThickEdgeFlux(k,iEdge) = -1.0e34_RKIND
         end do
         do k = kmin,kmax
            ! central differenced
            layerThickEdgeFlux(k,iEdge) = 0.5_RKIND * &
                                          ( layerThickness(k,cell1) &
                                          + layerThickness(k,cell2) )
         end do
      end do
   
      layerThicknessTend(:, :) = 0.0_RKIND

      if (computeOnInterface == 1) then
         do iRegion = 1, nRegions
            ! thickness tendency
            do ic = 1, nCellsInLTSRegion(iRegion,2)
               iCell = cellsInLTSRegion(iRegion,2,ic)
               do i = 1, nEdgesOnCell(iCell)
                  iEdge = edgesOnCell(i, iCell)
                  do k = minLevelEdgeBot(iEdge), maxLevelEdgeTop(iEdge)
                     flux = normalVelocity(k, iEdge) * dvEdge(iEdge) * layerThickEdgeFlux(k, iEdge)
                     layerThicknessTend(k, iCell) = layerThicknessTend(k, iCell) & 
                                                    + edgeSignOnCell(i, iCell) * flux * invAreaCell(iCell)
                  end do
               end do
            end do
         end do
      end if

      if (computeOnFineBig == 1) then
         ! thickness tendency
         do ic = 1, nCellsInLTSRegion(1,1)
            iCell = cellsInLTSRegion(1,1,ic)
            do i = 1, nEdgesOnCell(iCell)
               iEdge = edgesOnCell(i, iCell)
               do k = minLevelEdgeBot(iEdge), maxLevelEdgeTop(iEdge)
                  flux = normalVelocity(k, iEdge) * dvEdge(iEdge) * layerThickEdgeFlux(k, iEdge)
                  layerThicknessTend(k, iCell) = layerThicknessTend(k,iCell) &
                                                 + edgeSignOnCell(i, iCell) * flux * invAreaCell(iCell)
               end do
            end do
         end do
      end if

      if (computeOnFineSmall == 1) then
         ! thickness tendency
         do ic = 1, nCellsInLTSRegion(1,3)
            iCell = cellsInLTSRegion(1,3,ic)
            do i = 1, nEdgesOnCell(iCell)
               iEdge = edgesOnCell(i, iCell)
               do k = minLevelEdgeBot(iEdge), maxLevelEdgeTop(iEdge)
                  flux = normalVelocity(k, iEdge) * dvEdge(iEdge) * layerThickEdgeFlux(k, iEdge)
                  layerThicknessTend(k, iCell) = layerThicknessTend(k,iCell) &
                                                 + edgeSignOnCell(i, iCell) * flux * invAreaCell(iCell)
               end do
            end do
         end do
      end if

      if (computeOnCoarse == 1) then
         ! thickness tendency
         do ic = 1, nCellsInLTSRegion(2,1)
            iCell = cellsInLTSRegion(2,1,ic)
            do i = 1, nEdgesOnCell(iCell)
               iEdge = edgesOnCell(i, iCell)
               do k = minLevelEdgeBot(iEdge), maxLevelEdgeTop(iEdge)
                  flux = normalVelocity(k, iEdge) * dvEdge(iEdge) * layerThickEdgeFlux(k, iEdge)
                  layerThicknessTend(k, iCell) = layerThicknessTend(k,iCell) &
                                                 + edgeSignOnCell(i, iCell) * flux * invAreaCell(iCell)
               end do
            end do
         end do
      end if

   end subroutine ocn_lts_thick_tend!}}}


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_lts_vel_tend
!
!> \brief  Calculate velocity tendencies for FB_LTS
!> \author Giacomo Capodaglio
!> \date   October 2023 
!> \details
!>  This routine calculates the tendencies on different LTS regions
!
!-----------------------------------------------------------------------

   subroutine ocn_lts_vel_tend(statePool, LTSPool, tendPool, timeLevelIn, computeOnFineBig, &
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
         statePool  !< Input: state variables

      type (mpas_pool_type), intent(in) :: &
         LTSPool  !< Input: LTS data

      !-----------------------------------------------------------------
      ! Input/output variables
      !-----------------------------------------------------------------
      
      type (mpas_pool_type), intent(inout) :: &
         tendPool  !< Input/output: tendency variables
    
      !-----------------------------------------------------------------
      ! Local variables
      !-----------------------------------------------------------------
      
      integer, dimension(:,:), pointer :: &
         nEdgesInLTSRegion
      
      integer, dimension(:,:,:), pointer :: &
         edgesInLTSRegion
      
      real (kind=RKIND), dimension(:), pointer :: &
         ssh
      
      real (kind=RKIND), dimension(:,:), pointer :: &
         normalVelocityTend, &
         normalVelocity, & 
         layerThickness

      integer :: &
         iEdge, cell1, cell2, k, ie, iRegion, nRegions, &
         ic, i, iCell, kmin, kmax
      real (kind=RKIND) :: &
         invdcEdge, betaSelfAttrLoad, flux, ssh_sal_on
      
      !-----------------------------------------------------------------
      ! Begin routine
      !-----------------------------------------------------------------
      
      betaSelfAttrLoad = config_self_attraction_and_loading_beta
      nRegions = 2

      call mpas_pool_get_array(LTSPool, 'edgesInLTSRegion', edgesInLTSRegion)
      call mpas_pool_get_array(LTSPool, 'nEdgesInLTSRegion', nEdgesInLTSRegion)

      call mpas_pool_get_array(tendPool, 'normalVelocity', normalVelocityTend)

      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, timeLevelIn)
      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevelIn)
      call mpas_pool_get_array(statePool, 'ssh', ssh, timeLevelIn)

      if (config_use_self_attraction_loading) then
         ssh_sal_on = 1.0_RKIND
      else
         ssh_sal_on = 0.0_RKIND
      endif
      
      ! ssh
      do iCell = 1, nCellsAll
         k = maxLevelCell(iCell)
         zTop(k:nVertLevels,iCell) = -bottomDepth(iCell) + layerThickness(k,iCell)
         do k = maxLevelCell(iCell)-1, minLevelCell(iCell), -1
            zTop(k,iCell) = zTop(k+1,iCell) + layerThickness(k  ,iCell)
         end do
         ! copy zTop(1,iCell) into sea-surface height array
         ssh(iCell) = zTop(minLevelCell(iCell),iCell)
      end do

      normalVelocityTend(:,:) = 0.0_RKIND

      ! interface
      if (computeOnInterface == 1) then
         do iRegion = 1, nRegions
            ! velocity tendency
            do ie = 1, nEdgesInLTSRegion(iRegion,2)
               iEdge = edgesInLTSRegion(iRegion,2,ie)
               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)
               invdcEdge = 1.0_RKIND / dcEdge(iEdge)
               kMin = minLevelEdgeBot(iEdge)
               kMax = maxLevelEdgeTop(iEdge)
               do k=kMin,kMax
                  normalVelocityTend(k,iEdge) = normalVelocityTend(k,iEdge) &
                                                - edgeMask(k,iEdge) * invdcEdge &
                                                * ( gravity * ( (ssh(cell2) - ssh(cell1)) &
                                                - (1.0_RKIND - ssh_sal_on) * betaSelfAttrLoad &
                                                * (ssh(cell2) - ssh(cell1)) ) )
               end do
            end do
         end do
      end if

      if (computeOnFineBig == 1) then
         ! velocity tendency
         do ie = 1, nEdgesInLTSRegion(1,1)
            iEdge = edgesInLTSRegion(1,1,ie)
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            invdcEdge = 1.0_RKIND / dcEdge(iEdge)
            kMin = minLevelEdgeBot(iEdge)
            kMax = maxLevelEdgeTop(iEdge)
            do k=kMin,kMax
               normalVelocityTend(k,iEdge) = normalVelocityTend(k,iEdge) &
                                             - edgeMask(k,iEdge) * invdcEdge &
                                             * ( gravity * ( (ssh(cell2) - ssh(cell1)) &
                                             - (1.0_RKIND - ssh_sal_on) * betaSelfAttrLoad &
                                             * (ssh(cell2) - ssh(cell1)) ) )
            end do
         end do
      end if

      if (computeOnFineSmall == 1) then
         ! velocity tendency
         do ie = 1, nEdgesInLTSRegion(1,3)
            iEdge = edgesInLTSRegion(1,3,ie)
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            invdcEdge = 1.0_RKIND / dcEdge(iEdge)
            kMin = minLevelEdgeBot(iEdge)
            kMax = maxLevelEdgeTop(iEdge)
            do k=kMin,kMax
               normalVelocityTend(k,iEdge) = normalVelocityTend(k,iEdge) &
                                             - edgeMask(k,iEdge) * invdcEdge &
                                             * ( gravity * ( (ssh(cell2) - ssh(cell1)) &
                                             - (1.0_RKIND - ssh_sal_on) * betaSelfAttrLoad &
                                             * (ssh(cell2) - ssh(cell1)) ) )
            end do
         end do
      end if

      if (computeOnCoarse == 1) then
         ! velocity tendency
         do ie = 1, nEdgesInLTSRegion(2,1)
            iEdge = edgesInLTSRegion(2,1,ie)
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            invdcEdge = 1.0_RKIND / dcEdge(iEdge)
            kMin = minLevelEdgeBot(iEdge)
            kMax = maxLevelEdgeTop(iEdge)
            do k=kMin,kMax
               normalVelocityTend(k,iEdge) = normalVelocityTend(k,iEdge) &
                                             - edgeMask(k,iEdge) * invdcEdge &
                                             * ( gravity * ( (ssh(cell2) - ssh(cell1)) &
                                             - (1.0_RKIND - ssh_sal_on) * betaSelfAttrLoad &
                                             * (ssh(cell2) - ssh(cell1)) ) )
            end do
         end do
      end if

   end subroutine ocn_lts_vel_tend!}}}


end module ocn_time_integration_fblts

