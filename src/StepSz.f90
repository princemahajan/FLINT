!############################################################################################
!
! Copyright 2020 Bharat Mahajan
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!    http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!
!> \brief       StepSize Module
!! \details     This module provides procedures for step size computation.
!! \author      Bharat Mahajan
!! \date        Created: 01/28/2019    
!
!############################################################################################
    
!> Module containing all the algorithms for step size computations    
module StepSize

    use FLINT_base, only: WP, DiffEqSys
    
    ! Default constants
    
    !> Safety factor
    real(WP), parameter :: SF = 0.9_WP
    
    !> Minimum safety factor, so we dont decrease the step size too fast
    real(WP), parameter :: SFMIN = 0.333_WP
    
    !> Maximum safety factor, so we dont increase the step size too fast    
    real(WP), parameter :: SFMAX = 6.0_WP
    
    !> beta for Lund stabilization
    !> Positive values (< 0.04) make the step size control more stable
    real(WP), parameter :: BETA = 0.0_WP

    !> Also for Lund stabilization (see Hairer's code)
    real(WP), parameter :: hbyhoptOLD = 1.0e-4_WP
    
    

    contains
    
    
    !> \brief Subroutine for computing the starting step size using Hairer's algorithm.
    !! \details For details on the algorithm, see Hairer's book "Solving ODE I" or his code dop853.f
    !! at http://www.unige.ch/~hairer/prog/nonstiff/dop853.f.
    subroutine StepSz0Hairer(h, p, n, X0, Y0, F0, Sc, hSign, hMax, FCalls, pDiffEqSys, Params)
    
        implicit none
        
        intrinsic :: sqrt, max, min
    
        real(WP), intent(out)                 :: h     !< Computed step size
        integer, intent(in)                   :: p     !< Order of the method, which is used to advance the integration.
                                                       !! For example, for Dormand-Prince 8(7) ERK method, p=8.
        integer, intent(in)                   :: n     !< dimension of the DiffEq system 
        real(WP), intent(in)                  :: X0    !< Initial value of the independent variable         
        real(WP), dimension(n), intent(in)    :: Y0    !< Initial condition
        real(WP), dimension(n), intent(in)    :: F0    !< DiffEq Function evaluated at Y0
        real(WP), dimension(n), intent(in)    :: Sc    !< Error Scale factor
        real(WP), intent(in)                  :: hSign !< Must be either +1.0 or -1.0 depending 
                                                       !!on the forward or backward integration, respectively
        real(WP), intent(in)                  :: hMax  !< Maximum allowed step size 
        integer, intent(out)                  :: FCalls !< Number of Calls that this routine makes to DiffEq function "F"
        class(DiffEqSys), pointer, intent(in) :: pDiffEqSys
        real(WP), dimension(:), intent(in), optional :: Params !< Real parameter array to be passed to DEFunc
        
        real(WP) :: d0, d1, d2, dMax, h0, h1
        real(WP), dimension(n) :: F1
        
        ! Algorithm starts here
        
        d0 = norm2(Y0/Sc)
        d1 = norm2(F0/Sc)
    
        ! first guess
        if (d0 <= 1.0e-5 .or. d1 <= 1.0e-5) then
            h0 = 1.0e-6_WP
        else
            h0 = 0.01_WP*d0/d1
        end if
        
        h0 = hSign*min(h0, hMax)
        
        ! Explicit Euler step
        F1 = pDiffEqSys%F((X0 + h0), (Y0 + h0*F0), Params)
        
        ! 2nd derivative approximation
        d2 = norm2((F1 - F0)/Sc)/h0
        
        ! second guess based on 2nd drivative
        dMax = max(d1, abs(d2))
        
        !> \remark The equation to solve for the new guess for the starting step size is
        !! \f[ h_1^{p+1} max(d1, d2) = 0.01 \f]
        !! In Hairer's dop853 code, he actually uses \f$p\f$ instead of \f$p+1\f$
        !! as the power of \f$h_1\f$ for reasons unknown. He uses \f$p=8\f$ for his dop853.
        if (dMax <= 1.0e-15 ) then
            h1 = max(1.0e-6_WP, abs(h0)*1.0e-3_WP)
        else
            h1 = (0.01_WP/dMax)**(1.0/(p))
        end if
    
        ! choose the final value for the initial step size
        h = hSign*min(100.0_WP*abs(h0), h1, hMax)
        
        ! Total number of calls we made to F
        FCalls = 1
    
    end subroutine StepSz0Hairer

    
    
    !> \brief Pure function for computing the new step size using Hairer's algorithm.
    !! \details For details on the algorithm, see Hairer's book "Solving ODE I" or his code dop853.f
    !! at http://www.unige.ch/~hairer/prog/nonstiff/dop853.f.
    pure function StepSzHairer(h, IsStepRejected, q, e, StepSzParams)
 
        implicit none
        
        intrinsic :: sqrt, max, min
    
        real(WP), intent(in)                  :: h     !< current step size        
        logical, intent(in)                   :: IsStepRejected !< If true then the step size is decreased else increased.
        integer, intent(in)                   :: q     !< q = min(p,phat) but Hairer uses q=8 in DOP853.
        real(WP), intent(in)                  :: e     !< Error norm
        real(WP), dimension(5), intent(in) :: StepSzParams        
        
        real(WP) :: StepSzHairer
        
        real(WP) :: hbyhopt, hbyhoptOld, beta, SF, SFmin, SFmax
        
        ! Algorithm starts here    

        ! unpack the params
        SF = StepSzParams(1)   !< Safety factor
        SFmin = StepSzParams(2) !< Min safety factor so that we dont decrease the step size too fast
        SFmax = StepSzParams(3) !< Max safety factor so that we dont increase the step size too fast
        beta = StepSzParams(4) !< Lund stabilization parameter
        hbyhoptOld = StepSzParams(5) !< Needed for Lund stabilization                
        
        !> \remark Optimal step computed using the equations:
        !! \f[  err = C h^{(q+1)}  \f]
        !! \f[ 1 = C h_{opt}^{q+1} \f]
        hbyhopt = e**(1.0_WP/(q+1.0_WP) - beta*0.2_WP)
        
        if (IsStepRejected .EQV. .FALSE.) then
            
            ! The current step is accepted, so we can safely increase the step size
            
            ! \remark For Lund stabilization of the accepted step size, See Hairer's book.
            hbyhopt = hbyhopt/hbyhoptOld**beta
        
            !> \remark New step chosen such that
            !! \f[ SF_{min} <= SF\, h_{new}/h_{old} <= SF_{max}  \f]
            StepSzHairer = h*min(SFmax, max(SFmin, SF*1.0_WP/hbyhopt))
        
        else
            
            ! The current step is rejected, so we need to decrease the step size
            
            !> \remark New step chosen such that
            !! \f[ SF_{min} <= SF\, h_{new}/h_{old} <= SF_{max}  \f]
            StepSzHairer = h*max(SFmin, SF*1.0_WP/hbyhopt)
        
        end if
    
    end function StepSzHairer
    
    
    
end module StepSize
