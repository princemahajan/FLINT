!############################################################################################
!
! Copyright 2021 Bharat Mahajan
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
!> \brief       Explicit Runge-Kutta Module
!! \details     A module for Explicit Runge-Kutta type numerical solvers with dense-output,
!!              vector events detection, and stiffness detection.
!! \author      Bharat Mahajan
!! \date        Created: 01/25/2019    
!
!############################################################################################

module ERK
    
    use FLINT_base
    use ButcherTableaus
    use StepSize
    use, intrinsic :: IEEE_ARITHMETIC
    
    implicit none

    private

    real(WP), parameter :: DEFAULT_ABSTOL = 1.0e-9_WP !< default absolute tolerance
    real(WP), parameter :: DEFAULT_RELTOL = 1.0e-6_WP !< default relative tolerance

    !> If the final point is within the range of this expansion factor times
    !! the current step size, the current step size is adjusted to attain the
    !! the final point in this last step. It is not used for constant step size.
    real(WP), parameter :: LASTSTEP_EXPFAC = 1.01_WP


    !> Maximum number of events that can be detected during event-checking inside the 
    !! integrator's internal step using the event step size
    integer, parameter :: MAXNUMSTEPSZEVENTS = 4
    
    real(WP), parameter :: DEFAULT_EVENTTOL = 1.0e-6_WP !< this value is used to detect events
    
    !> Stiffness test default number of steps. See Hairer's DOP853 code.
    integer, parameter :: STIFFTEST_STEPS = 1000
    
    !> All the supported ERK Methods
    enum, bind(C)
        enumerator :: ERK_DOP853 = 17       !< Hairer's DOP853
        enumerator :: ERK_DOP54             !< Dormand-Prince 5(4)
        enumerator :: ERK_VERNER98R         !< Verner's 9th-order Robust coefficients
        !enumerator :: ERK_VERNER98E        !< Verner's 9th-order Efficient coefficients
        !enumerator :: ERK_VERNER87R        !< Verner's 8th-order Robust coefficients
        !enumerator :: ERK_VERNER87E        !< Verner's 8th-order Efficient coefficients
        !enumerator :: ERK_VERNER76R        !< Verner's 7th-order Robust coefficients
        !enumerator :: ERK_VERNER76E        !< Verner's 7th-order Efficient coefficients
        enumerator :: ERK_VERNER65E         !< Verner's 6th-order Efficient coefficients
    end enum

    public :: ERK_DOP853, ERK_DOP54, ERK_VERNER98R, ERK_VERNER65E  !, ERK_VERNER98E, ERK_VERNER87R, ERK_VERNER87E 
    !public :: ERK_VERNER76R, ERK_VERNER76E
    
    !> The ERK class for all the ERK methods available to the user. It inherits FLINT_class.
    type, extends(FLINT_class), public :: ERK_class
    
        private
        
        !> Type of ERK method to use for integration, DOP853 by default
        integer(kind(ERK_DOP853)), public :: method = ERK_DOP853
        
        ! Specs of the chosen method, DOP853 by default
        integer :: p        = DOP853_p
        integer :: phat     = DOP853_phat
        integer :: q        = DOP853_q
        integer :: pstar    = DOP853_pstar
        
        ! Butcher Tableau coefficients of the chosen method
        integer :: s        !< Total number of stages
        integer :: sint     !< Number of Integration stages (including the FSAL stage)
        ! real(WP), dimension(ERK_MAX_A) :: a             !< a_ij
        ! real(WP), dimension(1:ERK_MAXSTAGES) :: b       !< b_i
        ! real(WP), dimension(2:ERK_MAXSTAGES) :: c       !< c_i
        ! real(WP), dimension(1:ERK_MAXSTAGES) :: e       !< e_i        

        ! k values for each each stage, used in integration step routines
        real(WP), dimension(:, :), allocatable :: k

        !> Stages that gets multiplied by non-zero dij, in other words the stages specified by dinz
        !! are the ones that contribute to the interpolating polynomial coefficients.
        !! Specify the stages in dinz and terminate the list by -1.
        integer, dimension(1:ERK_MAXSTAGES+1) :: dinz 

        real(WP), dimension(1:ERK_MAXSTAGES,1:ERK_MAXIPDEGREE) :: d !< d_ij coefficients for interpolation
        
        logical :: IsFSALMethod !< If true then the last stage is reused as the first for the new step
        
    contains
        
        procedure, public :: Init => erk_init
        procedure, public :: Integrate => erk_int
        procedure, public :: Interpolate => erk_interp
        procedure, public :: Info => erk_info
        
        procedure, private :: StepInt => erk_stepint
        procedure, private :: InterpCoeff => erk_intpcoeff

        final :: erk_destroy !< destructor
        
    end type ERK_class
    
    
    
    interface
    
    module subroutine erk_init(me, DE, MaxSteps, Method, ATol, RTol, InterpOn, &
        InterpStates, MinStepSize, MaxStepSize, StepSzParams, &
        EventsOn, EventStepSz, EventOptions, EventTol)
    
        implicit none
        
        class(ERK_class), intent(inout) :: me
        class(DiffEqSys), target :: DE
        integer, intent(in) :: MaxSteps        
        integer, intent(in), optional :: Method
        real(WP), dimension(:), intent(in), optional :: ATol
        real(WP), dimension(:), intent(in), optional :: RTol
        logical, intent(in), optional :: InterpOn
        integer, dimension(:), intent(in), optional :: InterpStates
        real(WP), intent(in), optional :: MinStepSize, MaxStepSize        
        real(WP), dimension(6), intent(in), optional :: StepSzParams
        logical, intent(in), optional :: EventsOn
        real(WP), dimension(:), intent(in), optional :: EventStepSz
        integer(kind(FLINT_EVENTOPTION_ROOTFINDING)), dimension(:), &
                                intent(in), optional :: EventOptions
        real(WP), dimension(:), intent(in), optional :: EventTol    

    end subroutine erk_init

    
    module subroutine erk_int(me, X0, Y0, Xf, Yf, StepSz, UseConstStepSz, StepAdvance, IntStepsOn, &
                                Xint, Yint, EventMask, EventStates, StiffTest, params)
    
        implicit none
        
        class(ERK_class), intent(inout) :: me
        
        real(WP), intent(in) :: X0
        real(WP), dimension(me%pDiffEqSys%n), intent(in) :: Y0
        real(WP), intent(inout) :: Xf            
        real(WP), dimension(size(Y0)), intent(out) :: Yf
        real(WP), intent(inout) :: StepSz
        logical, intent(in), optional :: UseConstStepSz 
        logical, intent(in), optional :: StepAdvance
        logical, intent(in), optional :: IntStepsOn
        real(WP), allocatable, dimension(:), intent(out), optional :: Xint            
        real(WP), allocatable, dimension(:,:), intent(out), optional :: Yint
        logical, dimension(me%pDiffEqSys%m), intent(in), optional :: EventMask
        real(WP), allocatable, dimension(:,:), intent(out), optional :: EventStates
        integer, intent(inout), optional :: StiffTest        
        real(WP), dimension(:), intent(in), optional :: params        
        
    end subroutine erk_int

    
    
    module subroutine erk_interp(me, Xarr, Yarr, SaveInterpCoeffs)
    
        implicit none
        
        class(ERK_class), intent(inout) :: me
        
        real(WP), dimension(1:), intent(in) :: Xarr
        real(WP), dimension(size(me%InterpStates),size(Xarr)), intent(out) :: Yarr
        logical, intent(in), optional :: SaveInterpCoeffs
        
        
    end subroutine erk_interp
    
    
    module subroutine erk_info(me, LastStatus, StatusMsg, nSteps, nAccept, nReject, nFCalls, &
                                InterpReady, h0, X0, Y0, hf, Xf, Yf)

            class(ERK_class), intent(inout) :: me
            integer(kind(FLINT_SUCCESS)), intent(out) :: LastStatus
            character(len=:), allocatable, intent(out), optional :: StatusMsg
            integer, intent(out), optional :: nSteps
            integer, intent(out), optional :: nAccept
            integer, intent(out), optional :: nReject
            integer, intent(out), optional :: nFCalls
            logical, intent(out), optional :: InterpReady
            real(WP), intent(out), optional :: h0
            real(WP), intent(out), optional :: X0
            real(WP), dimension(:), intent(out), optional :: Y0
            real(WP), intent(out), optional :: hf
            real(WP), intent(out), optional :: Xf
            real(WP), dimension(:), intent(out), optional :: Yf
            
    end subroutine erk_info
    
    
    pure module function InterpY(np, X, X0, h, pstar, bCoeffs) result(Yip)
        integer, intent(in) :: np
        real(WP), intent(in) :: X, X0, h
        integer, intent(in) :: pstar
        real(WP), dimension(np,0:pstar), intent(in) :: bCoeffs
        real(WP), dimension(np) :: Yip
    end function InterpY


    module subroutine erk_destroy(me)
    
        type(ERK_class) :: me
        
    end subroutine erk_destroy    
    


    module subroutine erk_stepint(me, X0, Y0, F0, h, Y1, Yint12, FCalls, EstimateErr, Err, params)
        class(ERK_class), intent(inout) :: me
        real(WP), intent(in) :: X0
        real(WP), dimension(me%pDiffEqSys%n), intent(in) :: Y0
        real(WP), dimension(me%pDiffEqSys%n), intent(inout) :: F0
        real(WP), intent(in) :: h
        real(WP), dimension(size(Y0)), intent(out) :: Y1, Yint12
        integer, intent(out) :: FCalls
        logical, intent(in) :: EstimateErr
        real(WP), intent(out) :: Err
        real(WP), dimension(:), intent(in), optional :: params
    end subroutine erk_stepint

    module subroutine erk_intpcoeff(me, X0, Y0, h, Y1, Bip, Fcalls, params)
        class(ERK_class), intent(inout) :: me
        real(WP), intent(in) :: X0
        real(WP), dimension(me%pDiffEqSys%n), intent(in) :: Y0
        real(WP), intent(in) :: h
        real(WP), dimension(size(Y0)), intent(in) :: Y1
        real(WP), dimension(size(me%InterpStates),0:me%pstar), intent(out) :: Bip        
        integer, intent(out) :: FCalls
        real(WP), dimension(:), intent(in), optional :: params        
    end subroutine erk_intpcoeff


    end interface

    
    contains
    
    
    
end module ERK
        
