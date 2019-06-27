!############################################################################################
!
! Copyright 2019 Bharat Mahajan
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
!> \brief       FLINT Base Module
!! \details     It includes the abstract intefaces for every solver in FLINT in addition to 
!!              the intefaces for the differential equation and event functions provided by the
!!              user. FLINT error codes are defined here as well.
!! \author      Bharat Mahajan
!! \date        01/25/2019    
!
!############################################################################################    

    
    
module FLINT_base
    
    use FLINTUtils

    implicit none

    public
    
    integer, parameter :: WORK_MAXSIZE = 20    
    
    !> Default number of maximum steps that intgerator is allowed to take
    integer, parameter :: INT_MAXSTEPS = 1024
    
    !> Supported Error Codes
    enum, bind(C)
        enumerator :: FLINT_SUCCESS                 = 0   !< Successfully completed
        enumerator :: FLINT_EVENT_TERM              = 1   !< Termination event occured
        enumerator :: FLINT_STIFF_PROBLEM           = 2   !< Stiff problem
        enumerator :: FLINT_ERROR                   = 3   !< Some unknown error
        enumerator :: FLINT_ERROR_PARAMS            = 4   !< Wrong user-provided parameters
        enumerator :: FLINT_ERROR_MAXSTEPS          = 5   !< Termination due to maximum steps reached
        enumerator :: FLINT_ERROR_MEMALLOC          = 6   !< Memory allocation error
        enumerator :: FLINT_ERROR_MEMDEALLOC        = 7   !< Memory deallocation error
        enumerator :: FLINT_ERROR_DIMENSION         = 8   !< Wrong dimensions of the differential equations specified
        enumerator :: FLINT_ERROR_STEPSZ_TOOSMALL   = 9   !< Step size reduced to very small value
        enumerator :: FLINT_ERROR_INTPSTATES        = 10  !< Wrong interpolation states specified
        enumerator :: FLINT_ERROR_INTP_ARRAY        = 11  !< Array provided for interpolation is wrong
        enumerator :: FLINT_ERROR_INTP_OFF          = 12  !< Interpolation is not enabled
        enumerator :: FLINT_ERROR_INIT_REQD         = 13  !< FLINT Init is required
        enumerator :: FLINT_ERROR_EVENTPARAMS       = 14  !< Required event parameters are not specified
        enumerator :: FLINT_ERROR_EVENTROOT         = 15  !< Event root find is failed
    end enum

    !> FLINT_class
    !! Abstract base class must be inherited by all the numerical integrator classes
    type, abstract, public :: FLINT_class
        
        !> Status of the last operation executed
        integer(kind(FLINT_SUCCESS)), public :: status = FLINT_SUCCESS
        
        !> Pointer to the user-supplied Differential Equation object
        class(DiffEqSys), pointer, public :: pDiffEqSys => null()

        !> Tolerance related constants: rtol and atol
        !! Dimension can be either 1 or n (dimension of the initial 
        !! condition vector). These are used to compute the scale factor (SC) needed
        !! for error computations using the following formula: 
        !! \f[ SC_i = atol_i + rtol_i |Y_{i}| \quad i=1..n,\f]
        !! where \f$|Y|\f$ is the absolute value of the current solution vector. 
        logical :: IsScalarTol = .TRUE.
        real(WP), dimension(:), allocatable  :: RTol
        real(WP), dimension(:), allocatable  :: ATol
        
        integer :: TotalSteps = 0           !< Total number of integrator steps
        integer :: AcceptedSteps = 0        !< Number of accepted integrator steps
        integer :: RejectedSteps = 0        !< Number of rejected integrator steps
        
        integer :: FCalls = 0 !< Total number of Diff Eq. function calls made
        
        integer :: MaxSteps = INT_MAXSTEPS !< Maximum number of steps integrator is allowed to take

        logical :: InterpOn = .FALSE. !< True for dense output
        integer, dimension(:), allocatable :: InterpStates !< States between 1 and n for dense output
        
        real(WP), allocatable, dimension(:) :: Xint     !< Array of the natural integrator steps 
        !real(WP), allocatable, dimension(:,:) :: Yint   !< Solution array at the natural integrator steps
        
        real(WP), allocatable, dimension(:,:,:) :: Bip  !< Array to store interpolation coefficients at integrator steps
        
        !> Minimum step-size of the integrator. It is always a positive number.
        !! Default value is used if user do not specifify it.
        real(WP) :: MinStepSize
        
        !> Maximum step-size of the integrator. It is always a positive number.
        !! Default value is the difference between the initial and final value of the independent variable.
        real(WP) :: MaxStepSize
        
        real(WP) :: InitialStepSize     !< Initial step-size
        
        real(WP), dimension(5) :: StepSzParams  !< Step size method-specific tuning parameters
        
        logical :: EventsOn = .FALSE. !< Events checking off by default
        
        !> A 0 value of event step size means events will be checked
        !! only at the integrator natural steps, otherwise events will be checked 
        !! after an interval specified by EventStepSz. Howeever, if EventStepSz is greater
        !! than the current integrator step-size then event is checked at the latter step.
        real(WP) :: EventStepSz = 0.0
        
        real(WP) :: X0 !< Initial value of the independent variable
        real(WP) :: Xf !< Final value of the independent variable
        real(WP) :: h0 !< Initial step-size that was accepted
        real(WP) :: hf !< The accepted step-size that was used at the last step
        real(WP), allocatable, dimension(:) :: Y0 !< Initial condition
        real(WP), allocatable, dimension(:) :: Yf !< Final solution
        
    contains

        !> Init is called once and only to specify new differential equaions, events, 
        !! and a few integrator options
        procedure (Init), deferred :: Init
        
        !> Integrate the solution from the initial to final value of the independent variable
        procedure (Integrate), deferred :: Integrate
        
        !> If interpolation is enabled, then it outputs the solution at the user-specified grid points
        procedure (Interpolate), deferred :: Interpolate
        
        !> Current FLINT object query function
        procedure (Info), deferred :: Info
        
    end type FLINT_class
    
    
    !> User must extend DiffEqSys and provide the DiffEq and Events function
    !! + n must be set to the number of the differential equations.
    !! + m must be set equal to the number of events, if events are used.
    !! + 'TBD' p must be set to the number of the parameters for sensitivity analysis
    type, abstract, public :: DiffEqSys
        integer :: n = 1
        integer :: m = 1
        integer :: p = 0
    contains
        procedure(DEFunc), deferred :: F
        procedure(EventFunc), deferred :: G
    end type DiffEqSys
   
    
    abstract interface

        !> Interface of the user-supplied Differential Equation function. It must
        !! return a vector of double precision values with dimension same as that of
        !! the input vector Y.
        pure function DEFunc(me, X, Y, params)

            import :: WP, DiffEqSys       
            implicit none
            intrinsic :: size
            
            class(DiffEqSys), intent(in) :: me      !< Differential Equation object
            real(WP), intent(in) :: X   !< Independent variable value
            real(WP), intent(in), dimension(:) :: Y !< Initial condition
            
            !> 'TBD' Optional array of double-precision parameters for sensitivity analysis
            real(WP), dimension(:), intent(in), optional :: params 
            
            real(WP), dimension(size(Y)) :: DEFunc !< First-derivatives
        
        end function DEFunc

        
        
        !> Interface for a user-supplied procedure for computing events
        pure subroutine EventFunc(me, EventID, X, Y, Value, Direction, Terminal)
            
            import :: WP, DiffEqSys
            implicit none
            
            class(DiffEqSys), intent(in) :: me  !< Object of class type DiffEqSys            
            
            !> Value is either 0 or between 1 and m. If 0, the event function must compute
            !! values for all the events. If its value is between 1 and m, then it needs 
            !! only compute only one event with index EventID. Rest of the event are ignored.
            integer, intent(in) :: EventID      
            
            real(WP), intent(in) :: X               !< Current independent variable value
            real(WP), dimension(:), intent(in) :: Y !< Current solution value
            
            !> Size is m. 
            !! If Value(i) changes sign, then the event with EventID i may trigger based on Direction(i).
            !! If Value(i) changes from/to NaN, then the event is not triggered.
            real(WP), dimension(:), intent(out) :: Value

            !> Size must be same as Value. 
            !! + Direction=1: only sign changes from -ve to +ve are detected
            !! + Direction=-1: only sign changes from +ve to -ve are detected
            !! + Direction=0: all sign changes are detected
            integer, dimension(:), intent(out) :: Direction
            
            !> Size must be same as that of Value. If true for any of the detected event, 
            !! then all the events are detected at the current step and the integration 
            !! is stopped with appropriate status code.
            logical, dimension(:), intent(out) :: Terminal
        
        end subroutine EventFunc        
    
    end interface
        
    abstract interface
        
        !> Interface for the Initializaition method. It must be called before Integration.
        subroutine Init(me, DE, MaxSteps, Method, ATOl, RTol, InterpOn, InterpStates, &
            MinStepSize, MaxStepSize, StepSzParams, EventsOn, EventStepSz)

            import :: FLINT_class, DiffEqSys, WP
        
            class(FLINT_class), intent(inout) :: me !< Object
            
            class(DiffEqSys), target :: DE      !< Differential Equation system object
            
            
            !> Must specify max no. of steps
            !! The implementation can use this number to decide how and how much to 
            !! allocate internal storage for storing integrator natural steps 
            !! as well as the interpolation coefficients for delayed interpolation.
            integer, intent(in) :: MaxSteps
            
            !> Type of Intergrator to use. Possible values are:
            !! + ERK_DOP853 (default)
            !! + ERK_DOP54
            !! + ERK_VERNER98R
            !! + ERK_VERNER65E
            integer, intent(in), optional :: Method
            
            !> Absolute tolerance, size must be either 1 or n
            real(WP), dimension(:), intent(in), optional :: ATol 
            
            !> Relative tolerance, size must be either 1 or n
            real(WP), dimension(:), intent(in), optional :: RTol    
            
            !> If true, then the interpolation coefficients along with integrator natural
            !! steps are computed and stored internally. User can use Interpolate method
            !! to compute solutions at any value of the indepedent variable without
            !! doing the actual integrating again. Note that if InterpOn is true, then
            !! IntStepsOn option of Integrate must be either omitted or set to false.
            logical, intent(in), optional :: InterpOn

            !> Specify the interpolation states for which the interpolation coefficients
            !! will be computed and stored, e.g. [1,4,6,n].
            integer, dimension(:), intent(in), optional :: InterpStates
            
            !> User can specify minimum and maximum step-sizes allowed.
            real(WP), intent(in), optional :: MinStepSize, MaxStepSize
            
            !> Step-size computation specific parameters. User can specify these for a 
            !! soecific step-size computation method, otherwise defaults are used.
            !! For Hairer's DOP853 code:
            !! + StepSzParams(1): Safety-factor 
            !! + StepSzParams(2): Minimum safety-factor
            !! + StepSzParams(3): Maximum safety-factor
            !! + StepSzParams(4): Lund Stabilization parameter beta
            !! + StepSzParams(1): Lund Stabilization parameter
            real(WP), dimension(5), intent(in), optional :: StepSzParams
            
            logical, intent(in), optional :: EventsOn !< Set true for events detection
            
            !> If events detection is enabled, then user can specify the maximum step-size
            !! at which events are detected. If there are even number of sign changes between
            !! X and X+EventStepSz, that event will not be detected. If EventStepSz is greater
            !! than the current natural step-size of the integrator, then EventStepSz is ignored.
            !! If it is not specified, then events are always checked at the natural steps.
            real(WP), intent(in), optional :: EventStepSz
            end subroutine Init
    
        !> Interface for the main Integrate method. It must be called after initialization and
        !! can be called multiple times with different IC and options without calling the init
        !! routine first every time.
        subroutine Integrate(me, X0, Y0, Xf, Yf, StepSz, IntStepsOn, Xint, Yint, EventMask, EventStates, StiffTest, params)

            import :: FLINT_class, WP, WORK_MAXSIZE
        
            class(FLINT_class), intent(inout) :: me !< Object of the class type FLINT
            
            real(WP), intent(in) :: X0              !< Initial value of the independent variable
            real(WP), dimension(me%pDiffEqSys%n), intent(in) :: Y0 !< Initial conditions
            real(WP), intent(inout) :: Xf           !< Final value of the independent variable
            real(WP), dimension(me%pDiffEqSys%n), intent(out) :: Yf !< Final Solution
            
            !> User can specify the initial step-size to use. If it is 0, then FLINT will 
            !! compute the initial step-size internally. After the completion of the integration,
            !! the last accepted step-size will be returned back.
            real(WP), intent(inout) :: StepSz
            
            !> If true, then FLINT will return the solution at the integrator's natural step size.
            logical, intent(in), optional :: IntStepsOn
            
            !> If IntStepsOn=true, user must provide this parameter. FLINT will allocate the
            !! memory and return the integrator's natural steps in Xint. FLINT will not explicitly
            !! deallocate this memory.
            real(WP), allocatable, dimension(:), intent(out), optional :: Xint            

            !> If IntStepsOn=true, user must provide this parameter. FLINT will allocate the
            !! memory and return the solution at integrator's natural steps in Yint. The
            !! first column is the initial condition Y0 and the final column is Yf. FLINT will 
            !! not explicitly deallocate this memory.
            real(WP), allocatable, dimension(:,:), intent(out), optional :: Yint

            !> If EventMask(i)=false, then the event with Eventid 'i' will not be checked.
            !! If all the events are masked, then event checking will be turned off.
            logical, dimension(me%pDiffEqSys%m), intent(in), optional :: EventMask
            
            !> Solutions at the event locations for each of the unmasked event. If events
            !! are enabled, then user must specify this parameter. FLINT will internally
            !! allocate the storage and it will not explicitly deallocate this memory.
            !! EventStates(1:(n+2),i) = [Xevent, Yevent, EventID], where EventID is the ID
            !! of the event to which Xevent and Yevent corresponds.
            real(WP), allocatable, dimension(:,:), intent(out), optional :: EventStates

            !> FLINT will use Hairer's DOP853 alogorithm to test for stiffness
            !! of the differential equations if StffTest == 1 or StffTest == 2. 
            !! Test is disabled by default or if Stifftest == 0.
            !! + StiffTest == 1: If DiffEq are detected to be stiff then -1 is returned
            !!  in StiffTest and the integration halts at the current step with appropriate
            !! status code.
            !! + StiffTest == 2: If DiffEq are detected to be stiff then -1 is returned
            !!  in StiffTest but the integration continues.
            integer, intent(inout), optional :: StiffTest
            
            !> Parameters for sensitivity analysis. Not implemented yet.
            real(WP), dimension(:), intent(in), optional :: params
        end subroutine Integrate
        

        
        !> Method for delayed interpolation. It must be called after at least one invocation
        !! of the Init and Integrate method with InterpOn option enabled. It can be called by the 
        !! user any number of times with different points at which the solution is sought.
        !! FLINT computes the interpolation coefficients during integration and stores them
        !! internally at each natural step of the integrator. During interpolation, these
        !! coefficients are used to compute the solutions quickly. This is much less
        !! computationally intensive than the Integrate method.
        subroutine Interpolate(me, Xarr, Yarr, SaveInterpCoeffs)

            import :: FLINT_class, WP, WORK_MAXSIZE
        
            class(FLINT_class), intent(inout) :: me !< Object of class type FLINT
            
            !> Xarr is provided by the user with values of the independent variable at
            !! which the solution is sought. Xzrr must be in ascending order and the min
            !! and max values must not be smaller and larger than X0 and Xf, respectively.
            real(WP), dimension(1:), intent(in) :: Xarr

            !> Yarr provides the solution at each point specified in Xarr. 
            real(WP), dimension(size(me%InterpStates),size(Xarr)), intent(out) :: Yarr           
            
            !> If SaveInterpCoeffs=false, then the internal memory which stores interpolatio
            !! coefficients is freed. After that, Interpolate function can not be called again.
            !! By default, SaveInterpCoeffs=false. In order to be able to call Interpolate again,
            !! user must specify true for this parameter at each invocation of Interpolate.
            logical, intent(in), optional :: SaveInterpCoeffs
           
        end subroutine Interpolate

        !> Query method to get the current state information of FLINT. Should be called 
        !! after each Integration to get the statistics.
        subroutine Info(me, LastStatus, nSteps, nAccept, nReject, nFCalls, &
                                InterpReady, h0, X0, Y0, hf, Xf, Yf)

            import :: FLINT_class, WP
        
            class(FLINT_class), intent(inout) :: me !< Object of class type FLINT
            
            !> The status code of the last operation
            integer(kind(FLINT_SUCCESS)), intent(out) :: LastStatus
            
            integer, intent(out), optional :: nSteps !< Total number of steps taken
            integer, intent(out), optional :: nAccept !< Number of Accepted steps
            integer, intent(out), optional :: nReject !< Number of rejected steps
            integer, intent(out), optional :: nFCalls !< Function calls made not including event calls
            logical, intent(out), optional :: InterpReady !< Whether Interpolation is enabled
            
            real(WP), intent(out), optional :: h0 !< The initial accepted step-size
            real(WP), intent(out), optional :: X0 !< The initial value of the independent variable
            real(WP), dimension(:), intent(out), optional :: Y0 !< Initial condition
            real(WP), intent(out), optional :: hf !< Final accepted step-size
            real(WP), intent(out), optional :: Xf !< The final value of the independent variable
            real(WP), dimension(:), intent(out), optional :: Yf !< Final solution
            
        end subroutine Info
        
        
    end interface
    
    contains
    
    
end module FLINT_base
    
    
    