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
!> \brief       FLINT Base Module
!! \details     It includes the abstract intefaces for every solver in FLINT in addition to 
!!              the intefaces for the differential equation and event functions provided by
!!              the user. FLINT error codes are defined here as well.
!! \author      Bharat Mahajan
!! \date        Created: 01/25/2019    
!
!############################################################################################    

    
    
module FLINT_base
    
    ! All the basic definitions
    use FLINTUtils

    implicit none

    public
    
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
        enumerator :: FLINT_ERROR_INIT_REQD         = 13  !< FLINT Init is required first
        enumerator :: FLINT_ERROR_EVENTPARAMS       = 14  !< Event related parameters are wrong or missing
        enumerator :: FLINT_ERROR_EVENTROOT         = 15  !< Event root find failed
        enumerator :: FLINT_ERROR_CONSTSTEPSZ       = 16  !< Wrong value of StepSz for constant step size option
    end enum

    !> Error message strings for the supported error codes
    character(len=*), parameter, dimension(0:16) :: FLINT_ErrorStrings = [character(len=64) :: &
                            'Successfully Completed',                           &
                            'Termination Event Occurred',                       &
                            'Problem appears to be stiff',                      &
                            'Unknown error',                                    &
                            'Incorrect user-provided parameters',               &
                            'Maximum integration steps reached',                &
                            'Memory allocation error',                          &
                            'Memory deallocation error',                        &
                            'Incorrect dimensions',                             &
                            'Step size reduced to very small value',            &
                            'Incorrect interpolation states',                   &
                            'Incorrect interpolation array provided',           &
                            'Interpolation is not enabled',                     &
                            'Init is required',                            		&
                            'Incorrect event-related parameters ',              &
                            'Event root finding failed',                        &
                            'Incorrect step size for the non-adaptive option'   &
                            ]


    !> Supported FLINT Event Actions
    enum, bind(C)
        !> Continue the integration
        enumerator :: FLINT_EVENTACTION_CONTINUE    = 0   
        
        !> Terminate the integration
        enumerator :: FLINT_EVENTACTION_TERMINATE   = 1

        !> Change the solution
        enumerator :: FLINT_EVENTACTION_CHANGESOLY  = 2

        !> Restart the integration
        enumerator :: FLINT_EVENTACTION_RESTARTINT  = 4
        
        !> Mask the event for future detection. It can be combined with 
        !! any of the other event actions.
        enumerator :: FLINT_EVENTACTION_MASK        = 128   
    end enum



    !> Supported FLINT Event Options
    enum, bind(C)
        !> Locate the event as closely as possible using root-finding
        enumerator :: FLINT_EVENTOPTION_ROOTFINDING   = 0
        
        !> Return the solution in EventStates at the beginning of the integrator's step 
        !! (with natural or event step size) just after which the event function 
        !! changed sign.
        enumerator :: FLINT_EVENTOPTION_STEPBEGIN  = 1

        !> Return the solution in EventStates at the end of the integrator's step 
        !! (with natural or event step size) just before which the event function
        !! changed sign.        
        enumerator :: FLINT_EVENTOPTION_STEPEND    = 2
    end enum



    !> Abstract base class that must be inherited by all the numerical integrator classes
    type, abstract, public :: FLINT_class
        
        !> Status of the last operation completed
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
        
        real(WP), dimension(:), allocatable :: Xint     !< Array of the natural integrator steps 
        
        !> Array to store interpolation coefficients at integrator steps
        real(WP), dimension(:,:,:), allocatable :: Bip  
        
        !> Minimum step size of the integrator. It is always a positive number.
        !! Default value is used if user do not specifify it.
        real(WP) :: MinStepSize
        
        !> Maximum step-size of the integrator. It is always a positive number.
        !! Default value is the difference between the initial and final value of the
        !! independent variable.
        real(WP) :: MaxStepSize
        
        real(WP) :: InitialStepSize     !< Initial step size
        
        real(WP), dimension(6) :: StepSzParams  !< Step size method-specific tuning parameters
        
        logical :: EventsOn = .FALSE. !< Events checking off by default
        
        !> This specifies the maximum time after which the event function must be called 
        !! for each event detection. A 0 value in EventStepSz means that event will be checked 
        !! after the integrator takes a natural step always. Otherwise, events will be checked
        !! after the interval specified by EventStepSz if that interval size is smaller than the
        !! length of the current step size. In case it is bigger, then event will be checked at
        !! the integrator's natural step boundary. Setting is for each event spearately.
        real(WP), dimension(:), allocatable :: EventStepSz
        
        !> Event options for each event separately. Must be initialized
        !! by the Init Function. Specifies what event related info to return to the user. 
        integer, dimension(:), allocatable :: EventOptions

        !> Event tolerance value for each event in case of root-finding is used for event loation
        real(WP), dimension(:), allocatable  :: ETol

        real(WP) :: X0 !< Initial value of the independent variable
        real(WP) :: Xf !< Final value of the independent variable
        real(WP) :: h0 !< Initial step-size that was accepted
        real(WP) :: hf !< The accepted step-size that was used at the last step
        real(WP), dimension(:), allocatable :: Y0 !< Initial condition
        real(WP), dimension(:), allocatable :: Yf !< Final solution
        
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
        function DEFunc(me, X, Y, params)

            import :: WP, DiffEqSys       
            implicit none
            intrinsic :: size
            
            class(DiffEqSys), intent(inout) :: me      !< Differential Equation object
            real(WP), intent(in) :: X   !< Independent variable value
            real(WP), intent(in), dimension(me%n) :: Y !< Initial condition
            
            !> 'TBD' Optional array of double-precision parameters for sensitivity analysis
            real(WP), dimension(:), intent(in), optional :: params 
            
            real(WP), dimension(size(Y)) :: DEFunc !< First derivatives
        
        end function DEFunc

        
        
        !> Interface for a user-supplied procedure for computing events values during 
        !! the integration. It is also called after an event is detected in order for 
        !! the user to specify the event actions to be taken by the event handler. Note
        !! that if multiple events are triggered at the same location then only the event
        !! with the lowest index will be reported.
        subroutine EventFunc(me, X, Y, EvalEvents, Value, Direction, LocEvent, LocEventAction)
            
            import :: WP, DiffEqSys, FLINT_EVENTACTION_CONTINUE
            implicit none
            
            class(DiffEqSys), intent(inout) :: me  !< Object of class type DiffEqSys            
            
            !> The current independent variable value. When an event is detected, it contains
            !! the location of that event.
            real(WP), intent(in) :: X                  

            !> Current solution value. User can modify this value from inside the 
            !! event function, however it will only take effect if an event is located
            !! and the LocEvent contains the index of the located event. In that case, user
            !! user also must set the appropriate action in EventAction. 
            real(WP), dimension(:), intent(inout) :: Y 
            
            !> Size is m. 
            !! User should evaluate only those event functions for which the corresponding
            !! index entry is 1 in EvalEvents. 
            integer, dimension(:), intent(in) :: EvalEvents

            !> Size is m. 
            !! If Value(i) changes sign, then the event "i" may trigger based on Direction(i).
            !! If Value(i) changes from/to NaN, then the sign change is ignored.
            real(WP), dimension(:), intent(out) :: Value

            !> Size must be same as Value. 
            !! + Direction=1: only sign changes from -ve to +ve are detected
            !! + Direction=-1: only sign changes from +ve to -ve are detected
            !! + Direction=0: all sign changes are detected
            integer, dimension(:), intent(out) :: Direction
            
            !> FLINT will only call EventFunc with this parameter after an event is located.
            !! In that case, it will contain the index between 1 and m of the 
            !! located event. User must check whether this parameter is "present" before
            !! using it.
            integer, intent(in), optional :: LocEvent

            !> This parameter will be only present if an event has been located and
            !! LocEvent contains a valid index of the located event. In that case,
            !! user can specify a one of following actions that the event handler
            !! will take:
            !! + FLINT_EVENTACTION_CONTINUE: Continue the integration (default)
            !! + FLINT_EVENTACTION_TERMINATE: Terminate the integration at the current value of X
            !! + FLINT_EVENTACTION_CHANGESOLY: Change the solution value Y for the given X and use
            !! this value as the initial condition for computing the next step. The new
            !! value must be returned in "Y" by the user. EventStates will still contain the
            !! old value of Y that triggered the event.
            !! +  FLINT_EVENTACTION_RESTARTINT: Restart the integration from the current value of X
            integer(kind(FLINT_EVENTACTION_CONTINUE)), intent(out), optional :: LocEventAction

        end subroutine EventFunc
    end interface
        

    

    abstract interface
        
        !> Interface for the Initializaition method. It must be called before Integration.
        subroutine Init(me, DE, MaxSteps, Method, ATOl, RTol, InterpOn, InterpStates, &
            MinStepSize, MaxStepSize, StepSzParams, EventsOn, EventStepSz, EventOptions, EventTol)

            import :: FLINT_class, DiffEqSys, WP
        
            class(FLINT_class), intent(inout) :: me !< Object
            
            class(DiffEqSys), target :: DE      !< Differential Equation system object
            
            
            !> Must specify max no. of steps
            !! The implementation can use this number to decide how and how much to 
            !! allocate internal storage for storing integrator natural steps 
            !! as well as the interpolation coefficients for delayed interpolation.
            integer, intent(in) :: MaxSteps
            
            !> Type of Intergrator to use. Implemented methods are:
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
            !! will be computed and stored, e.g. [1,3,4,6,n].
            integer, dimension(:), intent(in), optional :: InterpStates
            
            !> User can specify minimum and maximum step sizes allowed.
            !! Note that these are positive quantities. Sign will be handled bu FLINT.
            real(WP), intent(in), optional :: MinStepSize, MaxStepSize
            
            !> Step-size computation specific parameters. User can specify these for a 
            !! soecific step-size computation method, otherwise defaults are used.
            !! For Hairer's DOPRI5/DOP853 codes:
            !! + StepSzParams(1): Safety-factor 
            !! + StepSzParams(2): Minimum safety-factor
            !! + StepSzParams(3): Maximum safety-factor
            !! + StepSzParams(4): Lund Stabilization parameter beta
            !! + StepSzParams(1): Lund Stabilization parameter beta multiplier
            real(WP), dimension(6), intent(in), optional :: StepSzParams
            
            logical, intent(in), optional :: EventsOn !< Set true for events detection
            
            !> If events detection is enabled, then user can specify the maximum interval
            !! after which events must be checked. If there are even number of sign changes between
            !! X and X+EventStepSz, that event cannot be detected. If EventStepSz is not specified,
            !! then events are always checked at the integrator's natural steps. Note that these 
            !! are positive quantities and FLINT will take care of the appropriate sign
            !! based on the direction of integration.
            real(WP), dimension(:), intent(in), optional :: EventStepSz


            !> Event Options that control the event information returned by FLINT
            !! Default is FLINT_EVENTOPTION_ROOTFINDING.
            integer(kind(FLINT_EVENTOPTION_ROOTFINDING)), dimension(:), &
                                intent(in), optional :: EventOptions

            !> Event tolerance values used for root-finding if enabled, size must be m
            real(WP), dimension(:), intent(in), optional :: EventTol    

        end subroutine Init
    


        
        !> Interface for the main Integrate method. It must be called after initialization and
        !! can be called multiple times with different IC and options without calling the init
        !! routine first every time.
        subroutine Integrate(me, X0, Y0, Xf, Yf, StepSz, UseConstStepSz, IntStepsOn, Xint, Yint, &
                                EventMask, EventStates, StiffTest, params)

            import :: FLINT_class, WP
        
            class(FLINT_class), intent(inout) :: me !< Object of the class type FLINT
            
            real(WP), intent(in) :: X0              !< Initial value of the independent variable
            real(WP), dimension(me%pDiffEqSys%n), intent(in) :: Y0 !< Initial conditions
            real(WP), intent(inout) :: Xf           !< Final value of the independent variable
            real(WP), dimension(me%pDiffEqSys%n), intent(out) :: Yf !< Final Solution
            
            !> User can specify the initial step-size to use. If it is 0, then FLINT will 
            !! compute the initial step-size internally. After the completion of the integration,
            !! the last accepted step size will be returned back. Note that user must provide
            !! the correct sign of the desired step size.
            real(WP), intent(inout) :: StepSz

            !> If true, then FLINT will not use adaptive step size algorithm, rather it will use 
            !! the initial step size value given in StepSz as the integrator's natural step. StepSz
            !! must take into account the sign of the step otherwise FLINT will return an error. 
            !! By default, this option is set to False.
            logical, intent(in), optional :: UseConstStepSz            

            !> If true, then FLINT will return the solution at the integrator's natural step size.
            logical, intent(in), optional :: IntStepsOn
            
            !> If IntStepsOn=true, user must provide this parameter. FLINT will allocate the
            !! memory and return the integrator's natural steps taken in Xint array. 
            !! FLINT will not explicitly deallocate this memory.
            real(WP), allocatable, dimension(:), intent(out), optional :: Xint            

            !> If IntStepsOn=true, user must provide this parameter. FLINT will allocate the
            !! memory and return the solution at the integrator's natural steps in Yint. The
            !! first column is the initial condition Y0 and the final column is Yf. FLINT will 
            !! not explicitly deallocate this memory. In case of an event action changing the
            !! solution at the location of the event, the new solution will be returned in Yint.
            real(WP), allocatable, dimension(:,:), intent(out), optional :: Yint

            !> If EventMask(i)=false, then the event "i" will not be detected.
            !! If all the events are masked, then the event detection will be turned off.
            logical, dimension(me%pDiffEqSys%m), intent(in), optional :: EventMask
            
            !> Solutions at the event locations for each of the unmasked event. If events
            !! are enabled, then user must specify this parameter. FLINT will internally
            !! allocate the storage and it will not explicitly deallocate this memory.
            !! EventStates(1:(n+2),i) = [Xevent, Yevent, EventID], where EventID is the
            !! index of the triggered event (between 1 and m) to which Xevent and Yevent
            !! corresponds. In case of an event action changing the solution at the 
            !! location of the event, EventStates will still return the orginal state.
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
        !! of the Init and Integrate method with the InterpOn option enabled. It can be called by the 
        !! user any number of times with different points at which the solution is sought.
        !! FLINT computes the interpolation coefficients during integration and stores them
        !! internally at each natural step of the integrator. During interpolation, these
        !! coefficients are used to compute the solutions quickly. This is much less
        !! computationally intensive than using the Integrate method.
        subroutine Interpolate(me, Xarr, Yarr, SaveInterpCoeffs)

            import :: FLINT_class, WP
        
            class(FLINT_class), intent(inout) :: me !< Object of class type FLINT
            
            !> Xarr is provided by the user with values of the independent variable at
            !! which the solution is sought. Xarr must be in ascending order and the min
            !! and max values must not be smaller and larger than X0 and Xf, respectively.
            real(WP), dimension(1:), intent(in) :: Xarr

            !> Yarr provides the solution at each point specified in Xarr. 
            real(WP), dimension(size(me%InterpStates),size(Xarr)), intent(out) :: Yarr           
            
            !> If SaveInterpCoeffs=False, then the internal memory which stores interpolation
            !! coefficients is freed. After that, Interpolate function can not be called again 
            !! without calling Integrate first. By default, SaveInterpCoeffs=False. In order
            !! to be able to call Interpolate again, user must set it to Trueat for every 
            !! call to Interpolate.
            logical, intent(in), optional :: SaveInterpCoeffs
           
        end subroutine Interpolate
        
        
        

        !> Query method to get the current state information of FLINT. Can be called 
        !! after each Integration to get the integration stats.
        subroutine Info(me, LastStatus, StatusMsg, nSteps, nAccept, nReject, nFCalls, &
                                InterpReady, h0, X0, Y0, hf, Xf, Yf)

            import :: FLINT_class, WP
        
            class(FLINT_class), intent(inout) :: me !< Object of class type FLINT
            
            !> The status code of the last operation completed
            integer(kind(FLINT_SUCCESS)), intent(out) :: LastStatus
            
            !> Status message corresponding to the status code in LastStatus
            character(len=:), allocatable, intent(out), optional :: StatusMsg

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
    
    
    