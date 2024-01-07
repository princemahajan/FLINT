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
!> \brief       test_TBP Main Program
!! \details     Main program to test FLINT performance for Two-Body Orbits
!! \author      Bharat Mahajan
!! \date        02/04/2023    
!
!############################################################################################    

module TwoBodyDynamics

    use FLINT
    !use ddeabm_module

    use, intrinsic :: IEEE_ARITHMETIC

    implicit none
    
    ! user diff eq system

    type, extends(DiffEqSys) :: TBSys
        real(WP) :: GM = 398600.436233_wp
    contains
        procedure :: F => TwoBodyDE
        ! procedure :: G => TBPEvent
        procedure :: Energy
    end type TBSys
        
    contains
    

    function TwoBodyDE(me, X, Y, Params)        
        intrinsic :: size  

        class(TBSys), intent(inout) :: me
        real(WP), intent(in) :: X
        real(WP), intent(in), dimension(me%n) :: Y
        real(WP), intent(in), dimension(:), optional :: Params
        
        real(WP), dimension(size(Y)) :: TwoBodyDE
        real(WP) :: MyParams
                
        TwoBodyDE(1:3) = Y(4:6)
        TwoBodyDE(4:6) = -me%GM/(norm2(Y(1:3))**3)*Y(1:3)                
    end function TwoBodyDE
           
    
    pure function Energy(me, X)
        class(TBSys), intent(in) :: me
        real(wp), dimension(:), intent(in) :: X
        real(wp) :: Energy
    
        ! Keplerian energy
        Energy = norm2(X(4:6))**2/2 - me%GM/norm2(X(1:3))
    end function Energy
    

    subroutine TBPEvent(me, X, Y, EvalEvents, Value, Direction, LocEvent, LocEventAction)
            
        implicit none
        class(DiffEqSys), intent(inout) :: me !< Differential Equation object            
        real(WP), intent(in) :: X
        real(WP), dimension(:), intent(inout) :: Y
        integer, dimension(:), intent(in) :: EvalEvents
        real(WP), dimension(:), intent(out) :: Value
        integer, dimension(:), intent(out) :: Direction
        integer, intent(in), optional :: LocEvent
        integer(kind(FLINT_EVENTACTION_CONTINUE)), intent(out), optional :: LocEventAction

        ! Event: detect all apsidal line crossings (no stopping)
        if (EvalEvents(1) == 1) then
            Value(1) = dot_product(Y(1:3), Y(4:6))        
            Direction(1) = 0
        end if
                
    end subroutine TBPEvent          

    
end module TwoBodyDynamics




program TestFLINT_TBP

    use iso_fortran_env, only: output_unit
    use FLINT
    use TwoBodyDynamics
        
    implicit none

    ! Simulation Parameters

    ! Number of periodic orbits to run and loops for benchmark
    real(WP), parameter :: norb = 10
    integer, parameter :: nloops = 10

    logical, parameter :: EventsEnabled = .FALSE.

    ! solvers to run from 1 to 4
    integer, parameter :: nmethod = 4

    ! scalar tolerance range for testing
    real(WP), parameter :: rtol   = 1.0e-12_WP
    real(WP), parameter :: atol   = rtol*1.0e-3_WP  
    real(WP), parameter, dimension(1) :: evtol = rtol

    ! max no. of steps
    integer, parameter :: MAX_STEPS = 10000
    logical, parameter :: CONST_STEPSZ = .FALSE.

    ! step size computation settings
    real(WP), dimension(6), parameter :: StepSzParams = [0.9_WP, 0.333_WP, 6.0_WP, &
                                                        0.0_WP, 0.0_WP, 0.0_WP]

    ! random dispersion of IC: this multiplies the random number
    real(WP), parameter :: randon = 1.0e-6_WP

    ! Interpolation points
    integer, parameter :: nIp = int(1000.0*norb)

    real(WP), parameter :: pi = 4*atan(1.0_WP)

    ! results file
    character(len=20) :: fname = 'perf_TBP.txt'
    character(len=20) :: sfname = 'states_TBP.txt'
    character(len=20) :: spfname = 'states_int_TBP.txt'
    character(len=20) :: efname = 'events_TBP.txt'

    integer, dimension(4), parameter :: un = [17, 18, 19, 20]

    type(ERK_class) erkvar
    type(TBSys) :: TBDobj

    real(WP) :: x0, xf, rnum, x0r
    real(WP), dimension(6) :: yf, y0r

    character(len=:), allocatable :: errstr

    real(WP) :: stepsz0, MinStepSz, MaxStepSz
    real(WP), dimension(:), allocatable, target :: Xint, Xarr
    real(WP), dimension(:,:), allocatable, target :: Yint, Yarr
    real(WP), allocatable, dimension(:,:) :: EventStates
    integer :: stiffstatus, stifftestval
    integer :: nsteps, naccpt, nrejct, fcalls, itr, nevents, ctr
    character(len=10) :: mname
    integer :: method

    real(WP) :: t0, tf, P, sma

    ! Two-Body Orbit - Molniya (26600 km, 0.74, 63.4 deg, 270 deg, 0, 0)
    real(WP), parameter, dimension(6) :: y0 = [0.0_WP, &
                                                -3096.70185149293_WP, &
                                                -6183.97070198107_WP, &
                                                10.0141943725295_WP, &
                                                0.0_WP, &
                                                0.0_WP]

    
    TBDobj%n = 6
    TBDobj%m = 0
                                            
    sma = - TBDobj%GM/(2*TBDobj%Energy(y0))
    P = 2.0_WP*pi*sqrt(sma**3/TBDobj%GM)
    x0 = 0
    xf = P*norb
                                            
    ! Integration settings
    stepsz0 = 0.0E-3_WP    ! let FLINT compute the initial step-size    
    stifftestval = 1  ! check for stiffness and stop integration if stiff
    MinStepSz = 1.0e-6_WP
    MaxStepSz = 1.0e+6_WP
    
    ! create results file
    open(unit=un(1), file= fname, status = 'replace')
    write(un(1), *) '--- FLINT TBP Stats (',norb, ' orbits, ', nloops, ' loops)  ---'
        
    write(un(1), *) ' '
    write(un(1), *) 'A. Internal Step Size'
    write(un(1), '(7A12)') 'Method', 'Time(s)', 'Closing Err','Jacobi Err', 'FCalls', 'Accepted', 'Rejected'

    ! states with natural step size
    open(unit=un(2), file= sfname, status = 'replace')
    write(un(2), '(7A12)') 'X','Y1','Y2','Y3','Y4','Y5','Y6'             

    ! ! states using interpolation
    ! open(unit=un(3), file= spfname, status = 'replace')
    ! write(un(3), '(7A12)') 'X','Y1','Y2','Y3','Y4','Y5','Y6'             

    ! event states
    if (EventsEnabled) then
        open(unit=un(4), file= efname, status = 'replace')
        write(un(4), '(8A12)') 'Event', 'X','Y1','Y2','Y3','Y4','Y5','Y6'             
    end if

    ! Solution at internal step size

    do itr = 1, nmethod
        
        select case (itr)
        case (1)
            mname = 'DOP54'
            method = ERK_DOP54
        case (2)
            mname = 'DOP853'
            method = ERK_DOP853
        case (3)
            mname = 'Verner65E'
            method = ERK_Verner65E
        case (4)
            mname = 'Verner98R'
            method = ERK_Verner98R
        end select            

        call erkvar%Init(TBDobj, MAX_STEPS, Method=method, ATol=[atol], RTol=[rtol],&
            MinStepSize = 1.0e-6_WP, MaxStepSize = 1.0e+6_WP,  &
            InterpOn=.FALSE., EventsOn=EventsEnabled, EventTol = evtol)

        if (erkvar%status == FLINT_SUCCESS) then
              
            call CPU_TIME(t0)
            do ctr = 1, nloops
                stiffstatus = stifftestval

                ! forward propagate
                call erkvar%Integrate(x0, y0, xf, yf, StepSz=stepsz0, UseConstStepSz=CONST_STEPSZ, &
                    IntStepsOn=.TRUE., Xint=Xint, Yint=Yint, EventStates=EventStates, StiffTest=stiffstatus)

                ! propagate back to IC
                call erkvar%Integrate(xf, yf, x0r, y0r, StepSz=stepsz0, UseConstStepSz=CONST_STEPSZ, &
                    IntStepsOn=.FALSE., EventStates=EventStates, StiffTest=stiffstatus)

            end do
            call CPU_TIME(tf)
        
            if (stiffstatus == -1) write(un(1), *) mname//': problem is stiff'
        
            if (erkvar%status == FLINT_SUCCESS .OR. erkvar%status == FLINT_EVENT_TERM &
                        .OR. erkvar%status == FLINT_ERROR_MAXSTEPS) then 
                call erkvar%Info(stiffstatus, errstr, nAccept=naccpt, nReject=nrejct, nFCalls=fcalls)
            
                write(un(1), '(A12,1F12.3,2E12.3,3I12.1)') mname, (tf-t0), norm2(y0r(1:3)-y0(1:3)), &
                        TBDobj%Energy(y0r)-TBDobj%Energy(y0), fcalls, naccpt, nrejct
            else
                call erkvar%Info(stiffstatus, errstr)
                write(un(1),*) mname//': Integrate failed: ', erkvar%status, ':', errstr
            end if
        
            ! write states
            if (itr == 1) then
                do ctr = 1, size(Xint)
                    write(un(2), '(7E18.9)') Xint(ctr), Yint(1:6,ctr)
                end do    
            end if

            ! write events
            if (EventsEnabled) then
                write(un(4), *) mname
                nevents = int(size(EventStates, 2))
                do ctr = 1,nevents
                    write(un(4), '(I3.1, 7E18.9)') int(EventStates(8,ctr)), EventStates(1,ctr), &
                        EventStates(2,ctr), EventStates(3,ctr), EventStates(4,ctr), &
                        EventStates(5,ctr), EventStates(6,ctr), EventStates(7,ctr)
                end do
            end if
        
            if (allocated(Xint)) deallocate(Xint)
            if (allocated(Yint)) deallocate(Yint)
        else
            write(un(1),*) mname//': Init failed: ', erkvar%status            
        end if
    end do

contains


end program TestFLINT_TBP

