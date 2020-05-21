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
!> \brief       TestFLINT1 Main Program
!! \details     Main program to run multiple test cases for FLINT.
!! \author      Bharat Mahajan
!! \date        01/25/2019    
!
!############################################################################################    


    program TestFLINT

    
    use iso_fortran_env, only: output_unit
    use FLINT
    use ButcherTableaus
    use MyDiffEq
      
    implicit none

    ! Simulation Parameters
    
    ! Number of periodic orbits to run and loops for benchmark
    real, parameter :: norb = 2.5
    integer, parameter :: nloops = 100
    
    ! Turn on the events
    logical, parameter :: EventsEnabled = .TRUE.
    logical, dimension(3), parameter :: EvMask = [.TRUE.,.TRUE.,.TRUE.]
    
    ! Event step size
    real, parameter :: evstepsz = 0.00001

    ! event tolerance
    real(WP), dimension(3) :: evtol = [0.001,1.0e-9,1.0e-9]

    ! scalar tolerance
    real(WP),parameter :: tol   = 1.0e-11_WP  
    
    ! max no. of steps
    integer, parameter :: MAX_STEPS = 100000
    logical, parameter :: CONST_STEPSZ = .FALSE.

    ! random dispersion of IC: this multiplies the random number
    real, parameter :: randon = 0.000000000001
    
    ! Interpolation points
    integer, parameter :: nIp = int(1000.0*norb)
    
    ! results file
    character(len=20) :: fname = 'results.txt'
    character(len=20) :: sfname = 'states.txt'
    character(len=20) :: spfname = 'states_intp.txt'
    character(len=20) :: efname = 'events.txt'

    type(ERK_class) erkvar
    type(CR3BPSys) :: CR3BPDEObject
    
    real(WP) :: x0, xf, xfval, rnum
    real(WP), dimension(6) :: y0, yf, y0r
 
    character(len=:), allocatable :: errstr

    real(WP) :: stepsz0, ipdx, stepsz
    real(WP), dimension(:), allocatable, target :: Xint, Xarr
    real(WP), dimension(:,:), allocatable, target :: Yint, Yarr
    real(WP), allocatable, dimension(:,:) :: EventStates
    integer :: stiffstatus, stifftestval
    integer :: nsteps, naccpt, nrejct, fcalls, itr, nevents, ctr
    character(len=10) :: mname
    integer :: method
    
    real(WP) :: t0, tf
 
    ! Arenstorf orbit: Earth-moon mass-ratio
    real(wp), parameter :: EMmu  = 0.012277471_WP
    x0    = 0.0
    xf    = 17.0652165601579625588917206249_WP
    y0    = [0.994_WP, 0.0_WP, 0.0_WP, 0.0_WP, -2.00158510637908252240537862224_WP, 0.0_WP]  !! initial `y` value
    ! final time
    xf = xf*norb
    CR3BPDEObject%mu = Emmu    
    CR3BPDEObject%n = 6
    CR3BPDEObject%m = 3

    ! arrays for interpolated data
    ipdx = (xf-x0)/(nIp-1)
    Xarr = [(x0+ipdx*itr,itr=0,(nIp-1))]
    allocate(Yarr(6,size(Xarr)))
    
    stepsz0 = 0.0E-3_WP    ! let FLINT compute the initial step-size    
    stifftestval = 1  ! check for stiffness and stop integration if stiff
    
    ! create results file
    open(unit=17,file= fname, status = 'replace')
    write(17, *) '--- FLINT Stat Results (',norb, ' orbits, ', nloops, ' loops)  ---'
    write(17, *) ' '
    write(17, *) 'A. Natural Step-size'
    write(17, '(7A12)') 'Method', 'Time(s)', 'Closing Err','Jacobi Err', 'FCalls', 'Accepted', 'Rejected'

    ! states with natural step size
    open(unit=18,file= sfname, status = 'replace')
    write(18, '(7A12)') 'X','Y1','Y2','Y3','Y4','Y5','Y6'             

    ! states using interpolation
    open(unit=19,file= spfname, status = 'replace')
    write(19, '(7A12)') 'X','Y1','Y2','Y3','Y4','Y5','Y6'             

    ! event states
    open(unit=20,file= efname, status = 'replace')
    write(20, '(8A12)') 'Event','X', 'Y1','Y2','Y3','Y4','Y5','Y6'             

    ! Solution at natural step size

    do itr = 1,4
        
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
          
        call erkvar%Init(CR3BPDEObject, MAX_STEPS, Method=method, ATol=[tol*1.0e-3], RTol=[tol],&
            InterpOn=.FALSE., EventsOn=EventsEnabled, EventStepSz=[real(WP)::evstepsz,0.0,0.0], &
            EventTol = evtol)
        if (erkvar%status == FLINT_SUCCESS) then
        
            y0r = y0
            call CPU_TIME(t0)
            do ctr = 1,nloops
                stiffstatus = stifftestval
                stepsz = stepsz0
                xfval = xf
                call RANDOM_NUMBER(rnum) ! randomize initial condition
                y0r(1) = y0r(1) + rnum*randon*y0r(1)
                call erkvar%Integrate(x0, y0r, xfval, yf, StepSz=stepsz, UseConstStepSz=CONST_STEPSZ, &
                    IntStepsOn=.TRUE.,Xint = Xint, Yint = Yint, &
                    EventStates=EventStates, EventMask = EvMask,StiffTest=stiffstatus)
            end do
            call CPU_TIME(tf)
        
            if (stiffstatus == -1) write(17, *) mname//': problem is stiff'
        
            if (erkvar%status == FLINT_SUCCESS .OR. erkvar%status == FLINT_EVENT_TERM &
                        .OR. erkvar%status == FLINT_ERROR_MAXSTEPS) then 
                call erkvar%Info(stiffstatus, errstr, nAccept=naccpt, nReject=nrejct, nFCalls=fcalls)
            
                write(17, '(A12,3E12.3,3I12.1)') mname, (tf-t0), norm2(y0(1:3)-yf(1:3)), &
                    JacobiC(Emmu, Yint(:,1:1))-JacobiC(Emmu, Yint(:,ubound(Yint,2):ubound(Yint,2))), &
                    fcalls, naccpt, nrejct
            else
                call erkvar%Info(stiffstatus, errstr)
                write(17,*) mname//': Integrate failed: ', erkvar%status, ':', errstr
            end if
        
            ! write states
            if (itr == 2) then
                do ctr = 1, size(Xint)
                    write(18, '(7E18.9)') Xint(ctr), Yint(1:6,ctr)
                end do    
            end if

            ! write events
            if (EventsEnabled) then
                write(20, *) mname
                nevents = int(size(EventStates,2))
                do ctr = 1,nevents
                 write(20, '(I3.1, 7E18.9)') int(EventStates(8,ctr)),EventStates(1,ctr),&
                       EventStates(2,ctr),EventStates(3,ctr),EventStates(4,ctr),EventStates(5,ctr),&
                       EventStates(6,ctr),EventStates(7,ctr)
                end do
            end if
        

            if (allocated(Xint)) deallocate(Xint)
            if (allocated(Yint)) deallocate(Yint)
        else
            write(17,*) mname//': Init failed: ', erkvar%status            
        end if
    end do


    if (EventsEnabled) then
        write(17, *) 'X-axis (Y1) and Y-axis (Y2) crossing events'
        write(17, '(5A12)') 'Event Index','X','Y1','Y2','Y3'             
             
        nevents = int(size(EventStates,2))
        do ctr = 1,nevents
         write(17, '(I12.1, 4E12.4)') int(EventStates(8,ctr)),EventStates(1,ctr),&
               EventStates(2,ctr),EventStates(3,ctr),EventStates(4,ctr)
        end do
    end if
    
    ! Solution at interpolated grid
    write(17, *) 'B. Interpolated Grid'
    write(17, '(7A12)') 'Method', 'Time(s)', 'Closing Err','Jacobi Err', 'FCalls', 'Accepted', 'Rejected'

    do itr = 1,4
        
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
          
        call erkvar%Init(CR3BPDEObject, MAX_STEPS, Method=method, ATol=[tol*1.0e-3], RTol=[tol],&
            InterpOn=.TRUE., EventsOn=EventsEnabled, EventStepSz=[real(WP)::evstepsz,0.0,0.0],&
            EventTol = evtol)
        

        if (erkvar%status == FLINT_SUCCESS) then
            
            call CPU_TIME(t0)
            do ctr = 1,nloops
                stiffstatus = stifftestval
                stepsz = stepsz0
                xfval = xf
                call RANDOM_NUMBER(rnum) ! randomize initial condition
                y0r(1) = y0r(1) + rnum*randon*y0r(1)                
                call erkvar%Integrate(x0, y0r, xfval, yf, StepSz=stepsz, UseConstStepSz=CONST_STEPSZ, &
                    EventStates=EventStates, EventMask = EvMask,StiffTest=stiffstatus)
                call erkvar%Interpolate(Xarr,Yarr,.TRUE.)            
            end do
            call CPU_TIME(tf)            
            
            if (stiffstatus == -1) write(17, *) mname//': problem is stiff'

            if (erkvar%status == FLINT_SUCCESS .OR. erkvar%status == FLINT_EVENT_TERM &
                        .OR. erkvar%status == FLINT_ERROR_MAXSTEPS) then 
                call erkvar%Info(stiffstatus, errstr, nAccept=naccpt, nReject=nrejct, nFCalls=fcalls)
            
                write(17, '(A12,3E12.3,3I12.1)') mname, (tf-t0), norm2(y0(1:3)-yf(1:3)), &
                    JacobiC(Emmu, Yarr(:,1:1))-JacobiC(Emmu, Yarr(:,ubound(Yarr,2):ubound(Yarr,2))), &
                    fcalls, naccpt, nrejct

            else
                call erkvar%Info(stiffstatus, errstr)
                write(17,*) mname//': Integrate failed: ', erkvar%status, ':',errstr
            end if
    

            ! write states
            if (itr == 2) then
                do ctr = 1, size(Xarr)
                    write(19, '(7E18.9)') Xarr(ctr), Yarr(1:6,ctr)
                end do    
            end if

            ! write events
            if (EventsEnabled) then
                write(20, *) mname
                nevents = int(size(EventStates,2))
                do ctr = 1,nevents
                 write(20, '(I3.1, 7E18.9)') int(EventStates(8,ctr)),EventStates(1,ctr),&
                       EventStates(2,ctr),EventStates(3,ctr),EventStates(4,ctr),EventStates(5,ctr),&
                       EventStates(6,ctr),EventStates(7,ctr)
                end do
            end if



        else
            write(17,*) mname//': Init failed: ', erkvar%status            
        end if
    end do


    if (EventsEnabled) then
        write(17, *) 'X-axis (Y1) and Y-axis (Y2) crossing events'
        write(17, '(5A12)') 'Event Index','X','Y1','Y2','Y3'             
        nevents = int(size(EventStates,2))
        do ctr = 1,nevents
        write(17, '(I12.1, 4E12.4)') int(EventStates(8,ctr)),EventStates(1,ctr),&
              EventStates(2,ctr),EventStates(3,ctr),EventStates(4,ctr)
        end do
    end if
    
    
    contains
    
    
    end program TestFLINT

