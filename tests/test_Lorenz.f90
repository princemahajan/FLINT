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
!> \brief       Lorenz Test Main Program
!! \details     Main program to run Lorenz test.
!! \author      Bharat Mahajan
!! \date        09/17/2021    
!
!############################################################################################    


module LorenzDiffEq

    use FLINT, WP => FLINT_WP

    use, intrinsic :: IEEE_ARITHMETIC

    implicit none
    
    ! user diff eq system

    type, extends(DiffEqSys) :: LorenzSys
        real(WP) :: sigma = 10.0_WP
        real(WP) :: rho = 28.0_WP
        real(WP) :: beta = 8.0_WP/3.0_WP
    contains
        procedure :: F => LorenzDE
    end type LorenzSys
    
    contains
    
    function LorenzDE(me, X, Y, Params)    
        intrinsic :: size  
        class(LorenzSys), intent(inout) :: me !< Differential Equation object
        real(WP), intent(in) :: X
        real(WP), intent(in), dimension(me%n) :: Y
        real(WP), intent(in), dimension(:), optional :: Params
        
        real(WP), dimension(size(Y)) :: LorenzDE

        LorenzDE(1) = me%sigma*(Y(2) - Y(1))
        LorenzDE(2) = Y(1)*(me%rho - Y(3)) - Y(2)
        LorenzDE(3) = Y(1)*Y(2) - me%beta*Y(3)
    end function LorenzDE
        
end module LorenzDiffEq




    program LorenzTest
    
    use iso_fortran_env, only: output_unit
    use FLINT
    use LorenzDiffEq
      
    implicit none

    ! Simulation Parameters
    
    ! Number of periodic orbits to run and loops for benchmark
    integer, parameter :: nloops = 100
    integer, parameter :: floops = 100000000
    
    ! solvers to run from 1 to 4
    integer, parameter :: nmethod = 4
    
    ! scalar tolerance
    real(WP), parameter :: rtol   = 1.0e-11_WP  
    real(WP), parameter :: atol   = rtol*1.0e-3_WP  
    
    ! max no. of steps
    integer, parameter :: MAX_STEPS = 100000
    logical, parameter :: CONST_STEPSZ = .FALSE.
    
    ! random dispersion of IC: this multiplies the random number
    real(WP), parameter :: randon = 0.000000000000_WP

    ! Interpolation points
    integer, parameter :: nIp = int(1000.0)
    
    ! results file
    character(len=20) :: fname = 'results_Lorenz.txt'
    character(len=20) :: spfname = 'states_ip_Lorenz.txt'

    type(ERK_class) :: erkvar
    type(LorenzSys) :: LorenzObj
    
    real(WP) :: x0, xf, xfval, rnum
    real(WP), dimension(3) :: y0, yf, y0r
 
    character(len=:), allocatable :: errstr

    real(WP) :: stepsz0, ipdx, stepsz
    real(WP), dimension(:), allocatable, target :: Xint, Xarr
    real(WP), dimension(:,:), allocatable, target :: Yint, Yarr
    integer :: stiffstatus, stifftestval
    integer :: nsteps, naccpt, nrejct, fcalls, steps, itr, nevents, ctr
    character(len=10) :: mname
    integer :: method
    
    real(WP) :: T0, Tf

    real(WP), dimension(3), parameter :: params = [10.0_WP, 28.0_WP, 8.0_WP/3.0_WP]


    y0 = [1.0_WP, 0.0_WP, 0.0_WP]

    x0 = 0.0_WP
    xF = 100.0_WP
    
    LorenzObj%n = 3
    LorenzObj%m = 0
    ! LorenzObj%G => SampleEvent
    
    stepsz0 = 0.0E-2_WP    ! let FLINT compute the initial step-size    
    stifftestval = 1*0  ! check for stiffness and stop integration if stiff

    ! arrays for interpolated data
    ipdx = (xf-x0)/(nIp-1)
    Xarr = [(x0+ipdx*itr,itr=0,(nIp-1))]
    Xarr(nIp) = xf
    allocate(Yarr(3,size(Xarr)))
    
    ! create results file
    open(unit=17,file= fname, status = 'replace')
    write(17, *) '--- FLINT Stats ---'
    write(17, *) ' '
    write(17, *) 'A. Internal Step-size'
    write(17, '(6A12)') 'Method', 'Time(s)', 'FCalls', 'Accepted', 'Rejected', 'Total Steps'

    ! states using interpolation
    open(unit=19,file= spfname, status = 'replace')
    write(19, '(4A12)') 'X','Y1','Y2','Y3'

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
          
        call erkvar%Init(LorenzObj, MAX_STEPS, Method=method, ATol=[atol], RTol=[rtol],&
            InterpOn=.FALSE.)
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
                    IntStepsOn=.FALSE.,StiffTest=stiffstatus)
            end do
            call CPU_TIME(tf)
        
            if (stiffstatus == -1) write(17, *) mname//': problem is stiff'
        
            if (erkvar%status == FLINT_SUCCESS .OR. erkvar%status == FLINT_EVENT_TERM &
                        .OR. erkvar%status == FLINT_ERROR_MAXSTEPS) then 
                call erkvar%Info(stiffstatus, errstr, nAccept=naccpt, nReject=nrejct, nFCalls=fcalls, nSteps = steps)
            
                write(17, '(A12,1F12.3,4I12.1)') mname, (tf-t0), fcalls, naccpt, nrejct, steps
            else
                call erkvar%Info(stiffstatus, errstr)
                write(17,*) mname//': Integrate failed: ', erkvar%status, ':', errstr
            end if
        
            ! ! write states
            ! if (itr == 2) then
            !     do ctr = 1, size(Xint)
            !         write(18, '(7E18.9)') Xint(ctr), Yint(1:3,ctr)
            !     end do    
            ! end if
        
            if (allocated(Xint)) deallocate(Xint)
            if (allocated(Yint)) deallocate(Yint)
        else
            write(17,*) mname//': Init failed: ', erkvar%status            
        end if
    end do
    
     ! Solution at interpolated grid
     write(17, *) 'B. Interpolated Grid'
     write(17, '(6A12)') 'Method', 'Time(s)', 'FCalls', 'Accepted', 'Rejected', 'Total Steps'
    
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
          
         call erkvar%Init(LorenzObj, MAX_STEPS, Method=method, ATol=[atol], RTol=[rtol],&
             InterpOn=.TRUE.)
        
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
                     StiffTest=stiffstatus, params=params)
             end do
             call CPU_TIME(tf)            
             call erkvar%Interpolate(Xarr,Yarr,.TRUE.)            
            
             if (stiffstatus == -1) write(17, *) mname//': problem is stiff'
    
             if (erkvar%status == FLINT_SUCCESS .OR. erkvar%status == FLINT_EVENT_TERM &
                         .OR. erkvar%status == FLINT_ERROR_MAXSTEPS) then 
                 call erkvar%Info(stiffstatus, errstr, nAccept=naccpt, nReject=nrejct, nFCalls=fcalls, nSteps=steps)
            
                 write(17, '(A12,1F12.3,4I12.1)') mname, (tf-t0), fcalls, naccpt, nrejct, steps
    
             else
                 call erkvar%Info(stiffstatus, errstr)
                 write(17,*) mname//': Integrate failed: ', erkvar%status, ':',errstr
             end if
    
    
             ! write states
             if (itr == 2) then
                 do ctr = 1, size(Xarr)
                     write(19, '(7E18.9)') Xarr(ctr), Yarr(1:3,ctr)
                 end do    
             end if    
    
         else
             write(17,*) mname//': Init failed: ', erkvar%status            
         end if
     end do
        
    
    ! function test
    y0r = y0*0.0_WP
    call CPU_TIME(t0)
    do ctr = 1,floops
        !y0 = y0 + 0.0001_WP*y0
        yf = LorenzObj%F(x0, y0, params)
        y0r = y0r + yf
    end do
    call CPU_TIME(tf)
    write(17, *) ' '    
    write(17, '(I12,A13,F12.2,A5)') floops, ' Func calls: ', (tf-t0)*1.0e3_WP, ' ms'    
    
    
    contains

    


    end program LorenzTest

