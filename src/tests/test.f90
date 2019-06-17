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
!! \version     0.9
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
    
    ! Number of periodic orbits to run
    real, parameter :: norb = 1
    
    ! scalar tolerance
    real(wp),parameter :: tol   = 1.0e-11_WP  
    
    ! Interpolation points
    integer, parameter :: nIp = int(1000.0*norb)
    
    ! results file
    character(len=20) :: fname = 'results.txt'
    
    type(ERK_class) erkvar
    type(CR3BPSys) :: CR3BPDEObject
    
    real(WP) :: x0, xf
    real(WP), dimension(6) :: y0, yf
 
    real(WP) :: stepsz, stepsz0, xfval, ipdx
    real(WP), dimension(:), allocatable :: Xint, Xarr
    real(WP), dimension(:,:), allocatable :: Yint, Yarr
    real(WP), allocatable, dimension(:,:) :: EventStates
    integer :: stiffstatus, stifftestval
    integer :: nsteps, naccpt, nrejct, fcalls, itr, nevents, ctr
    character(len=10) :: mname
    integer :: method
 
    ! Arenstorf orbit: Earth-moon mass-ratio
    real(wp), parameter :: EMmu  = 0.012277471_WP
    x0    = 0.0
    xf    = 17.0652165601579625588917206249_WP
    y0    = [0.994_WP, 0.0_WP, 0.0_WP, 0.0_WP, -2.00158510637908252240537862224_WP, 0.0_WP]  !! initial `y` value
    ! final time
    xf = xf*norb
    CR3BPDEObject%mu = Emmu    
    CR3BPDEObject%n = 6
    CR3BPDEObject%m = 2

    ! arrays for interpolated data
    ipdx = (xf-x0)/(nIp-1)
    Xarr = [(x0+ipdx*itr,itr=0,(nIp-1))]
    allocate(Yarr(6,size(Xarr)))
    
    stepsz0 = 0.0_WP    ! let FLINT compute the initial step-size    
    stifftestval = 1  ! check for stiffness and stop integration if stiff
    
    ! create results file
    open(unit=17,file= fname, status = 'replace')
    write(17, *) '--- FLINT Results ---'
    write(17, *) ' '
    write(17, *) 'A. Natural Step-size'
    write(17, '(6A12)') 'Method', 'Closing Err','Jacobi Err', 'FCalls', 'Accepted', 'Rejected'
    
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
          
        stiffstatus = stifftestval
        call erkvar%Init(CR3BPDEObject, 5000, Method=method, ATol=[tol*1.0e-3], RTol=[tol],&
            InterpOn=.FALSE., EventsOn=.TRUE.)
        if (erkvar%status == FLINT_SUCCESS) then
            call erkvar%Integrate(x0, y0, xf, yf, StepSz=stepsz0,  &
                IntStepsOn=.TRUE.,Xint = Xint, Yint = Yint, &
                EventStates=EventStates, EventMask = [.TRUE.,.TRUE.],StiffTest=stiffstatus)

            if (stiffstatus == -1) write(17, *) mname//': problem is stiff'

            if (erkvar%status == FLINT_SUCCESS) then 
                call erkvar%Info(stiffstatus, nAccept=naccpt, nReject=nrejct, nFCalls=fcalls)
            
                write(17, '(A12,2E12.3,3I12.1)') mname, norm2(y0(1:3)-yf(1:3)), &
                    maxval(JacobiC(Emmu, Yint))-minval(JacobiC(Emmu, Yint)), &
                    fcalls, naccpt, nrejct
            end if
            if (allocated(Xint)) deallocate(Xint)
            if (allocated(Yint)) deallocate(Yint)
            if (allocated(EventStates) .AND. itr < 4) deallocate(EventStates)
        else
            write(17,*) mname//': Init failed: ', erkvar%status            
        end if
    end do

    write(17, *) 'X-axis (Y1) and Y-axis (Y2) crossing events'
    write(17, '(5A12)') 'EventID','X','Y1','Y2','Y3'             
                
    nevents = int(size(EventStates,2))
    do ctr = 1,nevents
    write(17, '(I12.1, 4E12.3)') int(EventStates(8,ctr)),EventStates(1,ctr),&
              EventStates(2,ctr),EventStates(3,ctr),EventStates(4,ctr)
    end do
    

    ! Solution at interpolated grid
    write(17, *) 'B. Interpolated Grid'
    write(17, '(6A12)') 'Method', 'Closing Err','Jacobi Err', 'FCalls', 'Accepted', 'Rejected'
    
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
          
        stiffstatus = stifftestval
        call erkvar%Init(CR3BPDEObject, 10000, Method=method, ATol=[tol*1.0e-3], RTol=[tol],&
            InterpOn=.TRUE., EventsOn=.TRUE.)
        if (erkvar%status == FLINT_SUCCESS) then
            call erkvar%Integrate(x0, y0, xf, yf, StepSz=stepsz0,  &
                EventStates=EventStates, EventMask = [.TRUE.,.TRUE.],StiffTest=stiffstatus)
            
            if (stiffstatus == -1) write(17, *) mname//': problem is stiff'

            if (erkvar%status == FLINT_SUCCESS) then
                call erkvar%Interpolate(Xarr,Yarr,.TRUE.)
                call erkvar%Info(stiffstatus, nAccept=naccpt, nReject=nrejct, nFCalls=fcalls)
            
                write(17, '(A12,2E12.3,3I12.1)') mname, norm2(y0(1:3)-yf(1:3)), &
                    maxval(JacobiC(Emmu, Yarr))-minval(JacobiC(Emmu, Yarr)), &
                    fcalls, naccpt, nrejct
            end if
            if (allocated(EventStates) .AND. itr < 4) deallocate(EventStates)
        else
            write(17,*) mname//': Init failed: ', erkvar%status            
        end if
    end do

    write(17, *) 'X-axis (Y1) and Y-axis (Y2) crossing events'
    write(17, '(5A12)') 'EventID','X','Y1','Y2','Y3'             
    nevents = int(size(EventStates,2))
    do ctr = 1,nevents
    write(17, '(I12.1, 4E12.3)') int(EventStates(8,ctr)),EventStates(1,ctr),&
              EventStates(2,ctr),EventStates(3,ctr),EventStates(4,ctr)
    end do
    
    
    contains
    
    
    end program TestFLINT

