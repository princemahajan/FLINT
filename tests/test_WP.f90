!############################################################################################
!
! Copyright 2024 Bharat Mahajan
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
!> \brief       test_WorkPrecision Main Program
!! \details     Main program to generate FLINT solver's work-precision data
!! \author      Bharat Mahajan
!! \date        01/05/2024    
!
!############################################################################################    


module MyDiffEq

    use FLINT, WP => FLINT_WP
    use, intrinsic :: IEEE_ARITHMETIC

    implicit none

    ! two-body diffeq
    type, extends(DiffEqSys) :: TBSys
        real(WP) :: mu = 398600.436233_wp
        
        real(WP), dimension(6) :: Y0 = [0.0_WP, -3096.70185149293_WP, &
                                        -6183.97070198107_WP, &
                                        10.0141943725295_WP, &
                                        0.0_WP, 0.0_WP]
    contains
        procedure :: F => TwoBodyDE
        ! procedure :: G => TBPEvent
        procedure, nopass :: IOM => Energy
    end type TBSys
        
    ! CR3BP diff eq system    
    type, extends(DiffEqSys) :: CR3BPSys
        real(WP) :: mu = 1.0_WP/(81.30059_WP + 1.0_WP)

        real(WP), dimension(6) :: Y0 = [1.02202151273581740824714855590570360_WP, &   
                                        0.0_WP, &   
                                        -0.182096761524240501132977765539282777_WP, &   
                                        0.0_WP, &   
                                        -0.103256341062793815791764364248006121_WP, &   
                                        0.0_WP]

        real(WP) :: Period = 1.5111111111111111111111111111111111111111_WP
    contains
        procedure :: F => CR3BPDE
        procedure, nopass :: IOM => JacobiConstant
    end type CR3BPSys

    contains

    function TwoBodyDE(me, X, Y, Params)        
        intrinsic :: size  

        class(TBSys), intent(inout) :: me
        real(WP), intent(in) :: X
        real(WP), intent(in), dimension(me%n) :: Y
        real(WP), intent(in), dimension(:), optional :: Params
        
        real(WP), dimension(size(Y)) :: TwoBodyDE
                
        TwoBodyDE(1:3) = Y(4:6)
        TwoBodyDE(4:6) = -me%mu/(norm2(Y(1:3))**3)*Y(1:3)                
    end function TwoBodyDE
           
    pure function Energy(mu, Y)
        real(WP), intent(in) :: mu
        real(wp), dimension(:), intent(in) :: Y
        real(wp) :: Energy
    
        ! Keplerian energy
        Energy = norm2(Y(4:6))**2/2 - mu/norm2(Y(1:3))
    end function Energy            
        
    function CR3BPDE(me, X, Y, Params)    
        intrinsic :: size
        class(CR3BPSys), intent(inout) :: me !< Differential Equation object
        real(WP), intent(in) :: X
        real(WP), intent(in), dimension(me%n) :: Y
        real(WP), intent(in), dimension(:), optional :: Params
        
        real(WP), dimension(size(Y)) :: CR3BPDE
        real(WP) :: r1cube, r2cube, mu

        mu = me%mu

        r1cube =  norm2([(Y(1) + mu), Y(2), Y(3)])**3
        r2cube =  norm2([(Y(1) - (1.0_WP - mu)), Y(2), Y(3)])**3

        CR3BPDE(1:3) = Y(4:6)
        CR3BPDE(4) = Y(1) + 2.0_WP*Y(5) - (1.0_WP-mu)*(Y(1)+mu)/r1cube &
                                     - mu*(Y(1)-(1.0_WP-mu))/r2cube
        CR3BPDE(5) = Y(2) - 2.0_WP*Y(4) - (1.0_WP-mu)*Y(2)/r1cube &
                                        - mu*Y(2)/r2cube
        CR3BPDE(6) = - (1.0_WP-mu)*Y(3)/r1cube - mu*Y(3)/r2cube
    end function CR3BPDE    
            
    pure function JacobiConstant(mu, Y)
        real(WP), intent(in) :: mu
        real(WP), dimension(:), intent(in) :: Y
        real(WP) :: Omega, r1, r2

        real(WP) :: JacobiConstant
        
        r1 = sqrt((Y(1) + mu)**2 + Y(2)**2 + Y(3)**2)
        r2 = sqrt((Y(1) - 1.0 + mu)**2 + Y(2)**2+Y(3)**2)
        Omega = 1.0_WP/2.0_WP*sum(Y(1:2)**2) + (1.0-mu)/r1 + mu/r2
        JacobiConstant = sum(Y(4:6)**2)/2.0_WP - Omega
    end function JacobiConstant
    
end module MyDiffEq



program TestFLINT_WP

    use iso_fortran_env, only: output_unit
    use FLINT
    use MyDiffEq
        
    implicit none

    ! Simulation Parameters

    ! solvers to run from 1 to 4
    character(len=:), allocatable :: mname
    integer :: method

    ! Tolerances range
    integer :: i
    real(WP), parameter :: atol   = 1.0e-80_WP  
    real(WP), parameter :: rtol(*) = [ ( (1.0_WP/10.0_WP)**i, i=1,16)]

    ! timed loop
    integer, parameter :: tloops = 100

    ! max no. of steps
    integer, parameter :: MAX_STEPS = 30000

    real(WP), parameter :: stepsz00 = 0.0_WP
    integer :: stifftest = 1

    ! Interpolation points
    integer, parameter :: nsteps = 15

    real(WP), parameter :: pi = 4*atan(1.0_WP)

    real(WP) :: ts0, ts1, td0, td1

    ! results file
    character(len=:), allocatable :: fname
    character(len=:), allocatable :: errstr

    integer, dimension(4), parameter :: un = 17

    integer :: ctr, stiffstatus, fc, itr, mctr
    real, dimension(size(rtol)) :: IOMerr, IOMerr_d
    integer, dimension(size(rtol)) :: fcalls
    real(WP), dimension(size(rtol)) :: time_s, time_d

    real(WP) :: t, IOM0, stepsz0
    real(WP), dimension(6) :: Y

    real(WP), allocatable :: tin(:), Yin(:, :), IOMin(:)
    real(WP) :: tip(nsteps), Yip(6,nsteps), IOMip(nsteps)

    type(ERK_class) :: erk
    type(CR3BPSys) :: eom

    eom%n = 6
    t = eom%Period
    IOM0 = eom%IOM(eom%mu, eom%Y0)
              
    tip = [(eom%Period/(nsteps-1)*i, i = 0, nsteps-1)]

    do mctr = 1, 4

    select case (mctr)
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

    
    do ctr = 1, size(rtol)

        ! integrate sparse
        call erk%Init(eom, MAX_STEPS, Method=method, InterpOn=.FALSE., &
                                    RTol=[rtol(ctr)], ATol=[atol])
 
        call cpu_time(ts0)
        do itr = 1, tloops
            stiffstatus = stifftest            
            stepsz0 = stepsz00    
            call erk%Integrate(0.0_WP, eom%Y0, t, Y, IntStepsOn=.TRUE., &
                Xint=tin, Yint=Yin, StepSz=stepsz0, StiffTest=stiffstatus) 
        end do
        call cpu_time(ts1)
 
        if (erk%status /= FLINT_SUCCESS) then
            call erk%Info(stiffstatus, errstr)
            print *, rtol(ctr), mname//': Sparse int failed: '//errstr
        end if

        time_s(ctr) = ts1 - ts0

        ! integrate dense
        call erk%Init(eom, MAX_STEPS, Method=method, InterpOn=.TRUE., &
                                    RTol=[rtol(ctr)], ATol=[atol])
 
        call cpu_time(td0)
        do itr = 1, tloops
            stepsz0 = stepsz00
            stiffstatus = stifftest            
            call erk%Integrate(0.0_WP, eom%Y0, t, Y, &
                StepSz=stepsz0, StiffTest=stiffstatus) 
        end do
        call cpu_time(td1)
 
        if (erk%status /= FLINT_SUCCESS) then
            call erk%Info(stiffstatus, errstr)
            print *, rtol(ctr), mname//': Dense int failed: '//errstr
        end if

        time_d(ctr) = td1 - td0

        if (erk%status == FLINT_SUCCESS) then

            if (.NOT. allocated(IOMin)) allocate(IOMin(size(tin)))            

            ! interpolate       
            call erk%Interpolate(tip, Yip,.FALSE.)       
        
            ! compute stats
            call erk%Info(stiffstatus, errstr, nFCalls=fc)

            if (erk%status == FLINT_SUCCESS) then
                fcalls(ctr) = fc
                do i = 1, nsteps
                    IOMip(i) = eom%IOM(eom%mu, Yip(:,i))
                end do
                do i = 1, size(tin)
                    IOMin(i) = eom%IOM(eom%mu, Yin(:,i))
                end do
                IOMerr(ctr) = maxval(abs(IOMin-IOM0))
                IOMerr_d(ctr) = maxval(abs(IOMip-IOM0))
            else
                print *, mname//': Interpolation failed: ', erk%status, ',', errstr
            end if
        end if

        if (allocated(tin)) deallocate(tin)
        if (allocated(Yin)) deallocate(Yin)
        if (allocated(IOMin)) deallocate(IOMin)    
    
    end do


    ! create results file
    fname = "wp_"//mname//".txt"
    open(unit=un(1), file= fname, status = 'replace')
    write(un(1), '(2A10, 4A15)') 'Rel Tol', 'FCalls', 'Time-Sparse', 'Err-Sparse', 'Time-Dense', 'Err-Dense'
    do ctr = 1, size(rtol)
        write(un(1), '(1E10.1, 1I10, 1E15.6, 1E15.6, 1E15.6, 1E15.6)') rtol(ctr), fcalls(ctr), &
                time_s(ctr), IOMerr(ctr), time_d(ctr), IOMerr_d(ctr)
    end do

    end do

contains


end program TestFLINT_WP



