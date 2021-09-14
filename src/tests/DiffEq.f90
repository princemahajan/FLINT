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
!> \brief       DiffEq module
!! \details     It provides the differential equations and event functions
!!              for the test cases used for FLINT library testing.    
!! \author      Bharat Mahajan
!! \date        01/25/2019    
!
!############################################################################################    

    
module MyDiffEq

    use FLINT
    !use ddeabm_module

    use, intrinsic :: IEEE_ARITHMETIC

    implicit none
    
    real(WP), parameter :: mu = 0.012277471_WP

    ! user diff eq system

    type, extends(DiffEqSys) :: TBSys
        real(WP) :: mu = 0.2_WP
        real(WP) :: GM = 398600.436233_wp
    contains
        procedure :: F => TwoBodyDE
        procedure :: G => SampleEventTB
    end type TBSys
    
    type, extends(DiffEqSys) :: CR3BPSys
        real(WP) :: mu = mu
        real(WP) :: GM = 0.0_wp
    contains
        procedure :: F => CR3BPDE
        procedure :: G => SampleEventCR3BP
    end type CR3BPSys
    
    !type, extends(ddeabm_with_event_class) :: CR3BPSys_ddeabm
    !    real(WP) :: mu = mu         
    !    integer :: FEvals = 0       
    !    logical :: First = .TRUE. 
    !end type CR3BPSys_ddeabm



    contains
    
    function TwoBodyDE(me, X, Y, Params)
    
        implicit none
        
        intrinsic :: size  
        class(TBSys), intent(in) :: me !< Differential Equation object
        real(WP), intent(in) :: X
        real(WP), intent(in), dimension(:) :: Y
        real(WP), intent(in), dimension(:), optional :: Params
        
        real(WP), dimension(size(Y)) :: TwoBodyDE
        real(WP) :: MyParams
        
        if (present(Params)) then
            MyParams = Params(1)
        else
            MyParams = 1.0_WP        
        end if 
        
        TwoBodyDE(1:3) = Y(4:6)
        TwoBodyDE(4:6) = -me%GM/(norm2(Y(1:3))**3)*Y(1:3)
        
        ! add impulse with DV = 1 km/s
        !if(X > 100 .AND. X < 200) TwoBodyDE(4:6) = TwoBodyDE(4:6) + 0.006
        
        ! add finite maneuver
        !if(X > 50000 .AND. X < 100000) TwoBodyDE(4:6) = TwoBodyDE(4:6) + 0.006/100.0
        
    end function TwoBodyDE
        
        
    function CR3BPDE(me, X, Y, Params)    
        intrinsic :: size
        class(CR3BPSys), intent(in) :: me !< Differential Equation object
        real(WP), intent(in) :: X
        real(WP), intent(in), dimension(:) :: Y
        real(WP), intent(in), dimension(:), optional :: Params
        
        real(WP), dimension(size(Y)) :: CR3BPDE
        
        CR3BPDE = CR3BP_DiffEq(X, Y, me%mu)
    end function CR3BPDE    
    
    !subroutine CR3BPDE_DDEABM(me, X, Y, Ydot)    
    !    intrinsic :: size
    !    class(ddeabm_class), intent(inout) :: me !< Differential Equation object
    !    real(WP), intent(in) :: X
    !    real(WP), intent(in), dimension(:) :: Y
    !    real(WP), intent(out), dimension(:) :: Ydot
    !
    !    select type (me)
    !    class is (CR3BPSys_ddeabm)
    !        Ydot = CR3BP_DiffEq(X, Y, me%mu)
    !    end select
    !end subroutine CR3BPDE_DDEABM    
    !

    
    subroutine SampleEventTB(me, X, Y, EvalEvents, Value, Direction, LocEvent, LocEventAction)
            
        implicit none
        class(TBSys), intent(in) :: me !< Differential Equation object            
        real(WP), intent(in) :: X
        real(WP), dimension(:), intent(inout) :: Y
        integer, dimension(:), intent(in) :: EvalEvents
        real(WP), dimension(:), intent(out) :: Value
        integer, dimension(:), intent(out) :: Direction
        integer, intent(in), optional :: LocEvent
        integer(kind(FLINT_EVENTACTION_CONTINUE)), intent(out), optional :: LocEventAction
    
        if (EvalEvents(1)==1) Value(1) = norm2(Y(1:3)) - 20000.0
        
        if (EvalEvents(2)==1) Value(2) = 1 !Y(3)
        
        Direction = [0,0]
        
        if (present(LocEvent)) LocEventAction = FLINT_EVENTACTION_CONTINUE
        
    end subroutine SampleEventTB      
    

    subroutine SampleEventCR3BP(me, X, Y, EvalEvents, Value, Direction, LocEvent, LocEventAction)
            
        implicit none
        class(CR3BPSys), intent(in) :: me !< Differential Equation object            
        real(WP), intent(in) :: X
        real(WP), dimension(:), intent(inout) :: Y
        integer, dimension(:), intent(in) :: EvalEvents
        real(WP), dimension(:), intent(out) :: Value
        integer, dimension(:), intent(out) :: Direction
        integer, intent(in), optional :: LocEvent
        integer(kind(FLINT_EVENTACTION_CONTINUE)), intent(out), optional :: LocEventAction


        ! Event-1: detect a narrow square pulse or a discontinuous 
        ! event in both directions.
        ! This will need a small event step size to detect!
        ! We will also mask this event after the first detection to
        ! reduce unnecessary event checking after the first detection.
        if (EvalEvents(1) == 1) then
            if (X > 0.5_WP .AND. X < 0.5005_WP) then
                Value(1) = 1        
            else
                Value(1) = -1
            end if
            Direction(1) = 0
        end if

        ! Event-2: only detect y2-crossings in the region y1>0
        if (EvalEvents(2) == 1) then
            if (Y(1) < 0.0_WP) then
                Value(2) = IEEE_VALUE(1.0,IEEE_QUIET_NAN) 
            else
                Value(2) = Y(2)
            end if 
            Direction(2) = 1 ! detect in the increasing direction
        end if

        ! Event-3: Detect y2 crossings in region y1<0
        if (EvalEvents(3)==1) then
            if (Y(1) > 0.0_WP) then
                Value(3) = IEEE_VALUE(1.0,IEEE_QUIET_NAN) 
            else
                Value(3) = Y(2)
            end if 
            Direction(3) = 1 ! detect in the increasing direction
        end if
        
        ! Set actions for each event if located
        if (present(LocEvent) .AND. present(LocEventAction)) then

            if (LocEvent == 1) then
                LocEventAction = IOR(FLINT_EVENTACTION_CONTINUE, &
                                    FLINT_EVENTACTION_MASK)
            else if (LocEvent == 2) then
                ! Mask the event-2 after first trigger            
                LocEventAction = IOR(FLINT_EVENTACTION_CONTINUE, &
                                        FLINT_EVENTACTION_MASK)
            else if (LocEvent == 3) then
                ! Event-3 has been triggered
                ! now change the solution back to the initial condition
                ! and mask this event
                Y = [0.994_WP, 0.0_WP, 0.0_WP, 0.0_WP, -2.00158510637908252240537862224_WP, 0.0_WP]
                LocEventAction = IOR(FLINT_EVENTACTION_CHANGESOLY,&
                                        FLINT_EVENTACTION_MASK)
            end if
        end if
        
    end subroutine SampleEventCR3BP          

    
    pure function CR3BP_DiffEq(X, Y, mu)
        real(WP), intent(in) :: X
        real(WP), intent(in), dimension(:) :: Y
        real(WP), intent(in) :: mu
        
        real(WP), dimension(size(Y)) :: CR3BP_DiffEq
                
        real(WP) :: r1cube, r2cube
                
        r1cube =  norm2([(Y(1) + mu), Y(2)])**3
        r2cube =  norm2([(Y(1) - (1.0_WP - mu)), Y(2)])**3

        CR3BP_DiffEq(1:3) = Y(4:6)

        CR3BP_DiffEq(4) = Y(1) + 2.0_WP*Y(5) - (1.0_WP-mu)*(Y(1)+mu)/r1cube &
                                     - mu*(Y(1)-(1.0_WP-mu))/r2cube
        CR3BP_DiffEq(5) = Y(2) - 2.0_WP*Y(4) - (1.0_WP-mu)*Y(2)/r1cube &
                                     - mu*Y(2)/r2cube
        CR3BP_DiffEq(6) = 0.0_WP
    end function CR3BP_DiffEq
    
    function TBIOM(mu, X)
        real(wp), intent(in) :: mu
        real(wp), dimension(:), intent(in) :: X
        real(wp) :: TBIOM
    
        ! Keplerian energy
        TBIOM = norm2(X(4:6))**2/2 - mu/norm2(X(1:3))
    end function TBIOM
    
    function JacobiC(mu, X)
        real(wp), intent(in) :: mu
        real(wp), dimension(:), intent(in) :: X
        real(wp) :: JacobiC
        real(wp) :: Omega
        real(WP) :: r1, r2
        
        r1 = sqrt((X(1) + mu)**2 + X(2)**2 + X(3)**2)
        r2 = sqrt((X(1) - 1.0 + mu)**2 + X(2)**2+X(3)**2)
        Omega = 1.0_WP/2.0_WP*sum(X(1:3)**2) + (1.0-mu)/r1 + mu/r2
        ! Jacobian's Constant
        JacobiC = sum(X(4:6)**2)/2.0_WP - Omega
    end function JacobiC

    
end module MyDiffEq
