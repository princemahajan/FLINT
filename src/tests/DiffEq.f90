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
!> \brief       DiffEq module
!! \details     It provides the differential equations and event functions
!!              for the test cases used for FLINT library testing.    
!! \author      Bharat Mahajan
!! \date        01/25/2019    
!
!############################################################################################    

    
module MyDiffEq

    use FLINT

    implicit none
    
    ! user diff eq system

    type, extends(DiffEqSys) :: TBSys
        real(WP) :: mu = 0.2_WP
        real(WP) :: GM = 398600.436233_wp
    contains
        procedure :: F => TwoBodyDE
        procedure :: G => SampleEventTB
    end type TBSys
    
    type, extends(DiffEqSys) :: CR3BPSys
        real(WP) :: mu = 0.012277471_WP
        real(WP) :: GM = 0.0_wp
    contains
        procedure :: F => CR3BPDE
        procedure :: G => SampleEventCR3BP
    end type CR3BPSys
    
    contains
    
    pure function TwoBodyDE(me, X, Y, Params)
    
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
        
        
    pure function CR3BPDE(me, X, Y, Params)
    
        implicit none
        
        intrinsic :: size
        
        class(CR3BPSys), intent(in) :: me !< Differential Equation object
        real(WP), intent(in) :: X
        real(WP), intent(in), dimension(:) :: Y
        real(WP), intent(in), dimension(:), optional :: Params
        
        real(WP), dimension(size(Y)) :: CR3BPDE
        
        real(WP) :: MyParams
        
        real(WP), dimension(2) :: R1, R2
        
        if (present(Params)) then
            MyParams = Params(1)
        else
            MyParams = 1.0_WP        
        end if 
        
        R1 =  [(Y(1) + me%mu), Y(2)]
        R2 =  [(Y(1) - 1 + me%mu), Y(2)]

        CR3BPDE(1:3) = Y(4:6)
        CR3BPDE(4) = Y(1) + 2.0*Y(5) - (1.0-me%mu)*(Y(1) + me%mu)/norm2(R1)**3 &
                                - me%mu*(Y(1) - 1.0 + me%mu)/norm2(R2)**3

        CR3BPDE(5) = Y(2) - 2.0*Y(4) - (1.0-me%mu)*Y(2)/norm2(R1)**3 &
                                - me%mu*Y(2)/norm2(R2)**3
        CR3BPDE(6) = 0.0
        
        ! add impulse
        !if(X > 100 .AND. X < 200) TwoBodyDE(4:6) = TwoBodyDE(4:6) + 0.5*TwoBodyDE(4:6)
        
    end function CR3BPDE    
    
    

    
    pure subroutine SampleEventTB(me, EventID, X, Y, Value, Direction, Terminal)
            
        implicit none
        class(TBSys), intent(in) :: me !< Differential Equation object            
        integer, intent(in) :: EventID        
        real(WP), intent(in) :: X
        real(WP), dimension(:), intent(in) :: Y
        real(WP), dimension(:), intent(out) :: Value
        integer, dimension(:), intent(out) :: Direction
        logical, dimension(:), intent(out) :: Terminal
        
        Value = [1,1]
        
        if (EventID == 0 .OR. EventID == 1)   Value(1) = norm2(Y(1:3)) - 20000.0
        
        if (EventID == 0 .OR. EventID == 2)   Value(2) = 1 !Y(3)
        
        Direction = [0,0]
        Terminal = [.FALSE., .FALSE.]
        
    end subroutine SampleEventTB      
    
    pure subroutine SampleEventCR3BP(me, EventID, X, Y, Value, Direction, Terminal)
            
        implicit none
        class(CR3BPSys), intent(in) :: me !< Differential Equation object            
        integer, intent(in) :: EventID        
        real(WP), intent(in) :: X
        real(WP), dimension(:), intent(in) :: Y
        real(WP), dimension(:), intent(out) :: Value
        integer, dimension(:), intent(out) :: Direction
        logical, dimension(:), intent(out) :: Terminal
        
        Value = [Y(1),Y(2)]
        
        !if ((EventID == 0 .OR. EventID == 1))   Value(1) = Y(1)
        
        !if (EventID == 0 .OR. EventID == 2)   Value(2) = Y(2)
        
        Direction = [-1,1]
        Terminal = [.FALSE.,.FALSE.]
        
    end subroutine SampleEventCR3BP          
    
    function TBIOM(mu, X)
        real(wp), intent(in) :: mu
        real(wp), dimension(:), intent(in) :: X
        real(wp) :: TBIOM
    
        ! Keplerian energy
        TBIOM = norm2(X(4:6))**2/2 - mu/norm2(X(1:3))
    end function TBIOM
    
    function JacobiC(mu, X)
        real(wp), intent(in) :: mu
        real(wp), dimension(:,:), intent(in) :: X
        real(wp), dimension(size(X,2)) :: JacobiC
        real(wp), dimension(size(X,2)) :: Omega
        real(WP), dimension(size(X,2)) :: r1, r2
        
        r1 = sqrt((X(1,:) + mu)**2 + X(2,:)**2 + X(3,:)**2)
        r2 = sqrt((X(1,:) - 1.0 + mu)**2 + X(2,:)**2+X(3,:)**2)
        Omega = 1.0/2.0*sum(X(1:3,:)**2,1) + (1.0-mu)/r1 + mu/r2
        ! Jacobian's Constant
        JacobiC = sum(X(4:6,:)**2,1)/2 - Omega
    end function JacobiC

    
end module MyDiffEq