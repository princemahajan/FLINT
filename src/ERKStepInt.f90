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
!> \brief       ERK Step Integrate Submodule
!! \details     It provides the single integration step of the ERK methods.
!! \author      Bharat Mahajan
!! \date        Created: 09/21/2021    
!
!############################################################################################
    
submodule (ERK) ERKStepInt

    contains    

    module subroutine erk_stepint(me, X0, Y0, h, Y1, Yint12, FCalls, &
                                        EstimateErr, Err, params)
        class(ERK_class), intent(inout) :: me
        real(WP), intent(in) :: X0
        real(WP), dimension(me%pDiffEqSys%n), intent(in) :: Y0
        real(WP), intent(in) :: h
        real(WP), dimension(size(Y0)), intent(out) :: Y1, Yint12
        integer, intent(out) :: FCalls
        logical, intent(in) :: EstimateErr
        real(WP), intent(out) :: Err
        real(WP), dimension(:), intent(in), optional :: params

        select case (me%method)
        case (ERK_DOP853)
            if (EstimateErr) then
                call DOP853_stepint(me%pDiffEqSys, X0, Y0, h, me%k, Y1, Yint12, FCalls, Err, &
                    me%IsScalarTol, me%RTol, me%ATol, params)
            else
                call DOP853_stepint(me%pDiffEqSys, X0, Y0, h, me%k, Y1, Yint12, FCalls, Err, &
                    params=params)
            end if
        case (ERK_DOP54)
            if (EstimateErr) then
                call DOP54_stepint(me%pDiffEqSys, X0, Y0, h, me%k, Y1, Yint12, FCalls, Err, &
                    me%IsScalarTol, me%RTol, me%ATol, params)
            else
                call DOP54_stepint(me%pDiffEqSys, X0, Y0, h, me%k, Y1, Yint12, FCalls, Err, &
                    params=params)
            end if
        case (ERK_VERNER65E)
            if (EstimateErr) then
                call stepint(me%pDiffEqSys, X0, Y0, h, Verner65E_FSAL, Verner65E_sint, me%k, &
                        Verner65E_a, Verner65E_b, Verner65E_c, Y1, Yint12, FCalls, &
                        Err, me%IsScalarTol, me%RTol, me%ATol, Verner65E_e, params)        
            else
                call stepint(me%pDiffEqSys, X0, Y0, h, Verner65E_FSAL, Verner65E_sint, me%k, &
                        Verner65E_a, Verner65E_b, Verner65E_c, Y1, Yint12, FCalls, &
                        Err, params=params)        
            end if
        case (ERK_VERNER98R)
            if (EstimateErr) then
                call stepint(me%pDiffEqSys, X0, Y0, h, Verner98R_FSAL, Verner98R_sint, me%k, &
                        Verner98R_a, Verner98R_b, Verner98R_c, Y1, Yint12, FCalls, &
                        Err, me%IsScalarTol, me%RTol, me%ATol, Verner98R_e, params)        
            else
                call stepint(me%pDiffEqSys, X0, Y0, h, Verner98R_FSAL, Verner98R_sint, me%k, &
                        Verner98R_a, Verner98R_b, Verner98R_c, Y1, Yint12, FCalls, &
                        Err, params=params)        
            end if
        case default
            me%status = FLINT_ERROR
        end select
    end subroutine erk_stepint





    !> Procedure to compute the interpolation coefficients
    module subroutine erk_intpcoeff(me, X0, Y0, h, Y1, Bip, Fcalls, params)
        class(ERK_class), intent(inout) :: me
        real(WP), intent(in) :: X0
        real(WP), dimension(me%pDiffEqSys%n), intent(in) :: Y0
        real(WP), intent(in) :: h
        real(WP), dimension(size(Y0)), intent(in) :: Y1
        real(WP), dimension(size(me%InterpStates),0:me%pstar), intent(out) :: Bip        
        integer, intent(out) :: FCalls
        real(WP), dimension(:), intent(in), optional :: params        

        select case (me%method)
        case (ERK_DOP853)
            call DOP853_IntpCoeff(me%pDiffEqSys, X0, Y0, h, Y1, me%InterpStates, me%k, Bip, FCalls, params)
        case (ERK_DOP54)
            call DOP54_IntpCoeff(me%pDiffEqSys, X0, Y0, h, Y1, me%InterpStates, me%k, Bip, FCalls, params)
        case (ERK_VERNER65E)
            call IntpCoeff(me%pDiffEqSys, X0, Y0, h, me%InterpStates, Verner65E_sint, &
                    Verner65E_s, Verner65E_a, Verner65E_c(Verner65E_sint+1:Verner65E_s), Verner65E_d, &
                    [Verner65E_di_NZ,-1], Verner65E_pstar, me%k, Bip, FCalls, params)
        case (ERK_VERNER98R)
            call IntpCoeff(me%pDiffEqSys, X0, Y0, h, me%InterpStates, Verner98R_sint, &
                    Verner98R_s, Verner98R_a, Verner98R_c(Verner98R_sint+1:Verner98R_s), Verner98R_d, &
                    [Verner98R_di_NZ,-1], Verner98R_pstar, me%k, Bip, FCalls, params)
        case default
            me%status = FLINT_ERROR
        end select
    end subroutine erk_intpcoeff



    
    !> It advances integrator by 1 step by computing intermediate stages
    subroutine stepint(pDiffEqSys, X0, Y0, h, IsFSALMethod, sint, k, a, b, c, Y1, &
                        Yint12, FCalls, Err, IsScalarTol, RTol, ATol, e, params)
        class(DiffEqSys), intent(in)  :: pDiffEqSys
        real(WP), intent(in) :: X0
        real(WP), dimension(pDiffEqSys%n), intent(in) :: Y0
        real(WP), intent(in) :: h
        logical, intent(in) :: IsFSALMethod
        integer, intent(in) :: sint
        real(WP), dimension(size(Y0),1:sint), intent(inout) :: k
        real(WP), dimension(:), intent(in) :: a
        real(WP), dimension(1:sint), intent(in) :: b
        real(WP), dimension(2:sint), intent(in) :: c
        real(WP), dimension(size(Y0)), intent(out) :: Y1, Yint12
        integer, intent(out) :: FCalls
        real(WP), intent(out) :: Err
        logical, intent(in), optional :: IsScalarTol
        real(WP), dimension(size(Y0)), intent(in), optional :: RTol, ATol
        real(WP), dimension(1:sint), intent(in), optional :: e
        real(WP), dimension(:), intent(in), optional :: params
        
        real(WP), dimension(size(Y0)) :: Yint
        integer :: i, n
        logical :: EstimateErr
        
        EstimateErr = .FALSE.
        if (present(IsScalarTol) .AND. present(RTol) &
            .AND. present(ATol) .AND. present(e)) EstimateErr = .TRUE.

        n = pDiffEqSys%n

        ! first stage: we require k(:,1) to already contain F0
        !k(:,1) = F0
        
        ! compute rest of the stages
        ! This compact and beautiful code is slower than ugly hardcoded one!
        do i = 2,sint
            block
            integer :: astart, j
            real(WP), dimension(n) :: aijkj
            ! starting index of a_ij, where i=i, j=1:(i-1)
            ! It is a series sum: 1 + 2 + 3 + ...
            astart = int((i-1)/2.0*(i-2))
            ! In my testing, I am seeing matmul performing better than 
            ! MKL's gemv. However, do loop is faster than both but 
            ! but still slower than hard-coded computations
            ! call gemv(k(:,1:(i-1)), me%a((astart+1):(i-1)), aijkj)
            !aijkj = matmul(k(:,1:(i-1)), me%a((astart+1):(astart+i-1)))
            aijkj = 0.0_WP
            do j = 1,(i-1)
                aijkj = aijkj + k(:,j)*a(astart + j)
            end do
            
            Yint = Y0 + h*aijkj
            k(:,i) = pDiffEqSys%F(X0 + h*c(i), Yint, params)
            
            ! Yint for the last non-FSAL stage is needed for stiffness detection
            ! This is from Hairer's DOP853, so not sure its correct for other solvers
            if (i == (sint-1)) Yint12 = Yint 
            end block
        end do
        
        ! propagate the solution to the next step for error computation
        ! For FSAL methods, Yint already has the new solution
        
        ! For non-FSAL methods, compute the solution at the next step
        if (IsFSALMethod .EQV. .FALSE.) then
            ! for stiffness detection, we need last non-FSAL stage solution
            Yint12 = Yint
            Yint = Y0 + h*matmul(k(:,1:sint), b(1:sint))
        end if

        if (EstimateErr) then
            ! Estimate the error. We assume that the error coeffcients are precomputed,
            ! i.e., e_i = b_i - bhat_i
            block
            real(WP), dimension(n) :: temp, Sc
            integer :: j
            ! The scale factor for error computations
            if (IsScalarTol .EQV. .FALSE.) then
                Sc = ATol + RTol*max(abs(Y0), abs(Yint))
            else
                Sc = ATol(1) + RTol(1)*max(abs(Y0), abs(Yint))
            end if 
            temp = 0.0_WP
            do concurrent (j = 1:sint)
                temp = temp + k(:,j)*e(j)
            end do
            Err = sqrt(sum((h*temp/Sc)**2)/n)
            end block
        else
            Err = 0.0_WP
        end if
        
        ! If this step is accepted then return the new solution
        ! else just return the initial condition
        if (Err <= 1.0_WP) then 
            Y1 = Yint
            if (IsFSALMethod .EQV. .FALSE.) then
                ! Evaluate F0 for the next step
                k(:,1) = pDiffEqSys%F(X0 + h, Y1, params)
                ! Number of function calls made so far
                FCalls = sint            
            else
                ! In FSAL methods, the integration last stage is F0 for the next step
                k(:,1) = k(:,sint)
                ! Function calls made if the step is accepted
                FCalls = sint - 1                    
            end if
        else
            Y1 = Y0
            ! Function calls made if the step is rejected
            if (IsFSALMethod .EQV. .FALSE.) then
                FCalls = sint - 1            
            else
                FCalls = sint - 2            
            end if
        end if    
    end subroutine stepint

    
    subroutine IntpCoeff(pDiffEqSys, X0, Y0, h, InterpStates, sint, s, a, c, d, &
                                dinz, pstar, k, IpCoeff, FCalls, params)
        class(DiffEqSys), intent(in)  :: pDiffEqSys
        real(WP), intent(in) :: X0
        real(WP), dimension(pDiffEqSys%n), intent(in) :: Y0
        real(WP), intent(in) :: h
        !real(WP), dimension(size(Y0)), intent(in) :: Y1
        integer, dimension(:), intent(in) :: InterpStates
        integer, intent(in) :: sint, s
        real(WP), dimension(:), intent(in) :: a
        real(WP), dimension(sint+1:s), intent(in) :: c
        real(WP), dimension(:,:), intent(in) :: d
        integer, dimension(:), intent(in) :: dinz
        integer, intent(in) :: pstar
        real(WP), dimension(size(Y0),1:s), intent(inout) :: k
        real(WP), dimension(size(InterpStates),0:pstar), intent(out) :: IpCoeff
        integer, intent(out) :: Fcalls
        real(WP), dimension(:), intent(in), optional :: params
        
        integer :: i, j
        
        ! Compute extra stages needed for interpolation
    
        ! compute rest of the stages (the following compact code is slower!)   
        do i = (sint+1),s
            block
            integer :: astart, j
            real(WP), dimension(size(Y0)) :: aijkj
            ! starting index of a_ij, where i=i, j=1:(i-1)
            ! It is a series sum: 1 + 2 + 3 + ...
            astart = int((i-1)/2.0*(i-2))
            ! In my testing, I am seeing matmul performing better than 
            ! MKL's gemv. However, do concurrent is faster than both but 
            ! but still slower than hard-coded computations
            ! call gemv(me%k(:,1:(i-1)), me%a((astart+1):(i-1)), aijkj)
            ! aijkj = matmul(me%k(:,1:(i-1)), me%a((astart+1):(astart+i-1)))
            aijkj = 0.0_WP
            do concurrent (j = 1:(i-1))
                aijkj = aijkj + k(:,j)*a(astart + j)
            end do
            
            k(:,i) = pDiffEqSys%F(X0 + h*c(i), Y0 + h*aijkj, params)
            end block
        end do
        
        FCalls = s - sint

        !> \remark Dense output coefficients are coefficients of theta interpolating polynomial, where
        !! \f$0\,<\,\theta = (u_{desired} - x_0)/h\,<\,1\f$. These are computed using the contributions from 
        !! all the integration stages. An 8th order interpolating polynomial is
        !! \f[ u_{desired}(\theta) = Y0 + Bip(1)*\theta + Bip(2)*\theta^2 + ...+ Bip(8)*\theta^8 \f],
        !! where \f$ Bip(i) = d(i,1)*k_1 + d(i,2)*k_2 + ... + d(i,s)*k_s \f$.
        
        IpCoeff(:,0) = Y0(InterpStates)
        do concurrent (j = 1:pstar)
            block
                integer :: istage
                real(WP), dimension(size(InterpStates)) :: cont
                
                cont = 0.0_WP
        sloop: do i = 1,s
                    ! get the current stage that gets multiplied by non-zero dij
                    istage = dinz(i)
                    ! -1 means we are at the end of the list
                    if (istage == -1) exit sloop
                    
                    cont = cont + d(i,j)*k(InterpStates,istage)
                end do sloop
                IpCoeff(:,j) = h*cont

            end block
        end do
    end subroutine IntpCoeff    


    !> It advances integrator by 1 step by computing stages
    subroutine DOP853_stepint(pDiffEqSys, X0, Y0, h, k, Y1, Yint12, FCalls, Err, &
                                    IsScalarTol, RTol, ATol, params)
        class(DiffEqSys), intent(in)  :: pDiffEqSys
        real(WP), intent(in) :: X0
        real(WP), dimension(pDiffEqSys%n), intent(in) :: Y0
        real(WP), intent(in) :: h
        real(WP), dimension(size(Y0),1:DOP853_sint), intent(inout) :: k
        real(WP), dimension(size(Y0)), intent(out) :: Y1, Yint12
        integer, intent(out) :: FCalls
        real(WP), intent(out) :: Err
        logical, intent(in), optional :: IsScalarTol
        real(WP), dimension(size(Y0)), intent(in), optional :: RTol, ATol
        real(WP), dimension(:), intent(in), optional :: params
                
        real(WP), dimension(size(Y0)) :: Yint, aijkj, Sc
        integer :: i
        logical EstimateErr        
        real(WP) :: DenomErr, Err3
        
        EstimateErr = .FALSE.
        if (present(IsScalarTol) .AND. present(RTol) &
            .AND. present(ATol)) EstimateErr = .TRUE.


        ! first stage
        ! k(:,1) = F0  
        
        ! 2nd stage
        aijkj = k(:,1)*DOP853_a(1)           
        Yint = Y0 + h*aijkj
        k(:,2) = pDiffEqSys%F(X0 + h*DOP853_c(2), Yint, params)
        
        ! 3rd stage
        aijkj = k(:,1)*DOP853_a(2) + k(:,2)*DOP853_a(3)           
        Yint = Y0 + h*aijkj        
        k(:,3) = pDiffEqSys%F(X0 + h*DOP853_c(3), Yint, params)
        
        ! 4th stage
        aijkj = k(:,1)*DOP853_a(4) + k(:,3)*DOP853_a(6)    
        Yint = Y0 + h*aijkj        
        k(:,4) = pDiffEqSys%F(X0 + h*DOP853_c(4), Yint, params)
        
        ! 5th stage
        aijkj = k(:,1)*DOP853_a(7) + k(:,3)*DOP853_a(9) + k(:,4)*DOP853_a(10)           
        Yint = Y0 + h*aijkj
        k(:,5) = pDiffEqSys%F(X0 + h*DOP853_c(5), Yint, params)
        
        ! 6th stage
        aijkj = k(:,1)*DOP853_a(11) + k(:,4)*DOP853_a(14) &
                + k(:,5)*DOP853_a(15)
        Yint = Y0 + h*aijkj        
        k(:,6) = pDiffEqSys%F(X0 + h*DOP853_c(6), Yint, params)
        
        ! 7th stage
        aijkj = k(:,1)*DOP853_a(16) + k(:,4)*DOP853_a(19) &
                + k(:,5)*DOP853_a(20) + k(:,6)*DOP853_a(21)
        Yint = Y0 + h*aijkj        
        k(:,7) = pDiffEqSys%F(X0 + h*DOP853_c(7), Yint, params)
        
        ! 8th stage
        aijkj = k(:,1)*DOP853_a(22) + k(:,4)*DOP853_a(25) &
                + k(:,5)*DOP853_a(26) + k(:,6)*DOP853_a(27) +  k(:,7)*DOP853_a(28)
        Yint = Y0 + h*aijkj        
        k(:,8) = pDiffEqSys%F(X0 + h*DOP853_c(8), Yint, params)
        
        ! 9th stage
        aijkj = k(:,1)*DOP853_a(29) + k(:,4)*DOP853_a(32) &
                + k(:,5)*DOP853_a(33) + k(:,6)*DOP853_a(34) +  k(:,7)*DOP853_a(35) +  k(:,8)*DOP853_a(36)
        Yint = Y0 + h*aijkj        
        k(:,9) = pDiffEqSys%F(X0 + h*DOP853_c(9), Yint, params)
        
        ! 10th stage
        aijkj = k(:,1)*DOP853_a(37) + k(:,4)*DOP853_a(40) &
                + k(:,5)*DOP853_a(41) + k(:,6)*DOP853_a(42) +  k(:,7)*DOP853_a(43) +  k(:,8)*DOP853_a(44) &
                +  k(:,9)*DOP853_a(45)
        Yint = Y0 + h*aijkj        
        k(:,10) = pDiffEqSys%F(X0 + h*DOP853_c(10), Yint, params)
        
        ! 11th stage
        aijkj = k(:,1)*DOP853_a(46) + k(:,4)*DOP853_a(49) &
                + k(:,5)*DOP853_a(50) + k(:,6)*DOP853_a(51) +  k(:,7)*DOP853_a(52) +  k(:,8)*DOP853_a(53) &
                + k(:,9)*DOP853_a(54) + k(:,10)*DOP853_a(55)
        Yint = Y0 + h*aijkj        
        k(:,11) = pDiffEqSys%F(X0 + h*DOP853_c(11), Yint, params)
        
        ! 12th stage
        aijkj = k(:,1)*DOP853_a(56)  + k(:,4)*DOP853_a(59) &
                + k(:,5)*DOP853_a(60) + k(:,6)*DOP853_a(61) +  k(:,7)*DOP853_a(62) +  k(:,8)*DOP853_a(63) &
                + k(:,9)*DOP853_a(64) +  k(:,10)*DOP853_a(65) +  k(:,11)*DOP853_a(66)
        Yint12 = Y0 + h*aijkj        
        k(:,12) = pDiffEqSys%F(X0 + h*DOP853_c(12), Yint12, params)
        
        ! 13th stage
        aijkj = k(:,1)*DOP853_a(67)  &
                + k(:,6)*DOP853_a(72) +  k(:,7)*DOP853_a(73) +  k(:,8)*DOP853_a(74) &
                + k(:,9)*DOP853_a(75) +  k(:,10)*DOP853_a(76) +  k(:,11)*DOP853_a(77) +  k(:,12)*DOP853_a(78)
        Yint = Y0 + h*aijkj
        
        ! propagate the solution to the next step for error computation
        ! For FSAL, Yint is the new solution
        
        if (EstimateErr) then        
           ! The scale factor for error computations
           if (IsScalarTol .EQV. .FALSE.) then
               Sc = ATol + RTol*max(abs(Y0), abs(Yint))
           else
               Sc = ATol(1) + RTol(1)*max(abs(Y0), abs(Yint))
           end if 
        
           ! Estimate the error. We assume that the error coeffcients are precomputed,
           ! i.e., e_i = b_i - bhat_i
        
           Err = sum(((k(:,1)*DOP853_e(1) + k(:,6)*DOP853_e(6) + k(:,7)*DOP853_e(7) + k(:,8)*DOP853_e(8) &
               + k(:,9)*DOP853_e(9) + k(:,10)*DOP853_e(10) + k(:,11)*DOP853_e(11) + k(:,12)*DOP853_e(12))/Sc)**2)
        
           ! For DOP853, we need to apply a 3rd-order correction
           Err3 = sum(((aijkj-DOP853_bhh(1)*k(:,1)-DOP853_bhh(2)*k(:,9) &
                                - DOP853_bhh(3)*k(:,12))/Sc)**2)
        
           DenomErr = Err + 0.01_WP*Err3
           if (DenomErr == 0.0_WP) DenomErr = 1.0_WP
           Err = abs(h)*Err*sqrt(1.0/(pDiffEqSys%n*DenomErr))
        else
           Err = 0.0_WP
        end if
        
        
        ! If this step is accepted then return the new solution
        ! else just return the initial condition
        if (Err <= 1.0_WP) then 
            Y1 = Yint
            ! In FSAL methods, the intgeration last stage is F0 for the next step
            k(:,13) = pDiffEqSys%F(X0 + h, Yint, params)
            k(:,1) = k(:,13)
            ! Function calls made if the step is accepted
            FCalls = DOP853_sint - 1    
        else
            Y1 = Y0
            ! Function calls made if the step is rejected
            FCalls = DOP853_sint - 2            
        end if
    end subroutine DOP853_stepint

    
    subroutine DOP853_IntpCoeff(pDiffEqSys, X0, Y0, h, Y1, InterpStates, k, IntpCoeff, Fcalls, params)
        class(DiffEqSys), intent(in)  :: pDiffEqSys
        real(WP), intent(in) :: X0
        real(WP), dimension(pDiffEqSys%n), intent(in) :: Y0
        real(WP), intent(in) :: h
        real(WP), dimension(size(Y0)), intent(in) :: Y1
        integer, dimension(:), intent(in) :: InterpStates
        real(WP), dimension(:,:), intent(inout) :: k
        real(WP), dimension(size(InterpStates),0:7), intent(out) :: IntpCoeff
        integer, intent(out) :: Fcalls
        real(WP), dimension(:), intent(in), optional :: params

        real(WP), dimension(pDiffEqSys%n) :: Yint, aijkj
        real(WP), dimension(size(InterpStates)) :: cont0, cont1, cont2, cont3, cont4, cont5, cont6, cont7
        
        ! Compute extra stages needed for interpolation
    
        ! 14th stage
        aijkj = k(:,1)*DOP853_a(79) +  k(:,7)*DOP853_a(85) +  k(:,8)*DOP853_a(86) + k(:,9)*DOP853_a(87) &
            + k(:,10)*DOP853_a(88) +  k(:,11)*DOP853_a(89) +  k(:,12)*DOP853_a(90) +  k(:,13)*DOP853_a(91)
        Yint = Y0 + h*aijkj
        k(:,14) = pDiffEqSys%F(X0 + h*DOP853_c(14), Yint, params)
        
        ! 15th stage
        aijkj = k(:,1)*DOP853_a(92) +  k(:,6)*DOP853_a(97) +  k(:,7)*DOP853_a(98) + k(:,8)*DOP853_a(99) &
            + k(:,11)*DOP853_a(102) +  k(:,12)*DOP853_a(103) +  k(:,13)*DOP853_a(104) +  k(:,14)*DOP853_a(105)
        Yint = Y0 + h*aijkj
        k(:,15) = pDiffEqSys%F(X0 + h*DOP853_c(15), Yint, params)
        
        ! 16th stage
        aijkj = k(:,1)*DOP853_a(106) +  k(:,6)*DOP853_a(111) +  k(:,7)*DOP853_a(112) + k(:,8)*DOP853_a(113) &
            + k(:,9)*DOP853_a(114) +  k(:,13)*DOP853_a(118) +  k(:,14)*DOP853_a(119) +  k(:,15)*DOP853_a(120)
        Yint = Y0 + h*aijkj
        k(:,16) = pDiffEqSys%F(X0 + h*DOP853_c(16), Yint, params)
        
        FCalls = 3
        
        ! DOP853 dense output coefficients (see Hairer's code)
        
        !> \remark Dense output coefficients are coefficients of theta interpolating polynomial, where
        !! \f$0\,<\,\theta = (u_{desired} - x_0)/h\,<\,1\f$. These are computed using the contributions from 
        !! all the integration stages. A 7th order interpolating polynomial is
        !! \f[ u_{desired}(\theta) = Bip(0) + Bip(1)*\theta + Bip(2)*\theta^2 + ...+ Bip(7)*\theta^7 \f],
        !! where \f$ Bip(i) = d(i,1)*k_1 + d(i,2)*k_2 + ... + d(i,s)*k_s \f$.
        
        cont0(:) = Y0(InterpStates)
            
        cont1(:) = Y1(InterpStates) - Y0(InterpStates)
            
        cont2(:) = h*k(InterpStates,1) - cont1(:)
            
        cont3(:) = cont1(:) - h*k(InterpStates,13) - cont2(:)
            
        cont4(:) = h*(DOP853_d(1,4)*k(InterpStates,1) + DOP853_d(6-4,4)*k(InterpStates,6) &
                    + DOP853_d(7-4,4)*k(InterpStates,7) + DOP853_d(8-4,4)*k(InterpStates,8) &
                    + DOP853_d(9-4,4)*k(InterpStates,9) + DOP853_d(10-4,4)*k(InterpStates,10) &
                    + DOP853_d(11-4,4)*k(InterpStates,11) + DOP853_d(12-4,4)*k(InterpStates,12) &
                    + DOP853_d(13-4,4)*k(InterpStates,13) + DOP853_d(14-4,4)*k(InterpStates,14) &
                    + DOP853_d(15-4,4)*k(InterpStates,15) + DOP853_d(16-4,4)*k(InterpStates,16))
            
        cont5(:) = h*(DOP853_d(1,5)*k(InterpStates,1) + DOP853_d(6-4,5)*k(InterpStates,6) + DOP853_d(7-4,5)*k(InterpStates,7) &
                    + DOP853_d(8-4,5)*k(InterpStates,8) + DOP853_d(9-4,5)*k(InterpStates,9) + DOP853_d(10-4,5)*k(InterpStates,10) &
                    + DOP853_d(11-4,5)*k(InterpStates,11) + DOP853_d(12-4,5)*k(InterpStates,12) + DOP853_d(13-4,5)*k(InterpStates,13) &
                    + DOP853_d(14-4,5)*k(InterpStates,14) + DOP853_d(15-4,5)*k(InterpStates,15) + DOP853_d(16-4,5)*k(InterpStates,16))
            
        cont6(:) = h*(DOP853_d(1,6)*k(InterpStates,1) + DOP853_d(6-4,6)*k(InterpStates,6) + DOP853_d(7-4,6)*k(InterpStates,7) &
                    + DOP853_d(8-4,6)*k(InterpStates,8) + DOP853_d(9-4,6)*k(InterpStates,9) + DOP853_d(10-4,6)*k(InterpStates,10) &
                    + DOP853_d(11-4,6)*k(InterpStates,11) + DOP853_d(12-4,6)*k(InterpStates,12) + DOP853_d(13-4,6)*k(InterpStates,13) &
                    + DOP853_d(14-4,6)*k(InterpStates,14) + DOP853_d(15-4,6)*k(InterpStates,15) + DOP853_d(16-4,6)*k(InterpStates,16))
        
        cont7(:) = h*(DOP853_d(1,7)*k(InterpStates,1) + DOP853_d(6-4,7)*k(InterpStates,6) + DOP853_d(7-4,7)*k(InterpStates,7) &
                    + DOP853_d(8-4,7)*k(InterpStates,8) + DOP853_d(9-4,7)*k(InterpStates,9) + DOP853_d(10-4,7)*k(InterpStates,10) &
                    + DOP853_d(11-4,7)*k(InterpStates,11) + DOP853_d(12-4,7)*k(InterpStates,12) + DOP853_d(13-4,7)*k(InterpStates,13) &
                    + DOP853_d(14-4,7)*k(InterpStates,14) + DOP853_d(15-4,7)*k(InterpStates,15) + DOP853_d(16-4,7)*k(InterpStates,16))
            
        ! check Hairer's contd8 routine for verifying the following code
        IntpCoeff(:,0) = cont0(:)
        IntpCoeff(:,1) = cont1(:) + cont2(:)
        IntpCoeff(:,2) = cont3(:) + cont4(:) - cont2(:)
        IntpCoeff(:,3) = cont5(:) + cont6(:) - 2*cont4(:) - cont3(:)
        IntpCoeff(:,4) = cont7(:) - 3*cont6(:) - 2*cont5(:) + cont4(:)
        IntpCoeff(:,5) = -3*cont7(:) + 3*cont6(:) + cont5(:)
        IntpCoeff(:,6) = 3*cont7(:) - cont6(:)
        IntpCoeff(:,7) = -cont7(:)
    end subroutine DOP853_IntpCoeff
 


    !> It advances integrator by 1 step by computing stages
    subroutine DOP54_stepint(pDiffEqSys, X0, Y0, h, k, Y1, Yint12, FCalls, Err, &
                                    IsScalarTol, RTol, ATol, params) 
        class(DiffEqSys), intent(in)  :: pDiffEqSys
        real(WP), intent(in) :: X0
        real(WP), dimension(pDiffEqSys%n), intent(in) :: Y0
        real(WP), intent(in) :: h
        real(WP), dimension(size(Y0),1:DOP54_sint), intent(inout) :: k
        real(WP), dimension(size(Y0)), intent(out) :: Y1, Yint12
        integer, intent(out) :: FCalls
        real(WP), intent(out) :: Err
        logical, intent(in), optional :: IsScalarTol
        real(WP), dimension(size(Y0)), intent(in), optional :: RTol, ATol
        real(WP), dimension(:), intent(in), optional :: params


        real(WP), dimension(size(Y0)) :: Yint, Sc
        logical :: EstimateErr

        EstimateErr = .FALSE.
        if (present(IsScalarTol) .AND. present(RTol) &
            .AND. present(ATol)) EstimateErr = .TRUE.

        ! first stage: no computation needed
        ! k(:,1) is k(:,7) from the last step
        
        ! 2nd stage
        Yint = Y0 + h*(k(:,1)*DOP54_a(1))
        k(:,2) = pDiffEqSys%F(X0 + h*DOP54_c(2), Yint, params)
        
        ! 3rd stage
        Yint = Y0 + h*(k(:,1)*DOP54_a(2) + k(:,2)*DOP54_a(3))        
        k(:,3) = pDiffEqSys%F(X0 + h*DOP54_c(3), Yint, params)
        
        ! 4th stage
        Yint = Y0 + h*(k(:,1)*DOP54_a(4) + k(:,2)*DOP54_a(5) + k(:,3)*DOP54_a(6))        
        k(:,4) = pDiffEqSys%F(X0 + h*DOP54_c(4), Yint, params)
        
        ! 5th stage
        Yint = Y0 + h*(k(:,1)*DOP54_a(7) + k(:,2)*DOP54_a(8) + k(:,3)*DOP54_a(9) + k(:,4)*DOP54_a(10))
        k(:,5) = pDiffEqSys%F(X0 + h*DOP54_c(5), Yint, params)
        
        ! 6th stage
        Yint12 = Y0 + h*(k(:,1)*DOP54_a(11) + k(:,2)*DOP54_a(12) + k(:,3)*DOP54_a(13) &
        +  k(:,4)*DOP54_a(14) + k(:,5)*DOP54_a(15))                
        k(:,6) = pDiffEqSys%F(X0 + h*DOP54_c(6), Yint12, params)
        
        ! 7th stage
        Yint = Y0 + h*(k(:,1)*DOP54_a(16) + k(:,2)*DOP54_a(17) + k(:,3)*DOP54_a(18) &
        + k(:,4)*DOP54_a(19) + k(:,5)*DOP54_a(20) + k(:,6)*DOP54_a(21))  
        k(:,7) = pDiffEqSys%F(X0 + h, Yint, params)
        
        ! propagate the solution to the next step for error computation
        ! For FSAL, Yint is the new solution
        
        if (EstimateErr) then
            ! The scale factor for error computations
            if (IsScalarTol .EQV. .FALSE.) then
                Sc = ATol + RTol*max(abs(Y0), abs(Yint))
            else
                Sc = ATol(1) + RTol(1)*max(abs(Y0), abs(Yint))
            end if 
            
            ! Estimate the error. We assume that the error coeffcients are precomputed,
            ! i.e., e_i = b_i - bhat_i
            Err = sqrt(sum(((h*(k(:,1)*DOP54_e(1) + k(:,3)*DOP54_e(3) + k(:,4)*DOP54_e(4) &
                                + k(:,5)*DOP54_e(5) + k(:,6)*DOP54_e(6) + k(:,7)*DOP54_e(7)))/Sc)**2)/pDiffEqSys%n)
        else
            Err = 0.0_WP
        end if
        !
        ! If this step is accepted then return the new solution
        ! else just return the initial condition
        if (Err <= 1.0_WP) then 
            Y1 = Yint
            ! In FSAL methods, the intgeration first stage is same as the last stage
            k(:,1) = k(:,7)
            ! Function calls made if the step is accepted
            FCalls = DOP54_sint - 1    
        else
            Y1 = Y0
            ! Function calls made if the step is rejected
            FCalls = DOP54_sint - 2            
        end if
    end subroutine DOP54_stepint
    
    
    subroutine DOP54_IntpCoeff(pDiffEqSys, X0, Y0, h, Y1, InterpStates, k, IntpCoeff, FCalls, params)
        class(DiffEqSys), intent(in)  :: pDiffEqSys
        real(WP), intent(in) :: X0
        real(WP), dimension(pDiffEqSys%n), intent(in) :: Y0
        real(WP), intent(in) :: h
        real(WP), dimension(size(Y0)), intent(in) :: Y1
        integer, dimension(:), intent(in) :: InterpStates
        real(WP), dimension(size(Y0),1:DOP54_s), intent(inout) :: k
        real(WP), dimension(size(InterpStates),0:4), intent(out) :: IntpCoeff
        integer, intent(out) :: Fcalls
        real(WP), dimension(:), intent(in), optional :: params
        
        real(WP), dimension(size(Y0)) :: Yint
        real(WP), dimension(size(InterpStates)) :: cont0, cont1, cont2, cont3, cont4
        
        FCalls = 0
        
        ! DOPRI54 dense output coefficients (see Hairer's DOPRI5 code)

        !> \remark Dense output coefficients are coefficients of theta interpolating polynomial, where
        !! \f$0\,<\,\theta = (u_{desired} - x_0)/h\,<\,1\f$. These are computed using the contributions from 
        !! all the integration stages. A 7th order interpolating polynomial is
        !! \f[ u_{desired}(\theta) = Bip(0) + Bip(1)*\theta + Bip(2)*\theta^2 + ...+ Bip(7)*\theta^7 \f],
        !! where \f$ Bip(i) = d(i,1)*k_1 + d(i,2)*k_2 + ... + d(i,s)*k_s \f$.

        cont0(:) = Y0(InterpStates)
            
        cont1(:) = Y1(InterpStates) - Y0(InterpStates)
            
        cont2(:) = h*k(InterpStates,1) - cont1
            
        cont3(:) = cont1 - h*k(InterpStates,7) - cont2
            
        cont4(:) = h*(DOP54_d(1,1)*k(InterpStates,1) + DOP54_d(3-1,1)*k(InterpStates,3) & 
                    + DOP54_d(4-1,1)*k(InterpStates,4) + DOP54_d(5-1,1)*k(InterpStates,5) &
                    + DOP54_d(6-1,1)*k(InterpStates,6) + DOP54_d(7-1,1)*k(InterpStates,7))
            
        ! check Hairer's contd5 routine for verifying the following code
        IntpCoeff(:,0) = cont0
        IntpCoeff(:,1) = cont1 + cont2
        IntpCoeff(:,2) = cont3 + cont4 - cont2
        IntpCoeff(:,3) = -2.0_WP*cont4 - cont3
        IntpCoeff(:,4) = cont4

    end subroutine DOP54_IntpCoeff    
    
    
    
        
end submodule ERKStepInt


