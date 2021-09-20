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
!> \brief       ERK Interpolation Module
!! \details     Provides implementation for the dense output features.
!! \author      Bharat Mahajan
!! \date        Created: 02/04/2019    
!
!############################################################################################
    
submodule (ERK) ERKInterp

    contains    

    module subroutine erk_interp(me, Xarr, Yarr, SaveInterpCoeffs)
    
        implicit none
        
        class(ERK_class), intent(inout) :: me
        
        real(WP), dimension(1:), intent(in) :: Xarr
        real(WP), dimension(size(me%InterpStates),size(Xarr)), intent(out) :: Yarr
        logical, intent(in), optional :: SaveInterpCoeffs

        logical :: DeallocMem
        integer :: n, ctr, X0loc, X0start
        real(WP) :: h, X0
        real(WP) :: hSign
        
        ! if interpolation mode is not on then return with an error
        if (.NOT. me%InterpOn) then
            me%status = FLINT_ERROR_INTP_OFF
            return
        end if
        
        ! if user does not explicitly say to keep the internal memory for future interpolations
        ! then free this memory.
        DeallocMem = .TRUE.
        if (present(SaveInterpCoeffs)) DeallocMem = (.NOT. SaveInterpCoeffs)
        
        ! If empty array then just return without error
        n = size(Xarr)
        if (n < 1) return
        
        ! Xarr must be monotonically increasing or decreasing with no repetitions
        hSign = sign(1.0_WP, (me%Xint(me%AcceptedSteps+1)-me%Xint(1)))
        block
            real(WP), dimension(size(Xarr)-1) :: tXarr
            tXarr = Xarr(1:(n-1)) - Xarr(2:n)
            if (any((hSign*tXarr) >= 0.0_WP)) then
                me%status = FLINT_ERROR_INTP_ARRAY
                return            
            end if
        end block

        ! All Xarr elements must belong to [x0, xf] closed set but we need to
        ! check the first and last elements as this is a sorted array.
        if (hSign*(Xarr(1)-me%Xint(1)) < 0.0_WP &
                        .OR. hSign*(Xarr(n)-me%Xint(me%AcceptedSteps+1)) > 0.0_WP) then
            me%status = FLINT_ERROR_INTP_ARRAY
            return
        end if
                
        ! Binary search seems to be slower than the linear search. May be
        ! because Xarr is also a sorted array
        ! Xarrloc = BinarySearch(Xarr, me%Xint)
        
        ! Main Interpolation loop
        X0start = 1
        
        do concurrent (ctr = 1:n)
            
            ! Find X0loc s.t. Xarr(ctr) belongs to ( Xint(X0loc-1), Xint(X0loc) ]
            do while (X0start <= me%AcceptedSteps)
                X0start = X0start + 1                
                if (hSign*(me%Xint(X0start)-Xarr(ctr)) >= 0.0_WP) then
                    X0loc = X0start - 1
                    exit
                end if
            end do
        
            if (X0loc <= me%AcceptedSteps) then
                X0 = me%Xint(X0loc)
                h = me%Xint(X0loc+1)-me%Xint(X0loc)
            else
                ! There is a problem, we should never reach here!
                X0 = Xarr(ctr)
                h = 0.0_WP
            end if
            
            ! compute the interpolation polynomial  
            Yarr(:, ctr) = InterpY(size(me%InterpStates), Xarr(ctr),  X0, h, &
                                    me%pstar, me%Bip(:,:,X0loc))

            ! prepare for the next point
            X0start = X0loc      
        end do

        
        ! job done, deallocate memory
        if (DeallocMem .EQV. .TRUE.) then
            deallocate(me%Xint)
            !deallocate(me%Yint)
            deallocate(me%Bip)
        end if

        
    end subroutine erk_interp
    
    !> Function to interpolate the solution within integration steps
    pure module function InterpY(np, X, X0, h, pstar, bCoeffs) result(Yip)
        integer, intent(in) :: np
        real(WP), intent(in) :: X, X0, h
        integer, intent(in) :: pstar
        real(WP), dimension(np,0:pstar), intent(in) :: bCoeffs
        real(WP), dimension(np) :: Yip
            
        real(WP) :: theta
        integer :: k
            
        ! compute the interpolation polynomial           
        theta = (X - X0)/h
        Yip = 0.0_WP
        do k = 0,pstar
            Yip = Yip + bCoeffs(:,k)*(theta**k)
        end do
    end function InterpY
    
    
end submodule ERKInterp
