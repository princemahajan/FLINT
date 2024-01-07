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
!> \brief       FLINTUtils Module
!! \details     This module provides common utility functions.
!! \author      Bharat Mahajan
!! \date        Created: 04/25/2019    
!
!############################################################################################
 
    
!> Module containing all the utility functions
module FLINTUtils

    ! import the fortran environment module for precision-related constants
    use, intrinsic :: iso_fortran_env, only: real32, real64, real128

    !> By default, use 64-bit double precision (IEEE-754)
    integer, parameter, public :: FLINT_SP = real32   !> Single precision kind parameter
    integer, parameter, public :: FLINT_DP = real64   !> Double precision kind parameter
    integer, parameter, public :: FLINT_QP = real128  !> Quadruple precision kind parameter
    
    ! conditional compilation for choosing precision (default 64 bits)
#ifdef _REAL32
    integer, parameter, public :: FLINT_WP = FLINT_SP
#elif _REAL128
    integer, parameter, public :: FLINT_WP = FLINT_QP
#else
    integer, parameter, public :: FLINT_WP = FLINT_DP
#endif

    ! This is for use only in this module
    integer, parameter, private :: WP = FLINT_WP       !> choose the precision kind

    !> Smallest positive real satisfying 1.0_WP + eps > 1.0_WP
    real(WP), parameter, public :: FLINT_EPS = epsilon(1.0_WP)
    
    !> Max iterations for the root-finding
    integer, parameter :: FLINT_ROOTFIND_MAXITR = 50
    
    !> Interface of the function whose root needs to be computed
    abstract interface
        function func(x) result(y)
            import :: WP
            implicit none
            real(WP), intent(in) :: x !< independent variable
            real(WP) :: y
        end function func
    end interface
    
    contains
    
    !> \brief Subroutine for computing the root using Brent algorithm.
    !! \details For details on the algorithm, see Numerical receipes in C book
    !! at https://www2.units.it/ipl/students_area/imm2/files/Numerical_Recipes.pdf
    subroutine Root(a, b, ya, yb, tol, f, x, fx, niter, error)
    
        implicit none
        
        real(WP), intent(in)                 :: a     !< Lower bound on the root
        real(WP), intent(in)                 :: b     !< Upper bound on the root
        real(WP), intent(in)                 :: ya    !< f(a)
        real(WP), intent(in)                 :: yb    !< f(b)
        real(WP), intent(in)                 :: tol   !< tolerance
        procedure(func)                      :: f     !< function
        real(WP), intent(out)                :: x     !< root
        real(WP), intent(out)                :: fx    !< value at root
        integer, intent(out)                 :: niter !< no. of iterations
        integer, intent(out)                 :: error !< error status

        real(WP) :: aa, bb, cc, fa, fb, fc, xm, ss, dd, ee, pp, qq, rr
        
        real(WP) :: tolx
        
        real(WP), parameter :: NEARZERO = 1.0e-20_WP
        
        logical :: done
        integer :: itr
        
        aa = a
        bb = b
        fa = ya
        fb = yb

        error = 0
        niter = 0

        !! test
        !fa = f(a)
        !fb = f(b)
        
        ! Trivial cases
        if (abs(fa) <= NEARZERO) then
            x = aa
            fx = fa
            return
        else if (abs(fb) <= NEARZERO) then
            x = bb
            fx = fb
            return
        end if
        
        ! if root is not bracketed then return error
        if ((fa>0.0_WP .AND. fb>0.0_WP) .OR. (fa<0.0_WP .AND. fb<0.0_WP)) then
            error = 1
            niter = 0
            return
        end if
        
        ! algorithm start here
        fc = fb
        done = .FALSE.
        itr = 0
        do while (done .EQV. .FALSE. .AND. itr < FLINT_ROOTFIND_MAXITR)
            
            ! if root is NOT bracketed between C and B, then C = A
            if ((fc>0.0_WP .AND. fb>0.0_WP) .OR. (fc<0.0_WP .AND. fb<0.0_WP)) then
                cc = aa
                fc = fa
                dd = bb - aa
                ee = dd
            end if
            ! Make sure |f(c)| is bigger than |f(b)|
            if (abs(fc) < abs(fb)) then
                aa = bb
                bb = cc
                cc = aa
                fa = fb
                fb = fc
                fc = fa
            end if
            ! tolerance on the independent variable
            tolx = 2.0_WP*FLINT_EPS*abs(bb) + 0.5_WP*tol
            ! middle point
            xm = 0.5_WP*(cc - bb)
            
            if (abs(xm) <= tolx .OR. abs(fa) < NEARZERO) then
                ! root has been found
                x = bb
                done = .TRUE.
                fx = f(x)
            else
                if (abs(ee) >= tolx .AND. abs(fa) > abs(fb)) then
                    ss = fb/fa ! make sure |ss| < 1
                    if (abs(aa - cc) < NEARZERO) then
                        ! linear interpolation
                        pp = 2.0_WP*xm*ss
                        qq = 1.0_WP - ss
                    else
                        ! Inverse quadratic interpolation
                        qq = fa/fc
                        rr = fb/fc
                        pp = ss*(2.0_WP*xm*qq*(qq-rr)-(bb-aa)*(rr-1.0_WP))
                        qq = (qq - 1.0_WP)*(rr - 1.0_WP)*(ss - 1.0_WP)
                    end if
                    
                    if (pp > NEARZERO) qq = -qq
                    pp = abs(pp)
                    if ((2.0_WP*pp) < min(3.0_WP*xm*qq-abs(tolx*qq), abs(ee*qq))) then
                        ! use linear or inverse quadratic interpolation
                        ee = dd
                        dd = pp/qq
                    else
                        ! use bisection
                        dd = xm
                        ee = dd
                    end if
                else
                    ! use bisection
                    dd = xm
                    ee = dd
                end if

                aa = bb
                fa = fb
                if (abs(dd) > tolx) then
                    bb = bb + dd
                else
                    ! Correction term too small, advance by at least tolx
                    if (xm > 0.0_WP) then
                        bb = bb + abs(tolx)
                    else
                        bb = bb - abs(tolx)
                    end if
                end if
                fb = f(bb)
                itr = itr + 1
            end if
        end do
        
        if (itr >= FLINT_ROOTFIND_MAXITR) error = 2
        niter = itr

    end subroutine Root

    
    
    
    !> \brief Performs binary search on a provided sorted array
    !! \details TBD
    pure function BinarySearch(SortedX, SortedDataArray)

        !import :: WP
    
        implicit none    
    
        real(WP), dimension(1:), intent(in)  :: SortedX    !< Element to find
        real(WP), dimension(:), intent(in)  :: SortedDataArray   !< Ascending-order sorted vector
        integer, dimension(size(SortedX)) :: BinarySearch

        integer :: lb, ub, a, b, m, n, ctr, xloc
        real(WP) :: x
        
        ! initialize
        lb = 1
        ub = size(SortedDataArray)
        m = lb
        
        ! for each element
        do ctr = 1, size(SortedX)
        
            ! element to search for
            x = SortedX(ctr)
            
            ! current interval. Note we assume that the array of the elements
            ! to search for is also a sorted array.
            a = m
            b = ub
            
            ! loop for finding the element
            do
                ! size of the current interval
                n = b - a + 1
        
                ! mid-point
                m = floor((n+1)/2.0_WP) + (a-1)
        
                if (x < SortedDataArray(m)) then    
                    ! search the lower half
                    b = m
                    cycle
                elseif (x > SortedDataArray(m) .AND. (m+1) /= b) then
                    ! search the upper half
                    a = m
                    cycle
                else
                    ! We either have found the exact match 
                    ! Or x lies in between Data(m) and Data(m+1)
                    xloc = m
                    exit
                end if        
            end do

            ! save the location of the current element
            BinarySearch(ctr) = xloc
        end do
    
    end function BinarySearch    
    
    
    
end module FLINTUtils
