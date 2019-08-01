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
!> \brief       ERK Integration Submodule
!! \details     It provides the implementation for Integrate method of ERK_class.
!! \author      Bharat Mahajan
!! \date        02/06/2019    
!
!############################################################################################
    
submodule (ERK) ERKIntegrate

    contains    

    !> Integration main subroutine
    module subroutine erk_int(me, X0, Y0, Xf, Yf, StepSz, IntStepsOn, Xint, Yint, EventMask, EventStates, EventRootFindingOn, StiffTest, params)
    
        implicit none
        
        class(ERK_class), intent(inout) :: me
    
        real(WP), intent(in) :: X0
        real(WP), dimension(me%pDiffEqSys%n), intent(in) :: Y0
        real(WP), intent(inout) :: Xf            
        real(WP), dimension(me%pDiffEqSys%n), intent(out) :: Yf
        real(WP), intent(inout) :: StepSz
        logical, intent(in), optional :: IntStepsOn
        real(WP), allocatable, dimension(:), intent(out), optional :: Xint            
        real(WP), allocatable, dimension(:,:), intent(out), optional :: Yint
        logical, dimension(me%pDiffEqSys%m), intent(in), optional :: EventMask
        real(WP), allocatable, dimension(:,:), intent(out), optional :: EventStates
        logical, intent(in), optional :: EventRootFindingOn        
        integer, intent(inout), optional :: StiffTest        
        real(WP), dimension(:), intent(in), optional :: params                
    
        real(WP) :: h, hSign, hnew, hmax, X, Err, Eventh, EventX0, EventX1, EventStepSz
        real(WP), dimension(me%pDiffEqSys%n) :: Y1, Y2, F0, Sc0, Yint12

        !integer, dimension(:), allocatable :: IpStages
        real(WP), dimension(size(me%InterpStates), 0:me%pstar) :: Bip
        
        real(WP), dimension(me%pDiffEqSys%m) :: EV0, EV1
        integer, dimension(me%pDiffEqSys%m) :: EventDir
        logical, dimension(me%pDiffEqSys%m) :: EventTerm
        logical, dimension(me%pDiffEqSys%m) :: EvMask
        real(WP), allocatable, dimension(:) :: EventData
        real(WP), dimension(me%pDiffEqSys%n) :: EventY1, EventYm        
        real(WP), dimension(me%pDiffEqSys%n, me%s) :: k                
        real(WP), dimension(5) :: StepSzParams        

        integer :: n, m, FCalls, TotalFCalls, AcceptedSteps, RejectedSteps 
        integer :: TotalSteps, MaxSteps, s, sint, p, pstar, q, EventId
        integer :: StiffnessTest, StiffThreshold, NonStiffThreshold, StiffTestSteps
        
        logical :: IsProblemStiff, IntStepsNeeded, LastStepRejected, LAST_STEP
        logical :: IsScalarTol, IsFSALMethod, InterpOn, EventsOn, RootFindingOn, BipComputed
        
        class(DiffEqSys), pointer :: pDiffEqSys
        integer(kind(ERK_DOP853)) :: method

        integer(kind(FLINT_SUCCESS)) :: status, stat

        ! Init for this procedure
        IntStepsNeeded = .FALSE.
        LastStepRejected = .FALSE.
        LAST_STEP = .FALSE.
        status = FLINT_SUCCESS        
        
        ! Copy frequently needed scalars on the stack
        n = me%pDiffEqSys%n
        m = me%pDiffEqSys%m
        p = me%p
        pstar = me%pstar
        q = me%q
        pDiffEqSys => me%pDiffEqSys
        MaxSteps = me%MaxSteps
        hmax = me%MaxStepSize
        StepSzParams = me%StepSzParams        
        method = me%method
        IsScalarTol = me%IsScalarTol
        IsFSALMethod = me%IsFSALMethod
        InterpOn = me%InterpOn
        EventsOn = me%EventsOn
        EventStepSz = me%EventStepSz
        s = me%s
        sint = me%sint
        StiffTestSteps = STIFFTEST_STEPS
        
        ! If Interpolation is OFF and solution at integrator steps is needed,
        ! then allocate the storage for steps that will be returned to the user. 
        ! If Interpolation is ON, then no need as we already allocated the internal
        ! storage for the solution and we just use allocatable array assignment
        if (present(IntStepsOn)) then
            if (IntStepsOn .AND. InterpOn) then
                status = FLINT_ERROR_PARAMS
            else if (IntStepsOn .AND. (.NOT. InterpOn)) then
                IntStepsNeeded = .TRUE.
            end if            
        end if

        if (IntStepsNeeded) then
            allocate(Xint(MaxSteps), stat=stat)
            if (stat /= 0)  status = FLINT_ERROR_MEMALLOC
            allocate(Yint(pDiffEqSys%n,MaxSteps), stat=stat)
            if (stat /= 0)  status = FLINT_ERROR_MEMALLOC
        end if
        
        ! stiffness test options
        if (present(StiffTest)) then
            StiffnessTest = StiffTest
            select case (StiffnessTest)
                case (0)
                case (1)
                case (2)
                case default
                status = FLINT_ERROR_PARAMS
                end select
        else
            StiffnessTest = 0
        end if
        
        ! make sure memory is allocated if interpolation is on
        if (InterpOn) then
            if (allocated(me%Xint) .EQV. .FALSE. .OR. allocated(me%Bip) .EQV. .FALSE.) then
                me%status = FLINT_ERROR_INIT_REQD
            end if
            ! If size of Xint is less than the max size, then it means user is calling
            ! Intgerate again without doing Init first, then reallocate the memory
            if (allocated(me%Xint)) then
                if (size(me%Xint) < MaxSteps) then 
                    deallocate(me%Xint, stat = status)
                    deallocate(me%Bip, stat = status)
                    allocate(me%Xint(me%MaxSteps), stat=status)
                    allocate(me%Bip(size(me%InterpStates), 0:me%pstar, me%MaxSteps), stat=status)
                    if (status /= 0) me%status = FLINT_ERROR_MEMALLOC
                end if
            end if
        end if        

        ! Make sure optional event arguments are present if events are ON.        
        if (EventsOn) then
            ! If Mask is not specified, then all events are enabled by default
            if (present(EventMask)) then
                EvMask = EventMask
            else
                EvMask = .TRUE.
            end if

            if (.NOT. any(EvMask)) then
                ! If all the event masks are false then turn off events for this integration                
                EventsOn = .FALSE.
            elseif (.NOT. present(EventStates)) then
                ! If user has not provided the allocatable array for storage, retrun error
                me%status = FLINT_ERROR_EVENTPARAMS
            else
                ! Allocate EventData to 0 size to define it
                allocate(EventData(0))
            end if
            
            ! Root-finding by default is enabled
            if (present(EventRootFindingOn)) then
                RootFindingOn = EventRootFindingOn
            else
                RootFindingOn = .TRUE.
            end if
            
        end if

        ! Return if any error up to this point
        if (status /= FLINT_SUCCESS) then
            me%status = status
            return
        end if

        ! All the parameter checking should be before this point
        
        ! save the initial condition
        me%X0 = X0
        me%Y0 = Y0
        
        ! compute the derivative and scale factor at the initial condition
        F0 = pDiffEqSys%F(X0, Y0, params)
        TotalFCalls = 1
        
        !! \remark Scale factor is computed using
        !! \f[ SC_i = atol_i + rtol_i |Y_{i}| \quad i=1..n,\f]
        !! where \f$|Y|\f$ is the absolute value of the provided vector
        if (IsScalarTol .EQV. .FALSE.) then
            Sc0 = me%ATol + me%RTol*abs(Y0)
        else
            Sc0 = me%ATol(1) + me%RTol(1)*abs(Y0)
        end if        
        
        ! Determine the sign of the step: negative for backward integration
        hSign = sign(1.0_WP, (Xf-X0))

        ! Select the max step size if not provided by the user
        if (hmax == 0.0_WP) hmax = abs(Xf - X0)

        ! The initial step size: compute if user has not provided a non-zero value
        if (StepSz == 0.0_WP) then
            call StepSz0Hairer(h, p, n, X0, Y0, F0, Sc0, hSign, hmax, FCalls, pDiffEqSys, params)
            TotalFCalls = TotalFCalls + Fcalls
        else
            h = StepSz
        end if

        ! copy the initial condition as the first state
        if (IntStepsNeeded) then
            Xint(1) = X0
            Yint(:,1) = Y0
        end if
        if (InterpOn) then
            me%Xint(1) = X0
            !me%Yint(:,1) = Y0
        end if
        
        ! Initialize the number of steps taken and other counters
        X = X0
        Y1 = Y0
        TotalSteps = 0
        AcceptedSteps = 0
        RejectedSteps = 0
        StiffThreshold = 0
        NonStiffThreshold = 0
        status = FLINT_SUCCESS
        
        ! Main loop for computing the steps
        do

            if (status == FLINT_SUCCESS .AND. (.NOT. LAST_STEP)    &
                    .AND. (TotalSteps <= MaxSteps)) then

                ! Adjust the step size if (X+h) is bigger than Xf
                if (hSign*(X + 1.01_WP*h - Xf) > 0.0) then
                    h = Xf - X 
                    LAST_STEP = .TRUE.
                end if

                ! Advance one step using the chosen method
                if (method == ERK_DOP853) then
                    call DOP853_stepint(X, Y1, Y2)
                elseif (method == ERK_DOP54) then
                    call DOP54_stepint(X, Y1, Y2)
                else
                    call stepint(X, Y1, Y2)
                end if

                ! update steps taken and function calls made
                TotalSteps = TotalSteps + 1                
                TotalFCalls = TotalFCalls + FCalls                
                
                ! check error to see if the step needs to be accepted or rejected
                if (Err <= 1.0_WP) then

                    ! step is accepted
                    AcceptedSteps = AcceptedSteps + 1
                    
                    ! save the first accepted step-size
                    if (AcceptedSteps == 1) me%h0 = h

                    ! do we need stiffness detection?
                    if (StiffnessTest == 1 .OR. StiffnessTest == 2) then
                        call CheckStiffness(IsProblemStiff)
                        if (IsProblemStiff) then
                            StiffTest = -1 
                            if (StiffnessTest == 1) then
                                status = FLINT_STIFF_PROBLEM
                                LAST_STEP = .TRUE.
                            end if
                        end if
                    end if
                    
                    ! Events checking start here, ignore event checking at IC
                    BipComputed = .FALSE.
                    if (EventsOn .AND. AcceptedSteps > 1) then
                        ! generate X values at which event triggers are checked
                        if (EventStepSz == 0.0 .OR. EventStepSz >= abs(h)) then
                            Eventh = h
                        else
                            Eventh = sign(EventStepSz, h)
                            ! we need interpolation coefficinets to check signs at the intermediate points
                            if (.NOT. BipComputed) then
                                select case (method)
                                    case (ERK_DOP853)
                                    Bip = DOP853_IntpCoeff(X, Y1, X+h, Y2, me%InterpStates)
                                    case (ERK_DOP54)
                                    Bip = DOP54_IntpCoeff(X, Y1, X+h, Y2, me%InterpStates)                                
                                    case default
                                    Bip = IntpCoeff(X, Y1, X+h, Y2, me%InterpStates)
                                end select 
                                TotalFCalls = TotalFCalls + FCalls
                                BipComputed = .TRUE.
                            end if
                        end if 
                            
                        EventX0 = X
                        EventX1 = X + Eventh
                        ! The event value at the first point where event is checked
                        if (AcceptedSteps == 2) call pDiffEqSys%G(0, X, Y1, EV0, EventDir, EventTerm)

                        ! check for event triggers between X and X+h
             EventLoop: do
                            block
                                logical :: EventTriggered
                                real(WP) :: val0, val1, valm, EventXm
                                integer :: rootiter, rooterror

                                !! interpolate solution
                                if (Eventh == h) then
                                    EventY1 = Y2
                                else
                                    EventY1 = InterpY(size(me%InterpStates), EventX1, X, h, me%pstar, Bip)
                                end if
                                ! evaluate all event functions
                                call pDiffEqSys%G(0, EventX1, EventY1, EV1, EventDir, EventTerm)
                                    
                                ! check each event for trigger
                                do  EventId = 1,m
                                    ! if event is masked, no event check
                                    if (EvMask(EventId) .EQV. .False.) cycle
                                        
                                    val0 = EV0(EventId)
                                    val1 = EV1(EventId)
                                        
                                    ! If the previous or the new value is NaN, ignore the event
                                    if (IEEE_IS_NAN(val0) .OR. IEEE_IS_NAN(val1)) cycle 
                                    
                                    EventTriggered = .FALSE.
                                    select case (EventDir(EventId))
                                    case (0)
                                        if ((val0<=0 .AND. val1>=0) &
                                            .OR.(val0>=0 .AND. val1<=0)) EventTriggered = .TRUE.
                                    case (1)
                                        if ((val0<=0 .AND. val1>=0)) EventTriggered = .TRUE.
                                    case (-1)
                                        if ((val0>=0 .AND. val1<=0)) EventTriggered = .TRUE.
                                    end select

                                    if (EventTriggered) then
                                        ! Current event has indeed triggered
                                        if (RootFindingOn) then
                                            ! compute interpolation coefficients if not already computed
                                            if (.NOT. BipComputed) then
                                                select case (method)
                                                case (ERK_DOP853)
                                                Bip = DOP853_IntpCoeff(X, Y1, X+h, Y2, me%InterpStates)
                                                case (ERK_DOP54)
                                                Bip = DOP54_IntpCoeff(X, Y1, X+h, Y2, me%InterpStates)                                
                                                case default
                                                Bip = IntpCoeff(X, Y1, X+h, Y2, me%InterpStates)
                                                end select    
                                                TotalFCalls = TotalFCalls + FCalls
                                                BipComputed = .TRUE.
                                            end if
                                            ! Now find the exact event location using root-finding
                                            call Root(EventX0, EventX1, val0, val1, DEFAULT_EVENTTOL, &
                                                    SingleEvent, EventXm, valm, rootiter, rooterror)
                                            if (rooterror == 0 .OR. rooterror == 2) then
                                                ! save the event states
                                                EventYm = InterpY(size(me%InterpStates), EventXm, X, h, me%pstar, Bip)
                                                EventData = [EventData, EventXm, EventYm, real(EventId, WP)]
                                            else
                                                ! root finding failed, save the event states at the point "a"
                                                status = FLINT_ERROR_EVENTROOT
                                            end if
                                        else
                                            ! just return the state at which sign-change is detected
                                            EventXm = EventX1
                                            EventYm = val1
                                            EventData = [EventData, EventXm, EventYm, real(EventId, WP)]
                                        end if
                                        
                                        ! Do we need to stop at this event?
                                        if (EventTerm(EventId) .EQV. .TRUE.) then
                                            ! the last solution is the same as the event state
                                            h = EventXm - X
                                            Y2 = EventYm
                                            LAST_STEP = .TRUE.
                                            status = FLINT_EVENT_TERM
                                            exit EventLoop
                                        end if
                                    end if
                                end do
                 
                                ! new check point for event
                                EventX0 = EventX1
                                EV0 = EV1
                                if (EventX1 == X + h) exit EventLoop ! we are done
                                if (EventX1 + Eventh > X + h) then
                                    EventX1 = X + h
                                else
                                    EventX1 = EventX1 + Eventh
                                end if
                            end block
                        end do EventLoop
                    end if

                    ! return the solution at integrator's natural step size
                    if (IntStepsNeeded) then
                        Xint(AcceptedSteps + 1) = X + h
                        Yint(:,AcceptedSteps + 1) = Y2
                    end if

                    ! Do we need to store interpolation coefficients?
                    if (InterpOn) then
                        if (.NOT. BipComputed) then
                            ! Compute interpolation coefficients
                            select case (method)
                            case (ERK_DOP853)
                            Bip = DOP853_IntpCoeff(X, Y1, X+h, Y2, me%InterpStates)
                            case (ERK_DOP54)
                            Bip = DOP54_IntpCoeff(X, Y1, X+h, Y2, me%InterpStates)                                
                            case default
                            Bip = IntpCoeff(X, Y1, X+h, Y2, me%InterpStates)
                            end select
                            BipComputed = .TRUE.
                            TotalFCalls = TotalFCalls + FCalls
                        end if
                        ! store the solution internally at the natural step size
                        me%Xint(AcceptedSteps + 1) = X + h
                        !me%Yint(:,AcceptedSteps + 1) = Y2
                        me%Bip(:,:,AcceptedSteps) = Bip
                    end if
                    
                    ! get ready for the next iteration
                    X = X + h
                    Y1 = Y2

                    ! compute the new step size now
                    hnew = StepSzHairer(h, .FALSE., q, Err, StepSzParams)

                    ! make sure hnew is not greater than hmax
                    if (abs(hnew) > hmax) hnew = hsign*hmax
                    
                    ! if the last step was rejected, then dont increase step size
                    if (LastStepRejected .EQV. .TRUE.) hnew = hSign*min(abs(hnew),abs(h))

                    ! Update the Lund stabilization parameter hbyhoptOld
                    ! this should be after computing the new step size
                    StepSzParams(5) = max(Err, hbyhoptOLD)
                    
                    ! reset
                    LastStepRejected = .FALSE.
                else

                    ! step is rejected, now what
                    
                    ! increment the reject step counter, ignore the initial step rejections
                    if (AcceptedSteps >= 1) RejectedSteps = RejectedSteps + 1
                    
                    ! compute the new step size now
                    hnew = StepSzHairer(h, .TRUE., q, Err, StepSzParams)

                    ! Make note of your rejection. It is used in the next iteration
                    ! for computing the new step size.
                    LastStepRejected = .TRUE.
                    
                    ! Retry again if this was the last step
                    LAST_STEP = .FALSE.
                    
                end if
                
                ! update the step size for next iteration
                h = hnew
                
                ! Take another step only if step size is not too small
                ! Reference: Hairer's DOP853 code
                if (abs(h)*0.01_WP <= abs(X)*EPS) then
                    status = FLINT_ERROR_STEPSZ_TOOSMALL
                end if

            else
                    
                ! We have either accomplished our job or finished prematurely?!

                ! return the solution at the integrator's natural steps
                if (IntStepsNeeded) then
                    block
                        real(WP), allocatable, dimension(:) :: Xarr            
                        real(WP), allocatable, dimension(:,:) :: Yarr       
                    
                        allocate(Xarr(AcceptedSteps + 1))
                        allocate(Yarr(n,AcceptedSteps + 1))
                    
                        Xarr = Xint(1:AcceptedSteps + 1)
                        Yarr = Yint(:,1:AcceptedSteps + 1)
                    
                        call move_alloc(from=Xarr, to=Xint)
                        call move_alloc(from=Yarr, to=Yint)
                    end block
                end if
                
                ! Update internal storage size based on accepted steps
                if (InterpOn) then
                    block
                        real(WP), allocatable, dimension(:) :: Xarr            
                        !real(WP), allocatable, dimension(:,:) :: Yarr       
                        real(WP), allocatable, dimension(:,:,:) :: Biparr
                        
                        allocate(Xarr(AcceptedSteps + 1))
                        !allocate(Yarr(n,AcceptedSteps + 1))
                        allocate(Biparr(size(me%InterpStates),0:me%pstar,AcceptedSteps))
                    
                        Xarr = me%Xint(1:AcceptedSteps + 1)
                        !Yarr = me%Yint(:,1:AcceptedSteps + 1)
                        Biparr = me%Bip(:,:,1:AcceptedSteps)
                        
                        call move_alloc(from=Xarr, to=me%Xint)
                        !call move_alloc(from=Yarr, to=me%Yint)
                        call move_alloc(from=Biparr, to=me%Bip)
                    end block
                end if        
                
                ! Prepare event output data
                if (EventsOn .AND. allocated(EventData)) then                
                    EventStates = reshape(EventData, [n+2, int(size(EventData)/(n+2))])
                end if
                
                if (X == Xf .AND. status == FLINT_SUCCESS) then
                    ! we finished like we should
                elseif (EventsOn .AND. status == FLINT_EVENT_TERM) then
                    ! One of the terminal event has triggered, do nothing
                else if (TotalSteps > MaxSteps .AND. status == FLINT_SUCCESS) then
                    ! maximum steps reached
                    status = FLINT_ERROR_MAXSTEPS
                else if (status == FLINT_ERROR_STEPSZ_TOOSMALL) then
                    ! do nothing
                else if (EventsOn .AND. status == FLINT_ERROR_EVENTROOT) then
                    ! do nothing
                else if (status == FLINT_STIFF_PROBLEM) then
                    ! do nothing
                else
                    ! Some unknown error (We should never be here!)
                    status = FLINT_ERROR 
                end if

                ! save the stats
                Xf = X
                Yf = Y1
                StepSz = h                
                me%Xf = Xf
                me%hf = StepSz
                me%Yf = Yf
                me%AcceptedSteps = AcceptedSteps
                me%RejectedSteps = RejectedSteps
                me%FCalls = TotalFCalls                
                me%status = status
                
                ! exit the main integration loop
                exit
                
            end if
            
        end do
        
        
    contains
    
    
    
    !> It advances integrator by 1 step by computing stages
    subroutine stepint(X0, Y0, Y1)
    
        implicit none
        
        real(WP), intent(in) :: X0
        real(WP), dimension(n), intent(in) :: Y0
        real(WP), dimension(n), intent(out) :: Y1
        
        real(WP), dimension(n) :: Yint, Sc, aijkj
        integer :: i
    
        ! first stage
        k(:,1) = F0
        
        ! compute rest of the stages
        ! This compact and beautiful code is slower than ugly hardcoded one!
        do i = 2,sint
            block
            integer :: astart, j
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
                aijkj = aijkj + k(:,j)*me%a(astart + j)
            end do
            
            Yint = Y0 + h*aijkj
            k(:,i) = pDiffEqSys%F(X0 + h*me%c(i), Yint, params)
            
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
            Yint = Y0 + h*matmul(k(:,1:sint), me%b(1:sint))
        end if

        ! The scale factor for error computations
        if (IsScalarTol .EQV. .FALSE.) then
            Sc = me%ATol + me%RTol*max(abs(Y0), abs(Yint))
        else
            Sc = me%ATol(1) + me%RTol(1)*max(abs(Y0), abs(Yint))
        end if 

        ! Estimate the error. We assume that the error coeffcients are precomputed,
        ! i.e., e_i = b_i - bhat_i
        block
            real(WP), dimension(n) :: temp
            integer :: j
            temp = 0.0_WP
            do concurrent (j = 1:sint)
                temp = temp + k(:,j)*me%e(j)
            end do
            Err = abs(h)*sqrt(sum((temp/Sc)**2)/n)
        end block
        
        ! If this step is accepted then return the new solution
        ! else just return the initial condition
        if (Err <= 1.0_WP) then 
            Y1 = Yint
            if (IsFSALMethod .EQV. .FALSE.) then
                ! Evaluate F0 for the next step
                F0 = pDiffEqSys%F(X0 + h, Y1, params)
                ! Number of function calls made so far
                FCalls = sint            
            else
                ! In FSAL methods, the intgeration last stage is F0 for the next step
                F0 = k(:,sint)
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
    

    !> It advances integrator by 1 step by computing stages
    subroutine DOP853_stepint(X0, Y0, Y1)
    
        implicit none
        
        real(WP), intent(in) :: X0
        real(WP), dimension(n), intent(in) :: Y0
        real(WP), dimension(n), intent(out) :: Y1
        
        real(WP), dimension(n) :: Yint, aijkj, Sc
        integer :: i
        
        real(WP) :: DenomErr, Err3
        
        ! first stage
        k(:,1) = F0  
        
        ! 2nd stage
        aijkj = k(:,1)*me%a(1)           
        Yint = Y0 + h*aijkj
        k(:,2) = pDiffEqSys%F(X0 + h*me%c(2), Yint, params)

        ! 3rd stage
        aijkj = k(:,1)*me%a(2) + k(:,2)*me%a(3)           
        Yint = Y0 + h*aijkj        
        k(:,3) = pDiffEqSys%F(X0 + h*me%c(3), Yint, params)
        
        ! 4th stage
        aijkj = k(:,1)*me%a(4) + k(:,3)*me%a(6)    
        Yint = Y0 + h*aijkj        
        k(:,4) = pDiffEqSys%F(X0 + h*me%c(4), Yint, params)
        
        ! 5th stage
        aijkj = k(:,1)*me%a(7) + k(:,3)*me%a(9) + k(:,4)*me%a(10)           
        Yint = Y0 + h*aijkj
        k(:,5) = pDiffEqSys%F(X0 + h*me%c(5), Yint, params)
        
        ! 6th stage
        aijkj = k(:,1)*me%a(11) + k(:,4)*me%a(14) &
                + k(:,5)*me%a(15)
        Yint = Y0 + h*aijkj        
        k(:,6) = pDiffEqSys%F(X0 + h*me%c(6), Yint, params)
        
        ! 7th stage
        aijkj = k(:,1)*me%a(16) + k(:,4)*me%a(19) &
                + k(:,5)*me%a(20) + k(:,6)*me%a(21)
        Yint = Y0 + h*aijkj        
        k(:,7) = pDiffEqSys%F(X0 + h*me%c(7), Yint, params)
        
        ! 8th stage
        aijkj = k(:,1)*me%a(22) + k(:,4)*me%a(25) &
                + k(:,5)*me%a(26) + k(:,6)*me%a(27) +  k(:,7)*me%a(28)
        Yint = Y0 + h*aijkj        
        k(:,8) = pDiffEqSys%F(X0 + h*me%c(8), Yint, params)
        
        ! 9th stage
        aijkj = k(:,1)*me%a(29) + k(:,4)*me%a(32) &
                + k(:,5)*me%a(33) + k(:,6)*me%a(34) +  k(:,7)*me%a(35) +  k(:,8)*me%a(36)
        Yint = Y0 + h*aijkj        
        k(:,9) = pDiffEqSys%F(X0 + h*me%c(9), Yint, params)
        
        ! 10th stage
        aijkj = k(:,1)*me%a(37) + k(:,4)*me%a(40) &
                + k(:,5)*me%a(41) + k(:,6)*me%a(42) +  k(:,7)*me%a(43) +  k(:,8)*me%a(44) &
                +  k(:,9)*me%a(45)
        Yint = Y0 + h*aijkj        
        k(:,10) = pDiffEqSys%F(X0 + h*me%c(10), Yint, params)
        
        ! 11th stage
        aijkj = k(:,1)*me%a(46) + k(:,4)*me%a(49) &
                + k(:,5)*me%a(50) + k(:,6)*me%a(51) +  k(:,7)*me%a(52) +  k(:,8)*me%a(53) &
                + k(:,9)*me%a(54) + k(:,10)*me%a(55)
        Yint = Y0 + h*aijkj        
        k(:,11) = pDiffEqSys%F(X0 + h*me%c(11), Yint, params)
        
        ! 12th stage
        aijkj = k(:,1)*me%a(56)  + k(:,4)*me%a(59) &
                + k(:,5)*me%a(60) + k(:,6)*me%a(61) +  k(:,7)*me%a(62) +  k(:,8)*me%a(63) &
                + k(:,9)*me%a(64) +  k(:,10)*me%a(65) +  k(:,11)*me%a(66)
        Yint12 = Y0 + h*aijkj        
        k(:,12) = pDiffEqSys%F(X0 + h*me%c(12), Yint12, params)

        ! 13th stage
        aijkj = k(:,1)*me%a(67)  &
                + k(:,6)*me%a(72) +  k(:,7)*me%a(73) +  k(:,8)*me%a(74) &
                + k(:,9)*me%a(75) +  k(:,10)*me%a(76) +  k(:,11)*me%a(77) +  k(:,12)*me%a(78)
        Yint = Y0 + h*aijkj

        ! propagate the solution to the next step for error computation
        ! For FSAL, Yint is the new solution

        ! The scale factor for error computations
        if (IsScalarTol .EQV. .FALSE.) then
            Sc = me%ATol + me%RTol*max(abs(Y0), abs(Yint))
        else
            Sc = me%ATol(1) + me%RTol(1)*max(abs(Y0), abs(Yint))
        end if 

        ! Estimate the error. We assume that the error coeffcients are precomputed,
        ! i.e., e_i = b_i - bhat_i
        
        Err = sum(((k(:,1)*me%e(1) + k(:,6)*me%e(6) + k(:,7)*me%e(7) + k(:,8)*me%e(8) &
                + k(:,9)*me%e(9) + k(:,10)*me%e(10) + k(:,11)*me%e(11) + k(:,12)*me%e(12))/Sc)**2)

        ! For DOP853, we need to apply a 3rd-order correction
        Err3 = sum(((aijkj-DOP853_bhh(1)*k(:,1)-DOP853_bhh(2)*k(:,9) &
                                 - DOP853_bhh(3)*k(:,12))/Sc)**2)

        DenomErr = Err + 0.01_WP*Err3
        if (DenomErr == 0.0_WP) DenomErr = 1.0_WP
        Err = abs(h)*Err*sqrt(1.0/(n*DenomErr))
        
        ! If this step is accepted then return the new solution
        ! else just return the initial condition
        if (Err <= 1.0_WP) then 
            Y1 = Yint
            ! In FSAL methods, the intgeration last stage is F0 for the next step
            k(:,13) = pDiffEqSys%F(X0 + h, Yint, params)
            F0 = k(:,13)
            ! Function calls made if the step is accepted
            FCalls = sint - 1    
        else
            Y1 = Y0
            ! Function calls made if the step is rejected
            FCalls = sint - 2            
        end if
        
    end subroutine DOP853_stepint

    
    function DOP853_IntpCoeff(X0, Y0, X1, Y1, InterpStates) result(IntpCoeff)

        implicit none
    
        real(WP), intent(in) :: X0, X1
        real(WP), dimension(n), intent(in) :: Y0, Y1
        integer, dimension(:), intent(in) :: InterpStates
        
        real(WP), dimension(size(InterpStates),0:7) :: IntpCoeff
        
        real(WP), dimension(n) :: Yint, aijkj
        real(WP), dimension(size(InterpStates)) :: cont0, cont1, cont2, cont3, cont4, cont5, cont6, cont7
        integer :: i, j
        
        ! Compute extra stages needed for interpolation
    
        ! 14th stage
        aijkj = k(:,1)*me%a(79) +  k(:,7)*me%a(85) +  k(:,8)*me%a(86) + k(:,9)*me%a(87) &
            + k(:,10)*me%a(88) +  k(:,11)*me%a(89) +  k(:,12)*me%a(90) +  k(:,13)*me%a(91)
        Yint = Y0 + h*aijkj
        k(:,14) = pDiffEqSys%F(X0 + h*me%c(14), Yint, params)
        
        ! 15th stage
        aijkj = k(:,1)*me%a(92) +  k(:,6)*me%a(97) +  k(:,7)*me%a(98) + k(:,8)*me%a(99) &
            + k(:,11)*me%a(102) +  k(:,12)*me%a(103) +  k(:,13)*me%a(104) +  k(:,14)*me%a(105)
        Yint = Y0 + h*aijkj
        k(:,15) = pDiffEqSys%F(X0 + h*me%c(15), Yint, params)
        
        ! 16th stage
        aijkj = k(:,1)*me%a(106) +  k(:,6)*me%a(111) +  k(:,7)*me%a(112) + k(:,8)*me%a(113) &
            + k(:,9)*me%a(114) +  k(:,13)*me%a(118) +  k(:,14)*me%a(119) +  k(:,15)*me%a(120)
        Yint = Y0 + h*aijkj
        k(:,16) = pDiffEqSys%F(X0 + h*me%c(16), Yint, params)
        
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
            
        cont4(:) = h*(me%d(1,1)*k(InterpStates,1) + me%d(6-4,1)*k(InterpStates,6) + me%d(7-4,1)*k(InterpStates,7) &
                    + me%d(8-4,1)*k(InterpStates,8) + me%d(9-4,1)*k(InterpStates,9) + me%d(10-4,1)*k(InterpStates,10) &
                    + me%d(11-4,1)*k(InterpStates,11) + me%d(12-4,1)*k(InterpStates,12) + me%d(13-4,1)*k(InterpStates,13) &
                    + me%d(14-4,1)*k(InterpStates,14) + me%d(15-4,1)*k(InterpStates,15) + me%d(16-4,1)*k(InterpStates,16))
            
        cont5(:) = h*(me%d(1,2)*k(InterpStates,1) + me%d(6-4,2)*k(InterpStates,6) + me%d(7-4,2)*k(InterpStates,7) &
                    + me%d(8-4,2)*k(InterpStates,8) + me%d(9-4,2)*k(InterpStates,9) + me%d(10-4,2)*k(InterpStates,10) &
                    + me%d(11-4,2)*k(InterpStates,11) + me%d(12-4,2)*k(InterpStates,12) + me%d(13-4,2)*k(InterpStates,13) &
                    + me%d(14-4,2)*k(InterpStates,14) + me%d(15-4,2)*k(InterpStates,15) + me%d(16-4,2)*k(InterpStates,16))
            
        cont6(:) = h*(me%d(1,3)*k(InterpStates,1) + me%d(6-4,3)*k(InterpStates,6) + me%d(7-4,3)*k(InterpStates,7) &
                    + me%d(8-4,3)*k(InterpStates,8) + me%d(9-4,3)*k(InterpStates,9) + me%d(10-4,3)*k(InterpStates,10) &
                    + me%d(11-4,3)*k(InterpStates,11) + me%d(12-4,3)*k(InterpStates,12) + me%d(13-4,3)*k(InterpStates,13) &
                    + me%d(14-4,3)*k(InterpStates,14) + me%d(15-4,3)*k(InterpStates,15) + me%d(16-4,3)*k(InterpStates,16))
  
        cont7(:) = h*(me%d(1,4)*k(InterpStates,1) + me%d(6-4,4)*k(InterpStates,6) + me%d(7-4,4)*k(InterpStates,7) &
                    + me%d(8-4,4)*k(InterpStates,8) + me%d(9-4,4)*k(InterpStates,9) + me%d(10-4,4)*k(InterpStates,10) &
                    + me%d(11-4,4)*k(InterpStates,11) + me%d(12-4,4)*k(InterpStates,12) + me%d(13-4,4)*k(InterpStates,13) &
                    + me%d(14-4,4)*k(InterpStates,14) + me%d(15-4,4)*k(InterpStates,15) + me%d(16-4,4)*k(InterpStates,16))
            
        ! check Hairer's contd8 routine for verifying the following code
        IntpCoeff(:,0) = cont0(:)
        IntpCoeff(:,1) = cont1(:) + cont2(:)
        IntpCoeff(:,2) = cont3(:) + cont4(:) - cont2(:)
        IntpCoeff(:,3) = cont5(:) + cont6(:) - 2*cont4(:) - cont3(:)
        IntpCoeff(:,4) = cont7(:) - 3*cont6(:) - 2*cont5(:) + cont4(:)
        IntpCoeff(:,5) = -3*cont7(:) + 3*cont6(:) + cont5(:)
        IntpCoeff(:,6) = 3*cont7(:) - cont6(:)
        IntpCoeff(:,7) = -cont7(:)

    end function DOP853_IntpCoeff
 
    !> It advances integrator by 1 step by computing stages
    subroutine DOP54_stepint(X0, Y0, Y1)
    
        implicit none
        
        real(WP), intent(in) :: X0
        real(WP), dimension(n), intent(in) :: Y0
        real(WP), dimension(n), intent(out) :: Y1
        
        real(WP), dimension(n) :: Yint, aijkj, Sc
        integer :: i
        
        ! first stage
        k(:,1) = F0  
        
        ! 2nd stage
        aijkj = k(:,1)*me%a(1)           
        Yint = Y0 + h*aijkj
        k(:,2) = pDiffEqSys%F(X0 + h*me%c(2), Yint, params)

        ! 3rd stage
        aijkj = k(:,1)*me%a(2) + k(:,2)*me%a(3)           
        Yint = Y0 + h*aijkj        
        k(:,3) = pDiffEqSys%F(X0 + h*me%c(3), Yint, params)
        
        ! 4th stage
        aijkj = k(:,1)*me%a(4) + k(:,2)*me%a(5) + k(:,3)*me%a(6)    
        Yint = Y0 + h*aijkj        
        k(:,4) = pDiffEqSys%F(X0 + h*me%c(4), Yint, params)
        
        ! 5th stage
        aijkj = k(:,1)*me%a(7) + k(:,2)*me%a(8) + k(:,3)*me%a(9) + k(:,4)*me%a(10)           
        Yint = Y0 + h*aijkj
        k(:,5) = pDiffEqSys%F(X0 + h*me%c(5), Yint, params)
        
        ! 6th stage
        aijkj = k(:,1)*me%a(11) + k(:,2)*me%a(12) + k(:,3)*me%a(13) &
                +  k(:,4)*me%a(14) + k(:,5)*me%a(15)
        Yint12 = Y0 + h*aijkj                
        k(:,6) = pDiffEqSys%F(X0 + h*me%c(6), Yint12, params)
        
        ! 7th stage
        aijkj = k(:,1)*me%a(16) + k(:,2)*me%a(17) + k(:,3)*me%a(18) &
                + k(:,4)*me%a(19) + k(:,5)*me%a(20) + k(:,6)*me%a(21)
        Yint = Y0 + h*aijkj  
        k(:,7) = pDiffEqSys%F(X0 + h, Yint, params)
        
        ! propagate the solution to the next step for error computation
        ! For FSAL, Yint is the new solution

        ! The scale factor for error computations
        if (IsScalarTol .EQV. .FALSE.) then
            Sc = me%ATol + me%RTol*max(abs(Y0), abs(Yint))
        else
            Sc = me%ATol(1) + me%RTol(1)*max(abs(Y0), abs(Yint))
        end if 

        ! Estimate the error. We assume that the error coeffcients are precomputed,
        ! i.e., e_i = b_i - bhat_i
        
        Err = abs(h)*sqrt(sum(((k(:,1)*me%e(1) + k(:,3)*me%e(3) + k(:,4)*me%e(4) + k(:,5)*me%e(5) &
                + k(:,6)*me%e(6) + k(:,7)*me%e(7))/Sc)**2)/n)

        ! If this step is accepted then return the new solution
        ! else just return the initial condition
        if (Err <= 1.0_WP) then 
            Y1 = Yint
            ! In FSAL methods, the intgeration last stage is F0 for the next step
            F0 = k(:,7)
            ! Function calls made if the step is accepted
            FCalls = sint - 1    
        else
            Y1 = Y0
            ! Function calls made if the step is rejected
            FCalls = sint - 2            
        end if
        
    end subroutine DOP54_stepint
    
    
    function DOP54_IntpCoeff(X0, Y0, X1, Y1, InterpStates) result(IntpCoeff)

        implicit none
    
        real(WP), intent(in) :: X0, X1
        real(WP), dimension(n), intent(in) :: Y0, Y1
        integer, dimension(:), intent(in) :: InterpStates
        
        real(WP), dimension(size(InterpStates),0:4) :: IntpCoeff
        
        real(WP), dimension(n) :: Yint, aijkj
        real(WP), dimension(size(InterpStates)) :: cont0, cont1, cont2, cont3, cont4
        integer :: i, j
        
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
            
        cont4(:) = h*(me%d(1,1)*k(InterpStates,1) + me%d(3-1,1)*k(InterpStates,3) & 
                    + me%d(4-1,1)*k(InterpStates,4) + me%d(5-1,1)*k(InterpStates,5) &
                    + me%d(6-1,1)*k(InterpStates,6) + me%d(7-1,1)*k(InterpStates,7))
            
        ! check Hairer's contd5 routine for verifying the following code
        IntpCoeff(:,0) = cont0
        IntpCoeff(:,1) = cont1 + cont2
        IntpCoeff(:,2) = cont3 + cont4 - cont2
        IntpCoeff(:,3) = -2*cont4 - cont3
        IntpCoeff(:,4) = cont4

    end function DOP54_IntpCoeff    
    
    function IntpCoeff(X0, Y0, X1, Y1, InterpStates) result(IpCoeff)

        implicit none
    
        real(WP), intent(in) :: X0, X1
        real(WP), dimension(n), intent(in) :: Y0, Y1
        integer, dimension(:), intent(in) :: InterpStates
        
        real(WP), dimension(size(InterpStates),0:me%pstar) :: IpCoeff
        
        real(WP), dimension(n) :: Yint, aijkj
        integer :: i, j
        
        ! Compute extra stages needed for interpolation
    
        ! compute rest of the stages (the following compact code is slower!)   
        do i = (sint+1),s
            block
            integer :: astart, j
            ! starting index of a_ij, where i=i, j=1:(i-1)
            ! It is a series sum: 1 + 2 + 3 + ...
            astart = int((i-1)/2.0*(i-2))
            ! In my testing, I am seeing matmul performing better than 
            ! MKL's gemv. However, do concurrent is faster than both but 
            ! but still slower than hard-coded computations
            ! call gemv(k(:,1:(i-1)), me%a((astart+1):(i-1)), aijkj)
            ! aijkj = matmul(k(:,1:(i-1)), me%a((astart+1):(astart+i-1)))
            aijkj = 0.0_WP
            do concurrent (j = 1:(i-1))
                aijkj = aijkj + k(:,j)*me%a(astart + j)
            end do
            
            Yint = Y0 + h*aijkj
            k(:,i) = pDiffEqSys%F(X0 + h*me%c(i), Yint, params)
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
                
                cont = 0.0
        sloop: do i = 1,s
                    ! get the current stage that gets multiplied by non-zero dij
                    istage = me%dinz(i)
                    ! -1 means we are at the end of the list
                    if (istage == -1) exit sloop
                    
                    cont = cont + me%d(i,j)*k(InterpStates,istage)
                end do sloop
                IpCoeff(:,j) = h*cont

            end block
        end do
    end function IntpCoeff    
    
    
    !> This procedure is passed to Root function for finding its zeros.
    pure function SingleEvent(xx) result(gg)
        implicit none
        real(WP), intent(in) :: xx !< independent variable
        real(WP) :: gg
        
        real(WP), dimension(n) :: ff
        real(WP), dimension(m) :: eval
        integer, dimension(m) :: edir
        logical, dimension(m) :: eterm
        
        ! interpolate solution
        ff = InterpY(size(me%InterpStates), xx, X, h, me%pstar, Bip)
        ! evaluate the single event function
        call pDiffEqSys%G(EventId, xx, ff, eval, edir, eterm)
        ! return the event value
        gg = eval(EventId)
    end function SingleEvent    
    
    
    !> Stiffness detection algorithm. See Hairer' DOP853 code for reference.
    subroutine CheckStiffness(IsProblemStiff)
    
        implicit none
        logical, intent(out) :: IsProblemStiff
    
        real(WP) :: StiffN, StiffD, hLamb

        IsProblemStiff = .FALSE.
        
        if (mod(AcceptedSteps, StiffTestSteps) == 0 .OR. StiffThreshold > 0) then
            StiffN = norm2(F0 - k(:,sint-1))
            StiffD = norm2(Y2 - Yint12)

            if (StiffD > 0.0) hLamb = abs(h)*StiffN/StiffD
            
            if (hLamb > 6.1) then
                NonStiffThreshold = 0
                StiffThreshold = StiffThreshold + 1
                if (StiffThreshold == 15) then
                    IsProblemStiff = .TRUE.
                end if
            else
                NonStiffThreshold = NonStiffThreshold + 1
                if (NonStiffThreshold == 6) StiffThreshold = 0
            end if
        end if
        
    end subroutine CheckStiffness
    
    
    
    end subroutine erk_int

    
        
end submodule ERKIntegrate


