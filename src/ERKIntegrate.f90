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
!> \brief       ERK Integration Submodule
!! \details     It provides the implementation for Integrate method of the ERK class.
!! \author      Bharat Mahajan
!! \date        Created: 02/06/2019    
!
!############################################################################################
    
submodule (ERK) ERKIntegrate

    contains    

    !> Integration main subroutine
    module subroutine erk_int(me, X0, Y0, Xf, Yf, StepSz, UseConstStepSz, IntStepsOn, &
                                Xint, Yint, EventMask, EventStates, StiffTest, params)
    
        implicit none
        
        class(ERK_class), intent(inout) :: me
    
        real(WP), intent(in) :: X0
        real(WP), dimension(me%pDiffEqSys%n), intent(in) :: Y0
        real(WP), intent(inout) :: Xf            
        real(WP), dimension(me%pDiffEqSys%n), intent(out) :: Yf
        real(WP), intent(inout) :: StepSz
        logical, intent(in), optional :: UseConstStepSz            
        logical, intent(in), optional :: IntStepsOn
        real(WP), allocatable, dimension(:), intent(out), optional :: Xint            
        real(WP), allocatable, dimension(:,:), intent(out), optional :: Yint
        logical, dimension(me%pDiffEqSys%m), intent(in), optional :: EventMask
        real(WP), allocatable, dimension(:,:), intent(out), optional :: EventStates
        integer, intent(inout), optional :: StiffTest        
        real(WP), dimension(:), intent(in), optional :: params      


        class(DiffEqSys), pointer :: pDiffEqSys
        integer(kind(ERK_DOP853)) :: method

        integer(kind(FLINT_SUCCESS)) :: status, stat

        real(WP), dimension(me%pDiffEqSys%n, me%s) :: k                
        real(WP), dimension(size(me%InterpStates), 0:me%pstar) :: Bip
        real(WP), dimension(6) :: StepSzParams        

        real(WP) :: h, hSign, hnew, hmax, X, Err, LastStepSzFac
        real(WP), dimension(me%pDiffEqSys%n) :: Y1, Y2, F0, Sc0, Yint12

        integer :: n, m, FCalls, TotalFCalls, AcceptedSteps, RejectedSteps 
        integer :: TotalSteps, MaxSteps, s, sint, p, pstar, q, nInterpStates
        integer :: StiffnessTest, StiffThreshold, NonStiffThreshold, StiffTestSteps
        
        logical :: IntStepsNeeded, LastStepRejected, LAST_STEP, ConstStepSz, IsProblemStiff
        logical :: IsScalarTol, IsFSALMethod, InterpOn, EventsOn, BipComputed


        ! Events related variables
        integer :: EventId
        real(WP) :: EventX0, EventX1
        real(WP), dimension(me%pDiffEqSys%n) :: EventNewY
        integer, dimension(me%pDiffEqSys%m) :: EvalEvents, StepSzEvents, EventDir, nStepSzEvents
        logical, dimension(me%pDiffEqSys%m) :: EvMask, EventsDetected
        integer(kind(FLINT_EVENTACTION_CONTINUE)) :: EventAction
        real(WP), dimension(MAXNUMSTEPSZEVENTS, me%pDiffEqSys%m) :: EX0, EX1
        real(WP), dimension(MAXNUMSTEPSZEVENTS, me%pDiffEqSys%m) :: EV0, EV1
        real(WP), dimension(:), allocatable :: EventData     
        
        ! Init for this procedure
        IntStepsNeeded = .FALSE.
        LastStepRejected = .FALSE.
        LAST_STEP = .FALSE.
        status = FLINT_SUCCESS        
        
        ! Copy frequently needed data on the stack
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
        s = me%s
        sint = me%sint
        nInterpStates = size(me%InterpStates)
        StiffTestSteps = STIFFTEST_STEPS
        
        ! Sign of the step: negative for backward integration
        hSign = sign(1.0_WP, (Xf-X0))

        ! If UseConstStepSz option is present, then use it. Adaptive step size 
        ! algorithm is used by default. If constant step size is selected then
        ! StepSz must have proper value set by the user along with its sign.
        ConstStepSz = .FALSE.
        if (present(UseConstStepSz)) then
            ConstStepSz = UseConstStepSz
            if (ConstStepSz .AND. ((hSign*StepSz) <= 0)) status = FLINT_ERROR_CONSTSTEPSZ
        end if

        ! set the last step size expansion factor
        if (ConstStepSz) then
            LastStepSzFac = 1.0_WP
        else
            LastStepSzFac = LASTSTEP_EXPFAC
        end if

        ! If Interpolation is OFF and solution at integrator steps is needed,
        ! then allocate the storage for steps that will be returned to the user. 
        ! If Interpolation is ON, then no need as we already allocated the internal
        ! storage for the solution and we just use allocatable array assignment
        if (present(IntStepsOn)) then
            if (IntStepsOn .AND. InterpOn) then
                ! Both options mustn't be True at the same time
                status = FLINT_ERROR_PARAMS
            else if (IntStepsOn .AND. (.NOT. InterpOn)) then
                IntStepsNeeded = .TRUE.
            end if            
        end if

        if (IntStepsNeeded) then
            if (present(Xint) .AND. present(Yint)) then
                ! allocate 1 more than max steps to accommodate the initial condition
                allocate(Xint(MaxSteps+1), stat=stat)
                if (stat /= 0)  status = FLINT_ERROR_MEMALLOC
                allocate(Yint(pDiffEqSys%n,MaxSteps+1), stat=stat)
                if (stat /= 0)  status = FLINT_ERROR_MEMALLOC
            else
                ! With IntStepsNeeded, Xint and Yint must be provided
                status = FLINT_ERROR_PARAMS
            end if
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
            if ((.NOT. allocated(me%Xint)) .OR. (.NOT. allocated(me%Bip))) then
                me%status = FLINT_ERROR_INIT_REQD
            end if
            ! If size of Xint is less than the max size, then it means user is calling
            ! Intgerate again without doing Init first, so reallocate the memory
            if (allocated(me%Xint)) then
                if (size(me%Xint) < MaxSteps+1) then 
                    deallocate(me%Xint, stat = status)
                    deallocate(me%Bip, stat = status)
                    allocate(me%Xint(me%MaxSteps+1), stat=status)
                    allocate(me%Bip(nInterpStates, 0:me%pstar, me%MaxSteps+1), stat=status)
                    if (status /= 0) me%status = FLINT_ERROR_MEMALLOC
                end if
            end if
        end if        

        ! Make sure optional event arguments are present if events are enabled        
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
        end if

        ! Return if any error up to this point
        if (status /= FLINT_SUCCESS) then
            me%status = status
            return
        end if

        ! All the parameter checking should be done before this point

        ! Initialize the counters to prepare for the main integration loop
        X = X0
        Y1 = Y0
        TotalFCalls = 0
        TotalSteps = 0
        AcceptedSteps = 0
        RejectedSteps = 0
        StiffThreshold = 0
        NonStiffThreshold = 0
        status = FLINT_SUCCESS

        ! Evaluate the event function at the initial condition to check for
        ! events between X0 and the first integrated state
        if (EventsOn) then
            EventX0 = X
            EventNewY = Y1
            EvalEvents = 1
            call pDiffEqSys%G(EventX0, EventNewY, EvalEvents, EV0(1,:), EventDir)
        end if

        ! Compute the derivative and scale factor at the initial condition
        F0 = pDiffEqSys%F(X, Y1, params)
        TotalFCalls = TotalFCalls + 1
        
        !! \remark Scale factor is computed using
        !! \f[ SC_i = atol_i + rtol_i |Y_{i}| \quad i=1..n,\f]
        !! where \f$|Y|\f$ is the absolute value of the provided vector
        if (.NOT. IsScalarTol) then
            Sc0 = me%ATol + me%RTol*abs(Y1)
        else
            Sc0 = me%ATol(1) + me%RTol(1)*abs(Y1)
        end if        
        
        ! Select the max step size if not provided by the user
        if (hmax == 0.0_WP) hmax = abs(Xf - X)

        ! The initial step size: compute if user has not provided a non-zero value
        if (StepSz == 0.0_WP) then
            call StepSz0Hairer(h, p, n, X, Y1, F0, Sc0, hSign, hmax, FCalls, pDiffEqSys, params)
            TotalFCalls = TotalFCalls + Fcalls
        else
            h = StepSz
        end if

        ! save the initial condition
        me%X0 = X
        me%Y0 = Y1

        ! copy the initial condition as the first state

        if (IntStepsNeeded) then
            Xint(1) = X
            Yint(:,1) = Y1
        end if

        if (InterpOn) me%Xint(1) = X

        ! Main loop for computing the steps

        do
        
            if (status == FLINT_SUCCESS .AND. (.NOT. LAST_STEP)    &
                    .AND. (TotalSteps < MaxSteps)) then

                ! if we are only slightly away from the final state
                ! then adjust to step size to finish the integration
                if (hSign*(X + LastStepSzFac*h - Xf) > 0) then
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
                ! For constant step size, Err must be set to 0 by stepint function

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
                    
                    ! Event Detection start here after a step is accepted

                    ! reset the flag for computing interpolation coefficients
                    BipComputed = .FALSE.

                    if (EventsOn) then

                        EvalEvents = 1
                        EventsDetected = .FALSE.
                        EventX0 = X
                        EventX1 = X + h
                        EX0 = EventX0
                        EX1 = EventX1
                        nStepSzEvents = 1
      
                        ! First check for events sign change at the step boundaries
                        EventNewY = Y2
                        call pDiffEqSys%G(EventX1, EventNewY, EvalEvents, EV1(1,:), EventDir)
                        do EventId = 1,m
                            call DetectEvent(EventId, EV0(1,EventId), EV1(1,EventId), EventsDetected(EventId))
                        end do

                        ! Check do we need to check sign changes for any event at smaller intervals
                        ! than the current step length
                        StepSzEvents = 1
                        where (me%EventStepSz == 0 .OR. me%EventStepSz > hSign*h)
                            StepSzEvents = 0
                        end where

                        if (any(StepSzEvents == 1)) then

                            ! For computing solution between the steps
                            if (.NOT. BipComputed) then
                                call ComputeInterpCoefficients()
                                BipComputed = .TRUE.
                            end if

                            ! Check for sign changes one by one for these events
                            do EventId = 1, m
                            block
                                real(WP), dimension(m) :: EventV0, EventV1, EventVh
                                real(WP) :: Eventh
                                logical :: EvStepSzDetected, LastStepSzEvent
                            
                                ! Skip the event if statically masked
                                if ((EvMask(EventId)) .AND. &
                                    (StepSzEvents(EventId)==1)) then

                                    EvalEvents = 0
                                    EvalEvents(EventId) = 1
    
                                    Eventh = hSign*me%EventStepSz(EventId)
                                    EventX1 = X + Eventh
                                    EventV0 = EV0(1,:)
                                    EventV1 = EV1(1,:)
                                    
                                    ! Save the original event values at X+h
                                    EventVh(EventId) = EV1(1,EventId)

                                    EventsDetected(EventId) = .FALSE.
                                    LastStepSzEvent = .FALSE.

                                    ! detect events inside the integrator's step
                                    do while (EventX1 <= (X+h))
                                        if (EventX1 /= (X+h)) then
                                            ! interpolate solution
                                            EventNewY = InterpY(n, EventX1, X, h, me%pstar, Bip)
                                            ! evaluate this specific event function
                                            call pDiffEqSys%G(EventX1, EventNewY, EvalEvents, &
                                                                EventV1, EventDir)
                                        else
                                            EventV1(EventId) = EventVh(EventId)
                                        end if
                                        call DetectEvent(EventId, EventV0(EventId),EventV1(EventId),EvStepSzDetected)
                                        if (EvStepSzDetected) then
                                            ! store the this event details
                                            EventsDetected(EventId) = .TRUE.
                                            EX0(nStepSzEvents(EventId),EventId) = EventX0
                                            EX1(nStepSzEvents(EventId),EventId) = EventX1
                                            EV0(nStepSzEvents(EventId),EventId) = EventV0(EventId)
                                            EV1(nStepSzEvents(EventId),EventId) = EventV1(EventId)
                                            nStepSzEvents(EventId) = nStepSzEvents(EventId) + 1
                                            if (nStepSzEvents(EventId) > MAXNUMSTEPSZEVENTS) then
                                                ! max number of step size events has been detected
                                                exit
                                            end if
                                        end if
                                        ! prepare for the next step interval
                                        if (LastStepSzEvent) exit
                                        EventX0 = EventX1
                                        EventX1 = EventX1 + Eventh
                                        if (EventX1 > X+h) then
                                            EventX1 = X + h
                                            LastStepSzEvent = .TRUE.
                                        end if
                                        EventV0(EventId) = EventV1(EventId)                        
                                    end do
                                    ! nStepSzEvents always counts one more than actual no. of 
                                    ! step size events, so we need to subtract 1
                                    if (nStepSzEvents(EventId) > 1) nStepSzEvents(EventId) = nStepSzEvents(EventId) - 1
                                end if
                            end block
                            end do
                        end if

                        ! If any event is triggered?

                        if (any(EventsDetected)) then
                        block
                            
                            real(WP), dimension(MAXNUMSTEPSZEVENTS,m) :: EXm
                            real(WP), dimension(n) :: EYm
                            integer :: StepEvId, rootiter, rooterror, nEvents
                            integer, dimension(1) :: MinEvXIndex

                            ! Locate the events first
                            
                            do EventId = 1, m
                                do StepEvId = 1, nStepSzEvents(EventId)
                                    if (EventsDetected(EventId)) then

                                        select case (me%EventOptions(EventId))          
                                        ! Left bounrday of the integrator step is the event location
                                        case (FLINT_EVENTOPTION_STEPBEGIN)
                                            EXm(StepEvId,EventId) = EX0(StepEvId,EventId)

                                        ! Right bounrday of the integrator step is the event location
                                        case (FLINT_EVENTOPTION_STEPEND)
                                            EXm(StepEvId,EventId) = EX1(StepEvId,EventId)
                                    
                                        case default
                                            ! root-finding is enabled, so locate the event exactly
                                            if (.NOT. BipComputed) then
                                                call ComputeInterpCoefficients()
                                                BipComputed = .TRUE.
                                            end if

                                            call Root(EX0(StepEvId,EventId), EX1(StepEvId, EventId), &
                                                            EV0(StepEvId,EventId), EV1(StepEvId,EventId), &
                                                            me%ETol(EventId), SingleEvent, &
                                                            EXm(StepEvId,EventId), EV0(StepEvId,EventId), &
                                                            rootiter, rooterror)
                                            if (rooterror /= 0 .AND. rooterror /= 2) then
                                                ! root finding failed, dont exit. Try to go as far as you can!
                                                status = FLINT_ERROR_EVENTROOT
                                            end if
                                        end select

                                    end if
                                end do
                            end do

                            ! Handle the event action starting from the earliest event

                            ! iterate until all events are handled
                            do while (any(EventsDetected))
                                
                                ! Find the earliest event. Note that the events' location
                                ! increases with the 1st index of EXm.
                                MinEvXIndex = minloc(EXm(1,:), EventsDetected)

                                ! solution at the event
                                select case (me%EventOptions(MinEvXIndex(1)))
                                case (FLINT_EVENTOPTION_STEPBEGIN)
                                    EYm = Y1
                                case (FLINT_EVENTOPTION_STEPEND)
                                    EYm = Y2
                                case default
                                    EYm = InterpY(nInterpStates, EXm(1,MinEvXIndex(1)),&
                                                            X, h, me%pstar, Bip)
                                end select

                                ! Call the event function with the located event info

                                EvalEvents = 0
                                EvalEvents(MinEvXIndex(1)) = 1
                                EventNewY = EYm
                                call pDiffEqSys%G(EXm(1,MinEvXIndex(1)), EventNewY, EvalEvents, &
                                                    EV0(1,:), EventDir, MinEvXIndex(1), EventAction)

                                ! save the event state
                                ! Note that if an event action alters the solution Y, event data returned to
                                ! the user still contains the original state corresponding to the event location
                                EventData = [EventData, EXm(1,MinEvXIndex(1)), EYm, real(MinEvXIndex(1), WP)]

                                ! Handle the event actions. Note that we must recompute interpolation
                                ! coefficients in case we update the solution Y2 or h.

                                ! Handle the event mask action
                                if (btest(EventAction, 7)) then
                                    EvMask(MinEvXIndex(1)) = .FALSE.
                                end if

                                ! The following actions are mutually exclusive
                                select case (IAND(EventAction, b'01111111'))
                                case (FLINT_EVENTACTION_TERMINATE)
                                    ! Return the event solution as the last solution
                                    h = EXm(1,MinEvXIndex(1)) - X
                                    Y2 = EYm
                                    LAST_STEP = .TRUE.
                                    BipComputed = .FALSE.
                                    status = FLINT_EVENT_TERM
                                    exit
                                case (FLINT_EVENTACTION_RESTARTINT)
                                    ! Restart the integration from this point on
                                    h = EXm(1,MinEvXIndex(1)) - X
                                    Y2 = EYm
                                    ! If it was already a last step, then reset it
                                    LAST_STEP = .FALSE.
                                    BipComputed = .FALSE.
                                    exit
                                case (FLINT_EVENTACTION_CHANGESOLY)
                                    ! Change the solution Y and restart the integration
                                    h = EXm(1,MinEvXIndex(1)) - X
                                    Y2 = EYm
                                    ! If it was already a last step, then reset it
                                    LAST_STEP = .FALSE.
                                    BipComputed = .FALSE.
                                    exit
                                end select

                                ! This event has been handled, so remove it for the next event
                                ! If the handled event is a step-size event then we need to 
                                ! shift the contents of its columns to left
                                if (StepSzEvents(MinEvXIndex(1)) == 1) then
                                    nStepSzEvents(MinEvXIndex(1)) = nStepSzEvents(MinEvXIndex(1)) - 1
                                    if (nStepSzEvents(MinEvXIndex(1)) == 0) then
                                        EventsDetected(MinEvXIndex(1)) = .FALSE.
                                    else
                                        EXm(1:nStepSzEvents(MinEvXIndex(1)), MinEvXIndex(1)) &
                                        = EXm(2:(nStepSzEvents(MinEvXIndex(1))+1),MinEvXIndex(1))
                                    end if
                                else
                                    EventsDetected(MinEvXIndex(1)) = .FALSE.
                                end if
                            end do
                        end block
                        end if

                        ! save EV0 for the next integrator's step
                        EV0 = EV1

                    end if

                    ! Do we need to store interpolation coefficients?
                    if (InterpOn) then
                        ! If EventAction is performed to terminate, restart or change SolY
                        ! then we must compute interpolation coefficients again.
                        if (.NOT. BipComputed) then
                            call ComputeInterpCoefficients()
                            BipComputed = .TRUE.
                        end if
                        ! store the solution internally at the natural step size
                        me%Xint(AcceptedSteps + 1) = X + h
                        !me%Yint(:,AcceptedSteps + 1) = Y2
                        me%Bip(:,:,AcceptedSteps) = Bip
                    end if

                    
                    ! Restart the integration and update Y2 with the modified solution.
                    ! Note that the modified solution is returned back to the user in Yint.
                    ! The following code must be after the code for saving interpolation
                    ! coefficients as these coefficients need the original Y only.
                    if (EventsOn) then
                        select case (IAND(EventAction, b'01111111')) 
                        case(FLINT_EVENTACTION_RESTARTINT)
                            F0 = pDiffEqSys%F(X + h, Y2, params)
                            TotalFCalls = TotalFCalls + 1                
                        case(FLINT_EVENTACTION_CHANGESOLY)
                            Y2 = EventNewY
                            F0 = pDiffEqSys%F(X + h, Y2, params)
                            TotalFCalls = TotalFCalls + 1                
                        end select
                    end if

                    ! return the solution at integrator's natural step size
                    if (IntStepsNeeded) then
                        Xint(AcceptedSteps + 1) = X + h
                        Yint(:,AcceptedSteps + 1) = Y2
                    end if

                    ! get ready for the next iteration
                    X = X + h
                    Y1 = Y2
    
                    ! compute the new step size now
                    if (.NOT. ConstStepSz) then
                        hnew = StepSzHairer(h, .FALSE., q, Err, StepSzParams)
                        ! make sure hnew is not greater than hmax
                        if (abs(hnew) > hmax) hnew = hsign*hmax
                        ! if the last step was rejected, then dont increase step size
                        if (LastStepRejected .EQV. .TRUE.) hnew = hSign*min(abs(hnew),abs(h))
                        ! reset
                        LastStepRejected = .FALSE.
                    else
                        hnew = h
                    end if

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

                ! Take another step only if step size is not too small and it is
                ! not the last step
                ! Reference: Hairer's DOP853 code but added the last step check
                if ((.NOT. ConstStepSz) .AND. (.NOT. LAST_STEP) &
                                    .AND. abs(h)*0.01_WP <= abs(X)*EPS) then
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
                        allocate(Biparr(nInterpStates,0:me%pstar,AcceptedSteps))
                    
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

                if (LAST_STEP .AND. status == FLINT_SUCCESS) then
                    ! we finished like we should
                elseif (EventsOn .AND. LAST_STEP .AND. status == FLINT_EVENT_TERM) then
                    ! One of the terminal event has triggered, do nothing
                else if (TotalSteps >= MaxSteps .AND. status == FLINT_SUCCESS) then
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
                
                ! save the stats in any case
                Xf = X
                Yf = Y1
                if (.NOT. ConstStepSz) StepSz = h                
                me%Xf = Xf
                me%hf = StepSz
                me%Yf = Yf
                me%AcceptedSteps = AcceptedSteps
                me%RejectedSteps = RejectedSteps
                me%TotalSteps = TotalSteps
                me%FCalls = TotalFCalls                
                me%status = status
                
                ! exit the main integration loop
                exit
                
            end if
            
        end do
        
        
    contains
    


    !> Procedure to compute the interpolation coefficients
    subroutine ComputeInterpCoefficients()
        
        select case (method)
        case (ERK_DOP853)
            Bip = DOP853_IntpCoeff(X, Y1, X+h, Y2, me%InterpStates)
        case (ERK_DOP54)
            Bip = DOP54_IntpCoeff(X, Y1, X+h, Y2, me%InterpStates)                                
        case default
            Bip = IntpCoeff(X, Y1, X+h, Y2, me%InterpStates)
        end select
        
        TotalFCalls = TotalFCalls + FCalls
    
    end subroutine


    !> Procedure to detect events by checking sign-changes. It only checks for
    !! the events for which EvIndex has value "1".
    subroutine DetectEvent(EventId, val0, val1, EvDetected)
        implicit none

        integer, intent(in) :: EventId
        real(WP), intent(in) :: val0, val1
        logical, intent(out) :: EvDetected

        EvDetected = .FALSE.

        ! Skip the event check if statically masked
        if (EvMask(EventId) .AND. &
                (.NOT. IEEE_IS_NAN(val0)) .AND. (.NOT. IEEE_IS_NAN(val1))) then

            ! check for event trigger    
            select case (EventDir(EventId))
            case (0)
                if ((val0<=0.0_WP .AND. val1>=0.0_WP) &
                        .OR.(val0>=0.0_WP .AND. val1<=0.0_WP)) EvDetected = .TRUE.
            case (1)
                if ((val0<=0.0_WP .AND. val1>=0.0_WP)) EvDetected = .TRUE.
            case (-1)
                if ((val0>=0.0_WP .AND. val1<=0.0_WP)) EvDetected = .TRUE.
            end select

        end if             
        
    end subroutine DetectEvent


    
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

        if (.NOT. ConstStepSz) then
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
            !Err = abs(h)*sqrt(sum((temp/Sc)**2)/n)
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

        if (.NOT. ConstStepSz) then        
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
        else
            Err = 0.0_WP
        end if
        

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

        if (.NOT. ConstStepSz) then
            ! The scale factor for error computations
            if (IsScalarTol .EQV. .FALSE.) then
                Sc = me%ATol + me%RTol*max(abs(Y0), abs(Yint))
            else
                Sc = me%ATol(1) + me%RTol(1)*max(abs(Y0), abs(Yint))
            end if 

            ! Estimate the error. We assume that the error coeffcients are precomputed,
            ! i.e., e_i = b_i - bhat_i
            Err = sqrt(sum(((h*(k(:,1)*me%e(1) + k(:,3)*me%e(3) + k(:,4)*me%e(4) &
                                + k(:,5)*me%e(5) + k(:,6)*me%e(6) + k(:,7)*me%e(7)))/Sc)**2)/n)
        else
            Err = 0.0_WP
        end if
        
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
                
                cont = 0.0_WP
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
    function SingleEvent(xx) result(gg)
        implicit none
        real(WP), intent(in) :: xx !< independent variable
        real(WP) :: gg
        
        real(WP), dimension(n) :: ff
        real(WP), dimension(m) :: eval
        integer, dimension(m) :: EvalEvent, edir
        
        ! interpolate solution
        ff = InterpY(nInterpStates, xx, X, h, me%pstar, Bip)
        ! evaluate the single event function
        EvalEvent = 0
        EvalEvent(EventId) = 1
        call pDiffEqSys%G(xx, ff, EvalEvent, eval, edir)
        ! return the event value
        gg = eval(EventId)
    end function SingleEvent    
    
    
    !> Stiffness detection algorithm. See Hairer' DOP853 code for reference.
    subroutine CheckStiffness(IsProblemStiff)
    
        implicit none
        logical, intent(out) :: IsProblemStiff
    
        real(WP) :: StiffN, StiffD, hLamb

        IsProblemStiff = .FALSE.
        hLamb = 0.0_WP
        
        if (mod(AcceptedSteps, StiffTestSteps) == 0 .OR. StiffThreshold > 0) then
            StiffN = norm2(F0 - k(:,sint-1))
            StiffD = norm2(Y2 - Yint12)

            if (StiffD > 0.0_WP) hLamb = abs(h)*StiffN/StiffD
            
            if (hLamb > 6.1_WP) then
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


