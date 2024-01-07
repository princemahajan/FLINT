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
    module subroutine erk_int(me, X0, Y0, Xf, Yf, StepSz, UseConstStepSz, StepAdvance, IntStepsOn, &
                                Xint, Yint, EventMask, EventStates, StiffTest, params)
    
        implicit none
        
        class(ERK_class), intent(inout) :: me
    
        real(WP), intent(in) :: X0
        real(WP), dimension(me%pDiffEqSys%n), intent(in) :: Y0
        real(WP), intent(inout) :: Xf            
        real(WP), dimension(size(Y0)), intent(out) :: Yf
        real(WP), intent(inout) :: StepSz
        logical, intent(in), optional :: UseConstStepSz  
        logical, intent(in), optional :: StepAdvance
        logical, intent(in), optional :: IntStepsOn
        real(WP), allocatable, dimension(:), intent(out), optional :: Xint            
        real(WP), allocatable, dimension(:,:), intent(out), optional :: Yint
        logical, dimension(me%pDiffEqSys%m), intent(in), optional :: EventMask
        real(WP), allocatable, dimension(:,:), intent(out), optional :: EventStates
        integer, intent(inout), optional :: StiffTest        
        real(WP), dimension(:), intent(in), optional :: params      


        class(DiffEqSys), pointer :: pDiffEqSys

        integer(kind(FLINT_SUCCESS)) :: status, stat

        real(WP), dimension(size(me%InterpStates), 0:me%pstar) :: Bip
        real(WP), dimension(6) :: StepSzParams        

        real(WP) :: h, hSign, hnew, hmax, X, Err, LastStepSzFac
        real(WP), dimension(me%pDiffEqSys%n) :: Y, Y1, F0, Sc0, Yint12

        integer :: FCalls, nInterpStates
        integer :: StiffnessTest, StiffThreshold, NonStiffThreshold, StiffTestSteps
        
        logical :: IntStepsNeeded, LastStepRejected, LAST_STEP, ConstStepSz, IsProblemStiff
        logical :: StepOnce, InterpOn, EventsOn, BipComputed

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
        
        pDiffEqSys => me%pDiffEqSys
        hmax = me%MaxStepSize
        InterpOn = me%InterpOn
        EventsOn = me%EventsOn
        nInterpStates = size(me%InterpStates)
        StepSzParams = me%StepSzParams        
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

        ! The default choice is to integrate to the final time unless this is set true
        StepOnce = .FALSE.
        if (present(StepAdvance)) StepOnce = StepAdvance

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
                allocate(Xint(me%MaxSteps+1), stat=stat)
                if (stat /= 0)  status = FLINT_ERROR_MEMALLOC
                allocate(Yint(pDiffEqSys%n,me%MaxSteps+1), stat=stat)
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
                if (size(me%Xint) < me%MaxSteps+1) then 
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
        Y = Y0
        me%FCalls = 0
        me%TotalSteps = 0
        me%AcceptedSteps = 0
        me%RejectedSteps = 0
        StiffThreshold = 0
        NonStiffThreshold = 0
        status = FLINT_SUCCESS

        ! Evaluate the event function at the initial condition to check for
        ! events between X0 and the first integrated state
        if (EventsOn) then
            EventX0 = X
            EventNewY = Y
            EvalEvents = 1
            call pDiffEqSys%G(EventX0, EventNewY, EvalEvents, EV0(1,:), EventDir)
        end if

        ! Compute the derivative and scale factor at the initial condition
        F0 = pDiffEqSys%F(X, Y, params)
        me%FCalls = me%FCalls + 1

        !! \remark Scale factor is computed using
        !! \f[ SC_i = atol_i + rtol_i |Y_{i}| \quad i=1..n,\f]
        !! where \f$|Y|\f$ is the absolute value of the provided vector
        if (.NOT. me%IsScalarTol) then
            Sc0 = me%ATol + me%RTol*abs(Y)
        else
            Sc0 = me%ATol(1) + me%RTol(1)*abs(Y)
        end if        
        
        ! Select the max step size if not provided by the user
        if (hmax == 0.0_WP) hmax = abs(Xf - X)

        ! The initial step size: compute if user has not provided a non-zero value
        if (StepSz == 0.0_WP) then
            call StepSz0Hairer(h, me%p, pDiffEqSys%n, X, Y, F0, Sc0, hSign, hmax, FCalls, pDiffEqSys, params)
            me%FCalls = me%FCalls + Fcalls
        else
            h = StepSz
        end if

        ! save the initial conditions
        me%X0 = X
        me%Y0 = Y
        ! me%k(:,1) = F0

        ! copy the initial condition as the first state

        if (IntStepsNeeded) then
            Xint(1) = X
            Yint(:,1) = Y
        end if

        if (InterpOn) me%Xint(1) = X

        ! Main loop for computing the steps

        do while (status == FLINT_SUCCESS .AND. (.NOT. LAST_STEP)  &
                                .AND. (me%TotalSteps < me%MaxSteps))
        
            ! if we are only slightly away from the final state
            ! then adjust to step size to finish the integration
            if (hSign*(X + LastStepSzFac*h - Xf) > 0) then
                h = Xf - X 
                LAST_STEP = .TRUE.
            end if

            ! Advance one step using the chosen method
            call me%StepInt(X, Y, F0, h, Y1, Yint12, FCalls, (.NOT. ConstStepSz), Err, params)

            ! update steps taken and function calls made
            me%TotalSteps = me%TotalSteps + 1                
            me%FCalls = me%FCalls + FCalls                
            
            ! check error to see if the step needs to be accepted or rejected
            ! For constant step size, Err must be set to 0 by stepint function

            if (Err <= 1.0_WP) then

                ! step is accepted
                me%AcceptedSteps = me%AcceptedSteps + 1

                ! this is the last step if user wants only 1 step taken
                if (StepOnce) LAST_STEP = .TRUE.

                ! reset the flag for computing interpolation coefficients
                BipComputed = .FALSE.

                ! TBD: what is its use?
                ! save the first accepted step-size
                ! if (me%AcceptedSteps == 1) me%h0 = h

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


                if (EventsOn) then

                    EvalEvents = 1
                    EventsDetected = .FALSE.
                    EventX0 = X
                    EventX1 = X + h
                    EX0 = EventX0
                    EX1 = EventX1
                    nStepSzEvents = 1
    
                    ! First check for events sign change at the step boundaries
                    EventNewY = Y1
                    call pDiffEqSys%G(EventX1, EventNewY, EvalEvents, EV1(1,:), EventDir)
                    do EventId = 1,pDiffEqSys%m
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
                            call me%InterpCoeff(X, Y, h, Y1, Bip, FCalls, params)
                            me%FCalls = me%FCalls + FCalls
                            BipComputed = .TRUE.
                        end if

                        ! Check for sign changes one by one for these events
                        do EventId = 1, pDiffEqSys%m
                        block
                            real(WP), dimension(pDiffEqSys%m) :: EventV0, EventV1, EventVh
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
                                        EventNewY = InterpY(pDiffEqSys%n, EventX1, X, h, me%pstar, Bip)
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
                        
                        real(WP), dimension(MAXNUMSTEPSZEVENTS,pDiffEqSys%m) :: EXm
                        real(WP), dimension(pDiffEqSys%n) :: EYm
                        integer :: StepEvId, rootiter, rooterror, nEvents
                        integer, dimension(1) :: MinEvXIndex

                        ! Locate the events first
                        
                        do EventId = 1, pDiffEqSys%m
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
                                            call me%InterpCoeff(X, Y, h, Y1, Bip, FCalls, params)
                                            me%FCalls = me%FCalls + FCalls
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
                                EYm = Y
                            case (FLINT_EVENTOPTION_STEPEND)
                                EYm = Y1
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
                            ! coefficients in case we update the solution Y1 or h.

                            ! Handle the event mask action
                            if (btest(EventAction, 7)) then
                                EvMask(MinEvXIndex(1)) = .FALSE.
                            end if

                            ! The following actions are mutually exclusive
                            select case (IAND(EventAction, b'01111111'))
                            case (FLINT_EVENTACTION_TERMINATE)
                                ! Return the event solution as the last solution
                                h = EXm(1,MinEvXIndex(1)) - X
                                Y1 = EYm
                                LAST_STEP = .TRUE.
                                BipComputed = .FALSE.
                                status = FLINT_EVENT_TERM
                                exit
                            case (FLINT_EVENTACTION_RESTARTINT)
                                ! Restart the integration from this point on
                                h = EXm(1,MinEvXIndex(1)) - X
                                Y1 = EYm
                                ! If it was already a last step, then reset it
                                LAST_STEP = .FALSE.
                                BipComputed = .FALSE.
                                exit
                            case (FLINT_EVENTACTION_CHANGESOLY)
                                ! Change the solution Y and restart the integration
                                h = EXm(1,MinEvXIndex(1)) - X
                                Y1 = EYm
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
                        call me%InterpCoeff(X, Y, h, Y1, Bip, FCalls, params)
                        me%FCalls = me%FCalls + FCalls
                        BipComputed = .TRUE.
                    end if
                    ! store the solution internally at the natural step size
                    me%Xint(me%AcceptedSteps + 1) = X + h
                    me%Bip(:,:,me%AcceptedSteps) = Bip
                end if

                ! Restart the integration and update Y1 with the modified solution.
                ! Note that the modified solution is returned back to the user in Yint.
                ! The following code must be after the code for saving interpolation
                ! coefficients as these coefficients need the original Y only.
                if (EventsOn) then
                    select case (IAND(EventAction, b'01111111')) 
                    case(FLINT_EVENTACTION_RESTARTINT)
                        me%k(:,1) = pDiffEqSys%F(X + h, Y1, params)
                        me%FCalls = me%FCalls + 1                
                    case(FLINT_EVENTACTION_CHANGESOLY)
                        Y1 = EventNewY
                        me%k(:,1) = pDiffEqSys%F(X + h, Y1, params)
                        me%FCalls = me%FCalls + 1                
                    end select
                end if

                ! return the solution at integrator's natural step size
                if (IntStepsNeeded) then
                    Xint(me%AcceptedSteps + 1) = X + h
                    Yint(:,me%AcceptedSteps + 1) = Y1
                end if

                ! get ready for the next iteration
                X = X + h
                Y = Y1

                ! compute the new step size now
                if (.NOT. ConstStepSz) then
                    hnew = StepSzHairer(h, .FALSE., me%q, Err, StepSzParams)
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
                if (me%AcceptedSteps >= 1) me%RejectedSteps = me%RejectedSteps + 1
                
                ! compute the new step size now
                hnew = StepSzHairer(h, .TRUE., me%q, Err, StepSzParams)

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
                                .AND. abs(h)*0.01_WP <= abs(X)*FLINT_EPS) then
                status = FLINT_ERROR_STEPSZ_TOOSMALL
            end if

        end do
                    
        ! We have either accomplished our job or finished prematurely?!

        ! return the solution at the integrator's internal steps
        if (IntStepsNeeded) then
            block
                real(WP), allocatable, dimension(:) :: Xarr            
                real(WP), allocatable, dimension(:,:) :: Yarr       
            
                allocate(Xarr(me%AcceptedSteps + 1))
                allocate(Yarr(pDiffEqSys%n,me%AcceptedSteps + 1))
            
                Xarr = Xint(1:me%AcceptedSteps + 1)
                Yarr = Yint(:,1:me%AcceptedSteps + 1)
            
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
                
                allocate(Xarr(me%AcceptedSteps + 1))
                allocate(Biparr(nInterpStates,0:me%pstar,me%AcceptedSteps))
            
                Xarr = me%Xint(1:me%AcceptedSteps + 1)
                Biparr = me%Bip(:,:,1:me%AcceptedSteps)
                
                call move_alloc(from=Xarr, to=me%Xint)
                !call move_alloc(from=Yarr, to=me%Yint)
                call move_alloc(from=Biparr, to=me%Bip)
            end block
        end if        

        ! Prepare event output data
        if (EventsOn .AND. allocated(EventData)) then                
            EventStates = reshape(EventData, [pDiffEqSys%n+2, int(size(EventData)/(pDiffEqSys%n+2))])
        end if

        if (LAST_STEP .AND. status == FLINT_SUCCESS) then
            ! we finished like we should
        elseif (EventsOn .AND. LAST_STEP .AND. status == FLINT_EVENT_TERM) then
            ! One of the terminal event has triggered, do nothing
        else if (me%TotalSteps >= me%MaxSteps .AND. status == FLINT_SUCCESS) then
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
        Yf = Y
        if (.NOT. ConstStepSz) StepSz = h                
        me%Xf = Xf
        me%hf = StepSz
        me%Yf = Yf
        me%status = status
            
    contains
    

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
   
    
    
    !> This procedure is passed to Root function for finding its zeros.
    function SingleEvent(xx) result(gg)
        implicit none
        real(WP), intent(in) :: xx !< independent variable
        real(WP) :: gg
        
        real(WP), dimension(pDiffEqSys%n) :: ff
        real(WP), dimension(pDiffEqSys%m) :: eval
        integer, dimension(pDiffEqSys%m) :: EvalEvent, edir
        
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
        
        if (mod(me%AcceptedSteps, StiffTestSteps) == 0 .OR. StiffThreshold > 0) then
            StiffN = norm2(me%k(:,me%sint) - me%k(:,me%sint-1))
            StiffD = norm2(Y1 - Yint12)

            if (StiffD > 0.0_WP) hLamb = abs(h)*StiffN/StiffD
            
            ! TBD: this threshold for hLamb can be set specific to each
            ! integator method, check Hairer codes for DOP54 and DOP853 specific
            ! values
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


