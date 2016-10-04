
MODULE COMMONS
    IMPLICIT NONE

    INTEGER STEPS, indexlw
    DOUBLE PRECISION, ALLOCATABLE :: RATES(:,:)
    integer, ALLOCATABLE :: config1(:,:),SAVconfig(:,:,:),config1OLD(:,:), &
                            EMAT(:,:),EMATSAV(:,:,:,:),NEVENTS(:),&
                            checkmat(:,:,:),defectmat(:,:),indexmat(:,:)

    integer ::  posNx, posNy !number of surface atom in x and y !Read in
    integer :: ntriang ! number of triangles !Read in
    double precision :: Temp !temper ature in K !Eperimentally 400 to 450 K !read in
    integer :: isltype1,isl1sizex,isl1sizey,colorisl1, isl1posx,isl1posy,&
        isltype2,isl2sizex,isl2sizey,colorisl2, isl2posx,isl2posy !specifying islands on surface


    integer :: isl1dfct,isl1dfctx,isl1dfcty,isl1dfct2x,isl1dfct2y, &
        isl2dfct,isl2dfctx,isl2dfcty,isl2dfct2x,isl2dfct2y, &
        mdef1type,mdef2type,mdef1sizex,mdef1sizey,mdef2sizex,mdef2sizey,&
        mdef1posx,mdef1posy,mdef2posx,mdef2posy

    double precision ::   hbondconst,neibrepel !eV to lower barrier if a h bond is established
    integer,parameter :: color=4, molRot=3, osmov=4
    double precision, parameter :: kboltz = 8.617332478D0*10.D0**(-5) !eV/K
    logical :: initrand

    character(len=100) :: restart_file
    logical :: restart


END MODULE COMMONS

PROGRAM KMC
    USE COMMONS
    IMPLICIT NONE
    INTEGER NKMC,J1,J2,printlevel,statelw
    DOUBLE PRECISION :: LDUMMY, RANDOM,TOTOALSUMRATES,lowestwait,totalwait
    DOUBLE PRECISION, ALLOCATABLE :: SUMRATES(:),WAIT(:)

    character(len=100) :: buffer, label
    integer :: pos
    integer, parameter :: fh = 15 !unit=15
    integer :: ios = 0
    integer :: line = 0

    CALL init_random_seed()

    restart=.false.
    initrand=.true.

    print*,'At the moment possible two islands ', &
        &'with holes and two defects on the surface are possible.'

    open(fh, file='system.in')


    do while (ios == 0)
        read(fh, '(A)', iostat=ios) buffer
        if (ios == 0) then
            line = line + 1

            ! Find the first instance of whitespace.  Split label and data.
            pos = scan(buffer, '    ')
            label = buffer(1:pos)
            buffer = buffer(pos+1:)

            select case (label)
                case ('ntriangle')
                    read(buffer, *, iostat=ios) ntriang
                    print*,'Simulation with ntriangle',ntriang
                case ('posNxy')
                    read(buffer, *, iostat=ios) posNx, posNy
                    print *, 'Reading number of metal atoms posNx',posNx,'posNy',posNy
                case ('temp')
                    read(buffer, *, iostat=ios) temp
                    print*,'At',temp,'temperature (K)'
                case ('hbondstr')
                    read(buffer, *, iostat=ios) hbondconst
                    print*,'At',hbondconst,'hbondstrength (eV) standard 0.15D0'
                case ('nkmcsteps')
                    read(buffer, *, iostat=ios) nkmc
                    print*,'Number of KMC steps:',nkmc
                case ('printlevel')
                    read(buffer, *, iostat=ios) printlevel
                    print*,'Print after every',printlevel,'steps'
                case ('isl-type1')
                    read(buffer, *, iostat=ios) isltype1
                    print*,'Island 1 Type (0=off,1=different,2=same triangles):',isltype1
                case ('isl1-size-xy')
                    read(buffer, *, iostat=ios) isl1sizex,isl1sizey
                    print*,'Island 1 size xy:',isl1sizex,isl1sizey
                case ('color-isl1')
                    read(buffer, *, iostat=ios) colorisl1
                    print*,'color1 (if type 1: 12 34, if type 2: 1 2 3 4)',colorisl1
                case ('isl1-pos-xy')
                    read(buffer, *, iostat=ios) isl1posx,isl1posy
                    print*,'pos.1 on surf. xy',isl1posx,isl1posy
                case ('isl-type2')
                    read(buffer, *, iostat=ios) isltype2
                    print*,'Island 2 Type (0=off,1=different,2=same triangles):',isltype2
                case ('isl2-size-xy')
                    read(buffer, *, iostat=ios) isl2sizex, isl2sizey
                    print*,'Island 2 size xy:',isl2sizex,isl2sizey
                case ('color-isl2')
                    read(buffer, *, iostat=ios) colorisl2
                    print*,'color2 (if type 1: 12 34, if type 2: 1 2 3 4)',colorisl2
                case ('isl2-pos-xy')
                    read(buffer, *, iostat=ios) isl2posx,isl2posy
                    print*,'Pos.2 on surf. xy',isl2posx,isl2posy
                case ('isl1-defect-type')
                    read(buffer, *, iostat=ios) isl1dfct
                    print*,'Island 1 has defect(s) (0,1,2):',isl1dfct
                case ('isl1-defects-xy')
                    read(buffer, *, iostat=ios) isl1dfctx,isl1dfcty,isl1dfct2x,isl1dfct2y
                    print*,'Defect island 1 position',isl1dfctx,isl1dfcty,isl1dfct2x,isl1dfct2y
                case ('isl2-defect-type')
                    read(buffer, *, iostat=ios) isl2dfct
                    print*,'Island 2 has defect(s) (0,1,2):',isl2dfct
                case ('isl2-defects-xy')
                    read(buffer, *, iostat=ios) isl2dfctx,isl2dfcty,isl2dfct2x,isl2dfct2y
                    print*,'Defect island 2 position',isl2dfctx,isl2dfcty,isl2dfct2x,isl2dfct2y
                case ('mdefect1-type')
                    read(buffer, *, iostat=ios) mdef1type
                    print*,'Metal defect1 (0 no, 1 yes)',mdef1type
                case ('mdefect1-size-xy')
                    read(buffer, *, iostat=ios) mdef1sizex,mdef1sizey
                    print*,'Metal defect1 size-xy',mdef1sizex,mdef1sizey
                case ('mdefect1-pos-xy')
                    read(buffer, *, iostat=ios) mdef1posx,mdef1posy
                    print*,'Metal defect1 position-xy',mdef1posx,mdef1posy
                case ('mdefect2-type')
                    read(buffer, *, iostat=ios) mdef2type
                    print*,'Metal defect2 (0 no, 1 yes)',mdef2type
                case ('mdefect2-size-xy')
                    read(buffer, *, iostat=ios) mdef2sizex,mdef2sizey
                    print*,'Metal defect2 size-xy',mdef2sizex,mdef2sizey
                case ('mdefect2-pos-xy')
                    read(buffer, *, iostat=ios) mdef2posx,mdef2posy
                    print*,'Metal defect2 position-xy',mdef2posx,mdef2posy
                case ('restartTrue')
                    read(buffer, *, iostat=ios) restart_file
                    print*,'restart from file',restart_file
                    restart=.true.
                case default
                    print *, 'Skipping invalid label at line', line
            end select
        end if
    end do
    close(fh)


    if((isltype1.ne.0).or.(isltype2.ne.0)) then
    initrand=.false.!
    !sanity check
    if (ntriang.ne.isl2sizex*isl2sizey+isl1sizex*isl1sizey-isl1dfct-isl2dfct) then
        print*,'number of triangles',ntriang,'is not equal to',&
            isl2sizex*isl2sizey+isl1sizex*isl1sizey-isl1dfct-isl2dfct
        stop
    endif
    endif


    ALLOCATE(RATES(ntriang,10)) !total number of possibilities each trianlge can do 10 different moves
    ALLOCATE(SAVconfig(10,ntriang,5)) !total of possible configurations to be chosen from in KMC
    ALLOCATE(EMATSAV(ntriang,10,posNx,posNy)) !total of possible occupations to be chosen from in KMC
    ALLOCATE(config1(ntriang,5)) !ALLOCATE inital configuration !adjust number of triangles
    ALLOCATE(config1old(ntriang,5)) !ALLOCATE old configuration vector !adjust
    ALLOCATE(EMAT(posNx,posNy))! periodic boundary conditions
    ALLOCATE(NEVENTS(ntriang))! number of possible events per triangle
    ALLOCATE(SUMRATES(ntriang))
    ALLOCATE(WAIT(ntriang))
    ALLOCATE(checkmat(ntriang,2,3))!number of triangles, xy position, 1=N 2=O 3=O
    ALLOCATE(defectmat(posNx,posNy))

    ALLOCATE(indexmat(posNx,posNy))!triangle indices

    config1(:,:)=0
    EMAT(:,:)=0
    totalwait=0.D0
    config1OLD(:,:)=0
    WAIT(:)=0.0D0
    RATES(:,:)=0.D0
    SAVconfig(:,:,:)=0
    EMATSAV(:,:,:,:)=0
    checkmat(:,:,:)=0
    defectmat(:,:)=0

    indexmat(:,:)=0

    !neibrepel=0.D0
    neibrepel=2.D0*hbondconst


    DO STEPS=1,NKMC
        CALL GETEVENTS
        TOTOALSUMRATES=0.D0
        DO J1=1,ntriang
        SUMRATES(J1)=0.0D0 !change!!
        DO J2=1,NEVENTS(J1)
            SUMRATES(J1)=SUMRATES(J1)+RATES(J1,J2)
            !print*, RATES(J1,J2),'RATES(J1,J2)',J1,J2
        ENDDO

        TOTOALSUMRATES=TOTOALSUMRATES+SUMRATES(J1)
        ENDDO

        IF (dabs(TOTOALSUMRATES) < 1.0e-8) THEN
            PRINT '(A)','ERROR *** all rate sums are zero',TOTOALSUMRATES
            STOP
        ENDIF


        lowestwait=10.D0
        indexlw=0
        statelw=0

        DO J1=1,ntriang !triangle check
        ! Random choice of move the chosen triangle makes
        LDUMMY=0.D0
        CALL RANDOM_NUMBER(RANDOM)
        if(NEVENTS(J1).gt.0) then
        choose: DO J2=1,NEVENTS(j1)
            LDUMMY=LDUMMY+RATES(J1,J2)/SUMRATES(J1)
            IF (LDUMMY.GT.RANDOM) THEN
                CALL RANDOM_NUMBER(RANDOM)
                WAIT(J1)=dlog(1.0/RANDOM)/SUMRATES(J1)
                ! ln(1.0/RANDOM) gives a waiting time drawn from the
                ! Poisson distribution rather than the mean of that distribution.
                EXIT choose
            ENDIF
        ENDDO choose
        else
             WAIT(J1)=100.D0 !for no allowed events assign an "infinite" waiting time
        endif

        if(WAIT(J1).lt.lowestwait) then

        if(WAIT(J1).gt.0.D0) then
        lowestwait=WAIT(J1)
        indexlw=J1 !triangle
        statelw=J2 !event
        endif
        endif
        enddo !triangle check


        !Update configuration matrix
        config1(:,:)=config1OLD(:,:)
        config1(indexlw,:)=SAVconfig(statelw,indexlw,:)

        !Update waiting time
        totalwait=totalwait+lowestwait
        config1OLD(:,:)=config1(:,:)


        !Call routine for checkmat update of trianlge(indexlw)
        call checkmatgen(config1(indexlw,:),checkmat(indexlw,:,:))

        !Calc. new Emat from checkmat
        EMAT(:,:)=defectmat(:,:) !add metal defects to emat or reset to zero
        indexmat(:,:)=0
        do j1=1,ntriang
        EMAT(checkmat(j1,1,1),checkmat(j1,2,1))=1
        EMAT(checkmat(j1,1,2),checkmat(j1,2,2))=2
        EMAT(checkmat(j1,1,3),checkmat(j1,2,3))=2
        !update indexmat
        indexmat(checkmat(j1,1,1),checkmat(j1,2,1))=j1
        indexmat(checkmat(j1,1,2),checkmat(j1,2,2))=j1
        indexmat(checkmat(j1,1,3),checkmat(j1,2,3))=j1
        enddo

        !**********
        print*,printlevel,'printlevel',MODULO(STEPS, printlevel),'mod',steps
        if(MODULO(STEPS, printlevel).EQ.0) then !Print every printlevel-th step, modulo printlevel
            !Print chosen state
            call printstate(config1,EMAT,lowestwait,STEPS) !TODO adjust!
        endif
              !**********

        PRINT '(A,I10,A,G20.10,A,I6,A,I6,A,I6,A)','After ',STEPS,' KMC steps elapsed time=',totalwait,&
            & ' current triangle',indexlw,' state',statelw,' chosen from ', &
            & NEVENTS(indexlw),' possibilities'


    ENDDO

    DEALLOCATE(RATES)
    DEALLOCATE(SAVconfig)
    DEALLOCATE(EMATSAV)
    DEALLOCATE(config1)
    DEALLOCATE(config1old)
    DEALLOCATE(EMAT)
    DEALLOCATE(NEVENTS)
    DEALLOCATE(SUMRATES)
    DEALLOCATE(WAIT)
    DEALLOCATE(indexmat)
    STOP


END PROGRAM KMC

!
! This subroutine must return the possible events and rates for the current state.
!
SUBROUTINE GETEVENTS
    USE COMMONS
    implicit none
    double precision:: rateconst
    integer ::config2(ntriang,5)!!( (numbergly(1)),(color(4),posNx,posNy,molrot(3),osmov(4))) -- (1,5) entries
    !numberdim(0),numbertrim(0) for later use
    integer :: i,j,k,l,m,chosentriang,eventcount,triang
    logical :: eventOK

    character ( len = 24 ) :: data_filename
    double precision :: newinit(ntriang,3) ! ntriang, (Grid (posx, posy), color)
    logical :: checkaccept
    integer :: dx,dx1,dx2,dy,dy1,dy2,jcheck1,jcheck2
    double precision :: rand
    integer :: hbondcount,occupationcount,hbondmat(posNx+1,posNx+1)! Matrix with hbonds
    integer :: neighbourindex(ntriang)
    !****** Initialisation ****************
    hbondmat(:,:)=0
    hbondcount=0
    occupationcount=0
    rateconst=0.D0
    eventcount=0
    neighbourindex(:)=0

    IF (STEPS.EQ.1) then
        if(.not.restart) then !if restart is not true
             !Write points, system runs from (0,0) to (posNx+1,posNy+1)
            if((mdef1type.eq.0).and.(mdef2type.eq.0)) then
                call write_vtk_points(1,1,posNx,posNy)
            else
                call write_vtk_points_defect(1,1,posNx,posNy,&
                    mdef1type,mdef2type,mdef1sizex,mdef1sizey,mdef2sizex,mdef2sizey,&
                    mdef1posx,mdef1posy,mdef2posx,mdef2posy)

                !emat update
                if(mdef1type.eq.1)then
                    do i=1,mdef1sizex
                        do j=1,mdef1sizey
                            emat(mdef1posx-1+i,mdef1posy-1+j)=3 !3 for metal defect
                            defectmat(mdef1posx-1+i,mdef1posy-1+j)=3
                        enddo
                    enddo
!                 elseif(mdef1type.eq.2) then
!                    print*,'random init of defects T1 on metals'
!                    do i=1,mdef1sizex
!                        do j=1,mdef1sizey
!                            CALL RANDOM_NUMBER(rand)
!                            dx=floor(rand*posNx)+1 !random posX
!                            CALL RANDOM_NUMBER(rand)
!                            dy=floor(rand*posNy)+1 !random posY
!                            emat(dx,dy)=3 !3 for metal defect
!                            defectmat(dx,dy)=3
!                        enddo
!                    enddo

                print*,'Warning: Periodic boundaries not implemented for metal defects'
                elseif(mdef1type.ge.2) then
                print*,'Wrong defects type1 on metals .ne. 0,1'
                STOP
                elseif(mdef1type.eq.0) then
                print*,'NO defects type1 on metals'
                endif

                if(mdef2type.eq.1)then
                    do i=1,mdef2sizex
                        do j=1,mdef2sizey
                            emat(mdef2posx-1+i,mdef2posy-1+j)=3 !3 for metal defect
                            defectmat(mdef2posx-1+i,mdef2posy-1+j)=3
                        enddo
                    enddo
!                elseif(mdef2type.eq.2) then
!                    print*,'random init of defects T2 on metals'
!                    do i=1,mdef2sizex
!                        do j=1,mdef2sizey
!                            CALL RANDOM_NUMBER(rand)
!                            dx=floor(rand*posNx)+1 !random posX
!                            CALL RANDOM_NUMBER(rand)
!                            dy=floor(rand*posNy)+1 !random posY
!                            emat(dx,dy)=3 !3 for metal defect
!                            defectmat(dx,dy)=3
!                        enddo
!                    enddo
                print*,'Warning: Periodic boundaries not implemented for metal defects'
                elseif(mdef2type.ge.2) then
                print*,'Wrong defects type2 on metals .ne. 0,1'
                STOP
                elseif(mdef2type.eq.0) then
                print*,'NO defects type2 on metals'
                endif
            endif



            if(isl1dfct.eq.2) then
                if((isl1dfctx.lt.isl1posx).or.(isl1dfctx.gt.isl1posx+isl1sizex)) then
                    print*,'island defect1 pos x', isl1dfctx,'not in island range x',&
                        isl1posx,isl1posx+isl1sizex
                    STOP
                elseif((isl1dfcty.lt.isl1posy).or.(isl1dfcty.gt.isl1posy+isl1sizey)) then
                    print*,'island defect1 pos y', isl1dfcty,'not in island range y',&
                        isl1posy,isl1posy+isl1sizey
                    STOP
                elseif((isl1dfct2x.lt.isl1posx).or.(isl1dfct2x.gt.isl1posx+isl1sizex)) then
                    print*,'island defect1 pos x2', isl1dfct2x,'not in island range x',&
                        isl1posx,isl1posx+isl1sizex
                    STOP
                elseif((isl1dfct2y.lt.isl1posy).or.(isl1dfct2y.gt.isl1posy+isl1sizey)) then
                    print*,'island defect1 pos y2', isl1dfct2y,'not in island range y',&
                        isl1posy,isl1posy+isl1sizey
                    STOP
                else
                    print*,'island defect1 pos x', isl1dfctx,'is in island range x',&
                        isl1posx,isl1posx+isl1sizex
                    print*,'island defect1 pos y', isl1dfcty,'is in island range y',&
                        isl1posy,isl1posy+isl1sizey
                    print*,'island defect1 pos x2', isl1dfct2x,'is in island range x',&
                        isl1posx,isl1posx+isl1sizex
                    print*,'island defect1 pos y2', isl1dfct2y,'is in island range y',&
                        isl1posy,isl1posy+isl1sizey
                endif
            elseif(isl1dfct.eq.1) then
                if((isl1dfctx.lt.isl1posx).or.(isl1dfctx.gt.isl1posx+isl1sizex)) then
                    print*,'island defect1 pos x', isl1dfctx,'not in island range x',&
                        isl1posx,isl1posx+isl1sizex
                    STOP
                elseif((isl1dfcty.lt.isl1posy).or.(isl1dfcty.gt.isl1posy+isl1sizey)) then
                    print*,'island defect1 pos y', isl1dfcty,'not in island range y',&
                        isl1posy,isl1posy+isl1sizey
                    STOP
                else
                    print*,'island defect1 pos x', isl1dfctx,'is in island range x',&
                        isl1posx,isl1posx+isl1sizex
                    print*,'island defect1 pos y', isl1dfcty,'is in island range y',&
                        isl1posy,isl1posy+isl1sizey
                endif
            elseif( isl1dfct.gt.2) then
                print*,'wrong number of island 1 defects, only maximum 2 allowed'
                stop
            endif


            if(isl2dfct.eq.2) then
                if((isl2dfctx.lt.isl2posx).or.(isl2dfctx.gt.isl2posx+isl2sizex)) then
                    print*,'island defect2 pos x', isl2dfctx,'not in island range x',&
                        isl2posx,isl2posx+isl2sizex
                    STOP
                elseif((isl2dfcty.lt.isl2posy).or.(isl2dfcty.gt.isl2posy+isl2sizey)) then
                    print*,'island defect2 pos y', isl2dfcty,'not in island range y',&
                        isl2posy,isl2posy+isl2sizey
                    STOP
                elseif((isl2dfct2x.lt.isl2posx).or.(isl2dfct2x.gt.isl2posx+isl2sizex)) then
                    print*,'island defect2 pos x2', isl2dfct2x,'not in island range x',&
                        isl2posx,isl2posx+isl2sizex
                    STOP
                elseif((isl2dfct2y.lt.isl2posy).or.(isl2dfct2y.gt.isl2posy+isl2sizey)) then
                    print*,'island defect2 pos y2', isl2dfct2y,'not in island range y',&
                        isl2posy,isl2posy+isl2sizey
                    STOP
                else
                    print*,'island defect2 pos x', isl2dfctx,'is in island range x',&
                        isl2posx,isl2posx+isl2sizex
                    print*,'island defect2 pos y', isl2dfcty,'is in island range y',&
                        isl2posy,isl2posy+isl2sizey
                    print*,'island defect2 pos x2', isl2dfct2x,'is in island range x',&
                        isl2posx,isl2posx+isl2sizex
                    print*,'island defect2 pos y2', isl2dfct2y,'is in island range y',&
                        isl2posy,isl2posy+isl2sizey
                endif
            elseif(isl2dfct.eq.1) then
                if((isl2dfctx.lt.isl2posx).or.(isl2dfctx.gt.isl2posx+isl2sizex)) then
                    print*,'island defect2 pos x', isl2dfctx,'not in island range x',&
                        isl2posx,isl2posx+isl2sizex
                    STOP
                elseif((isl2dfcty.lt.isl2posy).or.(isl2dfcty.gt.isl2posy+isl2sizey)) then
                    print*,'island defect2 pos y', isl2dfcty,'not in island range y',&
                        isl2posy,isl2posy+isl2sizey
                    STOP
                else
                    print*,'island defect2 pos x', isl2dfctx,'is in island range x',&
                        isl2posx,isl2posx+isl2sizex
                    print*,'island defect2 pos y', isl2dfcty,'is in island range y',&
                        isl2posy,isl2posy+isl2sizey
                endif
            elseif( isl2dfct.gt.2) then
                print*,'wrong number of island 2 defects, only maximum 2 allowed'
                stop
            endif

            !emat update
            !generating islanddefects in emat for STEPS=1 only
            if(isl1dfct.eq.1) then
                emat(isl1dfctx,isl1dfcty)=4 !4 for island defect
            elseif(isl1dfct.eq.2) then
                emat(isl1dfctx,isl1dfcty)=4 !4 for island defect
                emat(isl1dfct2x,isl1dfct2y)=4 !4 for island defect
            endif

            !emat update
            if(isl2dfct.eq.1) then
                emat(isl2dfctx,isl2dfcty)=4 !4 for island defect
            elseif(isl2dfct.eq.2) then
                emat(isl2dfctx,isl2dfcty)=4 !4 for island defect
                emat(isl2dfct2x,isl2dfct2y)=4 !4 for island defect
            endif

            if(initrand) then
                do i=1,ntriang
                    checkaccept = .false.
                    do while(.not.checkaccept)

                        ! Random initialisation j = n + FLOOR((m+1-n)*rand) for interval [n,m]
100                     CALL RANDOM_NUMBER(rand)
                        config1(i,1)=floor(rand*color)+1 !random color 1 to 4
                        CALL RANDOM_NUMBER(rand)
                        config1(i,2)=floor(rand*posNx)+1 !random posNx 1 to 5 ,plot is from 0 to 6
                        CALL RANDOM_NUMBER(rand)
                        config1(i,3)=floor(rand*posNy)+1 !random posNy 1 to 5
                        CALL RANDOM_NUMBER(rand)
                        config1(i,4)=floor(rand*molRot) !random molRot 0 to 2
                        CALL RANDOM_NUMBER(rand)
                        config1(i,5)=floor(rand*osmov) !random molRot 0 to 3

                        print*,config1(i,:),'config1(i,:)'

                        !colors 1 green, 2 red, 3 yellow, 4 blu
                        if(config1(i,1).eq.1) then
                            !H O N green

                            dx = modulo(config1(i,2)-1,posNx)+1
                            dy = modulo(config1(i,3)-1,posNy)+1
                            dx1 = modulo(config1(i,2)-1,posNx)+1
                            dy1 = modulo(config1(i,3)-2,posNy)+1
                            dx2 = modulo(config1(i,2)-2,posNx)+1
                            dy2 = modulo(config1(i,3)-2,posNy)+1


                            if((emat(dx,dy).eq.0).and.(emat(dx1,dy1).eq.0).and.(emat(dx2,dy2).eq.0)) then
                                newinit(i,1)=dx2  !pos grid x
                                newinit(i,2)=dy2  !pos grid y
                                newinit(i,3)=1 !color

                                emat(dx,dy)=1
                                emat(dx1,dy1)=2
                                emat(dx2,dy2)=2

                               checkmat(i,1,1)=dx
                               checkmat(i,2,1)=dy
                               checkmat(i,1,2)=dx1
                               checkmat(i,2,2)=dy1
                               checkmat(i,1,3)=dx2
                               checkmat(i,2,3)=dy2

                            indexmat(dx,dy)=i !save triangle index
                            indexmat(dx1,dy1)=i !save triangle index
                            indexmat(dx2,dy2)=i !save triangle index


                                checkaccept = .true.
                            !print*,'green accepted'

                            else
                                checkaccept = .false.
                            endif

                        elseif(config1(i,1).eq.2) then
                            !H O N red

                            dx = modulo(config1(i,2)-1,posNx)+1
                            dy = modulo(config1(i,3)-1,posNy)+1
                            dx1 = modulo(config1(i,2)-1,posNx)+1
                            dy1 = modulo(config1(i,3)-2,posNy)+1
                            dx2 = modulo(config1(i,2),posNx)+1
                            dy2 = modulo(config1(i,3)-2,posNy)+1

                            if((emat(dx,dy).eq.0).and.(emat(dx1,dy1).eq.0).and.(emat(dx2,dy2).eq.0)) then
                                newinit(i,1)=dx  !pos grid x
                                newinit(i,2)=dy1  !pos grid y
                                newinit(i,3)=2 !color

                                emat(dx,dy)=1
                                emat(dx1,dy1)=2
                                emat(dx2,dy2)=2

                               checkmat(i,1,1)=dx
                               checkmat(i,2,1)=dy
                               checkmat(i,1,2)=dx1
                               checkmat(i,2,2)=dy1
                               checkmat(i,1,3)=dx2
                               checkmat(i,2,3)=dy2

                            indexmat(dx,dy)=i !save triangle index
                            indexmat(dx1,dy1)=i !save triangle index
                            indexmat(dx2,dy2)=i !save triangle index

                                checkaccept = .true.
                            ! print*,'red accepted'

                            else
                                checkaccept = .false.
                            endif


                        elseif(config1(i,1).eq.3) then
                             !H O N yelloW

                            dx = modulo(config1(i,2)-1,posNx)+1
                            dy = modulo(config1(i,3)-1,posNy)+1
                            dx1 = modulo(config1(i,2)-1,posNx)+1
                            dy1 = modulo(config1(i,3),posNy)+1
                            dx2 = modulo(config1(i,2)-2,posNx)+1
                            dy2 = modulo(config1(i,3),posNy)+1

                            if((emat(dx,dy).eq.0).and.(emat(dx1,dy1).eq.0).and.(emat(dx2,dy2).eq.0)) then
                                newinit(i,1)=dx2  !pos grid x
                                newinit(i,2)=dy  !pos grid y
                                newinit(i,3)=3 !yellow

                                emat(dx,dy)=1
                                emat(dx1,dy1)=2
                                emat(dx2,dy2)=2

                               checkmat(i,1,1)=dx
                               checkmat(i,2,1)=dy
                               checkmat(i,1,2)=dx1
                               checkmat(i,2,2)=dy1
                               checkmat(i,1,3)=dx2
                               checkmat(i,2,3)=dy2

                            indexmat(dx,dy)=i !save triangle index
                            indexmat(dx1,dy1)=i !save triangle index
                            indexmat(dx2,dy2)=i !save triangle index

                                checkaccept = .true.
                             ! print*,'yellow accepted',config1(i,2),config1(i,3)

                            else
                                checkaccept = .false.
                            endif


                        elseif(config1(i,1).eq.4) then
                            !H O N !BLUE
                            dx = modulo(config1(i,2)-1,posNx)+1
                            dy = modulo(config1(i,3)-1,posNy)+1
                            dx1 = modulo(config1(i,2)-1,posNx)+1
                            dy1 = modulo(config1(i,3),posNy)+1
                            dx2 = modulo(config1(i,2),posNx)+1
                            dy2 = modulo(config1(i,3),posNy)+1
                            if((emat(dx,dy).eq.0).and.(emat(dx1,dy1).eq.0).and.(emat(dx2,dy2).eq.0)) then
                                newinit(i,1)=dx  !pos grid x
                                newinit(i,2)=dy  !pos grid y
                                newinit(i,3)=4 !color
                                emat(dx,dy)=1
                                emat(dx1,dy1)=2
                                emat(dx2,dy2)=2

                               checkmat(i,1,1)=dx
                               checkmat(i,2,1)=dy
                               checkmat(i,1,2)=dx1
                               checkmat(i,2,2)=dy1
                               checkmat(i,1,3)=dx2
                               checkmat(i,2,3)=dy2

                            indexmat(dx,dy)=i !save triangle index
                            indexmat(dx1,dy1)=i !save triangle index
                            indexmat(dx2,dy2)=i !save triangle index

                                checkaccept = .true.
                               ! print*,'blue accepted'
                            else
                                checkaccept = .false.
                            endif
                        endif

                        if(.not.checkaccept)then
                            print*,'not accepted, trianlge',i
                            goto 100
                            stop
                        endif

                    enddo
                enddo

            else

                print*,'initrand false'

                i=1
                if((isltype1.eq.1).or.(isltype1.eq.2))then
                    !        checkaccept = .false.
                    !Number of triangles in row or column
                    do j=1,isl1sizex
                        do k=1,isl1sizey
                            if(i.gt.ntriang) then
                                print*, 'exceeded number of triangles in setup island 1'
                                stop
                            endif

                            if(isltype1.eq.1) then
                                if(modulo(k,2).eq.1) then
                                    if(colorisl1.eq.12) then
                                        config1(i,1)=1 !color 1 to 4
                                    elseif(colorisl1.eq.34) then
                                        config1(i,1)=3 !color 1 to 4
                                    else
                                        print*,'wrong color for cluster, 12 or 34'
                                        stop
                                    endif
                                    config1(i,2)=3*j-1+isl1posx-2
                                    config1(i,3)=k+isl1posy-1
                                else
                                    if(colorisl1.eq.12) then
                                        config1(i,1)=2 !color 1 to 4
                                    elseif(colorisl1.eq.34) then
                                        config1(i,1)=4 !color 1 to 4
                                    else
                                        print*,'wrong color for cluster, 12 or 34'
                                        stop
                                    endif
                                    config1(i,2)=3*j+isl1posx-2
                                    config1(i,3)=k+isl1posy-1
                                endif
                            elseif(isltype1.eq.2) then
                                config1(i,1)=colorisl1
                                config1(i,2)=2*(j-1)+isl1posx
                                config1(i,3)=2*(k-1)+isl1posy
                            endif


                            CALL RANDOM_NUMBER(rand)
                            config1(i,4)=int(rand*molRot) !random molRot 0 to 2
                            CALL RANDOM_NUMBER(rand)
                            config1(i,5)=int(rand*osmov) !random molRot 0 to 3



                            !colors 1 green, 2 red, 3 yellow, 4 blu
                            if(config1(i,1).eq.1) then
                                !H O N green

                                dx = modulo(config1(i,2)-1,posNx)+1
                                dy = modulo(config1(i,3)-1,posNy)+1
                                dx1 = modulo(config1(i,2)-1,posNx)+1
                                dy1 = modulo(config1(i,3)-2,posNy)+1
                                dx2 = modulo(config1(i,2)-2,posNx)+1
                                dy2 = modulo(config1(i,3)-2,posNy)+1



                                if((emat(dx,dy).eq.0).and.(emat(dx1,dy1).eq.0).and.(emat(dx2,dy2).eq.0)) then
                                    newinit(i,1)=dx2  !pos grid x
                                    newinit(i,2)=dy2  !pos grid y
                                    newinit(i,3)=1 !color

                                    emat(dx,dy)=1
                                    emat(dx1,dy1)=2
                                    emat(dx2,dy2)=2

                               checkmat(i,1,1)=dx
                               checkmat(i,2,1)=dy
                               checkmat(i,1,2)=dx1
                               checkmat(i,2,2)=dy1
                               checkmat(i,1,3)=dx2
                               checkmat(i,2,3)=dy2

                            indexmat(dx,dy)=i !save triangle index
                            indexmat(dx1,dy1)=i !save triangle index
                            indexmat(dx2,dy2)=i !save triangle index


                                else
                                    if((emat(dx,dy).eq.4).or.(emat(dx1,dy1).eq.4).or.(emat(dx2,dy2).eq.4)) then
                                        Print*,'defects in islands involved, reset i count'
                                        i=i-1
                                    else
                                        Print*,'warning: particles overlap or metal defect is in the way 1'
                                        !print*,emat(dx,dy),emat(dx1,dy1),emat(dx2,dy2),'emat dx dx1 dx2'
                                        stop
                                    endif
                                endif

                            elseif(config1(i,1).eq.2) then
                                !H O N red

                                dx = modulo(config1(i,2)-1,posNx)+1
                                dy = modulo(config1(i,3)-1,posNy)+1
                                dx1 = modulo(config1(i,2)-1,posNx)+1
                                dy1 = modulo(config1(i,3)-2,posNy)+1
                                dx2 = modulo(config1(i,2),posNx)+1
                                dy2 = modulo(config1(i,3)-2,posNy)+1

                                if((emat(dx,dy).eq.0).and.(emat(dx1,dy1).eq.0).and.(emat(dx2,dy2).eq.0)) then
                                    newinit(i,1)=dx  !pos grid x
                                    newinit(i,2)=dy1  !pos grid y
                                    newinit(i,3)=2 !color

                                    emat(dx,dy)=1
                                    emat(dx1,dy1)=2
                                    emat(dx2,dy2)=2

                               checkmat(i,1,1)=dx
                               checkmat(i,2,1)=dy
                               checkmat(i,1,2)=dx1
                               checkmat(i,2,2)=dy1
                               checkmat(i,1,3)=dx2
                               checkmat(i,2,3)=dy2

                            indexmat(dx,dy)=i !save triangle index
                            indexmat(dx1,dy1)=i !save triangle index
                            indexmat(dx2,dy2)=i !save triangle index

                                else
                                    if((emat(dx,dy).eq.4).or.(emat(dx1,dy1).eq.4).or.(emat(dx2,dy2).eq.4)) then
                                        Print*,'defects in islands involved, reset i count'
                                        i=i-1
                                    else
                                        Print*,'warning: particles overlap or metal defect is in the way 2'
                                        !print*,emat(dx,dy),emat(dx1,dy1),emat(dx2,dy2),'emat dx dx1 dx2'
                                        stop
                                    endif
                                endif

                            elseif(config1(i,1).eq.3) then
                                 !H O N yelloW

                                dx = modulo(config1(i,2)-1,posNx)+1
                                dy = modulo(config1(i,3)-1,posNy)+1
                                dx1 = modulo(config1(i,2)-1,posNx)+1
                                dy1 = modulo(config1(i,3),posNy)+1
                                dx2 = modulo(config1(i,2)-2,posNx)+1
                                dy2 = modulo(config1(i,3),posNy)+1

                                if((emat(dx,dy).eq.0).and.(emat(dx1,dy1).eq.0).and.(emat(dx2,dy2).eq.0)) then
                                    newinit(i,1)=dx2  !pos grid x
                                    newinit(i,2)=dy  !pos grid y
                                    newinit(i,3)=3 !yellow

                                    emat(dx,dy)=1
                                    emat(dx1,dy1)=2
                                    emat(dx2,dy2)=2

                               checkmat(i,1,1)=dx
                               checkmat(i,2,1)=dy
                               checkmat(i,1,2)=dx1
                               checkmat(i,2,2)=dy1
                               checkmat(i,1,3)=dx2
                               checkmat(i,2,3)=dy2



                            indexmat(dx,dy)=i !save triangle index
                            indexmat(dx1,dy1)=i !save triangle index
                            indexmat(dx2,dy2)=i !save triangle index


                                else
                                    if((emat(dx,dy).eq.4).or.(emat(dx1,dy1).eq.4).or.(emat(dx2,dy2).eq.4)) then
                                        Print*,'defects in islands involved, reset i count'
                                        i=i-1
                                    else
                                        Print*,'warning: particles overlap or metal defect is in the way 3'
                                        !print*,emat(dx,dy),dx,dy,emat(dx1,dy1),dx1,dy1,emat(dx2,dy2),dx2,dy2,'emat dx dx1 dx2',i
                                        stop
                                    endif
                                endif


                            elseif(config1(i,1).eq.4) then
                                !H O N !BLUE
                                dx = modulo(config1(i,2)-1,posNx)+1
                                dy = modulo(config1(i,3)-1,posNy)+1
                                dx1 = modulo(config1(i,2)-1,posNx)+1
                                dy1 = modulo(config1(i,3),posNy)+1
                                dx2 = modulo(config1(i,2),posNx)+1
                                dy2 = modulo(config1(i,3),posNy)+1
                                if((emat(dx,dy).eq.0).and.(emat(dx1,dy1).eq.0).and.(emat(dx2,dy2).eq.0)) then
                                    newinit(i,1)=dx  !pos grid x
                                    newinit(i,2)=dy  !pos grid y
                                    newinit(i,3)=4 !color
                                    emat(dx,dy)=1
                                    emat(dx1,dy1)=2
                                    emat(dx2,dy2)=2

                               checkmat(i,1,1)=dx
                               checkmat(i,2,1)=dy
                               checkmat(i,1,2)=dx1
                               checkmat(i,2,2)=dy1
                               checkmat(i,1,3)=dx2
                               checkmat(i,2,3)=dy2



                            indexmat(dx,dy)=i !save triangle index
                            indexmat(dx1,dy1)=i !save triangle index
                            indexmat(dx2,dy2)=i !save triangle index

                                else
                                    if((emat(dx,dy).eq.4).or.(emat(dx1,dy1).eq.4).or.(emat(dx2,dy2).eq.4)) then
                                        Print*,'defects in islands involved, reset i count'
                                        i=i-1
                                    else
                                        Print*,'warning: particles overlap or metal defect is in the way 4'
                                        !print*,emat(dx,dy),emat(dx1,dy1),emat(dx2,dy2),'emat dx dx1 dx2'
                                        stop
                                    endif
                                endif
                            endif

                            !print*,emat(dx,dy),dx,dy,emat(dx1,dy1),dx1,dy1,emat(dx2,dy2),dx2,dy2,'emat dx dx1 dx2',i

                            i=i+1 !update triangle count
                        enddo
                    enddo
                endif


                if(isl1dfct.ne.0) then
                    if(isl1dfct.eq.1) then
                        if((i-1).ne.isl1sizex*isl1sizey-1) then
                            print*,'defect in island 1 was not hit by a molecule: CHANGE',&
                                i-1,isl1sizex*isl1sizey-1
                            stop
                        endif
                    elseif(isl1dfct.eq.2) then
                        if((i-1).ne.isl1sizex*isl1sizey-2) then
                            print*,'defect 1 or 2 in island 1 was not hit by a molecule: CHANGE',&
                                i-1,isl1sizex*isl1sizey-2
                            stop
                        endif
                    endif
                endif


                if((isltype2.eq.1).or.(isltype2.eq.2)) then
                    !        checkaccept = .false.
                    !Number of triangles in row or column
                    do j=1,isl2sizex
                        do k=1,isl2sizey
                            if(i.gt.ntriang) then
                                print*, 'exceeded number of triangles in setup island2', &
                                    'check if defects are actually hit'
                                stop
                            endif

                                  ! Random initialisation int(rand(0)*(j+1-i))+i for interval [i,j]
                            if(isltype2.eq.1) then
                                if(modulo(k,2).eq.1) then
                                    if(colorisl2.eq.12) then
                                        config1(i,1)=1 !color 1 to 4
                                    elseif(colorisl2.eq.34) then
                                        config1(i,1)=3 !color 1 to 4
                                    else
                                        print*,'wrong color for cluster, 12 or 34'
                                        stop
                                    endif
                                    config1(i,2)=3*j-1+isl2posx-2
                                    config1(i,3)=k+isl2posy-1
                                else
                                    if(colorisl2.eq.12) then
                                        config1(i,1)=2 !color 1 to 4
                                    elseif(colorisl2.eq.34) then
                                        config1(i,1)=4 !color 1 to 4
                                    else
                                        print*,'wrong color for cluster, 12 or 34'
                                        stop
                                    endif
                                    config1(i,2)=3*j+isl2posx-2
                                    config1(i,3)=k+isl2posy-1
                                endif
                            elseif(isltype2.eq.2) then
                                config1(i,1)=colorisl2
                                config1(i,2)=2*(j-1)+isl2posx !random posNx 1 to 5 ,plot is from 0 to 6
                                config1(i,3)=2*(k-1)+isl2posy !random posNy 1 to 5
                            endif

                            CALL RANDOM_NUMBER(rand)
                            config1(i,4)=int(rand*molRot) !random molRot 0 to 2
                            CALL RANDOM_NUMBER(rand)
                            config1(i,5)=int(rand*osmov) !random molRot 0 to 3


                            !colors 1 green, 2 red, 3 yellow, 4 blu
                            if(config1(i,1).eq.1) then
                                !H O N green

                                dx = modulo(config1(i,2)-1,posNx)+1
                                dy = modulo(config1(i,3)-1,posNy)+1
                                dx1 = modulo(config1(i,2)-1,posNx)+1
                                dy1 = modulo(config1(i,3)-2,posNy)+1
                                dx2 = modulo(config1(i,2)-2,posNx)+1
                                dy2 = modulo(config1(i,3)-2,posNy)+1


                                if((emat(dx,dy).eq.0).and.(emat(dx1,dy1).eq.0).and.(emat(dx2,dy2).eq.0)) then
                                    newinit(i,1)=dx2  !pos grid x
                                    newinit(i,2)=dy2  !pos grid y
                                    newinit(i,3)=1 !color

                                    emat(dx,dy)=1
                                    emat(dx1,dy1)=2
                                    emat(dx2,dy2)=2

                               checkmat(i,1,1)=dx
                               checkmat(i,2,1)=dy
                               checkmat(i,1,2)=dx1
                               checkmat(i,2,2)=dy1
                               checkmat(i,1,3)=dx2
                               checkmat(i,2,3)=dy2

                            indexmat(dx,dy)=i !save triangle index
                            indexmat(dx1,dy1)=i !save triangle index
                            indexmat(dx2,dy2)=i !save triangle index


                                else
                                    if((emat(dx,dy).eq.4).or.(emat(dx1,dy1).eq.4).or.(emat(dx2,dy2).eq.4)) then
                                        Print*,'defects in islands involved, reset i count'
                                        i=i-1
                                    else
                                        Print*,'warning: particles overlap or metal defect is in the way 5'
                                        !print*,emat(dx,dy),emat(dx1,dy1),emat(dx2,dy2),'emat dx dx1 dx2'
                                        stop
                                    endif
                                endif

                            elseif(config1(i,1).eq.2) then
                                !H O N red

                                dx = modulo(config1(i,2)-1,posNx)+1
                                dy = modulo(config1(i,3)-1,posNy)+1
                                dx1 = modulo(config1(i,2)-1,posNx)+1
                                dy1 = modulo(config1(i,3)-2,posNy)+1
                                dx2 = modulo(config1(i,2),posNx)+1
                                dy2 = modulo(config1(i,3)-2,posNy)+1

                                if((emat(dx,dy).eq.0).and.(emat(dx1,dy1).eq.0).and.(emat(dx2,dy2).eq.0)) then
                                    newinit(i,1)=dx  !pos grid x
                                    newinit(i,2)=dy1  !pos grid y
                                    newinit(i,3)=2 !color

                                    emat(dx,dy)=1
                                    emat(dx1,dy1)=2
                                    emat(dx2,dy2)=2

                               checkmat(i,1,1)=dx
                               checkmat(i,2,1)=dy
                               checkmat(i,1,2)=dx1
                               checkmat(i,2,2)=dy1
                               checkmat(i,1,3)=dx2
                               checkmat(i,2,3)=dy2

                            indexmat(dx,dy)=i !save triangle index
                            indexmat(dx1,dy1)=i !save triangle index
                            indexmat(dx2,dy2)=i !save triangle index

                                    checkaccept = .true.
                                else
                                    if((emat(dx,dy).eq.4).or.(emat(dx1,dy1).eq.4).or.(emat(dx2,dy2).eq.4)) then
                                        Print*,'defects in islands involved, reset i count'
                                        i=i-1
                                    else
                                        Print*,'warning: particles overlap or metal defect is in the way 6'
                                        !print*,emat(dx,dy),emat(dx1,dy1),emat(dx2,dy2),'emat dx dx1 dx2'
                                        stop
                                    endif
                                endif


                            elseif(config1(i,1).eq.3) then
                                 !H O N yelloW

                                dx = modulo(config1(i,2)-1,posNx)+1
                                dy = modulo(config1(i,3)-1,posNy)+1
                                dx1 = modulo(config1(i,2)-1,posNx)+1
                                dy1 = modulo(config1(i,3),posNy)+1
                                dx2 = modulo(config1(i,2)-2,posNx)+1
                                dy2 = modulo(config1(i,3),posNy)+1

                                if((emat(dx,dy).eq.0).and.(emat(dx1,dy1).eq.0).and.(emat(dx2,dy2).eq.0)) then
                                    newinit(i,1)=dx2  !pos grid x
                                    newinit(i,2)=dy  !pos grid y
                                    newinit(i,3)=3 !yellow

                                    emat(dx,dy)=1
                                    emat(dx1,dy1)=2
                                    emat(dx2,dy2)=2

                               checkmat(i,1,1)=dx
                               checkmat(i,2,1)=dy
                               checkmat(i,1,2)=dx1
                               checkmat(i,2,2)=dy1
                               checkmat(i,1,3)=dx2
                               checkmat(i,2,3)=dy2

                            indexmat(dx,dy)=i !save triangle index
                            indexmat(dx1,dy1)=i !save triangle index
                            indexmat(dx2,dy2)=i !save triangle index

                                else
                                    if((emat(dx,dy).eq.4).or.(emat(dx1,dy1).eq.4).or.(emat(dx2,dy2).eq.4)) then
                                        Print*,'defects in islands involved, reset i count'
                                        i=i-1
                                    else
                                        Print*,'warning: particles overlap or metal defect is in the way 7'
                                        !print*,emat(dx,dy),emat(dx1,dy1),emat(dx2,dy2),'emat dx dx1 dx2'
                                        stop
                                    endif
                                endif


                            elseif(config1(i,1).eq.4) then
                                !H O N !BLUE
                                dx = modulo(config1(i,2)-1,posNx)+1
                                dy = modulo(config1(i,3)-1,posNy)+1
                                dx1 = modulo(config1(i,2)-1,posNx)+1
                                dy1 = modulo(config1(i,3),posNy)+1
                                dx2 = modulo(config1(i,2),posNx)+1
                                dy2 = modulo(config1(i,3),posNy)+1
                                if((emat(dx,dy).eq.0).and.(emat(dx1,dy1).eq.0).and.(emat(dx2,dy2).eq.0)) then
                                    newinit(i,1)=dx  !pos grid x
                                    newinit(i,2)=dy  !pos grid y
                                    newinit(i,3)=4 !color
                                    emat(dx,dy)=1
                                    emat(dx1,dy1)=2
                                    emat(dx2,dy2)=2

                               checkmat(i,1,1)=dx
                               checkmat(i,2,1)=dy
                               checkmat(i,1,2)=dx1
                               checkmat(i,2,2)=dy1
                               checkmat(i,1,3)=dx2
                               checkmat(i,2,3)=dy2

                            indexmat(dx,dy)=i !save triangle index
                            indexmat(dx1,dy1)=i !save triangle index
                            indexmat(dx2,dy2)=i !save triangle index

                                else
                                    if((emat(dx,dy).eq.4).or.(emat(dx1,dy1).eq.4).or.(emat(dx2,dy2).eq.4)) then
                                        Print*,'defects in islands involved, reset i count'
                                        i=i-1
                                    else
                                        Print*,'warning: particles overlap or metal defect is in the way 8'
                                        !print*,emat(dx,dy),emat(dx1,dy1),emat(dx2,dy2),'emat dx dx1 dx2'
                                        stop
                                    endif
                                endif
                            endif

                            i=i+1 !update triangle count
                        enddo
                    enddo
                endif
            endif !initrand
        elseif(restart) then

            print*,'reading file: ',restart_file,'for restart (do not change grid size).'
            call read_vtk_triangles(ntriang,config1,newinit,emat,indexmat,posNx,posNy,restart_file)

            !todo call routine for checkmat generation from config for ALL triangles
            print*,'todo restart checkmat'

            !Writing point surf
            if((mdef1type.eq.0).and.(mdef2type.eq.0)) then
                call write_vtk_points(1,1,posNx,posNy)
            else
                call write_vtk_points_defect(1,1,posNx,posNy,&
                    mdef1type,mdef2type,mdef1sizex,mdef1sizey,mdef2sizex,mdef2sizey,&
                    mdef1posx,mdef1posy,mdef2posx,mdef2posy)
                !emat update
                if(mdef1type.ne.0)then
                    do i=1,mdef1sizex
                        do j=1,mdef1sizey
                            emat(mdef1posx-1+i,mdef1posy-1+j)=3 !3 for metal defect
                            defectmat(mdef1posx-1+i,mdef1posy-1+j)=3
                        enddo
                    enddo
                endif
                if(mdef2type.ne.0)then
                    do i=1,mdef2sizex
                        do j=1,mdef2sizey
                            emat(mdef2posx-1+i,mdef2posy-1+j)=3 !3 for metal defect
                            defectmat(mdef2posx-1+i,mdef2posy-1+j)=3
                        enddo
                    enddo
                endif
            endif

        do j=1,ntriang
        call checkmatgen(config1(j,:),checkmat(j,:,:))
        enddo

        endif

        WRITE(data_filename,'(a)')'initialpos.vtk'
        call vtk_triangles(ntriang,config1,newinit,data_filename)

        WRITE(data_filename,'(a)')'COMinit.vtk'
        call centreofmasscompare(ntriang,newinit,data_filename)


        !initial H-bond check
        do j=1,posNx
            do k=1,posNy
                if(emat(j,k).ne.0) occupationcount=occupationcount+1
                if(emat(j,k).eq.2) then
                    jcheck1=modulo(j-1,posNx)
                    if(jcheck1.eq.0) jcheck1 = posNx
                    jcheck2=modulo(j-1,posNx)+2
                    if(jcheck2.eq.posNx+1) jcheck2 = 1
                    if(emat(jcheck1,k).eq.1) then
                        hbondmat(jcheck1,k)=1
                        hbondcount=hbondcount+1
                    elseif(emat(jcheck2,k).eq.1) then
                        hbondmat(j,k)=1
                        hbondcount=hbondcount+1
                    endif
                endif
            enddo
        enddo
        !periodic boundaries, double row and line 1
        hbondmat(posNx+1,1:posNy)=hbondmat(1,1:posNy)
        hbondmat(1:posNx,posNy+1)=hbondmat(1:posNx,1)
        hbondmat(posNx+1,posNy+1)=hbondmat(1,1)

        WRITE(data_filename,'(a,i9.9,a)')'hlinesinit.vtk'
        call write_vtk_lines( 1,1,posNx,posNy,hbondmat,data_filename)

    endif



    !****** END Initialisation ****************
    !for testing the number of possibilites


    IF (STEPS.EQ.1) then  ! create rate list with probablity to move for each triangle

    config1OLD(:,:)=config1(:,:)

    do triang=1,ntriang
    eventcount=0
    dx = config1(triang,2)
    dy = config1(triang,3)

    do i=1,color
        do j= dx-2,dx+2 ! Saving looping time
            !do j= -1,posNx+2 !start at 0
            do k=dy-3,dy+3! Saving looping time
                !do k=-1,posNy+2 !start at 0
                do l=0,molRot-1 !start at 0
                    do m=0,osmov-1 !start at 0
                        config2(triang,:)=(/ i,j,k,l,m /)
                        call glyevents(ntriang,triang,posNx,posNy,config1,config2, &
                            rateconst,eventOK)
                        if(eventOK) then
                            eventcount=eventcount+1
                            RATES(triang,eventcount)=(10.D0**(13))*dexp(-rateconst/(kboltz*Temp))
                            SAVconfig(eventcount,triang,:)=config2(triang,:)
                        endif
                    enddo
                enddo
            enddo
        enddo
    enddo


    NEVENTS(triang)=eventcount !preparation for David's KMC routine

    enddo ! End of triangle scans for list



    else ! Else for Step 2 and higher (updating the rate list)

    dx = config1(indexlw,2)
    dy = config1(indexlw,3)

    do chosentriang=1,ntriang
            if(config1(indexlw,1).eq.1) then !green
                if(modulo(config1(chosentriang,3)-(dy-4),posNy).le.6) then
                    if(modulo(config1(chosentriang,2)-(dx-4),posNx).le.7) then
                        neighbourindex(chosentriang)=1
                    endif
                endif
            elseif(config1(indexlw,1).eq.2) then !red
                if(modulo(config1(chosentriang,3)-(dy-4),posNy).le.6) then
                    if(modulo(config1(chosentriang,2)-(dx-3),posNx).le.7) then
                        neighbourindex(chosentriang)=1
                    endif
                endif
            elseif(config1(indexlw,1).eq.3) then !yellow
                if(modulo(config1(chosentriang,3)-(dy-2),posNy).le.6) then
                    if(modulo(config1(chosentriang,2)-(dx-4),posNx).le.7) then
                        neighbourindex(chosentriang)=1
                    endif
                endif
            elseif(config1(indexlw,1).eq.4) then !blue
                if(modulo(config1(chosentriang,3)-(dy-2),posNy).le.6) then
                    if(modulo(config1(chosentriang,2)-(dx-3),posNx).le.7) then
                        neighbourindex(chosentriang)=1
                    endif
                endif
            else
                print*,'wrong colour code',config1(chosentriang,1),'config1(chosentriang,1)',chosentriang,'chosentriang'
                stop
            endif

            if(neighbourindex(chosentriang).eq.1) then

                eventcount=0 !Reset eventcount for each triangle

                dx1 = config1(chosentriang,2) !modulo spater?
                dy1 = config1(chosentriang,3)

    do i=1,color
        do j= dx1-2,dx1+2 ! Saving looping time
            do k=dy1-3,dy1+3! Saving looping time
                do l=0,molRot-1 !start at 0
                    do m=0,osmov-1 !start at 0

                        config2(chosentriang,:)=(/ i,j,k,l,m /)

                        call glyevents(ntriang,chosentriang,posNx,posNy,config1,config2, &
                            rateconst,eventOK)

                        if (eventOK) then
                            eventcount=eventcount+1
                            RATES(chosentriang,eventcount)=(10.D0**(13))*dexp(-rateconst/(kboltz*Temp))
                            SAVconfig(eventcount,chosentriang,:)=config2(chosentriang,:)
                        endif
                    enddo
                enddo
            enddo
        enddo
    enddo

    NEVENTS(chosentriang)=eventcount !Number of possible events


    endif

    enddo

    endif


    RETURN

END SUBROUTINE GETEVENTS


subroutine glyevents(ntriang,chosentriang,posNx,posNy,config1,config2, &
    rateconst,eventOK)
    USE COMMONS, only: emat
    implicit none
    integer, intent(in) :: ntriang,config1(ntriang,5),config2(ntriang,5),&
        posNx,posNy,chosentriang
    logical, intent(out) :: eventOK
    double precision, intent(out) :: rateconst

    integer :: deltax, deltay,icount
    double precision::deltarconst
    integer :: dx,dy,dx1,dy1,dx2,dy2,hbondcount,hbondmatnew(posNx+1,posNy+1)



    !************************************
    rateconst=0.D0 !rateconst is a barrier here, changes in Getevents routine
    eventOK=.false.
    deltarconst=0.D0
    hbondcount=0
    hbondmatnew(:,:)=0
    !************************************

    do icount=chosentriang,chosentriang !1,ntriang !count over all triangles

        deltax=config2(icount,2)-config1(icount,2)
        deltay=config2(icount,3)-config1(icount,3)


        if(abs(deltax).gt.2) then
            !print*, 'step not allowed'
            !print*, 'abs(deltax) > 2',deltax,'deltax'
            return
        endif
        if(abs(deltay).gt.2) then
            !print*, 'step not allowed'
            !print*, 'abs(deltay) > 2',deltay,'deltay'
            return
        endif


            !colors 1 green, 2 red, 3 yellow, 4 blu
        if(config1(icount,1).eq.1) then !colorcheck : GREEN moves
            if(config2(icount,1).eq.1) then !from green to green
                !print*, 'step not allowed'
                return
            elseif((config2(icount,1).eq.2).AND.(config2(icount,4).eq.0)) then !from green to red

                if(deltay.eq.0) then !deltay checking
                    if((deltax.eq.0).AND.(config2(icount,5).ne.0))  then !deltax checking
                        !Periodic boundaries
                        dx = modulo(config1(icount,2),posNx)+1
                        dy = modulo(config1(icount,3)-2,posNy)+1
                        if(emat(dx,dy).eq.0) then
                            if(config2(icount,5).eq.1) then ! Os translation
                                rateconst = 0.372D0
                                eventOK=.true.
                                call ematandhbondsC(config1(icount,2),config1(icount,3), &
                                    config1(icount,1),posNx,posNy,deltarconst)!,ematnew)

                                rateconst = rateconst - deltarconst
                            elseif(config2(icount,5).eq.2) then !Os flipping
                                rateconst = 0.537D0
                                eventOK=.true.
                                call ematandhbondsC(config1(icount,2),config1(icount,3), &
                                    config1(icount,1),posNx,posNy,deltarconst)!,ematnew)

                                rateconst = rateconst - deltarconst
                            elseif(config2(icount,5).eq.3) then !Os back-flipping
                                rateconst = 0.724D0
                                eventOK=.true.
                                call ematandhbondsC(config1(icount,2),config1(icount,3), &
                                    config1(icount,1),posNx,posNy,deltarconst)!,ematnew)

                                rateconst = rateconst - deltarconst
                            endif
                        else
                            !print*, 'step not allowed emat not empty (position 1)'
                            return
                        endif

                    elseif(deltax.eq.-1 .AND.config2(icount,5).eq.0)  then !deltax checking

                        dx = modulo(config1(icount,2)-2,posNx)+1
                        dy = modulo(config1(icount,3)-1,posNy)+1
                        if(emat(dx,dy).eq.0) then

                            rateconst = 0.358D0
                            eventOK=.true.
                            call ematandhbondsD(config1(icount,2),config1(icount,3), &
                                config1(icount,1),posNx,posNy,deltarconst)!,ematnew)

                            rateconst = rateconst - deltarconst

                        else
                            !print*, 'step not allowed emat not empty (position 4)'
                            return
                        endif

                    elseif(deltax.eq.-2.AND.config2(icount,5).eq.0)  then !deltax checking


                        dx = modulo(config1(icount,2)-3,posNx)+1
                        dy = modulo(config1(icount,3)-1,posNy)+1
                        dx1 = modulo(config1(icount,2)-3,posNx)+1
                        dy1 = modulo(config1(icount,3)-2,posNy)+1

                        if((emat(dx,dy).eq.0).and.(emat(dx1,dy1).eq.0)) then

                            rateconst = 0.582D0
                            eventOK=.true.

                            call ematandhbondsF(config1(icount,2),config1(icount,3), &
                                config1(icount,1),posNx,posNy,deltarconst)!,ematnew)
                            rateconst = rateconst - deltarconst
                        else
                            !print*, 'step not allowed emat not empty (position 5 and 6)'
                            return
                        endif


                    else
                        !print*, 'step not allowed'
                        return
                    endif
                else
                    !print*, 'step not allowed'
                    return
                endif
            elseif(config2(icount,1).eq.3 .AND.config2(icount,5).eq.0) then !from green to yellow

                if(deltax.eq.0)  then !deltax checking
                    if(deltay.eq.-2 .AND. config2(icount,4).eq.0)  then !deltay checking

                        dx = modulo(config1(icount,2)-1,posNx)+1
                        dy = modulo(config1(icount,3)-3,posNy)+1

                        if((emat(dx,dy).eq.0)) then
                            rateconst = 0.384D0
                            eventOK=.true.
                            call ematandhbondsA(config1(icount,2),config1(icount,3), &
                                config1(icount,1),posNx,posNy,deltarconst)!,ematnew)
                            rateconst = rateconst - deltarconst

                        else
                            !print*, 'step not allowed emat not empty (position 2)'
                            return
                        endif

                    elseif(deltay.eq.-1) then
                        dx = modulo(config1(icount,2)-2,posNx)+1
                        dy = modulo(config1(icount,3)-1,posNy)+1
                        if(emat(dx,dy).eq.0) then
                            if(config2(icount,4).eq.1) then !rotation plus flipping
                                rateconst = 0.582D0
                                eventOK=.true.
                                call ematandhbondsE(config1(icount,2),config1(icount,3), &
                                    config1(icount,1),posNx,posNy,deltarconst)!,ematnew)
                                rateconst = rateconst - deltarconst
                            elseif(config2(icount,4).eq.2) then !pure rotation
                                rateconst = 0.610D0
                                eventOK=.true.
                                call ematandhbondsE(config1(icount,2),config1(icount,3), &
                                    config1(icount,1),posNx,posNy,deltarconst)!,ematnew)
                                rateconst = rateconst - deltarconst
                            else
                                ! print*, 'step not allowed'
                                return
                            endif
                        else
                            !print*, 'step not allowed emat not empty (position 4)'
                            return
                        endif

                    else
                        !print*, 'step not allowed'
                        return
                    endif
                else
                    !print*, 'step not allowed'
                    return
                endif
            elseif((config2(icount,1).eq.4).AND. (config2(icount,4).eq.0).AND. (config2(icount,5).eq.0)) then !from green to blue

                if(deltax.eq.-1)  then !deltax checking
                    if(deltay.eq.-2)  then !deltay checking

                        dx = modulo(config1(icount,2)-2,posNx)+1
                        dy = modulo(config1(icount,3)-3,posNy)+1

                        if((emat(dx,dy).eq.0)) then
                            rateconst = 0.384D0
                            eventOK=.true.
                            call ematandhbondsB(config1(icount,2),config1(icount,3), &
                                config1(icount,1),posNx,posNy,deltarconst)!,ematnew)
                            rateconst = rateconst - deltarconst

                        else
                            ! print*, 'step not allowed emat not empty (position 3)'
                            return
                        endif


                    else
                        !print*, 'step not allowed'
                        return
                    endif


                elseif(deltax.eq.-2)  then !deltax checking
                    if(deltay.eq.-1)  then !deltay checking

                        dx = modulo(config1(icount,2)-2,posNx)+1
                        dy = modulo(config1(icount,3)-1,posNy)+1
                        dx1 = modulo(config1(icount,2)-3,posNx)+1
                        dy1 = modulo(config1(icount,3)-1,posNy)+1
                        dx2 = modulo(config1(icount,2)-3,posNx)+1
                        dy2 = modulo(config1(icount,3)-2,posNy)+1
                        if((emat(dx,dy).eq.0).and.(emat(dx1,dy1).eq.0).and.(emat(dx2,dy2).eq.0)) then

                            rateconst = 0.582D0
                            eventOK=.true.
                            call ematandhbondsG(config1(icount,2),config1(icount,3), &
                                config1(icount,1),posNx,posNy,deltarconst)!,ematnew)
                            rateconst = rateconst - deltarconst


                        else
                            ! print*, 'step not allowed emat not empty (position 5 and 6)'
                            return
                        endif
                    else
                        ! print*, 'step not allowed'
                        return
                    endif
                else
                    ! print*, 'step not allowed'
                    return
                endif

            endif

        !*********************

        elseif(config1(icount,1).eq.2) then ! RED moves
            if((config2(icount,1).eq.1).AND.(config2(icount,4).eq.0)) then !from red to green
                if(deltay.eq.0) then !deltay checking
                    if((deltax.eq.0).AND.(config2(icount,5).ne.0))   then !deltax checking
                        !Periodic boundaries
                        dx = modulo(config1(icount,2)-2,posNx)+1
                        dy = modulo(config1(icount,3)-2,posNy)+1
                        if(emat(dx,dy).eq.0) then
                            if(config2(icount,5).eq.1) then ! Os translation
                                rateconst = 0.372D0
                                eventOK=.true.
                                call ematandhbondsC(config1(icount,2),config1(icount,3), &
                                    config1(icount,1),posNx,posNy,deltarconst)!,ematnew)
                                rateconst = rateconst - deltarconst
                            elseif(config2(icount,5).eq.2) then !Os flipping
                                rateconst = 0.537D0
                                eventOK=.true.
                                call ematandhbondsC(config1(icount,2),config1(icount,3), &
                                    config1(icount,1),posNx,posNy,deltarconst)!,ematnew)
                                rateconst = rateconst - deltarconst
                            elseif(config2(icount,5).eq.3) then !Os back-flipping
                                rateconst = 0.724D0
                                eventOK=.true.
                                call ematandhbondsC(config1(icount,2),config1(icount,3), &
                                    config1(icount,1),posNx,posNy,deltarconst)!,ematnew)
                                rateconst = rateconst - deltarconst
                            endif
                        else
                            !    print*, 'step not allowed emat not empty (position 1)'
                            return
                        endif

                    elseif(deltax.eq.1 .AND.config2(icount,5).eq.0)  then !deltax checking

                        dx = modulo(config1(icount,2),posNx)+1
                        dy = modulo(config1(icount,3)-1,posNy)+1
                        if(emat(dx,dy).eq.0) then
                            rateconst = 0.358D0
                            eventOK=.true.

                            call ematandhbondsD(config1(icount,2),config1(icount,3), &
                                config1(icount,1),posNx,posNy,deltarconst)!,ematnew)

                            rateconst = rateconst - deltarconst


                        else
                            !   print*, 'step not allowed emat not empty (position 4)'
                            return
                        endif

                    elseif(deltax.eq.2 .AND.config2(icount,5).eq.0)  then !deltax checking

                        dx = modulo(config1(icount,2)+1,posNx)+1
                        dy = modulo(config1(icount,3)-1,posNy)+1
                        dx1 = modulo(config1(icount,2)+1,posNx)+1
                        dy1 = modulo(config1(icount,3)-2,posNy)+1

                        if((emat(dx,dy).eq.0).and.(emat(dx1,dy1).eq.0)) then
                            rateconst = 0.582D0
                            eventOK=.true.
                            call ematandhbondsF(config1(icount,2),config1(icount,3), &
                                config1(icount,1),posNx,posNy,deltarconst)!,ematnew)
                            rateconst = rateconst - deltarconst


                        else
                            !  print*, 'step not allowed emat not empty (position 5 and 6)'
                            return
                        endif

                    else
                        ! print*, 'step not allowed'
                        return
                    endif
                else
                    !print*, 'step not allowed'
                    return
                endif
            elseif(config2(icount,1).eq.2) then !from red to red
                !print*, 'step not allowed'
                return
            elseif((config2(icount,1).eq.3) .AND. (config2(icount,4).eq.0).AND.(config2(icount,5).eq.0)) then
                !from red to yellow
                if(deltax.eq.1)  then !deltax checking
                    if(deltay.eq.-2)  then !deltay checking

                        dx = modulo(config1(icount,2),posNx)+1
                        dy = modulo(config1(icount,3)-3,posNy)+1

                        if((emat(dx,dy).eq.0)) then
                            rateconst = 0.384D0
                            eventOK=.true.
                            call ematandhbondsB(config1(icount,2),config1(icount,3), &
                                config1(icount,1),posNx,posNy,deltarconst)!,ematnew)
                            rateconst = rateconst - deltarconst

                        else
                            ! print*, 'step not allowed emat not empty (position 3)'
                            return
                        endif
                    else
                        !print*, 'step not allowed'
                        return
                    endif

                elseif(deltax.eq.2)  then !deltax checking
                    if(deltay.eq.-1)  then !deltay checking

                        dx = modulo(config1(icount,2),posNx)+1
                        dy = modulo(config1(icount,3)-1,posNy)+1
                        dx1 = modulo(config1(icount,2)+1,posNx)+1
                        dy1 = modulo(config1(icount,3)-1,posNy)+1
                        dx2 = modulo(config1(icount,2)+1,posNx)+1
                        dy2 = modulo(config1(icount,3)-2,posNy)+1
                        if((emat(dx,dy).eq.0).and.(emat(dx1,dy1).eq.0).and.(emat(dx2,dy2).eq.0)) then

                            rateconst = 0.582D0
                            eventOK=.true.

                            call ematandhbondsG(config1(icount,2),config1(icount,3), &
                                config1(icount,1),posNx,posNy,deltarconst)!,ematnew)
                            rateconst = rateconst - deltarconst

                        else
                            !   print*, 'step not allowed emat not empty (position 5 and 6)'
                            return
                        endif
                    else
                        !print*, 'step not allowed'
                        return
                    endif

                else
                    !print*, 'step not allowed'
                    return
                endif

            elseif(config2(icount,1).eq.4 .AND.config2(icount,5).eq.0) then !from red to blue

                if(deltax.eq.0)  then !deltax checking
                    if(deltay.eq.-2 .AND. config2(icount,4).eq.0)  then !deltay checking

                        dx = modulo(config1(icount,2)-1,posNx)+1
                        dy = modulo(config1(icount,3)-3,posNy)+1

                        if((emat(dx,dy).eq.0)) then
                            rateconst = 0.384D0
                            eventOK=.true.

                            call ematandhbondsA(config1(icount,2),config1(icount,3), &
                                config1(icount,1),posNx,posNy,deltarconst)!,ematnew)
                            rateconst = rateconst - deltarconst

                        else
                            !   print*, 'step not allowed emat not empty (position 2)'
                            return
                        endif

                    elseif(deltay.eq.-1) then

                        dx = modulo(config1(icount,2),posNx)+1
                        dy = modulo(config1(icount,3)-1,posNy)+1
                        if(emat(dx,dy).eq.0) then
                            if(config2(icount,4).eq.1) then !rotation plus flipping
                                rateconst = 0.582D0
                                eventOK=.true.

                                call ematandhbondsE(config1(icount,2),config1(icount,3), &
                                    config1(icount,1),posNx,posNy,deltarconst)!,ematnew)
                                rateconst = rateconst - deltarconst
                            elseif(config2(icount,4).eq.2) then !pure rotation
                                rateconst = 0.610D0
                                eventOK=.true.

                                call ematandhbondsE(config1(icount,2),config1(icount,3), &
                                    config1(icount,1),posNx,posNy,deltarconst)!,ematnew)
                                rateconst = rateconst - deltarconst
                            else
                                !   print*, 'step not allowed'
                                return
                            endif
                        else
                            ! print*, 'step not allowed emat not empty (position 4)'
                            return
                        endif
                    else
                        ! print*, 'step not allowed'
                        return
                    endif
                else
                    ! print*, 'step not allowed'
                    return
                endif
            endif


        !******************
        elseif(config1(icount,1).eq.3) then! YELLOW moves
            if(config2(icount,1).eq.1.AND.config2(icount,5).eq.0) then!from yellow to green
                if(deltax.eq.0)  then !deltax checking
                    if(deltay.eq.2.AND. config2(icount,4).eq.0)  then !deltay checking

                        dx = modulo(config1(icount,2)-1,posNx)+1
                        dy = modulo(config1(icount,3)+1,posNy)+1

                        if((emat(dx,dy).eq.0)) then
                            rateconst = 0.384D0
                            eventOK=.true.

                            call ematandhbondsA(config1(icount,2),config1(icount,3), &
                                config1(icount,1),posNx,posNy,deltarconst)!,ematnew)
                            rateconst = rateconst - deltarconst


                        else
                            ! print*, 'step not allowed emat not empty (position 2)'
                            return
                        endif

                    elseif(deltay.eq.1) then
                        dx = modulo(config1(icount,2)-2,posNx)+1
                        dy = modulo(config1(icount,3)-1,posNy)+1
                        if(emat(dx,dy).eq.0) then
                            if(config2(icount,4).eq.1) then !rotation plus flipping
                                rateconst = 0.582D0
                                eventOK=.true.

                                call ematandhbondsE(config1(icount,2),config1(icount,3), &
                                    config1(icount,1),posNx,posNy,deltarconst)!,ematnew)
                                rateconst = rateconst - deltarconst
                            elseif(config2(icount,4).eq.2) then !pure rotation
                                rateconst = 0.610D0
                                eventOK=.true.

                                call ematandhbondsE(config1(icount,2),config1(icount,3), &
                                    config1(icount,1),posNx,posNy,deltarconst)!,ematnew)
                                rateconst = rateconst - deltarconst
                                !else
                                !   print*, 'step not allowed'
                                return
                            endif
                        else
                            !print*, 'step not allowed emat not empty (position 4)'
                            return
                        endif

                    else
                        !      print*, 'step not allowed'
                        return
                    endif
                else
                    !  print*, 'step not allowed'
                    return
                endif

            elseif((config2(icount,1).eq.2).AND.(config2(icount,4).eq.0).AND.(config2(icount,5).eq.0)) then! from yellow to red

                if(deltax.eq.-1)  then !deltax checking
                    if(deltay.eq.2)  then !deltay checking

                        dx = modulo(config1(icount,2)-2,posNx)+1
                        dy = modulo(config1(icount,3)+1,posNy)+1

                        if((emat(dx,dy).eq.0)) then
                            rateconst = 0.384D0
                            eventOK=.true.

                            call ematandhbondsB(config1(icount,2),config1(icount,3), &
                                config1(icount,1),posNx,posNy,deltarconst)!,ematnew)
                            rateconst = rateconst - deltarconst

                        else
                            !   print*, 'step not allowed emat not empty (position 3)'
                            return
                        endif

                    else
                        ! print*, 'step not allowed'
                        return
                    endif
                elseif(deltax.eq.-2)  then !deltax checking
                    if(deltay.eq.1)  then !deltay checking

                        dx = modulo(config1(icount,2)-2,posNx)+1
                        dy = modulo(config1(icount,3)-1,posNy)+1
                        dx1 = modulo(config1(icount,2)-3,posNx)+1
                        dy1 = modulo(config1(icount,3)-1,posNy)+1
                        dx2 = modulo(config1(icount,2)-3,posNx)+1
                        dy2 = modulo(config1(icount,3),posNy)+1
                        if((emat(dx,dy).eq.0).and.(emat(dx1,dy1).eq.0).and.(emat(dx2,dy2).eq.0)) then

                            rateconst = 0.582D0
                            eventOK=.true.

                            call ematandhbondsG(config1(icount,2),config1(icount,3), &
                                config1(icount,1),posNx,posNy,deltarconst)!,ematnew)
                            rateconst = rateconst - deltarconst

                        else
                            !    print*, 'step not allowed emat not empty (position 5 and 6)'
                            return
                        endif
                    else
                        ! print*, 'step not allowed'
                        return
                    endif
                else
                    ! print*, 'step not allowed'
                    return
                endif

            elseif(config2(icount,1).eq.3) then!from yellow to yellow
                ! print*, 'step not allowed'
                return
            elseif((config2(icount,1).eq.4).AND.(config2(icount,4).eq.0)) then!from yellow to blue
                if(deltay.eq.0) then !deltay checking
                    if((deltax.eq.0).AND.(config2(icount,5).ne.0))  then !deltax checking
                        !Periodic boundaries
                        dx = modulo(config1(icount,2),posNx)+1
                        dy = modulo(config1(icount,3),posNy)+1
                        if(emat(dx,dy).eq.0) then
                            if(config2(icount,5).eq.1) then ! Os translation
                                rateconst = 0.372D0
                                eventOK=.true.

                                call ematandhbondsC(config1(icount,2),config1(icount,3), &
                                    config1(icount,1),posNx,posNy,deltarconst)!,ematnew)
                                rateconst = rateconst - deltarconst
                            elseif(config2(icount,5).eq.2) then !Os flipping
                                rateconst = 0.537D0
                                eventOK=.true.

                                call ematandhbondsC(config1(icount,2),config1(icount,3), &
                                    config1(icount,1),posNx,posNy,deltarconst)!,ematnew)
                                rateconst = rateconst - deltarconst
                            elseif(config2(icount,5).eq.3) then !Os back-flipping
                                rateconst = 0.724D0
                                eventOK=.true.

                                call ematandhbondsC(config1(icount,2),config1(icount,3), &
                                    config1(icount,1),posNx,posNy,deltarconst)!,ematnew)
                                rateconst = rateconst - deltarconst
                            endif
                        else
                            !    print*, 'step not allowed emat not empty (position 1)'
                            return
                        endif
                    elseif(deltax.eq.-1 .AND.config2(icount,5).eq.0)  then !deltax checking
                        dx = modulo(config1(icount,2)-2,posNx)+1
                        dy = modulo(config1(icount,3)-1,posNy)+1
                        if(emat(dx,dy).eq.0) then
                            rateconst = 0.358D0
                            eventOK=.true.

                            call ematandhbondsD(config1(icount,2),config1(icount,3), &
                                config1(icount,1),posNx,posNy,deltarconst)!,ematnew)

                            rateconst = rateconst - deltarconst
                        else
                            !   print*, 'step not allowed emat not empty (position 4)'
                            return
                        endif


                    elseif(deltax.eq.-2 .AND.config2(icount,5).eq.0)  then !deltax checking

                        dx = modulo(config1(icount,2)-3,posNx)+1
                        dy = modulo(config1(icount,3)-1,posNy)+1
                        dx1 = modulo(config1(icount,2)-3,posNx)+1
                        dy1 = modulo(config1(icount,3),posNy)+1

                        if((emat(dx,dy).eq.0).and.(emat(dx1,dy1).eq.0)) then
                            rateconst = 0.582D0
                            eventOK=.true.

                            call ematandhbondsF(config1(icount,2),config1(icount,3), &
                                config1(icount,1),posNx,posNy,deltarconst)!,ematnew)
                            rateconst = rateconst - deltarconst
                        else
                            !  print*, 'step not allowed emat not empty (position 5 and 6)'
                            return
                        endif

                    else
                        ! print*, 'step not allowed'
                        return
                    endif
                else
                    ! print*, 'step not allowed'
                    return
                endif
            endif

            !*****************************

        elseif(config1(icount,1).eq.4) then! BLUE moves
            if((config2(icount,1).eq.1) .AND.(config2(icount,4).eq.0).AND.(config2(icount,5).eq.0)) then!from blue to green
                if(deltax.eq.1)  then !deltax checking
                    if(deltay.eq.2)  then !deltay checking
                        dx = modulo(config1(icount,2),posNx)+1
                        dy = modulo(config1(icount,3)+1,posNy)+1
                        if((emat(dx,dy).eq.0)) then
                            rateconst = 0.384D0
                            eventOK=.true.

                            call ematandhbondsB(config1(icount,2),config1(icount,3), &
                                config1(icount,1),posNx,posNy,deltarconst)!,ematnew)
                            rateconst = rateconst - deltarconst
                        else
                            !        print*, 'step not allowed emat not empty (position 3)'
                            return
                        endif


                    else
                        ! print*, 'step not allowed'
                        return
                    endif
                elseif(deltax.eq.2)  then !deltax checking
                    if(deltay.eq.1)  then !deltay checking
                        dx = modulo(config1(icount,2),posNx)+1
                        dy = modulo(config1(icount,3)-1,posNy)+1
                        dx1 = modulo(config1(icount,2)+1,posNx)+1
                        dy1 = modulo(config1(icount,3)-1,posNy)+1
                        dx2 = modulo(config1(icount,2)+1,posNx)+1
                        dy2 = modulo(config1(icount,3),posNy)+1

                        if((emat(dx,dy).eq.0).and.(emat(dx1,dy1).eq.0).and.(emat(dx2,dy2).eq.0)) then

                            rateconst = 0.582D0
                            eventOK=.true.

                            call ematandhbondsG(config1(icount,2),config1(icount,3), &
                                config1(icount,1),posNx,posNy,deltarconst)!,ematnew)
                            rateconst = rateconst - deltarconst

                        else
                            !        print*, 'step not allowed emat not empty (position 5 and 6)'
                            return
                        endif
                    else
                        !  print*, 'step not allowed'
                        return
                    endif
                else
                    ! print*, 'step not allowed'
                    return
                endif
            elseif(config2(icount,1).eq.2.AND.config2(icount,5).eq.0) then !blue to red
                if(deltax.eq.0)  then !deltax checking
                    if(deltay.eq.2.AND. config2(icount,4).eq.0)  then !deltay checking
                        dx = modulo(config1(icount,2)-1,posNx)+1
                        dy = modulo(config1(icount,3)+1,posNy)+1
                        if((emat(dx,dy).eq.0)) then
                            rateconst = 0.384D0
                            eventOK=.true.

                            call ematandhbondsA(config1(icount,2),config1(icount,3), &
                                config1(icount,1),posNx,posNy,deltarconst)!,ematnew)
                            rateconst = rateconst - deltarconst

                        else
                            !   print*, 'step not allowed emat not empty (position 2)'
                            return
                        endif
                    elseif(deltay.eq.1) then
                        dx = modulo(config1(icount,2),posNx)+1
                        dy = modulo(config1(icount,3)-1,posNy)+1

                        if(emat(dx,dy).eq.0) then
                            if(config2(icount,4).eq.1) then !rotation plus flipping
                                rateconst = 0.582D0
                                eventOK=.true.

                                call ematandhbondsE(config1(icount,2),config1(icount,3), &
                                    config1(icount,1),posNx,posNy,deltarconst)!,ematnew)
                                rateconst = rateconst - deltarconst

                            elseif(config2(icount,4).eq.2) then !pure rotation
                                rateconst = 0.610D0
                                eventOK=.true.

                                call ematandhbondsE(config1(icount,2),config1(icount,3), &
                                    config1(icount,1),posNx,posNy,deltarconst)!,ematnew)
                                rateconst = rateconst - deltarconst

                            else
                                !  print*, 'step not allowed'
                                return
                            endif

                        else
                            !  print*, 'step not allowed emat not empty (position 4)'
                            return
                        endif
                    else
                        !  print*, 'step not allowed'
                        return
                    endif
                else
                    ! print*, 'step not allowed'
                    return
                endif
            elseif((config2(icount,1).eq.3).AND.(config2(icount,4).eq.0)) then!from blue to yellow
                if(deltay.eq.0) then !deltay checking
                    if((deltax.eq.0).AND.(config2(icount,5).ne.0))  then !deltax checking
                        !Periodic boundaries
                        dx = modulo(config1(icount,2)-2,posNx)+1
                        dy = modulo(config1(icount,3),posNy)+1

                        if(emat(dx,dy).eq.0) then
                            if(config2(icount,5).eq.1) then ! Os translation
                                rateconst = 0.372D0
                                eventOK=.true.

                                call ematandhbondsC(config1(icount,2),config1(icount,3), &
                                    config1(icount,1),posNx,posNy,deltarconst)!,ematnew)
                                rateconst = rateconst - deltarconst
                            elseif(config2(icount,5).eq.2) then !Os flipping
                                rateconst = 0.537D0
                                eventOK=.true.

                                call ematandhbondsC(config1(icount,2),config1(icount,3), &
                                    config1(icount,1),posNx,posNy,deltarconst)!,ematnew)
                                rateconst = rateconst - deltarconst
                            elseif(config2(icount,5).eq.3) then !Os back-flipping
                                rateconst = 0.724D0
                                eventOK=.true.

                                call ematandhbondsC(config1(icount,2),config1(icount,3), &
                                    config1(icount,1),posNx,posNy,deltarconst)!,ematnew)
                                rateconst = rateconst - deltarconst
                            endif

                        else
                            !  print*, 'step not allowed emat not empty (position 1)'
                            return
                        endif


                    elseif(deltax.eq.1.AND.config2(icount,5).eq.0)  then !deltax checking
                        dx = modulo(config1(icount,2),posNx)+1
                        dy = modulo(config1(icount,3)-1,posNy)+1

                        if(emat(dx,dy).eq.0) then
                            rateconst = 0.358D0
                            eventOK=.true.

                            call ematandhbondsD(config1(icount,2),config1(icount,3), &
                                config1(icount,1),posNx,posNy,deltarconst)!,ematnew)

                            rateconst = rateconst - deltarconst

                        else
                            !  print*, 'step not allowed emat not empty (position 4)'
                            return
                        endif

                    elseif(deltax.eq.2.AND.config2(icount,5).eq.0)  then !deltax checking
                        dx = modulo(config1(icount,2)+1,posNx)+1
                        dy = modulo(config1(icount,3)-1,posNy)+1
                        dx1 = modulo(config1(icount,2)+1,posNx)+1
                        dy1 = modulo(config1(icount,3),posNy)+1

                        if((emat(dx,dy).eq.0).and.(emat(dx1,dy1).eq.0)) then
                            rateconst = 0.582D0
                            eventOK=.true.

                            call ematandhbondsF(config1(icount,2),config1(icount,3), &
                                config1(icount,1),posNx,posNy,deltarconst)!,ematnew)
                            rateconst = rateconst - deltarconst

                        else
                            !  print*, 'step not allowed emat not empty (position 5 and 6)'
                            return
                        endif

                    else
                        !  print*, 'step not allowed'
                        return
                    endif
                else
                    ! print*, 'step not allowed'
                    return
                endif



            elseif(config2(icount,1).eq.4) then!from blue to blue
                !  print*, 'step not allowed'
                return
            endif

        else
            !  print*, 'step not allowed config1 has a wrong color'
            return
        endif

    enddo


end subroutine glyevents

          subroutine init_random_seed()
            use iso_fortran_env, only: int64
            implicit none
            integer, allocatable :: seed(:)
            integer :: n, un, istat

            !fixed unit (dangerous but works for now)
            un=55

            call random_seed(size = n)
            allocate(seed(n))
            ! First try if the OS provides a random number generator
            open(unit=un, file='/dev/urandom', access='stream', &
            form='unformatted', action='read', status='old', iostat=istat)
            if (istat == 0) then
               read(un) seed
               close(un)
            else
            print*,'file=/dev/urandom not readable'
            end if
            call random_seed(put=seed)

          contains
            ! This simple PRNG might not be good enough for real work, but is
            ! sufficient for seeding a better PRNG.
            function lcg(s)
              integer :: lcg
              integer(int64) :: s
              if (s == 0) then
                 s = 104729
              else
                 s = mod(s, 4294967296_int64)
              end if
              s = mod(s * 279470273_int64, 4294967291_int64)
              lcg = int(mod(s, int(huge(0), int64)), kind(0))
            end function lcg
          end subroutine init_random_seed


subroutine checkmatgen(config,checkmatg)
use commons, only:posNx,posNy
    implicit none
    integer, intent(in) :: config(5)
    integer, intent(out) :: checkmatg(2,3)

checkmatg(:,:)=0

if (config(1).eq.1) then !green

 checkmatg(1,1)=modulo(config(2)-1,posNx)+1 !x Nickel
 checkmatg(2,1)=modulo(config(3)-1,posNy)+1 !y Nickel
 checkmatg(1,2)=modulo(config(2)-1,posNx)+1 !x Oxygen
 checkmatg(2,2)=modulo(config(3)-2,posNy)+1 !y Oxygen
 checkmatg(1,3)=modulo(config(2)-2,posNx)+1 !x Oxygen
 checkmatg(2,3)=modulo(config(3)-2,posNy)+1 !y Oxygen

elseif (config(1).eq.2) then !red

 checkmatg(1,1)=modulo(config(2)-1,posNx)+1 !x Nickel
 checkmatg(2,1)=modulo(config(3)-1,posNy)+1 !y Nickel
 checkmatg(1,2)=modulo(config(2)-1,posNx)+1 !x Oxygen
 checkmatg(2,2)=modulo(config(3)-2,posNy)+1 !y Oxygen
 checkmatg(1,3)=modulo(config(2),posNx)+1 !x Oxygen
 checkmatg(2,3)=modulo(config(3)-2,posNy)+1 !y Oxygen

elseif (config(1).eq.3) then !yellow

 checkmatg(1,1)=modulo(config(2)-1,posNx)+1 !x Nickel
 checkmatg(2,1)=modulo(config(3)-1,posNy)+1 !y Nickel
 checkmatg(1,2)=modulo(config(2)-1,posNx)+1 !x Oxygen
 checkmatg(2,2)=modulo(config(3),posNy)+1 !y Oxygen
 checkmatg(1,3)=modulo(config(2)-2,posNx)+1 !x Oxygen
 checkmatg(2,3)=modulo(config(3),posNy)+1 !y Oxygen

elseif (config(1).eq.4) then !blue

 checkmatg(1,1)=modulo(config(2)-1,posNx)+1 !x Nickel
 checkmatg(2,1)=modulo(config(3)-1,posNy)+1 !y Nickel
 checkmatg(1,2)=modulo(config(2)-1,posNx)+1 !x Oxygen
 checkmatg(2,2)=modulo(config(3),posNy)+1 !y Oxygen
 checkmatg(1,3)=modulo(config(2),posNx)+1 !x Oxygen
 checkmatg(2,3)=modulo(config(3),posNy)+1 !y Oxygen

endif

end subroutine checkmatgen


subroutine ematandhbondsA(xPos,yPos,color,posNx,posNy, &
    deltarconst)!,ematnew)
    use commons, only:emat,hbondconst,indexmat,neibrepel
    implicit none
    integer, intent(in) :: xPos,yPos,posNx,posNy,color
    double precision, intent(out) :: deltarconst
!    integer, intent(out) :: ematnew(posNx,posNy)
    integer :: dx1,dy1,dy2,dx3,dy3,nbind1,nbind2

  nbind1=0
  nbind2=0

    if (color.eq.1) then
        !updating for green to yellow (case A)
        dx1 = modulo(xPos-1,posNx)+1
        dy1 = modulo(yPos-1,posNy)+1
        dy2 = modulo(yPos-3,posNy)+1

!        ematnew(dx1,dy1)=0
!        ematnew(dx1,dy2)=1

        !hbond check (2 broken 2 stay 2 build)
        dx3 = modulo(xPos-3,posNx)+1
        dy3 = modulo(yPos-2,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen   !already there
            deltarconst = deltarconst + hbondconst
            nbind1=indexmat(dx3,dy3)
        endif
        dx3 = modulo(xPos,posNx)+1
        dy3 = modulo(yPos-2,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen  !already there
            deltarconst = deltarconst + hbondconst
            nbind2=indexmat(dx3,dy3)
        endif
        !build
        dx3 = modulo(xPos,posNx)+1
        dy3 = modulo(yPos-3,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst + hbondconst
            if((indexmat(dx3,dy3).eq.nbind1).or.(indexmat(dx3,dy3).eq.nbind2)) then
            deltarconst = deltarconst - neibrepel
            endif
        endif
        dx3 = modulo(xPos-2,posNx)+1
        dy3 = modulo(yPos-3,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst + hbondconst
            if((indexmat(dx3,dy3).eq.nbind1).or.(indexmat(dx3,dy3).eq.nbind2)) then
            deltarconst = deltarconst - neibrepel
            endif
        endif

        !broken
        dx3 = modulo(xPos,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst - hbondconst
        endif
        dx3 = modulo(xPos-2,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst - hbondconst
        endif

    elseif (color.eq.2) then
        !updating for red to blue (case A)
        dx1 = modulo(xPos-1,posNx)+1
        dy1 = modulo(yPos-1,posNy)+1
        dy2 = modulo(yPos-3,posNy)+1

!        ematnew(dx1,dy1)=0
!        ematnew(dx1,dy2)=1

        !hbond check (2 broken 2 stay 2 build)
        dx3 = modulo(xPos-2,posNx)+1
        dy3 = modulo(yPos-2,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen   !already there
            deltarconst = deltarconst + hbondconst
            nbind1=indexmat(dx3,dy3)
        endif
        dx3 = modulo(xPos+1,posNx)+1
        dy3 = modulo(yPos-2,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen  !already there
            deltarconst = deltarconst + hbondconst
            nbind2=indexmat(dx3,dy3)
        endif
        !build

        dx3 = modulo(xPos,posNx)+1
        dy3 = modulo(yPos-3,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst + hbondconst
            if((indexmat(dx3,dy3).eq.nbind1).or.(indexmat(dx3,dy3).eq.nbind2)) then
            deltarconst = deltarconst - neibrepel
            endif
        endif
        dx3 = modulo(xPos-2,posNx)+1
        dy3 = modulo(yPos-3,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst + hbondconst
            if((indexmat(dx3,dy3).eq.nbind1).or.(indexmat(dx3,dy3).eq.nbind2)) then
            deltarconst = deltarconst - neibrepel
            endif
        endif
        !broken
        dx3 = modulo(xPos,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst - hbondconst
        endif
        dx3 = modulo(xPos-2,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst - hbondconst
        endif
    elseif (color.eq.3) then
        !updating for yellow to green (case A)
        dx1 = modulo(xPos-1,posNx)+1
        dy1 = modulo(yPos-1,posNy)+1
        dy2 = modulo(yPos+1,posNy)+1

!        ematnew(dx1,dy1)=0
!        ematnew(dx1,dy2)=1

        !hbond check (2 broken 2 stay 2 build)
        dx3 = modulo(xPos-3,posNx)+1
        dy3 = modulo(yPos,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen   !already there
            deltarconst = deltarconst + hbondconst
        nbind1=indexmat(dx3,dy3)
        endif
        dx3 = modulo(xPos,posNx)+1
        dy3 = modulo(yPos,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen  !already there
            deltarconst = deltarconst + hbondconst
        nbind2=indexmat(dx3,dy3)
        endif

        !build
        dx3 = modulo(xPos,posNx)+1
        dy3 = modulo(yPos+1,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst + hbondconst
            if((indexmat(dx3,dy3).eq.nbind1).or.(indexmat(dx3,dy3).eq.nbind2)) then
            deltarconst = deltarconst - neibrepel
            endif
        endif
        dx3 = modulo(xPos-2,posNx)+1
        dy3 = modulo(yPos+1,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst + hbondconst
            if((indexmat(dx3,dy3).eq.nbind1).or.(indexmat(dx3,dy3).eq.nbind2)) then
            deltarconst = deltarconst - neibrepel
            endif
        endif
        !broken
        dx3 = modulo(xPos,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst - hbondconst
        endif
        dx3 = modulo(xPos-2,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst - hbondconst
        endif
    elseif (color.eq.4) then
        !updating for blue to red (case A)
        dx1 = modulo(xPos-1,posNx)+1
        dy1 = modulo(yPos-1,posNy)+1
        dy2 = modulo(yPos+1,posNy)+1

!        ematnew(dx1,dy1)=0
!        ematnew(dx1,dy2)=1

        !hbond check (2 broken 2 stay 2 build)
        dx3 = modulo(xPos-2,posNx)+1
        dy3 = modulo(yPos,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen   !already there
            deltarconst = deltarconst + hbondconst
            nbind1=indexmat(dx3,dy3)
        endif
        dx3 = modulo(xPos+1,posNx)+1
        dy3 = modulo(yPos,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen  !already there
            deltarconst = deltarconst + hbondconst
            nbind2=indexmat(dx3,dy3)
        endif
        !build
        dx3 = modulo(xPos,posNx)+1
        dy3 = modulo(yPos+1,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst + hbondconst
           if((indexmat(dx3,dy3).eq.nbind1).or.(indexmat(dx3,dy3).eq.nbind2)) then
            deltarconst = deltarconst - neibrepel
            endif
        endif
        dx3 = modulo(xPos-2,posNx)+1
        dy3 = modulo(yPos+1,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst + hbondconst
            if((indexmat(dx3,dy3).eq.nbind1).or.(indexmat(dx3,dy3).eq.nbind2)) then
            deltarconst = deltarconst - neibrepel
            endif
        endif
        !broken
        dx3 = modulo(xPos,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst - hbondconst
        endif
        dx3 = modulo(xPos-2,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst - hbondconst
        endif

    else
        print*,'wrong color'
        stop
    endif

end subroutine ematandhbondsA


subroutine ematandhbondsB(xPos,yPos,color,posNx,posNy, &
    deltarconst)!,ematnew)
    use commons, only:emat,hbondconst,indexmat,neibrepel
    implicit none
    integer, intent(in) :: xPos,yPos,posNx,posNy,color
    double precision, intent(out) :: deltarconst
!    integer, intent(out) :: ematnew(posNx,posNy)
    integer :: dx1,dy1,dx2,dy2,dx3,dy3,nbind1,nbind2

    nbind1=0
    nbind2=0

    if (color.eq.1) then
        !updating for green to blue (case B)
        dx1 = modulo(xPos-1,posNx)+1
        dy1 = modulo(yPos-1,posNy)+1
        dx2 = modulo(xPos-2,posNx)+1
        dy2 = modulo(yPos-3,posNy)+1

!        ematnew(dx1,dy1)=0
!        ematnew(dx2,dy2)=1

        !hbond check (4 broken 4 build, 2 flip "stay")
        dx3 = modulo(xPos-3,posNx)+1
        dy3 = modulo(yPos-2,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen   !already there
            deltarconst = deltarconst + hbondconst
            nbind1=indexmat(dx3,dy3)
        endif
        dx3 = modulo(xPos,posNx)+1
        dy3 = modulo(yPos-2,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen  !already there
            deltarconst = deltarconst + hbondconst
            nbind2=indexmat(dx3,dy3)
        endif

        !build
        dx3 = modulo(xPos-3,posNx)+1
        dy3 = modulo(yPos-3,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst + hbondconst
            if((indexmat(dx3,dy3).eq.nbind1).or.(indexmat(dx3,dy3).eq.nbind2)) then
            deltarconst = deltarconst - neibrepel
            endif
        endif
        dx3 = modulo(xPos-1,posNx)+1
        dy3 = modulo(yPos-3,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst + hbondconst
            if((indexmat(dx3,dy3).eq.nbind1).or.(indexmat(dx3,dy3).eq.nbind2)) then
            deltarconst = deltarconst - neibrepel
            endif
        endif
        !broken
        dx3 = modulo(xPos-2,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst - hbondconst
        endif
        dx3 = modulo(xPos,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst - hbondconst
        endif

    elseif (color.eq.2) then
        !updating for red to yellow (case B)
        dx1 = modulo(xPos-1,posNx)+1
        dy1 = modulo(yPos-1,posNy)+1
        dx2 = modulo(xPos,posNx)+1
        dy2 = modulo(yPos-3,posNy)+1

!        ematnew(dx1,dy1)=0
!        ematnew(dx2,dy2)=1

        !hbond check (4 broken 4 build, 2 flip "stay")
        dx3 = modulo(xPos-2,posNx)+1
        dy3 = modulo(yPos-2,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen   !already there
            deltarconst = deltarconst + hbondconst
            nbind1=indexmat(dx3,dy3)
        endif
        dx3 = modulo(xPos+1,posNx)+1
        dy3 = modulo(yPos-2,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen  !already there
            deltarconst = deltarconst + hbondconst
            nbind2=indexmat(dx3,dy3)
        endif
        !build
        dx3 = modulo(xPos+1,posNx)+1
        dy3 = modulo(yPos-3,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst + hbondconst
            if((indexmat(dx3,dy3).eq.nbind1).or.(indexmat(dx3,dy3).eq.nbind2)) then
            deltarconst = deltarconst - neibrepel
            endif
        endif
        dx3 = modulo(xPos-1,posNx)+1
        dy3 = modulo(yPos-3,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst + hbondconst
            if((indexmat(dx3,dy3).eq.nbind1).or.(indexmat(dx3,dy3).eq.nbind2)) then
            deltarconst = deltarconst - neibrepel
            endif
        endif
        !broken
        dx3 = modulo(xPos-2,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst - hbondconst
        endif
        dx3 = modulo(xPos,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst - hbondconst
        endif

    elseif (color.eq.3) then
        !updating for yellow to red (case B)
        dx1 = modulo(xPos-1,posNx)+1
        dy1 = modulo(yPos-1,posNy)+1
        dx2 = modulo(xPos-2,posNx)+1
        dy2 = modulo(yPos+1,posNy)+1

!        ematnew(dx1,dy1)=0
!        ematnew(dx2,dy2)=1

        !hbond check (4 broken 4 build, 2 flip "stay")
        dx3 = modulo(xPos-3,posNx)+1
        dy3 = modulo(yPos,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen   !already there
            deltarconst = deltarconst + hbondconst
            nbind1=indexmat(dx3,dy3)
        endif
        dx3 = modulo(xPos,posNx)+1
        dy3 = modulo(yPos,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen  !already there
            deltarconst = deltarconst + hbondconst
            nbind2=indexmat(dx3,dy3)
        endif
        !build
        dx3 = modulo(xPos-3,posNx)+1
        dy3 = modulo(yPos+1,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst + hbondconst
            if((indexmat(dx3,dy3).eq.nbind1).or.(indexmat(dx3,dy3).eq.nbind2)) then
            deltarconst = deltarconst - neibrepel
            endif
        endif
        dx3 = modulo(xPos-1,posNx)+1
        dy3 = modulo(yPos+1,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst + hbondconst
            if((indexmat(dx3,dy3).eq.nbind1).or.(indexmat(dx3,dy3).eq.nbind2)) then
            deltarconst = deltarconst - neibrepel
            endif
        endif
        !broken
        dx3 = modulo(xPos-2,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst - hbondconst
        endif
        dx3 = modulo(xPos,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst - hbondconst
        endif
    elseif (color.eq.4) then
        !updating for blue to green (case B)
        dx1 = modulo(xPos-1,posNx)+1
        dy1 = modulo(yPos-1,posNy)+1
        dx2 = modulo(xPos,posNx)+1
        dy2 = modulo(yPos+1,posNy)+1

!        ematnew(dx1,dy1)=0
!        ematnew(dx2,dy2)=1

        !hbond check (4 broken 4 build, 2 flip "stay")
        dx3 = modulo(xPos-2,posNx)+1
        dy3 = modulo(yPos,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen   !already there
            deltarconst = deltarconst + hbondconst
            nbind1=indexmat(dx3,dy3)
        endif
        dx3 = modulo(xPos+1,posNx)+1
        dy3 = modulo(yPos,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen  !already there
            deltarconst = deltarconst + hbondconst
            nbind2=indexmat(dx3,dy3)
        endif
        !build
        dx3 = modulo(xPos+1,posNx)+1
        dy3 = modulo(yPos+1,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst + hbondconst
            if((indexmat(dx3,dy3).eq.nbind1).or.(indexmat(dx3,dy3).eq.nbind2)) then
            deltarconst = deltarconst - neibrepel
            endif
        endif
        dx3 = modulo(xPos-1,posNx)+1
        dy3 = modulo(yPos+1,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst + hbondconst
            if((indexmat(dx3,dy3).eq.nbind1).or.(indexmat(dx3,dy3).eq.nbind2)) then
            deltarconst = deltarconst - neibrepel
            endif
        endif
        !broken
        dx3 = modulo(xPos,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst - hbondconst
        endif
        dx3 = modulo(xPos-2,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst - hbondconst
        endif
    else
        print*,'wrong color'
        stop
    endif

end subroutine ematandhbondsB

subroutine ematandhbondsC(xPos,yPos,color,posNx,posNy, &
    deltarconst)!,ematnew)
    use commons, only:emat,hbondconst,indexmat,neibrepel
    implicit none
    integer, intent(in) :: xPos,yPos,posNx,posNy,color
    double precision, intent(out) :: deltarconst
!    integer, intent(out) :: ematnew(posNx,posNy)
    integer :: dx1,dy1,dx2,dy2,dx3,dy3,nbind1,nbind2

nbind1=0
nbind2=0

!    ematnew(:,:)=emat(:,:)
    if (color.eq.1) then
        !updating for green to red (case C)
        dx1 = modulo(xPos-2,posNx)+1
        dy1 = modulo(yPos-2,posNy)+1
        dx2 = modulo(xPos,posNx)+1
        dy2 = modulo(yPos-2,posNy)+1

!        ematnew(dx1,dy1)=0
!        ematnew(dx2,dy2)=2

        !hbond check (1 breaks 1 build 2 stay)
        !build
        dx3 = modulo(xPos,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen !already there
            deltarconst = deltarconst + hbondconst
            nbind1=indexmat(dx3,dy3)
        endif
        dx3 = modulo(xPos-2,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen !already there
            deltarconst = deltarconst + hbondconst
            nbind2=indexmat(dx3,dy3)
        endif
        dx3 = modulo(xPos+1,posNx)+1
        dy3 = modulo(yPos-2,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen
            deltarconst = deltarconst + hbondconst
            if((indexmat(dx3,dy3).eq.nbind1).or.(indexmat(dx3,dy3).eq.nbind2)) then
            deltarconst = deltarconst - neibrepel
            endif
        endif

        !broken
        dx3 = modulo(xPos-3,posNx)+1
        dy3 = modulo(yPos-2,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen
            deltarconst = deltarconst - hbondconst
        endif
    elseif (color.eq.2) then
        !updating for red to green (case C)
        dx1 = modulo(xPos,posNx)+1
        dy1 = modulo(yPos-2,posNy)+1
        dx2 = modulo(xPos-2,posNx)+1
        dy2 = modulo(yPos-2,posNy)+1

!        ematnew(dx1,dy1)=0
!        ematnew(dx2,dy2)=2

        !hbond check (1 breaks 1 build 2 stay)
        !build
        dx3 = modulo(xPos-2,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen !already there
            deltarconst = deltarconst + hbondconst
            nbind1=indexmat(dx3,dy3)
        endif
        dx3 = modulo(xPos,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen !already there
            deltarconst = deltarconst + hbondconst
            nbind2=indexmat(dx3,dy3)
        endif
        dx3 = modulo(xPos-3,posNx)+1
        dy3 = modulo(yPos-2,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen
            deltarconst = deltarconst + hbondconst
            if((indexmat(dx3,dy3).eq.nbind1).or.(indexmat(dx3,dy3).eq.nbind2)) then
            deltarconst = deltarconst - neibrepel
            endif
        endif

        !broken
        dx3 = modulo(xPos+1,posNx)+1
        dy3 = modulo(yPos-2,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen
            deltarconst = deltarconst - hbondconst
        endif
    elseif (color.eq.3) then
        !updating for yellow to blue (case C)
        dx1 = modulo(xPos-2,posNx)+1
        dy1 = modulo(yPos,posNy)+1
        dx2 = modulo(xPos,posNx)+1
        dy2 = modulo(yPos,posNy)+1

!        ematnew(dx1,dy1)=0
!        ematnew(dx2,dy2)=2

        !hbond check (1 breaks 1 build 2 stay)
        !build
        dx3 = modulo(xPos,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen !already there
            deltarconst = deltarconst + hbondconst
            nbind1=indexmat(dx3,dy3)
        endif
        dx3 = modulo(xPos-2,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen !already there
            deltarconst = deltarconst + hbondconst
            nbind2=indexmat(dx3,dy3)
        endif
        dx3 = modulo(xPos+1,posNx)+1
        dy3 = modulo(yPos,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen
            deltarconst = deltarconst + hbondconst
            if((indexmat(dx3,dy3).eq.nbind1).or.(indexmat(dx3,dy3).eq.nbind2)) then
            deltarconst = deltarconst - neibrepel
            endif
        endif
        !broken
        dx3 = modulo(xPos-3,posNx)+1
        dy3 = modulo(yPos,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen
            deltarconst = deltarconst - hbondconst
        endif
    elseif (color.eq.4) then
        !updating for blue to yellow (case C)
        dx1 = modulo(xPos,posNx)+1
        dy1 = modulo(yPos,posNy)+1
        dx2 = modulo(xPos-2,posNx)+1
        dy2 = modulo(yPos,posNy)+1

!        ematnew(dx1,dy1)=0
!        ematnew(dx2,dy2)=2

        !hbond check (1 breaks 1 build 2 stay)
        !build
        dx3 = modulo(xPos,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen !already there
            deltarconst = deltarconst + hbondconst
            nbind1=indexmat(dx3,dy3)
        endif
        dx3 = modulo(xPos-2,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen !already there
            deltarconst = deltarconst + hbondconst
            nbind2=indexmat(dx3,dy3)
        endif
        dx3 = modulo(xPos-3,posNx)+1
        dy3 = modulo(yPos,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen
            deltarconst = deltarconst + hbondconst
            if((indexmat(dx3,dy3).eq.nbind1).or.(indexmat(dx3,dy3).eq.nbind2)) then
            deltarconst = deltarconst - neibrepel
            endif
        endif
        !broken
        dx3 = modulo(xPos+1,posNx)+1
        dy3 = modulo(yPos,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen
            deltarconst = deltarconst - hbondconst
        endif
    else
        print*,'wrong color'
        stop
    endif

end subroutine ematandhbondsC



subroutine ematandhbondsD(xPos,yPos,color,posNx,posNy, &
    deltarconst)!,ematnew)
    use commons, only:emat,hbondconst,indexmat,neibrepel
    implicit none
    integer, intent(in) :: xPos,yPos,posNx,posNy,color
    double precision, intent(out) :: deltarconst
!    integer, intent(out) :: ematnew(posNx,posNy)
    integer :: dx1,dy1,dx2,dy2,dx3,dy3,nbind1,nbind2

nbind1=0
nbind2=0

    if (color.eq.1) then
        !updating for green to red (case D)
        dx1 = modulo(xPos-1,posNx)+1
        dy1 = modulo(yPos-1,posNy)+1
        dx2 = modulo(xPos-2,posNx)+1
        dy2 = modulo(yPos-1,posNy)+1

!        ematnew(dx1,dy1)=0
!        ematnew(dx2,dy2)=1

        !hbond check (3 break, 2 build, 2 flip "stay")
        !broken
        dx3 = modulo(xPos,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst - hbondconst
        endif
        dx3 = modulo(xPos,posNx)+1
        dy3 = modulo(yPos-2,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen  !already there
            deltarconst = deltarconst + hbondconst
            nbind2=indexmat(dx3,dy3)
        endif

        dx3 = modulo(xPos-3,posNx)+1
        dy3 = modulo(yPos-2,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen   !already there
            deltarconst = deltarconst + hbondconst
            nbind1=indexmat(dx3,dy3)
        endif
        !build
        dx3 = modulo(xPos-3,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst + hbondconst
            if((indexmat(dx3,dy3).eq.nbind1).or.(indexmat(dx3,dy3).eq.nbind2)) then
            deltarconst = deltarconst - neibrepel
            endif
        endif

    elseif (color.eq.2) then
        !updating for red to green (case D)
        dx1 = modulo(xPos-1,posNx)+1
        dy1 = modulo(yPos-1,posNy)+1
        dx2 = modulo(xPos,posNx)+1
        dy2 = modulo(yPos-1,posNy)+1

!        ematnew(dx1,dy1)=0
!        ematnew(dx2,dy2)=1

        !hbond check (3 break, 2 build, 2 flip "stay")
        !broken
        dx3 = modulo(xPos-2,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst - hbondconst
        endif
        dx3 = modulo(xPos+1,posNx)+1
        dy3 = modulo(yPos-2,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen  !already there
            deltarconst = deltarconst + hbondconst
            nbind2=indexmat(dx3,dy3)
        endif

        dx3 = modulo(xPos-2,posNx)+1
        dy3 = modulo(yPos-2,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen   !already there
            deltarconst = deltarconst + hbondconst
            nbind1=indexmat(dx3,dy3)
        endif
        !build
        dx3 = modulo(xPos+1,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst + hbondconst
            if((indexmat(dx3,dy3).eq.nbind1).or.(indexmat(dx3,dy3).eq.nbind2)) then
            deltarconst = deltarconst - neibrepel
            endif
        endif


    elseif (color.eq.3) then
        !updating for yellow to blue (case D)
        dx1 = modulo(xPos-1,posNx)+1
        dy1 = modulo(yPos-1,posNy)+1
        dx2 = modulo(xPos-2,posNx)+1
        dy2 = modulo(yPos-1,posNy)+1

!        ematnew(dx1,dy1)=0
!        ematnew(dx2,dy2)=1

        !hbond check (3 break, 2 build, 2 flip "stay")
        !broken
        dx3 = modulo(xPos,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst - hbondconst
        endif
        dx3 = modulo(xPos-3,posNx)+1
        dy3 = modulo(yPos,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen  !already there
            deltarconst = deltarconst + hbondconst
            nbind2=indexmat(dx3,dy3)
        endif

        dx3 = modulo(xPos,posNx)+1
        dy3 = modulo(yPos,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen   !already there
            deltarconst = deltarconst + hbondconst
             nbind1=indexmat(dx3,dy3)
        endif
        !build
        dx3 = modulo(xPos-3,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst + hbondconst
            if((indexmat(dx3,dy3).eq.nbind1).or.(indexmat(dx3,dy3).eq.nbind2)) then
            deltarconst = deltarconst - neibrepel
            endif
        endif


    elseif (color.eq.4) then
        !updating for blue to yellow (case D)
        dx1 = modulo(xPos-1,posNx)+1
        dy1 = modulo(yPos-1,posNy)+1
        dx2 = modulo(xPos,posNx)+1
        dy2 = modulo(yPos-1,posNy)+1

!        ematnew(dx1,dy1)=0
!        ematnew(dx2,dy2)=1

        !hbond check (3 break, 2 build, 2 flip "stay")

        !broken
        dx3 = modulo(xPos-2,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst - hbondconst
        endif
        dx3 = modulo(xPos+1,posNx)+1
        dy3 = modulo(yPos,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen  !already there
            deltarconst = deltarconst + hbondconst
            nbind2=indexmat(dx3,dy3)
        endif

        dx3 = modulo(xPos-2,posNx)+1
        dy3 = modulo(yPos,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen   !already there
            deltarconst = deltarconst + hbondconst
            nbind1=indexmat(dx3,dy3)
        endif
        !build
        dx3 = modulo(xPos+1,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst + hbondconst
            if((indexmat(dx3,dy3).eq.nbind1).or.(indexmat(dx3,dy3).eq.nbind2)) then
            deltarconst = deltarconst - neibrepel
            endif
        endif

    else
        print*,'wrong color'
        stop
    endif

end subroutine ematandhbondsD

subroutine ematandhbondsE(xPos,yPos,color,posNx,posNy, &
    deltarconst)!,ematnew)
    use commons, only:emat,hbondconst,indexmat,neibrepel
    implicit none
    integer, intent(in) :: xPos,yPos,posNx,posNy,color
    double precision, intent(out) :: deltarconst
!    integer, intent(out) :: ematnew(posNx,posNy)
    integer :: dx1,dy1,dx2,dy2,dx3,dy3,nbind1,nbind2

    nbind1=0
    nbind2=0

    if (color.eq.1) then
        !updating for green to yellow (case E)
        dx1 = modulo(xPos-1,posNx)+1
        dy1 = modulo(yPos-1,posNy)+1
        dx2 = modulo(xPos-2,posNx)+1
        dy2 = modulo(yPos-1,posNy)+1
!        ematnew(dx1,dy1)=2
!        ematnew(dx2,dy2)=2

        dx2 = modulo(xPos-1,posNx)+1
        dy2 = modulo(yPos-2,posNy)+1
!        ematnew(dx2,dy2)=1
        dx2 = modulo(xPos-2,posNx)+1
        dy2 = modulo(yPos-2,posNy)+1
!        ematnew(dx2,dy2)=0


        !hbond check (3 broken, 3 build)
        !broken
        dx3 = modulo(xPos-3,posNx)+1
        dy3 = modulo(yPos-2,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen
            deltarconst = deltarconst - hbondconst
            nbind1=indexmat(dx3,dy3)
        endif
        dx3 = modulo(xPos,posNx)+1
        dy3 = modulo(yPos-2,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen
            deltarconst = deltarconst - hbondconst
            nbind2=indexmat(dx3,dy3)
        endif
        dx3 = modulo(xPos,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst - hbondconst
        endif

        !build
        dx3 = modulo(xPos-3,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen
            deltarconst = deltarconst + hbondconst
        endif
        dx3 = modulo(xPos,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen
            deltarconst = deltarconst + hbondconst
        endif
        dx3 = modulo(xPos,posNx)+1
        dy3 = modulo(yPos-2,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst + hbondconst
            if((indexmat(dx3,dy3).eq.nbind1).or.(indexmat(dx3,dy3).eq.nbind2)) then
            deltarconst = deltarconst - neibrepel
            endif
        endif

    elseif (color.eq.2) then
        !updating for red to blue (case E)
        dx1 = modulo(xPos-1,posNx)+1
        dy1 = modulo(yPos-1,posNy)+1
        dx2 = modulo(xPos,posNx)+1
        dy2 = modulo(yPos-1,posNy)+1
!        ematnew(dx1,dy1)=2
!        ematnew(dx2,dy2)=2

        dx2 = modulo(xPos-1,posNx)+1
        dy2 = modulo(yPos-2,posNy)+1
!        ematnew(dx2,dy2)=1
        dx2 = modulo(xPos,posNx)+1
        dy2 = modulo(yPos-2,posNy)+1
!        ematnew(dx2,dy2)=0

        !hbond check (3 broken, 3 build)
        !broken
        dx3 = modulo(xPos-2,posNx)+1
        dy3 = modulo(yPos-2,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen
            deltarconst = deltarconst - hbondconst
            nbind1=indexmat(dx3,dy3)
        endif
        dx3 = modulo(xPos+1,posNx)+1
        dy3 = modulo(yPos-2,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen
            deltarconst = deltarconst - hbondconst
            nbind2=indexmat(dx3,dy3)
        endif
        dx3 = modulo(xPos-2,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst - hbondconst
        endif

        !build
        dx3 = modulo(xPos-2,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen
            deltarconst = deltarconst + hbondconst
        endif
        dx3 = modulo(xPos+1,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen
            deltarconst = deltarconst + hbondconst
        endif
        dx3 = modulo(xPos-2,posNx)+1
        dy3 = modulo(yPos-2,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst + hbondconst
            if((indexmat(dx3,dy3).eq.nbind1).or.(indexmat(dx3,dy3).eq.nbind2)) then
            deltarconst = deltarconst - neibrepel
            endif
        endif

    elseif (color.eq.3) then
        !updating for yellow to green (case E)
        dx1 = modulo(xPos-1,posNx)+1
        dy1 = modulo(yPos-1,posNy)+1
        dx2 = modulo(xPos-2,posNx)+1
        dy2 = modulo(yPos-1,posNy)+1
!        ematnew(dx1,dy1)=2
!        ematnew(dx2,dy2)=2

        dx2 = modulo(xPos-1,posNx)+1
        dy2 = modulo(yPos,posNy)+1
!        ematnew(dx2,dy2)=1
        dx2 = modulo(xPos-2,posNx)+1
        dy2 = modulo(yPos,posNy)+1
!        ematnew(dx2,dy2)=0

        !hbond check (3 broken, 3 build)
        !broken
        dx3 = modulo(xPos-3,posNx)+1
        dy3 = modulo(yPos,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen
            deltarconst = deltarconst - hbondconst
            nbind1=indexmat(dx3,dy3)
        endif
        dx3 = modulo(xPos,posNx)+1
        dy3 = modulo(yPos,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen
            deltarconst = deltarconst - hbondconst
            nbind2=indexmat(dx3,dy3)
        endif
        dx3 = modulo(xPos,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst - hbondconst
        endif

        !build
        dx3 = modulo(xPos-3,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen
            deltarconst = deltarconst + hbondconst
        endif
        dx3 = modulo(xPos,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen
            deltarconst = deltarconst + hbondconst
        endif
        dx3 = modulo(xPos,posNx)+1
        dy3 = modulo(yPos,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst + hbondconst
            if((indexmat(dx3,dy3).eq.nbind1).or.(indexmat(dx3,dy3).eq.nbind2)) then
            deltarconst = deltarconst - neibrepel
            endif
        endif

    elseif (color.eq.4) then
        !updating for blue to red (case E)
        dx1 = modulo(xPos-1,posNx)+1
        dy1 = modulo(yPos-1,posNy)+1
        dx2 = modulo(xPos,posNx)+1
        dy2 = modulo(yPos-1,posNy)+1
!        ematnew(dx1,dy1)=2
!        ematnew(dx2,dy2)=2

        dx2 = modulo(xPos-1,posNx)+1
        dy2 = modulo(yPos,posNy)+1
!        ematnew(dx2,dy2)=1
        dx2 = modulo(xPos,posNx)+1
        dy2 = modulo(yPos,posNy)+1
!        ematnew(dx2,dy2)=0

        !hbond check (3 broken, 3 build)
        !broken
        dx3 = modulo(xPos-2,posNx)+1
        dy3 = modulo(yPos,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen
            deltarconst = deltarconst - hbondconst
            nbind1=indexmat(dx3,dy3)
        endif
        dx3 = modulo(xPos+1,posNx)+1
        dy3 = modulo(yPos,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen
            deltarconst = deltarconst - hbondconst
            nbind2=indexmat(dx3,dy3)
        endif
        dx3 = modulo(xPos-2,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst - hbondconst
        endif

        !build
        dx3 = modulo(xPos-2,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen
            deltarconst = deltarconst + hbondconst
        endif
        dx3 = modulo(xPos+1,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen
            deltarconst = deltarconst + hbondconst
        endif
        dx3 = modulo(xPos-2,posNx)+1
        dy3 = modulo(yPos,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst + hbondconst
            if((indexmat(dx3,dy3).eq.nbind1).or.(indexmat(dx3,dy3).eq.nbind2)) then
            deltarconst = deltarconst - neibrepel
            endif
        endif

    else
        print*,'wrong color'
        stop
    endif

end subroutine ematandhbondsE


subroutine ematandhbondsF(xPos,yPos,color,posNx,posNy, &
    deltarconst)!,ematnew)
    use commons, only:emat,hbondconst,indexmat,neibrepel
    implicit none
    integer, intent(in) :: xPos,yPos,posNx,posNy,color
    double precision, intent(out) :: deltarconst
!    integer, intent(out) :: ematnew(posNx,posNy)
    integer :: dx1,dy1,dx2,dy2,dx3,dy3,nbind1,nbind2

nbind1=0
nbind2=0
!    ematnew(:,:)=emat(:,:)

    if (color.eq.1) then
        !updating for green to red (case F)
        !2 positions change
        dx1 = modulo(xPos-1,posNx)+1
        dy1 = modulo(yPos-1,posNy)+1
        dx2 = modulo(xPos-3,posNx)+1
        dy2 = modulo(yPos-1,posNy)+1
!        ematnew(dx1,dy1)=0
!        ematnew(dx2,dy2)=1

        dx1 = modulo(xPos-1,posNx)+1
        dy1 = modulo(yPos-2,posNy)+1
        dx2 = modulo(xPos-3,posNx)+1
        dy2 = modulo(yPos-2,posNy)+1
!        ematnew(dx1,dy1)=0
!        ematnew(dx2,dy2)=2

        !hbond check (2 broken 2 build)
        dx3 = modulo(xPos-4,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        nbind2=indexmat(dx3,dy3)
        !build
        dx3 = modulo(xPos-4,posNx)+1
        dy3 = modulo(yPos-2,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen
            deltarconst = deltarconst + hbondconst
            nbind1=indexmat(dx3,dy3)
            if((indexmat(dx3,dy3).eq.nbind2)) then
            deltarconst = deltarconst - neibrepel
            endif
        endif
        dx3 = modulo(xPos-4,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst + hbondconst
            if((indexmat(dx3,dy3).eq.nbind1)) then
            deltarconst = deltarconst - neibrepel
            endif
        endif
        !broken
        dx3 = modulo(xPos,posNx)+1
        dy3 = modulo(yPos-2,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen
            deltarconst = deltarconst - hbondconst
        endif
        dx3 = modulo(xPos,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst - hbondconst
        endif
    elseif (color.eq.2) then
        !updating for red to green (case F)
        !2 positions change
        dx1 = modulo(xPos-1,posNx)+1
        dy1 = modulo(yPos-1,posNy)+1
        dx2 = modulo(xPos+1,posNx)+1
        dy2 = modulo(yPos-1,posNy)+1
!        ematnew(dx1,dy1)=0
!        ematnew(dx2,dy2)=1

        dx1 = modulo(xPos-1,posNx)+1
        dy1 = modulo(yPos-2,posNy)+1
        dx2 = modulo(xPos+1,posNx)+1
        dy2 = modulo(yPos-2,posNy)+1
!        ematnew(dx1,dy1)=0
!        ematnew(dx2,dy2)=2

        !hbond check (2 broken 2 build)
        dx3 = modulo(xPos+2,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        nbind2=indexmat(dx3,dy3)
        !build
        dx3 = modulo(xPos+2,posNx)+1
        dy3 = modulo(yPos-2,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen
            deltarconst = deltarconst + hbondconst
            nbind1=indexmat(dx3,dy3)
            if((indexmat(dx3,dy3).eq.nbind2)) then
            deltarconst = deltarconst - neibrepel
            endif
        endif
        dx3 = modulo(xPos+2,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst + hbondconst
            if((indexmat(dx3,dy3).eq.nbind1)) then
            deltarconst = deltarconst - neibrepel
            endif
        endif
        !build
        dx3 = modulo(xPos-2,posNx)+1
        dy3 = modulo(yPos-2,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen
            deltarconst = deltarconst - hbondconst
        endif
        dx3 = modulo(xPos-2,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst - hbondconst
        endif
    elseif (color.eq.3) then
        !updating for yellow to blue (case F)
        !2 positions change
        dx1 = modulo(xPos-1,posNx)+1
        dy1 = modulo(yPos-1,posNy)+1
        dx2 = modulo(xPos-3,posNx)+1
        dy2 = modulo(yPos-1,posNy)+1
!        ematnew(dx1,dy1)=0
!        ematnew(dx2,dy2)=1

        dx1 = modulo(xPos-1,posNx)+1
        dy1 = modulo(yPos,posNy)+1
        dx2 = modulo(xPos-3,posNx)+1
        dy2 = modulo(yPos,posNy)+1
!        ematnew(dx1,dy1)=0
!        ematnew(dx2,dy2)=2

        !hbond check (2 broken 2 build)
        dx3 = modulo(xPos-4,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        nbind2=indexmat(dx3,dy3)
        !build
        dx3 = modulo(xPos-4,posNx)+1
        dy3 = modulo(yPos,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen
            deltarconst = deltarconst + hbondconst
            nbind1=indexmat(dx3,dy3)
            if((indexmat(dx3,dy3).eq.nbind2)) then
            deltarconst = deltarconst - neibrepel
            endif
        endif
        dx3 = modulo(xPos-4,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst + hbondconst
            if((indexmat(dx3,dy3).eq.nbind1)) then
            deltarconst = deltarconst - neibrepel
            endif
        endif
        !broken
        dx3 = modulo(xPos,posNx)+1
        dy3 = modulo(yPos,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen
            deltarconst = deltarconst - hbondconst
        endif
        dx3 = modulo(xPos,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst - hbondconst
        endif
    elseif (color.eq.4) then
        !updating for blue to yellow (case F)
        !2 positions change
        dx1 = modulo(xPos-1,posNx)+1
        dy1 = modulo(yPos-1,posNy)+1
        dx2 = modulo(xPos+1,posNx)+1
        dy2 = modulo(yPos-1,posNy)+1
!        ematnew(dx1,dy1)=0
!        ematnew(dx2,dy2)=1

        dx1 = modulo(xPos-1,posNx)+1
        dy1 = modulo(yPos,posNy)+1
        dx2 = modulo(xPos+1,posNx)+1
        dy2 = modulo(yPos,posNy)+1
!        ematnew(dx1,dy1)=0
!        ematnew(dx2,dy2)=2

        !hbond check (2 broken 2 build)
        dx3 = modulo(xPos+2,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        nbind2=indexmat(dx3,dy3)
        !build
        dx3 = modulo(xPos+2,posNx)+1
        dy3 = modulo(yPos,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen
            deltarconst = deltarconst + hbondconst
            nbind1=indexmat(dx3,dy3)
            if((indexmat(dx3,dy3).eq.nbind2)) then
            deltarconst = deltarconst - neibrepel
            endif
        endif
        dx3 = modulo(xPos+2,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst + hbondconst
            if((indexmat(dx3,dy3).eq.nbind1)) then
            deltarconst = deltarconst - neibrepel
            endif
        endif
        !broken
        dx3 = modulo(xPos-2,posNx)+1
        dy3 = modulo(yPos,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen
            deltarconst = deltarconst - hbondconst
        endif
        dx3 = modulo(xPos-2,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst - hbondconst
        endif
    else
        print*,'wrong color'
        stop
    endif

end subroutine ematandhbondsF


subroutine ematandhbondsG(xPos,yPos,color,posNx,posNy, &
    deltarconst)!,ematnew)
    use commons, only:emat,hbondconst,indexmat,neibrepel
    implicit none
    integer, intent(in) :: xPos,yPos,posNx,posNy,color
    double precision, intent(out) :: deltarconst
!    integer, intent(out) :: ematnew(posNx,posNy)
    integer :: dx1,dy1,dx2,dy2,dx3,dy3,nbind1,nbind2

nbind1=0
nbind2=0


    if (color.eq.1) then
        !updating for green to blue (case G)
        !3 positions change
        dx1 = modulo(xPos-1,posNx)+1
        dy1 = modulo(yPos-2,posNy)+1
        dx2 = modulo(xPos-2,posNx)+1
        dy2 = modulo(yPos-1,posNy)+1
!        ematnew(dx1,dy1)=0
!        ematnew(dx2,dy2)=2

        dx1 = modulo(xPos-2,posNx)+1
        dy1 = modulo(yPos-2,posNy)+1
        dx2 = modulo(xPos-3,posNx)+1
        dy2 = modulo(yPos-1,posNy)+1
!        ematnew(dx1,dy1)=0
!        ematnew(dx2,dy2)=2

        dx1 = modulo(xPos-1,posNx)+1
        dy1 = modulo(yPos-1,posNy)+1
        dx2 = modulo(xPos-3,posNx)+1
        dy2 = modulo(yPos-2,posNy)+1
!        ematnew(dx1,dy1)=0
!        ematnew(dx2,dy2)=1

        !hbond check (2 broken 2 build)
        dx3 = modulo(xPos-4,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        nbind2=indexmat(dx3,dy3)
        !build
        dx3 = modulo(xPos-4,posNx)+1
        dy3 = modulo(yPos-2,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst + hbondconst
            nbind1=indexmat(dx3,dy3)
            if((indexmat(dx3,dy3).eq.nbind2)) then
            deltarconst = deltarconst - neibrepel
            endif
        endif
        dx3 = modulo(xPos-4,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen
            deltarconst = deltarconst + hbondconst
            if((indexmat(dx3,dy3).eq.nbind1)) then
            deltarconst = deltarconst - neibrepel
            endif
        endif
        !broken
        dx3 = modulo(xPos,posNx)+1
        dy3 = modulo(yPos-2,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen
            deltarconst = deltarconst - hbondconst
        endif
        dx3 = modulo(xPos,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst - hbondconst
        endif
    elseif (color.eq.2) then
        !updating for red to yellow (case G)
        !3 positions change
        dx1 = modulo(xPos,posNx)+1
        dy1 = modulo(yPos-2,posNy)+1
        dx2 = modulo(xPos,posNx)+1
        dy2 = modulo(yPos-1,posNy)+1
!        ematnew(dx1,dy1)=0
!        ematnew(dx2,dy2)=2

        dx1 = modulo(xPos-1,posNx)+1
        dy1 = modulo(yPos-2,posNy)+1
        dx2 = modulo(xPos+1,posNx)+1
        dy2 = modulo(yPos-1,posNy)+1
!        ematnew(dx1,dy1)=0
!        ematnew(dx2,dy2)=2

        dx1 = modulo(xPos-1,posNx)+1
        dy1 = modulo(yPos-1,posNy)+1
        dx2 = modulo(xPos+1,posNx)+1
        dy2 = modulo(yPos-2,posNy)+1
!        ematnew(dx1,dy1)=0
!        ematnew(dx2,dy2)=1

        !hbond check (2 broken 2 build)
        dx3 = modulo(xPos+2,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        nbind2=indexmat(dx3,dy3)
        !build
        dx3 = modulo(xPos+2,posNx)+1
        dy3 = modulo(yPos-2,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst + hbondconst
            nbind1=indexmat(dx3,dy3)
            if((indexmat(dx3,dy3).eq.nbind2)) then
            deltarconst = deltarconst - neibrepel
            endif
        endif
        dx3 = modulo(xPos+2,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen
            deltarconst = deltarconst + hbondconst
            if((indexmat(dx3,dy3).eq.nbind1)) then
            deltarconst = deltarconst - neibrepel
            endif
        endif
        !broken
        dx3 = modulo(xPos-2,posNx)+1
        dy3 = modulo(yPos-2,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen
            deltarconst = deltarconst - hbondconst
        endif
        dx3 = modulo(xPos-2,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst - hbondconst
        endif
    elseif (color.eq.3) then
        !updating for yellow to red (case G)
        !3 positions change
        dx1 = modulo(xPos-2,posNx)+1
        dy1 = modulo(yPos,posNy)+1
        dx2 = modulo(xPos-2,posNx)+1
        dy2 = modulo(yPos-1,posNy)+1
!        ematnew(dx1,dy1)=0
!        ematnew(dx2,dy2)=2
        !
        dx1 = modulo(xPos-1,posNx)+1
        dy1 = modulo(yPos,posNy)+1
        dx2 = modulo(xPos-3,posNx)+1
        dy2 = modulo(yPos-1,posNy)+1
!        ematnew(dx1,dy1)=0
!        ematnew(dx2,dy2)=2

        dx1 = modulo(xPos-1,posNx)+1
        dy1 = modulo(yPos-1,posNy)+1
        dx2 = modulo(xPos-3,posNx)+1
        dy2 = modulo(yPos,posNy)+1
!        ematnew(dx1,dy1)=0
!        ematnew(dx2,dy2)=1

        !hbond check (2 broken 2 build)
        dx3 = modulo(xPos-4,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        nbind2=indexmat(dx3,dy3)
        !build
        dx3 = modulo(xPos-4,posNx)+1
        dy3 = modulo(yPos,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst + hbondconst
            nbind1=indexmat(dx3,dy3)
            if((indexmat(dx3,dy3).eq.nbind2)) then
            deltarconst = deltarconst - neibrepel
            endif
        endif
        dx3 = modulo(xPos-4,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen
            deltarconst = deltarconst + hbondconst
            if((indexmat(dx3,dy3).eq.nbind1)) then
            deltarconst = deltarconst - neibrepel
            endif
        endif
        !broken
        dx3 = modulo(xPos,posNx)+1
        dy3 = modulo(yPos,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen
            deltarconst = deltarconst - hbondconst
        endif
        dx3 = modulo(xPos,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst - hbondconst
        endif
    elseif (color.eq.4) then
        !updating for blue to green (case G)
        !3 positions change
        dx1 = modulo(xPos-1,posNx)+1
        dy1 = modulo(yPos,posNy)+1
        dx2 = modulo(xPos,posNx)+1
        dy2 = modulo(yPos-1,posNy)+1
!        ematnew(dx1,dy1)=0
!        ematnew(dx2,dy2)=2

        dx1 = modulo(xPos,posNx)+1
        dy1 = modulo(yPos,posNy)+1
        dx2 = modulo(xPos+1,posNx)+1
        dy2 = modulo(yPos-1,posNy)+1
!        ematnew(dx1,dy1)=0
!        ematnew(dx2,dy2)=2

        dx1 = modulo(xPos-1,posNx)+1
        dy1 = modulo(yPos-1,posNy)+1
        dx2 = modulo(xPos+1,posNx)+1
        dy2 = modulo(yPos,posNy)+1
!        ematnew(dx1,dy1)=0
!        ematnew(dx2,dy2)=1

        !hbond check (2 broken 2 build)
        dx3 = modulo(xPos+2,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        nbind2=indexmat(dx3,dy3)
        !build
        dx3 = modulo(xPos+2,posNx)+1
        dy3 = modulo(yPos,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst + hbondconst
            nbind1=indexmat(dx3,dy3)
            if((indexmat(dx3,dy3).eq.nbind2)) then
            deltarconst = deltarconst - neibrepel
            endif
        endif
        dx3 = modulo(xPos+2,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen
            deltarconst = deltarconst + hbondconst
            if((indexmat(dx3,dy3).eq.nbind1)) then
            deltarconst = deltarconst - neibrepel
            endif
        endif
        !broken
        dx3 = modulo(xPos-2,posNx)+1
        dy3 = modulo(yPos,posNy)+1
        if(emat(dx3,dy3).eq.1) then !Neighbouring Nitrogen
            deltarconst = deltarconst - hbondconst
        endif
        dx3 = modulo(xPos-2,posNx)+1
        dy3 = modulo(yPos-1,posNy)+1
        if(emat(dx3,dy3).eq.2) then !Neighbouring Oxygen
            deltarconst = deltarconst - hbondconst
        endif
    else
        print*,'wrong color'
        stop
    endif

end subroutine ematandhbondsG



subroutine printstate(config,EMAT,rate,STEPS)
    use commons, only: ntriang,posNx,posNy
    implicit none
    integer, intent(in) :: config(ntriang,5),STEPS
    integer, intent(in) :: emat(posNx,posNx)! changed emat matrix, posNx, posNy
    integer :: dx,dy
    double precision :: newinit(ntriang,3),rate ! ntriang, (Grid (posx, posy), color)
    integer :: jcheck1,jcheck2
    integer :: hbondcount,occupationcount,hbondmat(posNx+1,posNx+1)! Matrix with hbonds
    character ( len = 24 ) :: data_filename
    integer :: i,k,j
    logical :: exist
    !**********************
    hbondmat(:,:)=0
    hbondcount=0
    occupationcount=0

    do i=1,ntriang
        !colors 1 green, 2 red, 3 yellow, 4 blu
        if(config(i,1).eq.1) then
            !H O N green
            dx = modulo(config(i,2)-2,posNx)+1
            dy = modulo(config(i,3)-2,posNy)+1
            newinit(i,1)=dx  !pos grid x
            newinit(i,2)=dy  !pos grid y
            newinit(i,3)=1 !color
        elseif(config(i,1).eq.2) then
            !H O N red
            dx = modulo(config(i,2)-1,posNx)+1
            dy = modulo(config(i,3)-2,posNy)+1
            newinit(i,1)=dx  !pos grid x
            newinit(i,2)=dy  !pos grid y
            newinit(i,3)=2 !color
        elseif(config(i,1).eq.3) then
             !H O N yelloW
            dy = modulo(config(i,3)-1,posNy)+1
            dx = modulo(config(i,2)-2,posNx)+1
            newinit(i,1)=dx  !pos grid x
            newinit(i,2)=dy  !pos grid y
            newinit(i,3)=3 !yellow
        elseif(config(i,1).eq.4) then
            !H O N !BLUE
            dx = modulo(config(i,2)-1,posNx)+1
            dy = modulo(config(i,3)-1,posNy)+1
            newinit(i,1)=dx  !pos grid x
            newinit(i,2)=dy  !pos grid y
            newinit(i,3)=4 !color
        endif
    enddo

    !Printing triangle configuration
    WRITE(data_filename,'(a,i9.9,a)')'triangles',steps,'.vtk'
    call vtk_triangles(ntriang,config,newinit,data_filename)

    !Printing H bonds
    do j=1,posNx
        do k=1,posNy
            if(emat(j,k).ne.0) occupationcount=occupationcount+1
            if(emat(j,k).eq.2) then
                jcheck1=modulo(j-1,posNx)
                if(jcheck1.eq.0) jcheck1 = posNx
                jcheck2=modulo(j-1,posNx)+2
                if(jcheck2.eq.posNx+1) jcheck2 = 1
                if(emat(jcheck1,k).eq.1) then
                    hbondmat(jcheck1,k)=1
                    hbondcount=hbondcount+1
                elseif(emat(jcheck2,k).eq.1) then
                    hbondmat(j,k)=1
                    hbondcount=hbondcount+1
                endif
            endif
        enddo
    enddo
    !periodic boundaries, double row and line 1
    hbondmat(posNx+1,1:posNy)=hbondmat(1,1:posNy)
    hbondmat(1:posNx,posNy+1)=hbondmat(1:posNx,1)
    hbondmat(posNx+1,posNy+1)=hbondmat(1,1)
    !Write H bonds to file
    WRITE(data_filename,'(a,i9.9,a)')'hlines',STEPS,'.vtk'
    call write_vtk_lines( 1,1,posNx,posNy,hbondmat,data_filename)

    !Append the statistics to a file (number of H bonds and coverage)
    inquire(file="hbondcount.txt", exist=exist)
    if (exist) then
        open(12, file="hbondcount.txt", status="old", position="append", action="write")
    else
        open(12, file="hbondcount.txt", status="new", action="write")
    end if
    write(12, *) steps,hbondcount
    close(12)

    inquire(file="hbond_GLY.txt", exist=exist)
    if (exist) then
        open(12, file="hbond_GLY.txt", status="old", position="append", action="write")
    else
        open(12, file="hbond_GLY.txt", status="new", action="write")
    end if
    write(12, *) steps,hbondcount/ntriang
    close(12)

    inquire(file="rate.txt", exist=exist)
    if (exist) then
        open(12, file="rate.txt", status="old", position="append", action="write")
    else
        open(12, file="rate.txt", status="new", action="write")
    end if
    write(12, *) steps,rate
    close(12)

    inquire(file="coverage.txt", exist=exist)
    if (exist) then
        open(12, file="coverage.txt", status="old", position="append", action="write")
    else
        open(12, file="coverage.txt", status="new", action="write")
    end if
    write(12, *) steps,DBLE(occupationcount)/(DBLE(posNx*posNy))
    close(12)

end subroutine printstate

subroutine vtk_triangles(ntriang,config,colorcoords,data_filename)


    !*****************************************************************************
    !! VTK_PUVW_WRITE write VTK triangle

    implicit none
    integer, intent(in)  :: ntriang,config(ntriang,5)
    double precision, intent(in) :: colorcoords(ntriang,3)
    character ( len = 24 ) :: data_filename
    integer :: k,ios
    integer ( kind = 4 ) output_unit
    double precision :: sroot

     !*****************************************************************************
    sroot= dsqrt(2.D0)
     !*****************************************************************************

    output_unit=50

    open ( output_unit, file = data_filename, status = 'replace', iostat = ios )

    write ( output_unit, '(a)' ) '# vtk DataFile Version 3.1'
    write ( output_unit, '(a)' ) 'Gradient and points of surface'
    write ( output_unit, '(a)' ) 'ASCII'
    write ( output_unit, '(a)' ) 'DATASET POLYDATA'
    write ( output_unit, '(a)' )


    write ( output_unit,* ) 'POINTS ',3*ntriang,' double'
    do k=1,ntriang

        if(INT(colorcoords(k,3)).eq.1) then !green
            write ( output_unit,*) colorcoords(k,1)+1, sroot*(colorcoords(k,2)+1), 0.D0 !N
            write ( output_unit,*) colorcoords(k,1), sroot*colorcoords(k,2), 0.D0 !O
            write ( output_unit,*) colorcoords(k,1)+1, sroot*colorcoords(k,2), 0.D0 !O

        elseif(INT(colorcoords(k,3)).eq.2) then !red
            write ( output_unit,*) colorcoords(k,1), sroot*(colorcoords(k,2)+1), 0.D0 !N
            write ( output_unit,*) colorcoords(k,1), sroot*colorcoords(k,2), 0.D0 !O
            write ( output_unit,*) colorcoords(k,1)+1, sroot*colorcoords(k,2), 0.D0 !O

        elseif(INT(colorcoords(k,3)).eq.3) then !yellow
            write ( output_unit,*) colorcoords(k,1)+1, sroot*colorcoords(k,2), 0.D0 !N
            write ( output_unit,*) colorcoords(k,1), sroot*(colorcoords(k,2)+1), 0.D0 !O
            write ( output_unit,*) colorcoords(k,1)+1, sroot*(colorcoords(k,2)+1), 0.D0 !O

        elseif(INT(colorcoords(k,3)).eq.4) then !blue
            write ( output_unit,*) colorcoords(k,1), sroot*colorcoords(k,2), 0.D0 !N
            write ( output_unit,*) colorcoords(k,1), sroot*(colorcoords(k,2)+1), 0.D0 !O
            write ( output_unit,*) colorcoords(k,1)+1, sroot*(colorcoords(k,2)+1), 0.D0 !O

        else
            print*,'wrong colour code',colorcoords(k,3),'colorcoords(k,3)',k
            stop
        endif

    enddo

    write ( output_unit,* ) 'POLYGONS',ntriang,ntriang*5
    do k=1,ntriang
        write ( output_unit,*) 4, (3*(k-1)),(3*(k-1))+1 ,(3*(k-1))+2 , (3*(k-1))
    enddo


    write ( output_unit, * ) 'POINT_DATA', 3*ntriang
    write ( output_unit, '(a)' ) 'SCALARS AtomTriangle INT'
    write ( output_unit, '(a)' ) 'LOOKUP_TABLE default'
    do k=1,ntriang
        write ( output_unit,*) 3 !Nickel
        write ( output_unit,*) 2 !Oxygen
        write ( output_unit,*) 2 !Oxygen
    enddo

    write ( output_unit, * ) 'CELL_DATA',ntriang

    write ( output_unit, '(a)' ) 'SCALARS color int 1'
    write ( output_unit, '(a)' ) 'LOOKUP_TABLE default'
    do k=1,ntriang
        write ( output_unit, * ) config(k,1)
    enddo

    write ( output_unit, '(a)' ) 'SCALARS molrot int 1'
    write ( output_unit, '(a)' ) 'LOOKUP_TABLE default'
    do k=1,ntriang
        write ( output_unit,*) config(k,4)
    enddo

    write ( output_unit, '(a)' ) 'SCALARS triangle_index int 1'
    write ( output_unit, '(a)' ) 'LOOKUP_TABLE default'
    do k=1,ntriang
        write ( output_unit,*) k
    enddo

    write ( output_unit, '(a)' ) 'SCALARS osmov int 1'
    write ( output_unit, '(a)' ) 'LOOKUP_TABLE default'
    do k=1,ntriang
        write ( output_unit,*) config(k,5)
    enddo

    close(output_unit)
    return
end subroutine vtk_triangles

subroutine centreofmasscompare(ntriang,colorcoords,data_filename)


    !*****************************************************************************
    !! VTK_PUVW_WRITE write VTK triangle

    implicit none
    integer, intent(in)  :: ntriang
    double precision, intent(in) :: colorcoords(ntriang,3)
    character ( len = 24 ) :: data_filename
    integer :: k,ios
    integer ( kind = 4 ) output_unit
    double precision :: sroot

     !*****************************************************************************
    sroot= dsqrt(2.D0)
     !*****************************************************************************

    output_unit=50

    open ( output_unit, file = data_filename, status = 'replace', iostat = ios )

    write ( output_unit, '(a)' ) '# vtk DataFile Version 3.1'
    write ( output_unit, '(a)' ) 'Gradient and points of surface'
    write ( output_unit, '(a)' ) 'ASCII'
    write ( output_unit, '(a)' ) 'DATASET POLYDATA'
    write ( output_unit, '(a)' )


    write ( output_unit,* ) 'POINTS ',3*ntriang,' double'
    do k=1,ntriang

        if(INT(colorcoords(k,3)).eq.1) then !green
            write ( output_unit,*) colorcoords(k,1)+1, sroot*(colorcoords(k,2)+1), 0.D0 !N
            write ( output_unit,*) colorcoords(k,1), sroot*colorcoords(k,2), 0.D0 !O
            write ( output_unit,*) colorcoords(k,1)+1, sroot*colorcoords(k,2), 0.D0 !O

        elseif(INT(colorcoords(k,3)).eq.2) then !red
            write ( output_unit,*) colorcoords(k,1), sroot*(colorcoords(k,2)+1), 0.D0 !N
            write ( output_unit,*) colorcoords(k,1), sroot*colorcoords(k,2), 0.D0 !O
            write ( output_unit,*) colorcoords(k,1)+1, sroot*colorcoords(k,2), 0.D0 !O

        elseif(INT(colorcoords(k,3)).eq.3) then !yellow
            write ( output_unit,*) colorcoords(k,1)+1, sroot*colorcoords(k,2), 0.D0 !N
            write ( output_unit,*) colorcoords(k,1), sroot*(colorcoords(k,2)+1), 0.D0 !O
            write ( output_unit,*) colorcoords(k,1)+1, sroot*(colorcoords(k,2)+1), 0.D0 !O

        elseif(INT(colorcoords(k,3)).eq.4) then !blue
            write ( output_unit,*) colorcoords(k,1), sroot*colorcoords(k,2), 0.D0 !N
            write ( output_unit,*) colorcoords(k,1), sroot*(colorcoords(k,2)+1), 0.D0 !O
            write ( output_unit,*) colorcoords(k,1)+1, sroot*(colorcoords(k,2)+1), 0.D0 !O

        else
            print*,'wrong colour code',colorcoords(k,3),'colorcoords(k,3)',k
            stop
        endif
     enddo

end subroutine centreofmasscompare


subroutine write_vtk_points(posxmin,posymin,posxmax,posymax)


    !*****************************************************************************
    !! VTK_PUVW_WRITE write VTK file for screening and density test

    implicit none
    integer, intent(in)  ::posxmin,posymin,posxmax,posymax
    character ( len = 24 ) :: data_filename
    integer :: k,m,ios
    integer ( kind = 4 ) output_unit

     !*****************************************************************************

    output_unit=50
    WRITE(data_filename,'(a,i9.9,a,i9.9,a,i9.9,a)')'pointssurf.vtk'

    open ( output_unit, file = data_filename, status = 'replace', iostat = ios )

    write ( output_unit, '(a)' ) '# vtk DataFile Version 3.1'
    write ( output_unit, '(a)' ) 'Gradient and points of surface'
    write ( output_unit, '(a)' ) 'ASCII'
    write ( output_unit, '(a)' ) 'DATASET POLYDATA'
    write ( output_unit, '(a)' )
    write ( output_unit,* ) 'POINTS ',(posxmax+1) * (posymax+1),' double'

    do M = 1,posxmax+1
        do K = 1,posymax+1
            write ( output_unit,*) posxmin+M-1,(posymin+K-1)*sqrt(2.D0),-0.5D0
        end do
    end do


    write ( output_unit, * ) 'POINT_DATA', (posxmax+1) * (posymax+1)
    write ( output_unit, '(a)' ) 'SCALARS virtualPos INT'
    write ( output_unit, '(a)' ) 'LOOKUP_TABLE default'
    do M = 1,posxmax
        do K = 1,posymax
            write ( output_unit,*) 0
        end do
        write ( output_unit,*) 1
    end do
    do M = posxmax+1,posxmax+1
        do K = 1,posymax+1
            write ( output_unit,*) 1
        end do
    end do


    close(output_unit)

    return
end subroutine write_vtk_points


subroutine write_vtk_points_defect(posxmin,posymin,posxmax,posymax,&
    mdef1type,mdef2type,mdef1sizex,mdef1sizey,mdef2sizex,mdef2sizey,&
    mdef1posx,mdef1posy,mdef2posx,mdef2posy)


    !*****************************************************************************
    !! VTK_PUVW_WRITE write VTK file for screening and density test

    implicit none
    integer, intent(in)  ::posxmin,posymin,posxmax,posymax,&
        mdef1type,mdef2type,mdef1sizex,mdef1sizey,mdef2sizex,mdef2sizey,&
        mdef1posx,mdef1posy,mdef2posx,mdef2posy
    character ( len = 24 ) :: data_filename
    integer :: k,m,ios
    integer ( kind = 4 ) output_unit
    integer :: virtualpos((posxmax+1) * (posymax+1))

     !*****************************************************************************

    output_unit=50
    WRITE(data_filename,'(a,i9.9,a,i9.9,a,i9.9,a)')'pointssurf.vtk'

    open ( output_unit, file = data_filename, status = 'replace', iostat = ios )

    write ( output_unit, '(a)' ) '# vtk DataFile Version 3.1'
    write ( output_unit, '(a)' ) 'Gradient and points of surface'
    write ( output_unit, '(a)' ) 'ASCII'
    write ( output_unit, '(a)' ) 'DATASET POLYDATA'
    write ( output_unit, '(a)' )
    write ( output_unit,* ) 'POINTS ',(posxmax+1) * (posymax+1),' double'

    do M = 1,posxmax+1
        do K = 1,posymax+1
            write ( output_unit,*) posxmin+M-1,(posymin+K-1)*sqrt(2.D0),-0.5D0
        end do
    end do



    write ( output_unit, * ) 'POINT_DATA', (posxmax+1) * (posymax+1)
    write ( output_unit, '(a)' ) 'SCALARS virtualPos INT'
    write ( output_unit, '(a)' ) 'LOOKUP_TABLE default'

    virtualpos(:)=0 !init virtual positions
    !Add virtual boundary marks
    do M = 1,posxmax
        virtualpos(M*(posymax+1))=1
    end do
    do K = 1,posymax+1
        virtualpos((posxmax+1)*posymax+K)=1
    end do
    !Add defects
    if(mdef1type.eq.1) then
        do M=1,mdef1sizex
            do K=1,mdef1sizey
                virtualpos((posymax+1)*(mdef1posx-2+M)+ mdef1posy-1+K)=3 !3 for defect
            enddo
        enddo
    endif
    if(mdef2type.eq.1) then
        do M=1,mdef2sizex
            do K=1,mdef2sizey
                virtualpos((posymax+1)*(mdef2posx-2+M)+ mdef2posy-1+K)=3 !3 for defect
            enddo
        enddo
    endif

    do K = 1,(posymax+1)*(posxmax+1)
        write ( output_unit,*) virtualpos(K)
    end do



    close(output_unit)

    return
end subroutine write_vtk_points_defect

subroutine writematrix(stotal,natx,naty,data_filename)
    implicit none
    integer, intent(in) :: natx,naty
    integer, intent(in) :: stotal(natx,naty)


    character ( len = 24 ),intent(in) :: data_filename
    integer ( kind = 4 ) output_unit
    integer ::j,i,ios

    !*************************************************************



    output_unit=80
    !    WRITE(data_filename,'(a)')'matrixtest.txt'

    open ( output_unit, file = data_filename, status = 'replace', iostat = ios )


    do i=1, naty!+1
        do j=1, natx!+1
            write(output_unit, *) i,j, stotal(j,i)
        end do
        write(output_unit, *) ''
    end do

    close(output_unit)


end subroutine writematrix


subroutine write_vtk_lines(posxmin,posymin,posxmax,posymax,&
    hbondmat,data_filename)

    !*****************************************************************************
    !! VTK_PUVW_WRITE write VTK file for lines
    implicit none
    integer, intent(in)  :: posxmin,posymin,posxmax,posymax,hbondmat(posxmax+1,posymax+1)
    character ( len = 24 ), intent(in) :: data_filename
    !integer, intent(in)  :: i4,j4,k4
    integer :: k,n,m,ios,NATOMS,indexcount
    integer ( kind = 4 ) output_unit
    integer :: indexmat((posxmax+1) * (posymax+1))
     !*****************************************************************************
    NATOMS=(posxmax+1) * (posymax+1)
    indexcount=0
    indexmat(:)=0

    output_unit=50
    open ( output_unit, file = data_filename, status = 'replace', iostat = ios )

    write ( output_unit, '(a)' ) '# vtk DataFile Version 3.1'
    write ( output_unit, '(a)' ) 'Gradient and points of surface'
    write ( output_unit, '(a)' ) 'ASCII'
    write ( output_unit, '(a)' ) 'DATASET POLYDATA'
    write ( output_unit, '(a)' )

    write ( output_unit,* ) 'POINTS ',natoms,' double'
    do M = 1,posxmax+1
        do K = 1,posymax+1
            write ( output_unit,*) posxmin+M-1,(posymin+K-1)*sqrt(2.D0),0.D0
            indexcount=indexcount+1
            IF((hbondmat(M,K).eq.1).or.(hbondmat(M,K).eq.-1)) then
                indexmat(indexcount)=hbondmat(M,K)
            else
                indexmat(indexcount)=0
            endif
        end do
    end do

    !PRINT LINES between atoms
    write ( output_unit,* ) 'LINES',natoms,natoms*3
    do N = 1, natoms-posymax-1
        write ( output_unit,*) 2,N-1,N+posymax
    end do
    do N = natoms-posymax,natoms
        write ( output_unit,*) 2,N-1,N-1
    end do

    write ( output_unit, * ) 'CELL_DATA',NATOMS
    write ( output_unit, '(a)' ) 'SCALARS hbonds int 1'
    write ( output_unit, '(a)' ) 'LOOKUP_TABLE default'
    do N = 1, natoms
        write ( output_unit,*) indexmat(N)
    end do


    close(output_unit)

    return
end subroutine write_vtk_lines

subroutine read_vtk_triangles(ntriang,config,colorcoords,emat,&
            indexmat,posNx,posNy,data_filename)


    !*****************************************************************************
    !! VTK_PUVW_WRITE write VTK triangle

    implicit none
    integer, intent(in) :: ntriang,posNx,posNy
    integer, intent(out)  :: config(ntriang,5),emat(posNx,posNy),indexmat(posNx,posNy)
    double precision, intent(out) :: colorcoords(ntriang,3)
    double precision :: coordsN00xyz(ntriang,3,3)
    character ( len = 24 ), intent(in) :: data_filename
    integer :: k,ios
    integer ( kind = 4 ) output_unit
    double precision :: sroot
    character ( len = 24 ) :: test
    integer :: triang

     !*****************************************************************************
    sroot= dsqrt(2.D0)
         !*****************************************************************************

    config(:,:)=0
    indexmat(:,:)=0

    output_unit=50
    open ( output_unit, file = data_filename, status = 'old', iostat = ios )

    read ( output_unit, * )
    read ( output_unit, * )
    read ( output_unit, * )
    read ( output_unit, * )
    read ( output_unit, * )


    read ( output_unit,* )

    !reading colorcoords
    do k=1,ntriang
        read ( output_unit,*) coordsN00xyz(k,1,1), coordsN00xyz(k,1,2), coordsN00xyz(k,1,3) !N
        read ( output_unit,*) coordsN00xyz(k,2,1), coordsN00xyz(k,2,2), coordsN00xyz(k,2,3) !O
        read ( output_unit,*) coordsN00xyz(k,3,1), coordsN00xyz(k,3,2), coordsN00xyz(k,3,3) !O
    enddo

    !former polygons
    read ( output_unit,* ) test,triang
    if(triang.ne.ntriang) then
        print*,'number of triangles in file does not equal number in system.in'
        stop
    endif
    do k=1,ntriang
        read ( output_unit,*)
    enddo

    !atom triangle scalars
    read ( output_unit,*)
    read ( output_unit,*)
    read ( output_unit,*)

    do k=1,ntriang
        read ( output_unit,*) !3 Nickel
        read ( output_unit,*) !2 Oxygen
        read ( output_unit,*) !2 Oxygen
    enddo

    !color
    read ( output_unit, * )
    read ( output_unit, * )
    read ( output_unit, * )
    do k=1,ntriang
        read ( output_unit, * ) config(k,1)
    enddo

    !mol rot
    read ( output_unit, * )
    read ( output_unit, * )
    do k=1,ntriang
        read ( output_unit,*) config(k,4)
    enddo
    !triangle index
    read ( output_unit, * )
    read ( output_unit, * )
    do k=1,ntriang
        read ( output_unit, * )
    enddo

    !osmov
    read ( output_unit, * )
    read ( output_unit, * )
    do k=1,ntriang
        read ( output_unit, * ) config(k,5)
    enddo

    close(output_unit)

    colorcoords(:,3)=config(:,1)


    do k=1,ntriang
        colorcoords(k,1)=nint(coordsN00xyz(k,2,1))
        config(k,2)=nint(coordsN00xyz(k,1,1))

        if(INT(colorcoords(k,3)).eq.1) then !green
            colorcoords(k,2)=nint(coordsN00xyz(k,2,2)/sroot)
             config(k,3)=nint(coordsN00xyz(k,1,2)/sroot)
        elseif(INT(colorcoords(k,3)).eq.2) then !red
            colorcoords(k,2)=nint(coordsN00xyz(k,2,2)/sroot)
             config(k,3)=nint(coordsN00xyz(k,1,2)/sroot)
        elseif(INT(colorcoords(k,3)).eq.3) then !yellow
            colorcoords(k,2)=nint(coordsN00xyz(k,1,2)/sroot)
            config(k,3)=nint(coordsN00xyz(k,1,2)/sroot)
        elseif(INT(colorcoords(k,3)).eq.4) then !blue
            colorcoords(k,2)=nint(coordsN00xyz(k,1,2)/sroot)
            config(k,3)=nint(coordsN00xyz(k,1,2)/sroot)
        else
            print*,'wrong colour code',colorcoords(k,3),'colorcoords(k,3)',k
            stop
        endif


        emat(MOD(nint(coordsN00xyz(k,1,1))-1,posNx)+1,MOD(nint(coordsN00xyz(k,1,2)/sroot)-1,posNy)+1)=1
        emat(MOD(nint(coordsN00xyz(k,2,1))-1,posNx)+1,MOD(nint(coordsN00xyz(k,2,2)/sroot)-1,posNy)+1)=2
        emat(MOD(nint(coordsN00xyz(k,3,1))-1,posNx)+1,MOD(nint(coordsN00xyz(k,3,2)/sroot)-1,posNy)+1)=2

        indexmat(MOD(nint(coordsN00xyz(k,1,1))-1,posNx)+1,MOD(nint(coordsN00xyz(k,1,2)/sroot)-1,posNy)+1)=k
        indexmat(MOD(nint(coordsN00xyz(k,2,1))-1,posNx)+1,MOD(nint(coordsN00xyz(k,2,2)/sroot)-1,posNy)+1)=k
        indexmat(MOD(nint(coordsN00xyz(k,3,1))-1,posNx)+1,MOD(nint(coordsN00xyz(k,3,2)/sroot)-1,posNy)+1)=k


    enddo

    return
end subroutine read_vtk_triangles

