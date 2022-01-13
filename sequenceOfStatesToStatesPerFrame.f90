program sequenceOfStatesToStatesPerFrame
implicit none

! The maximum possible number of frames
integer,parameter :: n_measure_max = 50000

! The number of frames to spit out
integer :: n_measure != 9999

! The number of lipids in the system
integer :: Nlipids != 1687

integer :: Nstart, Nlipid
integer,allocatable :: states(:,:)
character(2*n_measure_max+100) :: aline
integer :: iostate
integer :: i, j, k, n
integer :: i1, i2
integer :: j1

!
! Example usage:
! gfortran sequenceOfStatesToStatesPerFrame.f90 -o sequenceOfStatesToStatesPerFrame.out; ./sequenceOfStatesToStatesPerFrame.out 5 1687 < <(sed 1,6d getHMMmostProbableStates.out | sed '$d') > hmm-hiddenstates.txt
!

if (iargc() /= 2) then
  write(0,FMT="(A)") "ERROR: Need exactly two arguments:"
  write(0,FMT="(A)") "      (1) the number of frames"
  write(0,FMT="(A)") "      (2) the number of lipids"
  stop
end if

call get_command_argument(1,aline)
read(aline,FMT=*,iostat=iostate) n_measure

if (iostate /= 0) then
  write(0,FMT="(A)") "ERROR: the first argument needs to be an integer"
  stop
end if

call get_command_argument(2,aline)
read(aline,FMT=*,iostat=iostate) Nlipids

if (iostate /= 0) then
  write(0,FMT="(A)") "ERROR: the second argument needs to be an integer"
  stop
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(states(n_measure,Nlipids))

! All lipids which are not included in a sequence
! are assumed to not be in the membrane of
! interest
states = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do
  read(5,FMT="(A)",iostat=iostate) aline
  if (is_iostat_end(iostate)) then
    exit
  elseif (iostate /= 0) then
    write(0,FMT="(A)") "Error. Bad formatting in at least one sequence"
    stop
  end if

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! The first number should be the frame number
  ! where the sequence starts
  
  Nstart = 0

  i1 = 0 !1
  do
    i2 = index(aline(i1+1:),' ')
    if (i2 == 0) then
      write(0,FMT="(A)") "Error. Bad formatting in sequence start"
      stop
    end if
    i2 = i1 + i2

    !if ((i2 > i1+1).or.((i1==1).and.(aline(i1:i1)/=' '))) then
    if (i2 > i1+1) then
      read(aline(i1+1:i2),FMT=*,iostat=iostate) Nstart
      if (iostate /= 0) then
        write(0,FMT="(A)") "Error. Bad formatting in sequence start"
        stop
      end if

      i1 = i2
      exit
    end if

    i1 = i2
  end do

  if (Nstart < 1) then
    write(0,FMT="(A)") "Error. Missing the sequence start"
    stop
  end if

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! The second number (which has a colon appended)
  ! should be the lipid number

  Nlipid = 0

  do
    i2 =index(aline(i1+1:),' ')
    if (i2 == 0) then
      write(0,FMT="(A)") "Error. Bad formatting in lipid number"
      stop
    end if
    i2 = i1 + i2

    if (i2 > i1+1) then
      read(aline(i1:i2-2),FMT=*,iostat=iostate) Nlipid
      if (iostate /= 0) then
        write(0,FMT="(A)") "Error. Bad formatting in lipid number"
        stop
      end if

      i1 = i2
      exit
    end if

    i1 = i2
  end do

  if (Nlipid < 1) then
    write(0,FMT="(A)") "Error. Missing the lipid number"
    stop
  end if

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! All following numbers are the hidden states in
  ! sequence

  i = index(aline(i1+1:),'1')
  j = index(aline(i1+1:),'2')
  do n = Nstart, n_measure
    if (i < j) then
      if (i > 0) then
        states(n,Nlipid) = 1
        i1 = i1 + i
        j = j - i

        i = index(aline(i1+1:),'1')
      else
        if (j == i) exit ! Both must be equal to zero
        states(n,Nlipid) = 2
        i1 = i1 + j
        j = index(aline(i1+1:),'2')
      end if
    else
      if (j > 0) then
        states(n,Nlipid) = 2
        i1 = i1 + j
        i = i - j

        j = index(aline(i1+1:),'2')
      else
        if (j == i) exit ! Both must be equal to zero
        states(n,Nlipid) = 1
        i1 = i1 + i
        i = index(aline(i1+1:),'1')
      end if
    end if
  end do

  !print *, "the line: ", trim(adjustl(aline))

end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Now, just print to the screen

do n = 1, n_measure
  ! Just make this number larger than the potential
  ! number of lipids in the system
  write(6,FMT="(10000(1x,I1))") states(n,1:Nlipids)
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

deallocate(states)

return
end program sequenceOfStatesToStatesPerFrame
