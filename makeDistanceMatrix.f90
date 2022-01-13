program makeDistanceMatrix
implicit none

integer :: n_observations, n_models
double precision, allocatable :: observations(:,:), DM(:,:)
integer :: i, j, k

character(100) :: aline
integer :: iostate

!
! Example usage:
!     (see script)
!

if (iargc() /= 2) then
  write(0,FMT="(A)") "ERROR: Need exactly two arguments:"
  write(0,FMT="(A)") "      (1) the number of observations for each model"
  write(0,FMT="(A)") "      (2) the number of models"
  stop
end if

call get_command_argument(1,aline)
read(aline,FMT=*,iostat=iostate) n_observations

if (iostate /= 0) then
  write(0,FMT="(A)") "ERROR: the first argument (n_observations) needs to be an integer"
  stop
end if

call get_command_argument(2,aline)
read(aline,FMT=*,iostat=iostate) n_models

if (iostate /= 0) then
  write(0,FMT="(A)") "ERROR: the second argument (n_models) needs to be an integer"
  stop
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(observations(n_models,n_observations),DM(n_models,n_models))

! Read in the models
do i = 1, n_models
  read(5,FMT=*,iostat=iostate) observations(i,1:n_observations)
  if (iostate /= 0) then
    write(0,FMT="(A)") "ERROR: incorrect formatting in one model (must be numbers)"
    stop
  end if
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Calculate the distance matrix
DM = 0.0d0
do i = 1, n_models-1
  do j = i+1, n_models
    ! Here, one may want to define a custom metric
    DM(i,j) = sqrt(sum((observations(i,1:n_observations)-&
                        observations(j,1:n_observations))**2) / n_observations)
    DM(j,i) = DM(i,j)
  end do
end do

! Write this matrix to the screen
do i = 1, n_models
  write(6,FMT="(999(ES15.6,1x))") DM(i,1:n_models)
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

deallocate(observations,DM)

return
end program makeDistanceMatrix
