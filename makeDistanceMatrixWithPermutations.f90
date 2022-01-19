program makeDistanceMatrixWithPermutations
implicit none

integer :: n_observations, n_models, n_states
double precision, allocatable :: observations(:,:,:), DM(:,:)

integer, dimension(24,4) :: permutations
double precision :: dist1
integer :: n_state
integer :: i, j, k

character(100) :: aline
integer :: iostate

!
! Example usage:
!

if (iargc() /= 3) then
  write(0,FMT="(A)") "ERROR: Need exactly three arguments:"
  write(0,FMT="(A)") "      (1) the number of observations for each model"
  write(0,FMT="(A)") "      (2) the number of hidden states for each model"
  write(0,FMT="(A)") "      (3) the number of models"
  stop
end if

call get_command_argument(1,aline)
read(aline,FMT=*,iostat=iostate) n_observations

if (iostate /= 0) then
  write(0,FMT="(A)") "ERROR: the first argument (n_observations) needs to be an integer"
  stop
end if

call get_command_argument(2,aline)
read(aline,FMT=*,iostat=iostate) n_states

if (iostate /= 0) then
  write(0,FMT="(A)") "ERROR: the second argument (n_states) needs to be an integer"
  stop
end if

!if (modulo(i,n_states) /= 0) then
!  write(0,FMT="(A)") "ERROR: the first argument (n_observations) must be"
!  write(0,FMT="(A)") "       divisible by the second argument (n_states)"
!  stop
!end if
!n_observations = i / n_states

call get_command_argument(3,aline)
read(aline,FMT=*,iostat=iostate) n_models

if (iostate /= 0) then
  write(0,FMT="(A)") "ERROR: the third argument (n_models) needs to be an integer"
  stop
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(observations(n_models,n_observations,n_states),DM(n_models,n_models))

! Read in the models
do i = 1, n_models
  read(5,FMT=*,iostat=iostate) observations(i,1:n_observations,1:n_states)
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

    DM(i,j) = 1.0e9
    select case (n_states)
      case(1)
        ! Here, one may want to define a custom metric
        dist1=sum((observations(i,1:n_observations,1)-&
                   observations(j,1:n_observations,1))**2)
        if (dist1 < DM(i,j)) DM(i,j) = dist1

      case(2)
        permutations(1,1:n_states) = (/ 1, 2 /)
        permutations(2,1:n_states) = (/ 2, 1 /)
        do k = 1, 2
          dist1 = 0.0
          do n_state = 1, n_states
            ! Here, one may want to define a custom metric
            dist1=dist1+&
             sum((observations(i,1:n_observations,permutations(1,n_state))-&
                  observations(j,1:n_observations,permutations(k,n_state)))**2)
          end do
          if (dist1 < DM(i,j)) DM(i,j) = dist1
        end do

      case(3)
        permutations(1,1:n_states) = (/ 1, 2, 3 /)
        permutations(2,1:n_states) = (/ 1, 3, 2 /)
        permutations(3,1:n_states) = (/ 3, 1, 2 /)
        permutations(4,1:n_states) = (/ 2, 1, 3 /)
        permutations(5,1:n_states) = (/ 2, 3, 1 /)
        permutations(6,1:n_states) = (/ 3, 2, 1 /)
        do k = 1, 6
          dist1 = 0.0
          do n_state = 1, n_states
            ! Here, one may want to define a custom metric
            dist1=dist1+&
             sum((observations(i,1:n_observations,permutations(1,n_state))-&
                  observations(j,1:n_observations,permutations(k,n_state)))**2)
          end do
          if (dist1 < DM(i,j)) DM(i,j) = dist1
        end do

      case default

        dist1=sum((observations(i,1:n_observations,1:n_states)-&
                   observations(j,1:n_observations,1:n_states))**2)
        if (dist1 < DM(i,j)) DM(i,j) = dist1

    end select

    DM(i,j) = sqrt(DM(i,j)/(n_states*n_observations))
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
end program makeDistanceMatrixWithPermutations
