program ViterbiMostProbableSequence_KF
implicit none

integer :: n_states, n_emis, n_measure
integer :: i, j, k, t, n
integer, allocatable :: obvs(:), MPsequence(:)
double precision, allocatable :: aij(:,:), bik(:,:), pi_i(:)

double precision :: probability_precision = 1.0e-6

double precision :: score

character(100) :: aline
integer :: iostate

!
! Example usage:
! line=$(head -1 4lipid-system-analysis/seventotenthousand/tmpsequences); gfortran ViterbiMostProbableSequence_KF.f95; ./a.out 2 84 $(echo "$line" | awk '{print NF}') < <(tac 4lipid-system-analysis/hmm-10000length-training1/4lipid-system-analysis_tenthousandlength.out | sed '/best model/q' | tac | sed 1,3d | awk '{if (NR==1) {print $1; print $2} else {print $0}}' | sed 3,3d | sed 5,5d; echo "$line" | xargs -n1) | xargs
!

if (iargc() /= 3) then
  write(0,FMT="(A)") "ERROR: Need exactly three arguments:"
  write(0,FMT="(A)") "      (1) the number of hidden states"
  write(0,FMT="(A)") "      (2) the number of observable states"
  write(0,FMT="(A)") "      (3) the number of observations"
  stop
end if

call get_command_argument(1,aline)
read(aline,FMT=*,iostat=iostate) n_states

if (iostate /= 0) then
  write(0,FMT="(A)") "ERROR: the first argument (n_states) needs to be an integer"
  stop
end if

call get_command_argument(2,aline)
read(aline,FMT=*,iostat=iostate) n_emis

if (iostate /= 0) then
  write(0,FMT="(A)") "ERROR: the second argument (n_emis) needs to be an integer"
  stop
end if

call get_command_argument(3,aline)
read(aline,FMT=*,iostat=iostate) n_measure

if (iostate /= 0) then
  write(0,FMT="(A)") "ERROR: the third argument (n_measure) needs to be an integer"
  stop
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! The initial, transition, and emission probabilities:
!  +> aij transition probability between state i and j, 1 <= i,j <=
!  n_states
!  +> bik emission probability of state i, 1<= k <= n_emis
!  +> pi_i initial probability of state i
allocate(aij(n_states,n_states),bik(n_states,n_emis),&
         pi_i(n_states))
allocate(obvs(n_measure),MPsequence(n_measure))

! Read in the initial probabilities pi_i
do i = 1, n_states
  read(5,FMT=*,iostat=iostate) pi_i(i)
  if (iostate /= 0) then
    write(0,FMT="(A)") "ERROR: incorrect formatting in the initial probabilities"
    stop
  else if ((pi_i(i) < 0.0).or.(pi_i(i) > 1.0)) then
    write(0,FMT="(A)") "ERROR: at least one pi_i probability is out of bounds [0,1]"
    stop
  end if
end do

if (abs(sum(pi_i(1:n_states))-1.0) > probability_precision) then
  write(0,FMT="(A)") "ERROR: initial probabilities do not sum to 1"
  stop
end if

! Read in the state transition probabilities aij
do i = 1, n_states
  read(5,FMT=*,iostat=iostate) aij(i,1:n_states)
  if (iostate /= 0) then
    write(0,FMT="(A)") "ERROR: incorrect formatting in the transition probabilities"
    stop
  else if (any(aij(i,1:n_states) < 0.0).or.any(aij(i,1:n_states) > 1.0)) then
    write(0,FMT="(A)") "ERROR: at least one aij probability is out of bounds [0,1]"
    stop
  end if

  if (abs(sum(aij(i,1:n_states))-1.0) > probability_precision) then
    write(0,FMT="(A)") "ERROR: at least one set of transition probabilities do not sum to 1"
    stop
  end if
end do

! Read in the state emission probabilities bik
do i = 1, n_states
  read(5,FMT=*,iostat=iostate) bik(i,1:n_emis)
  if (iostate /= 0) then
    write(0,FMT="(A)") "ERROR: incorrect formatting in the emission probabilities"
    stop
  else if (any(bik(i,1:n_emis) < 0.0).or.any(bik(i,1:n_emis) > 1.0)) then
    write(0,FMT="(A)") "ERROR: at least one bik probability is out of bounds [0,1]"
    stop
  end if

  if (abs(sum(bik(i,1:n_emis))-1.0) > probability_precision) then
    write(0,FMT="(A)") "ERROR: at least one set of emission probabilities do not sum to 1"
    stop
  end if
end do

! Read in the observations obvs
do t = 1, n_measure
  read(5,FMT=*,iostat=iostate) obvs(t)
  if (iostate /= 0) then
    write(0,FMT="(A)") "ERROR: incorrect formatting in the observations"
    stop
  else if ((obvs(t) < 1).or.(obvs(t) > n_emis)) then
    write(0,FMT="(A)") "ERROR: at least one obervation is out of bounds [1,n_emis]"
    stop
  end if
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (.false.) then
call makeMostProbableSequence(obvs(1:n_measure),&
            aij(1:n_states,1:n_states),&
            bik(1:n_states,1:n_emis),pi_i(1:n_states),&
            n_emis,n_states,n_measure,&
            MPsequence(1:n_measure),score)
write(6,FMT="(A,ES32.20)") "Score of the most probable path P: ", score

else
call makeMostProbableSequence_SCALED(obvs(1:n_measure),&
            aij(1:n_states,1:n_states),&
            bik(1:n_states,1:n_emis),pi_i(1:n_states),&
            n_emis,n_states,n_measure,&
            MPsequence(1:n_measure),score)
write(6,FMT="(A,ES32.20)") "Score of the most probable path log(P): ", score
end if



do t = 1, n_measure
  write(6,FMT="(I3)") MPsequence(t)
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

deallocate(obvs,MPsequence)
deallocate(aij,bik,pi_i)

return
end program ViterbiMostProbableSequence_KF

subroutine makeMostProbableSequence(obvs,aij,bik,pi_i,&
                                    n_emis,n_states,n_measure,&
                                    MPsequence,score)
      implicit none
      integer, intent(in) :: n_states, n_emis, n_measure
      integer :: i, j, k, t
      double precision :: dum
      double precision :: delta_ti(n_measure,n_states)
      integer :: psi_ti(n_measure,n_states)
      integer,intent(out) :: MPsequence(n_measure)
      double precision,intent(in) :: aij(n_states,n_states),&
                          bik(n_states,n_emis),&
                          pi_i(n_states)
      double precision,intent(out) :: score
      integer, intent(in) :: obvs(n_measure)

! ---- Forward Procedure ----
! 1. Initialization
! Eq (32) in Rabiner, 1989
      do i = 1, n_states 
         delta_ti(1,i) = pi_i(i)*bik(i,obvs(1))
         psi_ti(1,i) = 0
      enddo
! 2. Recursion
! Eq (33) in Rabiner, 1989
      do t = 2, n_measure
         k = t - 1
         do j = 1, n_states
            psi_ti(t,j) = 0
            delta_ti(t,j) = 0.0
            do i = 1, n_states
               dum = delta_ti(k,i)*aij(i,j)
               if (dum > delta_ti(t,j)) then
                 psi_ti(t,j) = i
                 delta_ti(t,j) = dum
               end if
            enddo
            delta_ti(t,j) = delta_ti(t,j)*bik(j,obvs(t))
         enddo
      enddo
! 3. Termination
! Eq (34) in Rabiner, 1989
      score = 0.0
      do i = 1, n_states
         if (delta_ti(n_measure,i) > score) then
           score = delta_ti(n_measure,i)
           MPsequence(n_measure) = i
         end if
      enddo
! 4. Path backtracking
! Eq (35) in Rabiner, 1989
      do t = n_measure-1, 1, -1
         k = t + 1
         MPsequence(t) = psi_ti(k,MPsequence(k))
      end do

      return
end subroutine makeMostProbableSequence

subroutine makeMostProbableSequence_SCALED(obvs,aij,bik,pi_i,&
                                           n_emis,n_states,n_measure,&
                                           MPsequence,score)
      implicit none
      integer, intent(in) :: n_states, n_emis, n_measure
      integer :: i, j, k, t
      double precision :: dum
      double precision :: delta_ti(n_measure,n_states)
      integer :: psi_ti(n_measure,n_states)
      integer,intent(out) :: MPsequence(n_measure)
      double precision,intent(in) :: aij(n_states,n_states),&
                          bik(n_states,n_emis),&
                          pi_i(n_states)
      double precision :: aij_SCALED(n_states,n_states),&
                          bik_SCALED(n_states,n_emis),&
                          pi_i_SCALED(n_states)
      double precision,intent(out) :: score
      integer, intent(in) :: obvs(n_measure)

      do i = 1, n_states
        if (pi_i(i) < tiny(1.0d0)) then
          pi_i_scaled(i) = -huge(1.0d0)
        else
          pi_i_scaled(i) = log10(pi_i(i))
        endif
      end do

      do i = 1, n_states
      do j = 1, n_states
        if (aij(i,j) < tiny(1.0d0)) then
          aij_scaled(i,j) = -huge(1.0d0)
        else
          aij_scaled(i,j) = log10(aij(i,j))
        endif
      end do
      end do

      do i = 1, n_states
      do k = 1, n_emis
        if (bik(i,k) < tiny(1.0d0)) then
          bik_scaled(i,k) = -huge(1.0d0)
        else
          bik_scaled(i,k) = log10(bik(i,k))
        endif
      end do
      end do

! ---- Forward Procedure ----
! 1. Initialization
! Eq (32 and 105a) in Rabiner, 1989
      do i = 1, n_states 
         delta_ti(1,i) = pi_i_scaled(i) + bik_scaled(i,obvs(1))
         psi_ti(1,i) = 0
      enddo
! 2. Recursion
! Eq (33 and 105b) in Rabiner, 1989
      do t = 2, n_measure
         k = t - 1
         do j = 1, n_states
            psi_ti(t,j) = 0
            delta_ti(t,j) = -huge(1.0d0)
            do i = 1, n_states
               dum = delta_ti(k,i) + aij_scaled(i,j)
! kazuumi's debugging:
! if (t > 9000) print *, t, dum, delta_ti(t,j), bik_scaled(j,obvs(t))
               if (dum > delta_ti(t,j)) then
                 psi_ti(t,j) = i
                 delta_ti(t,j) = dum
               end if
            enddo
            delta_ti(t,j) = delta_ti(t,j) + bik_scaled(j,obvs(t))
         enddo
      enddo
! 3. Termination
! Eq (34 and 105c) in Rabiner, 1989
      score = -huge(1.0d0)
      do i = 1, n_states
         if (delta_ti(n_measure,i) > score) then
           score = delta_ti(n_measure,i)
           MPsequence(n_measure) = i
         end if
      enddo
! 4. Path backtracking
! Eq (35) in Rabiner, 1989
      do t = n_measure-1, 1, -1
         k = t + 1
         MPsequence(t) = psi_ti(k,MPsequence(k))
      end do

      return
end subroutine makeMostProbableSequence_SCALED
