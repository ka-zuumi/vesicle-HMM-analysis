program HMMsolverMultipleSequences_KF
implicit none

integer,parameter :: n_measure_max = 50000
integer :: n_states, n_emis, n_measure, n_sequence, n_iter, q_viterbi, iter
integer :: i, j, k, t, n, dum_array(1)
integer, allocatable :: n_measures(:)
integer, allocatable :: obvs(:,:), phi_ti(:,:), state_track(:)
double precision, allocatable :: aij(:,:), bik(:,:), pi_i(:)
double precision, allocatable :: aij_numerator(:,:,:), &
                                 aij_denominator(:,:,:), &
                                 bik_numerator(:,:,:), &
                                 bik_denominator(:,:,:), &
                                 pi_numerator(:,:)
double precision, allocatable :: aij_p(:,:), bik_p(:,:), pi_i_p(:)
double precision, allocatable :: best_aij(:,:), best_bik(:,:), best_pi_i(:)
double precision, allocatable :: forw(:,:), bckw(:,:), &
                                 prob_ti(:,:)
double precision, allocatable :: delta_ti(:,:), delta_aij(:), &
                                 sai_tij(:,:,:)
double precision :: factor, prob, dum, prob_f, prob_b, p_viterbi
double precision :: aij_diff, bik_diff, pi_diff
double precision :: dum1, dum2

integer :: nO, nP, nN, nC

integer :: n_initial_guesses, initial_guess
double precision, allocatable :: scores(:), probabilities(:)
double precision :: probabilities_min
double precision :: max_score

character(100) :: aline
integer :: iostate

! The number of iterations for the algorithm (to get the
! most probable model)
n_iter = 100

! The number of initial conditions to test to potentially
! get a better model
n_initial_guesses = 500

if (iargc() /= 3) then
  write(0,FMT="(A)") "ERROR: Need exactly three arguments:"
  write(0,FMT="(A)") "      (1) the number of hidden states"
  write(0,FMT="(A)") "      (2) the number of observable states"
  write(0,FMT="(A)") "      (3) the number of sequences"
  stop
end if

call get_command_argument(1,aline)
read(aline,FMT=*,iostat=iostate) n_states

if (iostate /= 0) then
  write(0,FMT="(A)") "ERROR: the first argument needs to be an integer"
  stop
end if

call get_command_argument(2,aline)
read(aline,FMT=*,iostat=iostate) n_emis

if (iostate /= 0) then
  write(0,FMT="(A)") "ERROR: the second argument needs to be an integer"
  stop
end if

call get_command_argument(3,aline)
read(aline,FMT=*,iostat=iostate) n_sequence

if (iostate /= 0) then
  write(0,FMT="(A)") "ERROR: the third argument needs to be an integer"
  stop
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(obvs(n_sequence,n_measure_max),&
         n_measures(n_sequence))
!allocate(obvs(n_measure),prob_ti(n_measure,n_states),&
!        sai_tij(n_measure,n_states,n_states))

! Read in sequences, with each sequence
! starting with the total number of observables
! in it
do n = 1, n_sequence
  read(5,FMT=*,iostat=iostate) n_measures(n)
  if (iostate /= 0) then
    write(0,FMT="(A)") "ERROR: incorrect formatting in at least one sequence length"
    stop
  else if ((n_measures(n) < 1).or.(n_measures(n) > n_measure_max)) then
    write(0,FMT="(A)") "ERROR: at least one sequence length is out-of-bounds"
    stop
  end if

  ! After retrieving the sequence length,
  ! continue reading in observables for that
  ! sequence
  do i = 1, n_measures(n)
    read(5,FMT=*,iostat=iostate) obvs(n,i)
    if (iostate /= 0) then
      write(0,FMT="(A)") "ERROR: incorrect formatting in observables"
      stop
    else if ((obvs(n,i) < 1).or.(obvs(n,i) > n_emis)) then
      write(0,FMT="(A)") "ERROR: observed state is out-of-bounds"
      stop
    end if
  end do
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Initial estimate of the initial, transition, and emission possibility
! uniformly
!  +> aij transition probability between state i and j, 1 <= i,j <=
!  n_states
!  +> bik emission probability of state i, 1<= k <= n_emis
!  +> pi_i initial probability of state i
allocate(aij_numerator(n_sequence,n_states,n_states),&
         aij_denominator(n_sequence,n_states,n_states),&
         bik_numerator(n_sequence,n_states,n_emis),&
         bik_denominator(n_sequence,n_states,n_emis),&
         pi_numerator(n_sequence,n_states))
allocate(aij(n_states,n_states),bik(n_states,n_emis),&
         pi_i(n_states))
allocate(aij_p(n_states,n_states),bik_p(n_states,n_emis),&
         pi_i_p(n_states))
allocate(best_aij(n_states,n_states),best_bik(n_states,n_emis),&
         best_pi_i(n_states))

n_measure = maxval(n_measures(1:n_sequence))
allocate(prob_ti(n_measure,n_states),&
         sai_tij(n_measure,n_states,n_states))

prob_ti = 0.0; sai_tij = 0.0



call system_clock(i)
call random_seed(i)

! We will keep persistent "scores" to keep track
! of different models obtained from different
! initial conditions
max_score = -HUGE(1.0d0)
allocate(scores(n_initial_guesses),&
         probabilities(n_sequence))

! A giant do loop iterates over each initial
! condition
do initial_guess = 1, n_initial_guesses

write(6,FMT="(A,I3)") "Starting iteration No. ", initial_guess

! There are many ways to same the initial
! conditions
select case (2)

! This is only good for tests
case(1)

  ! A uniform probability
  do i = 1, n_states
    aij(i,1:n_states) = 1.0 / n_states
    bik(i,1:n_emis) = 1.0 / n_emis
  end do
  pi_i(1:n_states) = 1.0 / n_states

! This is good if you have no idea and want to
! sample a lot of different initial conditions
case(2)
  ! A random probability

  do i = 1, n_states
    do j = 1, n_states
      call random_number(aij(i,j))
    end do
    aij(i,1:n_states) = aij(i,1:n_states) / sum(aij(i,1:n_states))
  
    do j = 1, n_emis
      call random_number(bik(i,j))
    end do
    bik(i,1:n_emis) = bik(i,1:n_emis) / sum(bik(i,1:n_emis))
  
    call random_number(pi_i(i))
  end do
  pi_i(1:n_states) = pi_i(1:n_states) / sum(pi_i(1:n_states))

! This is good if we want to 'guide' the algorithm
! into picking a specific model that it may not have
! guessed
case(3)
  ! A Gaussian centered at one specific state + noise
  ! (one for each hidden state)

  ! The noise:
  do i = 1, n_states
    do j = 1, n_states
      call random_number(aij(i,j))
    end do

    do j = 1, n_emis
      call random_number(bik(i,j))
    end do

    call random_number(pi_i(i))
  end do
  pi_i(1:n_states) = pi_i(1:n_states) / sum(pi_i(1:n_states))


  ! For a 3-lipid system:
  i = 0
  do nP = 0, 6
  do nN = 0, 6
    nC = 6 - (nP+nN)
    if (nC < 0) cycle
    i = i + 1

    ! State 1: (2,3,1)
    aij(1,i) = 0.5*aij(1,i) + 1.0/(1 + ((2-nP)**2) + ((3-nN)**2) + ((1-nC)**2))
    ! State 2: (0,1,5)
    aij(2,i) = 0.5*aij(2,i) + 1.0/(1 + ((0-nP)**2) + ((1-nN)**2) + ((5-nC)**2))
  end do
  end do


  do i = 1, n_states
    aij(i,1:n_states) = aij(i,1:n_states) / sum(aij(i,1:n_states))
    bik(i,1:n_emis) = bik(i,1:n_emis) / sum(bik(i,1:n_emis))
  end do

end select


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     write(6,FMT="(A4,1x,3(A12))") "iter", "delta_pi",&
!                "delta_aij", "delta_bik"
!     write(6,FMT="(10x,1x,11A)") "log(P)"

! This is the main loop for the HMM solver;
! this is an iterative algorithm
do iter = 1, n_iter
  aij_p = aij
  bik_p = bik
  pi_i_p = pi_i

  ! Iterate over each sequence of observables
  do n = 1, n_sequence
    ! ---- Evaluate prob_ti and sai_tij with for/backward procedure ----
    call alpha_beta_gamma_sai_SCALED(obvs(n,1:n_measures(n)),&
              aij,bik,pi_i,n_emis,n_states,n_measures(n),&
              prob_ti(1:n_measures(n),1:n_states),&
              sai_tij(1:n_measures(n),1:n_states,1:n_states),&
              probabilities(n))
    !write(6,FMT=*) "i, probability(i): ", n, probabilities(n)

    ! ---- Reevaluate aij, bik, and pi_i with prob_ti and sai_tij ----
    ! Eq (40a) in Rabiner, 1989
    do i = 1, n_states
      !pi_i(i) = prob_ti(1,i)
      pi_numerator(n,i) = prob_ti(1,i)
    enddo
    ! Eq (40b) in Rabiner, 1989
    do i = 1, n_states
      do j = 1, n_states
        dum1 = 0.0
        dum2 = 0.0
        do t = 1, n_measures(n)-1
          dum1 = dum1 + sai_tij(t,i,j)
          dum2 = dum2 + prob_ti(t,i)
        enddo
        !aij(i,j) = dum1/dum2
        aij_numerator(n,i,j) = dum1
        aij_denominator(n,i,j) = dum2
      enddo
    enddo
    ! Eq (40c) in Rabiner, 1989
    do i = 1, n_states
      do k = 1, n_emis
        dum1 = 0.0
        dum2 = 0.0
        do t = 1, n_measures(n)
          if (obvs(n,t) .eq. k) then 
            dum1 = dum1 + prob_ti(t,i)
          endif
          dum2 = dum2 + prob_ti(t,i)
        enddo
        !bik(i,k) = dum1/dum2
        bik_numerator(n,i,k) = dum1
        bik_denominator(n,i,k) = dum2
      enddo
    enddo
  end do
  scores(initial_guess) = sum(probabilities(1:n_sequence))
  !write(6,FMT=*) "score: ", scores(initial_guess)

  ! If this set of initial conditions leads to a
  ! model with a better score, then keep it and
  ! print to the standard output
  if (scores(initial_guess) > max_score) then
    max_score = scores(initial_guess)
    write(6,FMT="(A,ES32.20)") "new best model with log(P): ", max_score
    write(6,FMT="(A)") ""
    write(6,FMT="(A)") "pi_i:"
    write(6,FMT="(999(1x,F11.8))")pi_i
    write(6,FMT="(A)") "aij:"
    do i = 1, n_states
      write(6,FMT="(999(1x,F11.8))") aij(i,1:n_states)
    enddo
    write(6,FMT="(A)") "bik:"
    do i = 1, n_states
      write(6,FMT="(999(1x,F11.8))") bik(i,1:n_emis)
    enddo

    best_pi_i = pi_i
    best_aij = aij
    best_bik = bik
  end if


  ! To get the reestimated parameters we must
  ! weight the sequences by their inverse probabilities;
  ! here I divide out by the smallest probability
  ! (which has the max score)
  !
  ! To get a feel on whether this precision is okay:
  ! https://stackoverflow.com/questions/872544/what-range-of-numbers-can-be-represented-in-a-16-32-and-64-bit-ieee-754-syste
  !
  probabilities_min = minval(probabilities(1:n_sequence))
  do n = 1, n_sequence
    ! Probablistic weighting used by Rabiner, 1989
    !probabilities(n) = 10.0**(probabilities_min-probabilities(n))
    ! Unit weighting used by  Lovel 2002
    probabilities(n) = 1.0d0

    ! Now use the modified reestimation formula to get the
    ! next set of parameters
    ! Eq (109) in Rabiner, 1989
    aij_numerator(n,1:n_states,1:n_states) = &
          probabilities(n) * aij_numerator(n,1:n_states,1:n_states)
    aij_denominator(n,1:n_states,1:n_states) = &
          probabilities(n) * aij_denominator(n,1:n_states,1:n_states)
    bik_numerator(n,1:n_states,1:n_emis) = &
          probabilities(n) * bik_numerator(n,1:n_states,1:n_emis)
    bik_denominator(n,1:n_states,1:n_emis) = &
          probabilities(n) * bik_denominator(n,1:n_states,1:n_emis)

    pi_numerator(n,1:n_states) = &
          probabilities(n) * pi_numerator(n,1:n_states)
  end do

  ! ---- Reevaluate aij, bik, and pi_i with prob_ti and sai_tij ----
  ! Eq (40a) in Rabiner, 1989
  do i = 1, n_states
    !pi_i(i) = prob_ti(1,i)
    pi_i(i) = sum(pi_numerator(1:n_sequence,i)) / sum(probabilities(1:n_sequence))
    ! Set the first state to be the one which has probability
    ! pi = 1 (?)
    ! Remark near Eq (110) in Rabiner, 1989
    !pi_i(i) = 0.0                !pi_numerator(n,i)
    !if (i == 1) pi_i(i) = 1.0    !pi_numerator(n,i)
    !write(6,FMT="(999(F4.2,1x))") pi_numerator(1:n_sequence,i)

    ! Remark by Kazuumi:
    ! I still have no idea why this works; if this is supposed
    ! to converge for all sequences, I could believe this but
    ! I cannot prove this
  enddo
  do i = 1, n_states
    ! Eq (109) in Rabiner, 1989
    do j = 1, n_states
      aij(i,j) = sum(aij_numerator(1:n_sequence,i,j)) / &
                 sum(aij_denominator(1:n_sequence,i,j))
    enddo

    ! Eq (110) in Rabiner, 1989
    do k = 1, n_emis
      bik(i,k) = sum(bik_numerator(1:n_sequence,i,k)) / &
                 sum(bik_denominator(1:n_sequence,i,k))
    enddo
  enddo

  ! ---- Check the difference between the updated aij, bik, and pi_i with
  ! the original ones ----
  pi_diff = sum((pi_i - pi_i_p)**2.0)
  pi_diff = sqrt(pi_diff)/n_states
  dum = 0.0
  do i = 1, n_states
    do j = 1, n_states
      dum = dum + (aij(i,j)-aij_p(i,j))**2.0
    enddo
  enddo
  aij_diff = dum/(n_states*n_states)
  dum = 0.0
  do i = 1, n_states
    do k = 1, n_emis
      dum = dum + (bik(i,k)-bik_p(i,k))**2.0
    enddo
  enddo
  bik_diff = dum/(n_states*n_emis)
  !     write(6,FMT="(I3,2x,3(1x,F11.8))") iter, pi_diff,&
  !              aij_diff, bik_diff

enddo

write(6,FMT="(A,ES32.20)") "new best model with log(P): ", scores(initial_guess)
write(6,FMT="(A)") ""
write(6,FMT="(A)") "pi_i:"
write(6,FMT="(999(1x,F11.8))")pi_i
write(6,FMT="(A)") "aij:"
do i = 1, n_states
  write(6,FMT="(999(1x,F11.8))") aij(i,1:n_states)
enddo
write(6,FMT="(A)") "bik:"
do i = 1, n_states
  write(6,FMT="(999(1x,F11.8))") bik(i,1:n_emis)
enddo

end do

write(6,FMT="(A,ES32.20)") "overall best model with log(P): ", max_score
write(6,FMT="(A)") ""
write(6,FMT="(A)") "pi_i:"
write(6,FMT="(999(1x,F11.8))") best_pi_i
write(6,FMT="(A)") "aij:"
do i = 1, n_states
  write(6,FMT="(999(1x,F11.8))") best_aij(i,1:n_states)
enddo
write(6,FMT="(A)") "bik:"
do i = 1, n_states
  write(6,FMT="(999(1x,F11.8))") best_bik(i,1:n_emis)
enddo

return
end program HMMsolverMultipleSequences_KF

subroutine alpha_beta_gamma_sai(obvs,aij,bik,pi_i,&
                                n_emis,n_states,n_measure,&
                                prob_ti,sai_tij,score)
      implicit none
      integer, intent(in) :: n_states, n_emis, n_measure
      integer :: i, j, k, t, q_viterbi, dum_array(1)
      double precision :: factor, prob, dum, prob_f, prob_b, p_viterbi
      double precision,intent(out) :: prob_ti(n_measure,n_states),&
                          sai_tij(n_measure,n_states,n_states)
      double precision,intent(in) :: aij(n_states,n_states),bik(n_states,n_emis),&
                          pi_i(n_states)
      double precision :: forw(n_measure,n_states),&
                          bckw(n_measure,n_states)
      double precision :: prob_ti_tmp(n_measure,n_states)
      double precision,intent(out) :: score
      integer :: obvs(n_measure)

! ---- Forward Procedure ----
! 1. Initialization
! Eq (19) in Rabiner, 1989
      do i = 1, n_states 
         forw(1,i) = pi_i(i)*bik(i,obvs(1))
      enddo
! 2. Induction
! Eq (20) in Rabiner, 1989
      do t = 2, n_measure
         k = t - 1
         do j = 1, n_states
            dum = 0.0
            do i = 1, n_states
               dum = dum + forw(k,i)*aij(i,j)
            enddo
            forw(t,j) = dum*bik(j,obvs(t))
         enddo
      enddo
! 3. Termination
! Eq (21) in Rabiner, 1989
      prob_f = 0.0
      do i = 1, n_states
         prob_f = prob_f + forw(n_measure,i)
      enddo

! ---- Backward Procedure ----
! 1. Initialization      
! Eq (24) in Rabiner, 1989
      do i = 1, n_states
         bckw(n_measure,i) = 1
      enddo 
! 2. Induction
! Eq (25) in Rabiner, 1989
      do t = n_measure-1, 1, -1
         k = t + 1
         do i = 1, n_states 
            dum = 0.0
            do j = 1, n_states
               dum = dum + aij(i,j)*bik(j,obvs(k))*bckw(k,j)
            enddo
            bckw(t,i) = dum
         enddo
      enddo
! 3. Termination
! Eq (23) in Rabiner, 1989 (by definition, although not obvious)
      prob_b = 0.0
      do i = 1, n_states
         prob_b = prob_b + bckw(1,i)*pi_i(i)*bik(i,obvs(1))
      enddo

! ---- Verify the forwrad prob and backward prob equals ----
      dum = sqrt((prob_f - prob_b)**2.0)
      if (dum .ge. 0.00001) then
         write(0,*) "Forward probability is not "//&
                    "equal to backward probability. Abort"
         stop
      endif

      ! The forward and backward probabilities
      ! should be the same, so just print one;
      ! take a logarithm to be consistent
      score = log10(prob_f)
      !write(6,FMT="(10x,1x,F11.4)") score

! ---- Compute gamma (prob_ti) ----
! prob_ti is the probability of being in state i at time t
! Eq (27) in Rabiner, 1989
      do t = 1, n_measure
         do i = 1, n_states
            dum = 0.0
           do j = 1, n_states
               dum = dum + forw(t,j)*bckw(t,j)
           enddo
           prob_ti(t,i) = forw(t,i)*bckw(t,i)/dum
         enddo
      enddo

! ---- Verify the legitmency of prob_ti ----
! sum of prob_ti at any given time should be 1
! Eq (28) in Rabiner, 1989
      do t = 1, n_measure
         dum = 0.0
         do i = 1, n_states
            dum = dum + prob_ti(t,i)
         enddo
         dum = sqrt((dum - 1.0)**2.0)
         if (dum .ge. 0.00001) then 
            write(0,*) "sum of state at time", t, "is not 1. Abort"
            stop
         endif
      enddo

! ---- Compute sai_tij ---
! sai_tij is the probability of being in state i at time t, and state j
! at time t + 1
! Eq (37) in Rabiner, 1989
      do t = 1, n_measure-1
         k = t + 1
         do i = 1, n_states
            do j = 1, n_states
               dum = forw(t,i)*aij(i,j)*bik(j,obvs(k))*bckw(k,j)
               sai_tij(t,i,j) = dum/prob_b 
            enddo
         enddo
      enddo
! this is in place to take care of sai_tij(n_measure,i,:) which won't be
! evaluated in the previous loop
      do i = 1, n_states
         sai_tij(n_measure,i,:) = prob_ti(n_measure,i)/n_states
      enddo

! ---- Compute gamma (prob_ti) via sai_tij and Verify ----
! Eq (38) in Rabiner, 1989
      do t = 1, n_measure
         do i = 1, n_states
            dum = 0.0
            do j = 1, n_states
               dum = dum + sai_tij(t,i,j)
            enddo
            prob_ti_tmp(t,i) = dum
            dum = prob_ti_tmp(t,i) - prob_ti(t,i) 
            dum = sqrt(dum **2.0)
            if (dum .ge. 0.00001) then
               write(0,*) prob_ti_tmp(t,i), prob_ti(t,i), t, i
               write(0,*) "sai_tij", sai_tij(t,i,:) 
               write(0,*) "sum of state from Baum at time", t, &
                          "is not 1. Abort"
               stop
            endif
         enddo
      enddo
end subroutine alpha_beta_gamma_sai

subroutine alpha_beta_gamma_sai_SCALED(obvs,aij,bik,pi_i,&
                                n_emis,n_states,n_measure,&
                                prob_ti,sai_tij,score)
      implicit none
      integer, intent(in) :: n_states, n_emis, n_measure
      integer :: i, j, k, t, q_viterbi, dum_array(1)
      double precision :: factor, prob, dum, prob_f, prob_b, p_viterbi
      double precision,intent(out) :: prob_ti(n_measure,n_states),&
                          sai_tij(n_measure,n_states,n_states)
      double precision,intent(in) :: aij(n_states,n_states),bik(n_states,n_emis),&
                          pi_i(n_states)
      double precision :: forw(n_measure,n_states),&
                          bckw(n_measure,n_states)
      double precision :: prob_ti_tmp(n_measure,n_states)
      double precision :: scaling(n_measure)
      double precision,intent(out) :: score
      integer :: obvs(n_measure)

      double precision :: zero_threshold = 1.0d-20

! ---- Forward Procedure ----
! 1. Initialization
! Eq (19) in Rabiner, 1989
      do i = 1, n_states 
         forw(1,i) = pi_i(i)*bik(i,obvs(1))
      enddo
! 2. Induction
! Eq (20) in Rabiner, 1989
      do t = 2, n_measure
         k = t - 1

         ! 2.5. Scaling
         ! Eq (91) in Rabiner, 1989
         scaling(k) = 1.0d0 / sum(forw(k,1:n_states))
         !write(6,FMT="(10x,999(1x,ES11.4))") forw(k,1:n_states)
         !write(6,FMT="(10x,999(1x,ES11.4))") scaling(k)
         !if (isnan(scaling(k))) stop
         !if ((scaling(k) > HUGE(scaling(k)))) stop
         if (isnan(scaling(k)).or.(scaling(k) > HUGE(scaling(k)))) then
           write(6,FMT="(10x,A,999(1x,ES11.4))") "pi_i: ", pi_i(1)
           write(6,FMT="(10x,A,999(1x,ES11.4))") "      ", pi_i(2)
           write(6,FMT="(10x,A,999(1x,ES11.4))") " bik: ", bik(1,1:n_emis)
           write(6,FMT="(10x,A,999(1x,ES11.4))") "      ", bik(2,1:n_emis)
           write(6,FMT="(10x,A,999(1x,I2))") "obvs: ", obvs(1:k)
           stop
         end if
         forw(k,1:n_states) = forw(k,1:n_states) * scaling(k)

         do j = 1, n_states
            dum = 0.0
            do i = 1, n_states
               dum = dum + forw(k,i)*aij(i,j)
            enddo
            forw(t,j) = dum*bik(j,obvs(t))
         enddo
      enddo

      k = n_measure
      scaling(k) = 1.0d0 / sum(forw(k,1:n_states))
      forw(k,1:n_states) = forw(k,1:n_states) * scaling(k)

! 3. Termination
! Eq (21) in Rabiner, 1989
! NOTE: this is now a SCALED probability (not the true probability)
      prob_f = 0.0
      do i = 1, n_states
         prob_f = prob_f + forw(n_measure,i)
      enddo

! ---- Backward Procedure ----
! 1. Initialization      
! Eq (24) in Rabiner, 1989
      do i = 1, n_states
         bckw(n_measure,i) = 1
      enddo 
! 2. Induction
! Eq (25) in Rabiner, 1989
      do t = n_measure-1, 1, -1
         k = t + 1

         ! 2.5. Scaling
         ! (Use the same constants as in the forward probabilities)
         ! Eq (94) in Rabiner, 1989
         bckw(k,1:n_states) = bckw(k,1:n_states) * scaling(k)

         do i = 1, n_states 
            dum = 0.0
            do j = 1, n_states
               dum = dum + aij(i,j)*bik(j,obvs(k))*bckw(k,j)
            enddo
            bckw(t,i) = dum
         enddo
      enddo

      k = 1
      bckw(k,1:n_states) = bckw(k,1:n_states) * scaling(k)

! 3. Termination
! Eq (23) in Rabiner, 1989 (by definition, although not obvious)
! NOTE: this is now a SCALED probability (not the true probability)
      prob_b = 0.0
      do i = 1, n_states
         prob_b = prob_b + bckw(1,i)*pi_i(i)*bik(i,obvs(1))
      enddo

! ---- Verify the forwrad prob and backward prob equals ----
      dum = sqrt((prob_f - prob_b)**2.0)
      if (dum .ge. 0.00001) then
         write(0,*) "Forward probability is not "//&
                    "equal to backward probability. Abort"
         stop
      endif

      ! The probability may be too small to be
      ! expressed on the machine but their logarithms
      ! can be; that is what is printed in the following
      ! Eq (103) in Rabiner, 1989
      score = -sum(log10(scaling(1:n_measure)))
      !write(6,FMT="(10x,999(1x,F11.4))") scaling(1:n_measure)
      !write(6,FMT="(10x,1x,F11.4)") score

! ---- Compute gamma (prob_ti) ----
! prob_ti is the probability of being in state i at time t
! Eq (27) in Rabiner, 1989
      do t = 1, n_measure
         do i = 1, n_states
           dum = 0.0
           do j = 1, n_states
               dum = dum + forw(t,j)*bckw(t,j)
           enddo
           !write(6,FMT="(A,ES11.3)") "dum: ", dum
           !if (dum < zero_threshold) write(6,FMT="(A)") "prob_ti will overflow"
           !if (dum < zero_threshold) dum = zero_threshold
           prob_ti(t,i) = forw(t,i)*bckw(t,i)/dum
           !if (isnan(prob_ti(t,i))) stop
           !if (prob_ti(t,i) > HUGE(prob_ti(t,i))) stop
         enddo
      enddo

! ---- Verify the legitmency of prob_ti ----
! sum of prob_ti at any given time should be 1
! Eq (28) in Rabiner, 1989
      do t = 1, n_measure
         dum = 0.0
         do i = 1, n_states
            dum = dum + prob_ti(t,i)
         enddo
         dum = sqrt((dum - 1.0)**2.0)
         if (dum .ge. 0.00001) then 
            write(0,*) "sum of state at time", t, "is not 1. Abort"
            stop
         endif
      enddo

! ---- Compute sai_tij ---
! sai_tij is the probability of being in state i at time t, and state j
! at time t + 1
! Eq (37) in Rabiner, 1989
      do t = 1, n_measure-1
         k = t + 1
         do i = 1, n_states
            do j = 1, n_states
               dum = forw(t,i)*aij(i,j)*bik(j,obvs(k))*bckw(k,j)
               sai_tij(t,i,j) = dum/prob_b 
            enddo
         enddo
      enddo
! this is in place to take care of sai_tij(n_measure,i,:) which won't be
! evaluated in the previous loop
      do i = 1, n_states
         sai_tij(n_measure,i,:) = prob_ti(n_measure,i)/n_states
      enddo

! ---- Compute gamma (prob_ti) via sai_tij and Verify ----
! Eq (38) in Rabiner, 1989
      do t = 1, n_measure
         do i = 1, n_states
            dum = 0.0
            do j = 1, n_states
               dum = dum + sai_tij(t,i,j)
            enddo
            prob_ti_tmp(t,i) = dum
            dum = prob_ti_tmp(t,i) - prob_ti(t,i) 
            dum = sqrt(dum **2.0)
            if (dum .ge. 0.00001) then
               write(0,*) prob_ti_tmp(t,i), prob_ti(t,i), t, i
               write(0,*) "sai_tij", sai_tij(t,i,:) 
               write(0,*) "sum of state from Baum at time", t, &
                          "is not 1. Abort"
               stop
            endif
         enddo
      enddo
end subroutine alpha_beta_gamma_sai_SCALED
