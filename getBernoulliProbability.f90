!
! WARNING: cannot handle n > 30 (approximately)
!          where n = sum of occupancies
!

program getBernoulliProbability
implicit none

integer,parameter :: Nstates = 4 !3
double precision,dimension(Nstates) :: &
       Pstates = (/ 0.24978759558d0, 0.0841121495d0, 0.33305012744d0, 0.24978759558d0 /)
! 392 CHL1, 99 POPC, 392 POPE, 294 POPS
!      Pstates = (/ 0.1724137931d0, 0.32758620689d0, 0.50d0 /)
! 60 CHOL, 114 DPPC, 174 DOPC
!      Pstates = (/ 0.0947192d0, 0.328583d0, 0.246438d0, 0.33026d0 /)
! 113 POPC, 392 POPE, 294 POPS, 394 CHL1

integer,dimension(Nstates) :: occupancies
double precision :: probability

integer :: factorial

character(100) :: aline
integer :: i, j, k

if (iargc() /= Nstates) then
  write(6,FMT="(A,I2)") "Error. Wrong number "//&
          "of arguments; need exactly:", Nstates
  write(6,FMT="(A,I2)") "  --each is the number "//&
          "of lipids of type 1,2,...,",Nstates
  stop
end if

do i = 1, Nstates
  call getarg(i,aline)
  read(aline,FMT=*) occupancies(i)
end do

probability = 1.0d0
do i = 1, Nstates
  j = occupancies(i)
  probability = probability * &
     (Pstates(i)**j) / factorial(j)
end do

probability = probability * &
    factorial(sum(occupancies(1:Nstates)))

write(6,FMT=*) probability

return
end program getBernoulliProbability

integer function factorial(n)
implicit none
integer,intent(in) :: n
integer :: i

if (n < 0) then
  write(6,FMT="(A)") "Error. Negative occupancy "//&
          "for Bernoulli probability calculation"
  stop

else if (n == 0) then
  factorial = 1
  return

else
  factorial = n
  do i = n-1,1,-1
    factorial = factorial * i
  end do

end if

return
end function factorial
