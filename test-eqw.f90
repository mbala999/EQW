      
      INCLUDE 'eqw.f90'

      program main
          use constants
          use eqw
          implicit none

          integer :: i,j,k
          real(8) :: tauw,qw,ue,te,pe,re,he

          re = 0.1d0
          ue = 1400.0d0
          te = 500.0d0
          pe = te*rgas*re
          he = 0.001d0

          call solve_eqw(ue,te,pe,he,tauw,qw)

          print*, tauw,qw

          return
      end program main
