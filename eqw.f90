      !==============================================================!
      !              Equilibrium Wall Model (EQWM) Module            !
      !==============================================================!
      !            Balachandra R. Mettu <brmettu@ncsu.edu>           !
      !            Pramod Subbareddy    <pksubbar@ncsu.edu>          !
      !                    Last modified: 12 Dec 2020                !
      !==============================================================!
      ! Cite the below the references when using this code:          !
      ! - Kawai & Larsson (2012), Physics of Fluids                  !
      ! - Mettu & Subbareddy (2021), AIAA Journal                    !
      !==============================================================!

      module constants

          ! Boundary conditions and free-stream values
          integer, parameter :: iwall = 2
          real(8), parameter :: uin   = 500.0d0
          real(8), parameter :: tin   = 300.0d0
          real(8), parameter :: twall = 300.0d0

          ! Gas constants
          real(8), parameter :: ugcon   = 8.3144622d3 
          real(8), parameter :: gamgas  = 1.4d0 
          real(8), parameter :: rgas    = ugcon/28.97d0  
          real(8), parameter :: cvgas   = rgas/(gamgas-1.0d0)  
          real(8), parameter :: cvgasi  = 1.0d0/cvgas
          real(8), parameter :: prandtl = 0.7d0  
          real(8), parameter :: prli    = 1.0d0/prandtl  
          real(8), parameter :: prti    = 1.0d0/0.9d0
          real(8), parameter :: cpav    = cvgas+rgas    
          real(8), parameter :: gamm1   = gamgas-1.0d0
          real(8), parameter :: gamm1i  = 1.0d0/gamm1

          real(8), parameter :: thd   = 1.0d0/3.0d0
          real(8), parameter :: tthd  = 2.0d0/3.0d0

      contains

          !---------------------------------------------------------------
          ! Compute gas viscosity - Sutherland's Law
          !---------------------------------------------------------------
          function fsuth(t11)
              implicit none

              real(8) :: fsuth,t11,tt11
              real(8) :: ub11,muw11
              real(8) :: ssuth,tsuth
              real(8), parameter :: asuth = 1.458d-6 
              real(8), parameter :: bsuth = 110.3d0

              fsuth = asuth*sqrt(t11*t11*t11)/(t11+bsuth)

              return
          end function fsuth

      end module constants

      module eqw
          use constants

          integer, parameter :: idebug   = 1
          integer, parameter :: iturb_wm = 1
          integer, parameter :: kmax_wm  = 99
          integer, parameter :: nwm      = 99
          real(8), parameter :: res_wm   = 1.0d-6

          ! Constants mixing length model
          real(8), parameter :: kt_wm   = 0.41d0
          real(8), parameter :: apl1_wm = 17.0d0

          ! Grid stretching factor
          real(8), parameter :: stch_wm = 2.0d0

      contains

          !---------------------------------------------------------------------
          ! Equilibrium model
          !---------------------------------------------------------------------
          ! INPUT: ue,te,pe,dwe
          ! wall-parallel velocity, temperature, pressure, wall distance
          !---------------------------------------------------------------------
          ! OUTOUT: tauw,qw
          !---------------------------------------------------------------------
          subroutine solve_eqw(ue,te,pe,he,tauw,qw)
              use constants
              implicit none

              real(8), intent(in)  :: ue,te,pe,he
              real(8), intent(out) :: tauw,qw

              integer :: k,ncr_wm

              real(8) :: re,tw,rw,muw
              real(8) :: tauwm,rtauw,vdamp,ystar

              real(8) :: urms,trms
              real(8) :: uwmf,twmf,rwmf,muwmf,mutwmf,dywm,duwm

              real(8) :: a_low(nwm),b_dia(nwm),c_upp(nwm),d_rhs(nwm)
              real(8) :: xuwm(nwm),xtwm(nwm)

              real(8) :: ycwml(0:nwm+1),yfwml(0:nwm)

              real(8) :: mutwml(0:nwm+1)
              real(8) :: uwml(0:nwm+1),twml(0:nwm+1),rwml(0:nwm+1)
              real(8) :: uwml_n(0:nwm+1),twml_n(0:nwm+1)

              real(8) :: dnk,sf(0:nwm+1),sc(0:nwm+1)

              !
              !- Generate 1D stetched grid
              !
              !- Note that you can do this as a pre-processing step instead
              !  of generating a 1D grid at every iteration.
              !
              ! Face nodes
              do k = nwm+1,0,-1
              dnk   = dble(nwm+1-k)/dble(nwm+1)
              sf(k) = 1.0d0 - tanh(stch_wm*dnk)/tanh(stch_wm)
              enddo

              ! Cell centers
              do k = 1,nwm+1
              sc(k) = (sf(k)+sf(k-1))*0.5d0
              enddo
              sc(0) = -sc(1)

              ! Adjust nodes at top so that the top ghost cell center
              ! matches the LES exchange location cell center
              sc(nwm-1) = sc(nwm-2) + (sf(nwm+1)-sc(nwm-2))*thd
              sc(nwm  ) = sc(nwm-1) + (sf(nwm+1)-sc(nwm-2))*thd
              sc(nwm+1) = sf(nwm+1)

              !
              !- Get cell centers and face centers
              !
              ycwml(:) = he*sc(:) ! cell centers
              do k = 0,nwm
              yfwml(k) = 0.5d0*(ycwml(k)+ycwml(k+1)) ! face centers 'j+1/2'
              enddo

              !
              !- Initial guess for U and T profiles.
              !
              !- Note that you can also store these profiles for the next LES
              !  iteration for faster convergence. Instead of intializing
              !  at every iteration with a linear profile.
              !
              re = pe/(rgas*te)
              uwml(1:nwm) = ue*sc(1:nwm)
              twml(1:nwm) = te*sc(1:nwm) + twall*(1.0d0 - sc(1:nwm))
              rwml(1:nwm) = pe/(rgas*twml(1:nwm))

              ! Isothermal wall
              if (iwall==2) then 
                  tw = twall
              elseif (iwall==1) then
                  tw = twml(1)
              endif
              muw = fsuth(tw)

              ! Explicit BCs - U
              uwml(0)     = -uwml(1)    ! no-slip wall
              uwml(nwm+1) = ue          ! dirichlet

              ! Explicit BCs - RHO & T 
              if (iwall==1) then                    ! adiabatic
                  rw  = rwml(1)
                  tw  = twml(1)
                  muw = fsuth(tw)
                  rwml(0) = rwml(1)
                  twml(0) = twml(1)
              elseif (iwall==2 .or. iwall==0) then  ! isothermal
                  rw = rwml(1)*twml(1)/tw
                  rwml(0) = 2.0d0*rw - rwml(1)
                  twml(0) = 2.0d0*tw - twml(1)
              endif
              rwml(nwm+1) = re          ! dirichlet
              twml(nwm+1) = te 

              ! Set old
              uwml_n(:) = uwml(:)
              twml_n(:) = twml(:)

              !
              !- Iterate until convergence
              !
              DO ncr_wm = 1,kmax_wm

              !---------------------------------------------------
              !
              !- Solve momentum equation
              !
              !---------------------------------------------------

              !
              !- Compute mut
              !
              tauwm = dabs(muw*uwml(1)/ycwml(1))
              select case (iturb_wm)
              case(1) ! Mixing length model, JK
                  do k = 1,nwm+1
                  muwmf = fsuth(twml(k))
                  rtauw = sqrt(rwml(k)*tauwm)
                  ystar = ycwml(k)*rtauw/muwmf       ! semi-local scaling
                  vdamp = (1.0d0-exp(-ystar/apl1_wm))**2

                  mutwml(k) = kt_wm*ycwml(k)*rtauw*vdamp
                  enddo
                  mutwml(0) = -mutwml(1)
              case default
                  mutwml = 0
              end select

              !
              !- Form matrix system
              !
              d_rhs = 0.0d0
              do  k = 1,nwm

              twmf = 0.5d0*( twml(k) + twml(k-1) )
              rwmf = 0.5d0*( rwml(k) + rwml(k-1) )

              dywm   = ycwml(k) - ycwml(k-1)
              muwmf  = fsuth(twmf)
              mutwmf = 0.5d0*( mutwml(k) + mutwml(k-1) )

              !- Lower
              a_low(k) = (muwmf + mutwmf)/dywm

              twmf  = 0.5d0*( twml(k+1) + twml(k) )
              rwmf  = 0.5d0*( rwml(k+1) + rwml(k) )

              dywm   = ycwml(k+1) - ycwml(k)
              muwmf  = fsuth(twmf)
              mutwmf = 0.5d0*( mutwml(k+1) + mutwml(k) )

              !- Upper
              c_upp(k) = (muwmf + mutwmf)/dywm

              !- Diagonal
              b_dia(k) = -(a_low(k) + c_upp(k))

              enddo ! nwm

              !- Implicit BCs: no-slip wall
              b_dia(1) = b_dia(1) - a_low(1)
              d_rhs(1) = d_rhs(1) + 0.0d0

              !- Implicit BCs: dirichlet top
              d_rhs(nwm) = d_rhs(nwm) - ue*c_upp(nwm)

              !- Solve tridiagonal system for U
              call tdma(nwm,a_low,b_dia,c_upp,d_rhs,xuwm)

              ! Get Urms
              urms = sum((uwml(1:nwm) - xuwm(1:nwm))**2)
              urms = sqrt(urms/dble(nwm))/uin

              ! Update - U
              uwml(1:nwm) = xuwm(1:nwm)

              ! Explicit BCs - U
              uwml(0)     = -uwml(1)    ! no-slip wall
              uwml(nwm+1) = ue          ! dirichlet

              !---------------------------------------------------
              !
              !- Solve energy equation
              !
              !---------------------------------------------------

              !
              !- Compute mut
              !
              tauwm = dabs(muw*uwml(1)/ycwml(1))
              select case (iturb_wm)
              case(1) ! Mixing length model, JK
                  do k = 1,nwm+1
                  muwmf = fsuth(twml(k))
                  rtauw = sqrt(rwml(k)*tauwm)
                  ystar = ycwml(k)*rtauw/muwmf       ! semi-local scaling
                  vdamp = (1.0d0-exp(-ystar/apl1_wm))**2

                  mutwml(k) = kt_wm*ycwml(k)*rtauw*vdamp
                  enddo
                  mutwml(0) = -mutwml(1)
              case default
                  mutwml = 0
              end select

              !
              !- Form matrix system
              !
              d_rhs = 0.0d0
              do  k = 1,nwm

              uwmf = 0.5d0*( uwml(k) + uwml(k-1) )
              twmf = 0.5d0*( twml(k) + twml(k-1) )
              rwmf = 0.5d0*( rwml(k) + rwml(k-1) )

              duwm   = uwml (k) - uwml (k-1)
              dywm   = ycwml(k) - ycwml(k-1)
              muwmf  = fsuth(twmf)
              mutwmf = 0.5d0*( mutwml(k) + mutwml(k-1) )

              !- Lower
              a_low(k) = cpav*(muwmf*prli + mutwmf*prti)/dywm 

              !- RHS
              d_rhs(k) = d_rhs(k) + (muwmf+mutwmf)*uwmf*duwm/dywm

              uwmf = 0.5d0*( uwml(k+1) + uwml(k) )
              twmf = 0.5d0*( twml(k+1) + twml(k) )
              rwmf = 0.5d0*( rwml(k+1) + rwml(k) )

              duwm   = uwml (k+1) - uwml (k)
              dywm   = ycwml(k+1) - ycwml(k)
              muwmf  = fsuth(twmf)
              mutwmf = 0.5d0*( mutwml(k+1) + mutwml(k) )

              !- Upper
              c_upp(k) = cpav*(muwmf*prli + mutwmf*prti)/dywm 

              !- Diagonal
              b_dia(k) = -(a_low(k) + c_upp(k))

              !- RHS: viscous work 
              d_rhs(k) = d_rhs(k) - (muwmf+mutwmf)*uwmf*duwm/dywm

              enddo ! nwm

              !- Implicit BCs: wall
              if (iwall==1) then 
                  ! Adiabatic
                  b_dia(1) = b_dia(1) + a_low(1)
                  d_rhs(1) = d_rhs(1) + 0.0d0
              elseif (iwall==2 .or. iwall==0) then      
                  ! Isothermal
                  b_dia(1) = b_dia(1) - a_low(1)
                  d_rhs(1) = d_rhs(1) - 2.0d0*tw*a_low(1)
              endif

              !- Implicit BCs: dirichlet top
              d_rhs(nwm) = d_rhs(nwm) - te*c_upp(nwm)

              !- Solve tridiagonal system for T
              call tdma(nwm,a_low,b_dia,c_upp,d_rhs,xtwm)

              ! Get Trms
              trms = sum((twml(1:nwm) - xtwm(1:nwm))**2)
              trms = sqrt(trms/dble(nwm))/tin

              ! Update - RHO & T
              twml(1:nwm) = xtwm(1:nwm)           ! T
              rwml(1:nwm) = pe/(rgas*twml(1:nwm)) ! Rho

              ! Explicit BCs - RHO & T 
              if (iwall==1) then      ! adiabatic
                  rw  = rwml(1)
                  tw  = twml(1)
                  muw = fsuth(tw)
                  rwml(0) = rwml(1)
                  twml(0) = twml(1)
              elseif (iwall==2) then  ! isothermal
                  rw = rwml(1)*twml(1)/tw
                  rwml(0) = 2.0d0*rw - rwml(1)
                  twml(0) = 2.0d0*tw - twml(1)
              endif
              rwml(nwm+1) = re        ! dirichlet
              twml(nwm+1) = te 

              !- Check if urms/trms < 0.0001% free stream values
              if ((urms<res_wm .and. trms<res_wm) .or. ncr_wm==kmax_wm) then

                  ! Print profiles
                  if (idebug==1) then
                      write(100,*) '# ', ncr_wm 
                      do k = 0,nwm+1
                      write(100,'(99es15.6)') ycwml(k),rwml(k),uwml(k),twml(k),mutwml(k)
                      enddo
                  endif

                  tauw = muw*uwml(1)/ycwml(1)
                  qw   = cpav*prli*muw*(twml(1)-tw)/ycwml(1)
                  EXIT
              endif

              ENDDO ! kmax_wm

              return
          end subroutine solve_eqw


          !-----------------------------------------------------------
          ! Solves a tri-diagonal matrix system                         
          !-----------------------------------------------------------
          ! n - No. of columns/rows                          
          ! a - Lower diagonal                                        
          ! b - Main  diagonal                                        
          ! c - Upper diagonal                                          
          ! d - Right hand side                                         
          ! x - Solution vector                                      
          !-----------------------------------------------------------
          subroutine tdma(n,a,b,c,d,x)
              implicit none
              integer, intent(in)  :: n
              real(8), intent(in)  :: a(n),c(n),b(n),d(n)
              real(8), intent(out) :: x(n)

              integer :: i
              real(8) :: m,cp(n),dp(n)

              cp(1) = c(1)/b(1)
              dp(1) = d(1)/b(1)

              !- Forward elimination
              do i = 2,n
              m = b(i)-cp(i-1)*a(i)
              cp(i) = c(i)/m
              dp(i) = (d(i)-dp(i-1)*a(i))/m
              enddo

              !- Back substitution
              x(n) = dp(n)
              do i = n-1,1,-1
              x(i) = dp(i)-cp(i)*x(i+1)
              enddo

              return
          end subroutine tdma


      end module eqw
      !==============================================================!
      !==============================================================!
