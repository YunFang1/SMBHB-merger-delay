subroutine mergerrate(M_BH, z, theta, params_galaxy, params_stellar, num_samples, integration) 
    ! merger rate of SMBH binary
    real(8), intent(in)  :: M_BH, z
    integer, intent(in)  :: num_samples
    real(8), dimension(3), intent(in)    :: theta
    real(8), dimension(5,14), intent(in) :: params_stellar
    real(8), dimension(10), intent(in)   :: params_galaxy
    real(8), intent(out) :: integration
  
    integration = bhmergerrate(M_BH, z, theta, params_galaxy, params_stellar, num_samples)
  
    contains
      
      !##############################
      ! calculate merger event number  
      real(8) function merger_number(theta, params_galaxy, params_stellar, points_inner, points_outer)
      real(8), dimension(3) :: theta
      real(8), dimension(5,14) :: params_stellar
      real(8), dimension(10) :: params_galaxy
      integer :: points_inner, points_outer, ii
      real(8) :: z_min, z_max, z_i, M_bh_min, M_bh_max, M_bh_i, pp, qq, summ

      z_min = 0.0001
      z_max = 10
      M_bh_min = 5
      M_bh_max = 8
      summ = 0
      do ii = 1, points_outer
        pp = rand()
        qq = rand()
        z_i    = z_min + pp * (z_max - z_min)
        M_bh_i = M_bh_min + qq * (M_bh_max - M_bh_min)
        summ = summ + bhmergerrate(M_bh_i, z_i, theta, params_galaxy, params_stellar, points_inner)
      end do
      merger_number = (z_max - z_min) * (M_bh_max - M_bh_min) * summ/points_outer
      end function merger_number     
      
      !*******************************************************
      ! merger rate of SMBH binary as a function of (M_bh, z),
      ! for given alpha, Beta, sigma, and params_Rg, params_phi
      real(8) function bhmergerrate(M_BH, z, theta, params_galaxy, params_stellar, num_samples) 
      ! merger rate of SMBH binary
      real(8) :: M_BH, z
      integer :: num_samples
      real(8), dimension(3)    :: theta
      real(8), dimension(5,14) :: params_stellar
      real(8), dimension(10)   :: params_galaxy
      real(8) :: aa, bb, sigma_M, M_g_mean, t_min, t_max, M_g_min, M_g_max, dMgddMbh
      real(8) :: p, q, t_i, M_g_i, f, summ
      integer :: i 
  
      aa=8.690000
      bb=1.170000
      sigma_M = 0.280000/bb
      M_g_mean = M_BH/bb + 11 - aa/bb
      dMgddMbh = 1/bb

      t_min = 0.000100
      t_max = 13.500000 - t_z(z) 
      ! redshit of t_delay + t_z(z_event) is limited up to 15, since galaxy merger rate is ignorable 
      M_g_min = M_g_mean - 2*sigma_M
      M_g_max = M_g_mean + 2*sigma_M

      summ = 0
      do i = 1, num_samples
        p = rand()
        q = rand()
        t_i = t_min + p * (t_max - t_min)
        M_g_i = M_g_min + q * (M_g_max - M_g_min)
        f = pdfdelay(M_BH, t_i, theta) &
            * rategalaxy(M_g_i, z_t(t_i + t_z(z)), params_galaxy) * dMgddMbh &
            * phi(M_g_i, z_t(t_i + t_z(z)), params_stellar) * massrelation(M_g_i, M_BH) &
            * dvddz( z_t(t_i + t_z(z)) ) &
            * dzddtl( z_t(t_i + t_z(z)) ) * 1 / dzddtl(z)
        summ = summ + f
      end do 
      bhmergerrate = (t_max - t_min) * (M_g_max - M_g_min) * summ/num_samples

      end function bhmergerrate

      ! stellar mass function Phi
      real(8) function phi(M_g, z, params_stellar)
      real(8)     :: M_g, z
      real(8), dimension(5,14) :: params_stellar
      real(8), dimension(14)   :: param_Mstar, param_phi1star, param_phi2star, param_alpha1, param_alpha2, phi_i
      integer     :: i
      real(8)     :: x
  
      param_Mstar    = params_stellar(1,:)
      param_phi1star = params_stellar(2,:)
      param_phi2star = params_stellar(3,:)
      param_alpha1   = params_stellar(4,:)
      param_alpha2   = params_stellar(5,:)
      x = 10
    
      do i = 1, 14
        phi_i(i) = log(x) * exp( - 10**(M_g - param_Mstar(i))) * 10**(M_g - param_Mstar(i)) &
                   * ( param_phi1star(i) * 10**((M_g - param_Mstar(i)) * param_alpha1(i)) &
                   + param_phi2star(i) * 10**((M_g - param_Mstar(i)) * param_alpha2(i) ) )
      end do
      
      if (z >= 0 .AND. z < 0.125) then
        phi = phi_i(1)
      else if ( z >= 0.125.AND. z < 0.5) then
        phi = ( 1 - (z - 0.125) / (0.5 - 0.125) ) * phi_i(1) + (z-0.125)/(0.5-0.125) * phi_i(2)
      else if ( z >= 0.5  .AND. z < 1  ) then
        phi = ( 1 - (z - 0.5  ) / (1   - 0.5  ) ) * phi_i(2) + (z - 0.5)/(1 - 0.5  ) * phi_i(3)
      else if ( z >= 1    .AND. z < 1.5) then 
        phi = ( 1 - (z - 1    ) / (1.5 - 1    ) ) * phi_i(3) + (z - 1  )/(1.5 - 1  ) * phi_i(4)
      else if ( z >= 1.5  .AND. z < 2  ) then
        phi = ( 1 - (z - 1.5  ) / (2   - 1.5  ) ) * phi_i(4) + (z - 1.5)/(2   - 1.5) * phi_i(5)
      else if ( z >= 2    .AND. z < 2.5) then
        phi = ( 1 - (z - 2    ) / (2.5 - 2    ) ) * phi_i(5) + (z - 2  )/(2.5 - 2  ) * phi_i(6)
      else if ( z >= 2.5  .AND. z <3.25) then
        phi = ( 1 - (z - 2.5  ) / (3.25- 2.5  ) ) * phi_i(6) + (z - 2.5)/(3.25- 2.5) * phi_i(7)
      else if ( z >= 3.25 .AND. z <4.125) then
        phi = ( 1 - (z - 3.25 ) / (4.125-3.25 ) ) * phi_i(7) + (z -3.25)/(4.125-3.25)* phi_i(8)
      else if ( z >= 4.125.AND. z < 5  ) then
        phi = ( 1 - (z - 4.125) / (5 - 4.125  ) ) * phi_i(8) + (z-4.125)/(5 - 4.125) * phi_i(9)
      else if ( z >= 5    .AND. z < 6  ) then 
        phi = ( 1 - (z - 5    ) / (6 - 5      ) ) * phi_i(9) + ( z - 5 )/( 6 - 5   ) * phi_i(10)
      else if ( z >= 6    .AND. z < 7  ) then
        phi = ( 1 - (z - 6    ) / (7 - 6      ) ) * phi_i(10) + ( z - 6 )/( 7 - 6   ) * phi_i(11)
      else if ( z >= 7    .AND. z < 8  ) then
        phi = ( 1 - (z - 7    ) / (8 - 7      ) ) * phi_i(11) + ( z - 7 )/( 8 - 7   ) * phi_i(12)
      else if ( z >= 8    .AND. z < 9  ) then
        phi = ( 1 - (z - 8    ) / (9 - 8      ) ) * phi_i(12) + ( z - 8 )/( 9 - 8   ) * phi_i(13)
      else if ( z >= 9    .AND. z < 10  ) then
        phi = ( 1 - (z - 9    ) / (10 - 9     ) ) * phi_i(13) + ( z - 9 )/( 10 - 9  ) * phi_i(14)
      else if ( z >= 10    .AND. z < 10.5) then
        phi = ( 1 - (z - 10   ) / (10.5 - 10  ) ) * phi_i(14)
      else if ( z >= 10.5 )               then 
        phi = 0
      end if
      end function phi

      ! galaxy merger rate per galaxy R_g
      real(8) function rategalaxy(M_g, z, params_galaxy)

      real(8)  :: M_g, z
      real(8), dimension(10) :: params_galaxy
      real(8)  :: param_M0, param_A0, param_eta, param_alpha0_g, param_alpha1_g, param_beta0, param_beta1
      real(8)  :: param_gamma, param_delta0, param_delta1
      real(8)  :: Az, alphaz, betaz, deltaz
    
      param_M0       = params_galaxy(1)
      param_A0       = params_galaxy(2)
      param_eta      = params_galaxy(3)
      param_alpha0_g = params_galaxy(4)
      param_alpha1_g = params_galaxy(5)
      param_beta0    = params_galaxy(6)
      param_beta1    = params_galaxy(7)
      param_gamma    = params_galaxy(8)
      param_delta0   = params_galaxy(9)
      param_delta1   = params_galaxy(10)
    
      Az     = param_A0 * (1 + z)**param_eta
      alphaz = param_alpha0_g * (1 + z)**param_alpha1_g
      betaz  = param_beta0 * (1 + z)**param_beta1
      deltaz = param_delta0 * (1 + z)**param_delta1
    
      rategalaxy = Az * 10**( (M_g - 10) * alphaz ) * ( 1 + 10**( (M_g - log10(param_M0)) * deltaz ) ) &
            * (1 - 0.25**( 1 + betaz + param_gamma * (M_g - 10) ) ) / ( 1 + betaz + param_gamma * (M_g - 10) )
      end function rategalaxy
  
      ! relationship between the mass of SMBHB and host galaxy
      real(8) function massrelation(M_g, M_bh)
      real(8) :: M_g, M_bh 
      real(8) :: aa, bb, sigma_M, M_g_mean, pi
      pi = 3.14
      aa = 8.690000
      bb = 1.170000
      sigma_M  = 0.280000/bb
      M_g_mean = M_bh/bb + 11.000000 - aa/bb
      massrelation = 1 / (sqrt(2*pi) * sigma_M) * exp( - (M_g - M_g_mean)**2 / (2 *sigma_M**2) ) 
      end function massrelation

      ! delay time distribution
      real(8) function pdfdelay(Log10Mbh, t, theta)
      real(8) :: Log10Mbh, t, alpha, beta, sigma
      real(8), dimension(3) :: theta
      real(8) :: pi

      alpha = theta(1)
      beta  = theta(2)
      sigma = theta(3)

      pi = 3.141593
      pdfdelay = (sqrt(2/pi) * exp(-0.5*(t - (10**(-6 + Log10Mbh))**alpha * beta)**2 / sigma**2)) &
                 /(sigma * (1 + erf(((10**(-6 + Log10Mbh))**alpha * beta) / (sqrt(2.000000) * sigma ))))
      end function pdfdelay

      !!!!!!!!!!!!!!Cosmological quantities!!!!!!!!!!!!!!!!!
      ! d z/d t_L
      real(8) function dzddtl(z)
      real(8) :: z

      dzddtl = 0.069444 * (1 + z) * (0.700 + 0.300 * (1 + z)**3 )**0.500
      end function dzddtl
      
      ! t_L(z)
      real(8) function t_z(z)
      real(8) :: z
        
      t_z = 10.429376 - 5.737097*log((0.836660 + 1*Sqrt(1 + z*(0.900000 + (0.900000 + 0.300000*z)*z))) &
            /(-1.527525 + 1.*Sqrt(3.333333 + z*(3 + z*(3 + 1*z)))))
      end function t_z
      
      ! z(t_L)
      real(8) function z_t(x) 
      real(8) :: x

      z_t = -1.000000 + (0.02744296535135245*(-0.000030517578124999997*exp(-1.0458250331675938*x) &
            + 8.117551404190137e10*exp(-0.8715208609729945*x) - 2.887677987263166e10*exp(-0.6972166887783953*x) &
            + 3.8521549215967855e9*exp(-0.5229125165837964*x) - 2.283895396683414e8*exp(-0.3486083443891972*x) &
            + 5.077850861229274e6*exp(-0.17430417219459837*x))**0.3333333333333333*exp(0.5229125165837971*x))&
            /(126.43652557190416*exp(0.17430417219459904*x) - 22.488799485246357*exp(0.3486083443891981*x) &
            + 1.000000*exp(0.5229125165837971*x))
      end function z_t


    !*****************************************************************************
    !   Calculate function dVddz(z), which is a function contains hyp2f1
    !*****************************************************************************
      real(8) function dvddz(z)
      real(8) :: z
      real(8) :: a_para, b_para, c_para, x_para
      a_para = 0.333333333
      b_para = 0.500000000
      c_para = 1.333333333
      x_para = -0.42857143*(1 + z)**3
      dvddz = ( -0.7195300 + 0.75394744*(1 + z) * hyg2f1(a_para, b_para, c_para, x_para) )**2 &
              * 2485.98970031 / Sqrt(0.7 + 0.3*(1 + z)**3)
      end function dvddz


    !*****************************************************************************
      real(8) function hyg2f1(a_param, b_param, c_param, x_param)
    !*****************************************************************************
        !
        !! hyg2f1 evaluates the hypergeometric function 2F1(A,B,C,X).
        !
        !  Parameters:
        !
        !    Input, real ( kind = rk ) A, B, C, X, the arguments of the function.
        !    C must not be equal to a nonpositive integer.
              
      real (8)  :: a_param, b_param, c_param, x_param
    
      if ( abs(x_param) <= 1 ) then 
        hyg2f1 = hygfx ( a_param, b_param, c_param, x_param)
      else 
        hyg2f1 =   gamma ( c_param ) * gamma ( b_param - a_param ) / (gamma ( b_param ) * gamma ( c_param - a_param )) &
                 * (- x_param)**(- a_param) * hygfx( a_param, 1 - c_param + a_param, 1 - b_param + a_param, 1/x_param) &
                 + gamma ( c_param ) * gamma ( a_param - b_param ) / (gamma ( a_param ) * gamma ( c_param - b_param )) &
                 * (- x_param)**(- b_param) * hygfx( b_param, 1 - c_param + b_param, 1 - a_param + b_param, 1/x_param) 
      end if
      end function hyg2f1
          
      real (8) function hygfx ( a, b, c, x)
      real (8) :: a, b, c, x, hf
      real (8) :: a0, aa, bb, c0, c1, eps, f0, f1
      real (8) :: g0, g1, g2, g3, ga, gabc, gam, gb
      real (8) :: gbm, gc, gca, gcab, gcb, gm, hw
      integer  :: j, k, m,nm
      logical  :: l0, l1, l2, l3, l4, l5
      real (8) :: pa, pb, r, r0, r1, rm, rp, sm, sp, sp0, x1
      real (8), parameter :: el = 0.5772156649015329D+00
      real (8), parameter :: pi = 3.141592653589793D+00
    
      l0 = ( c == aint ( c ) ) .and. ( c < 0.0D+00 )
      l1 = ( 1.0D+00 - x < 1.0D-15 ) .and. ( c - a - b <= 0.0D+00 )
      l2 = ( a == aint ( a ) ) .and. ( a < 0.0D+00 )
      l3 = ( b == aint ( b ) ) .and. ( b < 0.0D+00 )
      l4 = ( c - a == aint ( c - a ) ) .and. ( c - a <= 0.0D+00 )
      l5 = ( c - b == aint ( c - b ) ) .and. ( c - b <= 0.0D+00 )
    
      ! if ( l0 .or. l1 ) then
      !   write ( *, '(a)' ) ' '
      !   write ( *, '(a)' ) 'HYGFX - Fatal error!'
      !   write ( *, '(a)' ) '  The hypergeometric series is divergent.'
      !   return
      ! end if
    
      if ( 0.95D+00 < x ) then
        eps = 1.0D-08
      else
        eps = 1.0D-15
      end if
    
      if ( x == 0.0D+00 .or. a == 0.0D+00 .or. b == 0.0D+00 ) then
    
        hf = 1.0D+00
    
      else if ( 1.0D+00 - x == eps .and. 0.0D+00 < c - a - b ) then
    
        gc = gamma ( c )
        gcab = gamma ( c - a - b )
        gca = gamma ( c - a )
        gcb = gamma ( c - b )
        hf = gc * gcab /( gca *gcb )
    
      else if ( 1.0D+00 + x <= eps .and. abs ( c - a + b - 1.0D+00 ) <= eps ) then
    
        g0 = sqrt ( pi ) * 2.0D+00**( - a )
        g1 = gamma ( c )
        g2 = gamma ( 1.0D+00 + a / 2.0D+00 - b )
        g3 = gamma ( 0.5D+00 + 0.5D+00 * a )
        hf = g0 * g1 / ( g2 * g3 )
    
      else if ( l2 .or. l3 ) then
    
        if ( l2 ) then
          nm = int ( abs ( a ) )
        end if
    
        if ( l3 ) then
          nm = int ( abs ( b ) )
        end if
    
        hf = 1.0D+00
        r = 1.0D+00
    
        do k = 1, nm
          r = r * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) &
            / ( k * ( c + k - 1.0D+00 ) ) * x
          hf = hf + r
        end do
    
      else if ( l4 .or. l5 ) then
    
        if ( l4 ) then
          nm = int ( abs ( c - a ) )
        end if
    
        if ( l5 ) then
          nm = int ( abs ( c - b ) )
        end if
    
        hf = 1.0D+00
        r  = 1.0D+00
        do k = 1, nm
          r = r * ( c - a + k - 1.0D+00 ) * ( c - b + k - 1.0D+00 ) &
            / ( k * ( c + k - 1.0D+00 ) ) * x
          hf = hf + r
        end do
        hf = ( 1.0D+00 - x )**( c - a - b ) * hf
    
      end if
    
      aa = a
      bb = b
      x1 = x
    !
    
      if ( 0.75D+00 <= x ) then
    
        gm = 0.0D+00
    
        if ( abs ( c - a - b - aint ( c - a - b ) ) < 1.0D-15 ) then
    
          m = int ( c - a - b )
          ga = gamma (a)
          gb = gamma ( b )
          gc = gamma ( c )
          gam = gamma ( a + m )
          gbm = gamma ( b + m )
          pa = psi ( a )
          pb = psi ( b )
    
          if ( m /= 0 ) then
            gm = 1.0D+00
          end if
    
          do j = 1, abs ( m ) - 1
            gm = gm * j
          end do
    
          rm = 1.0D+00
          do j = 1, abs ( m )
            rm = rm * j
          end do
    
          f0 = 1.0D+00
          r0 = 1.0D+00
          r1 = 1.0D+00
          sp0 = 0.0D+00
          sp = 0.0D+00
    
          if ( 0 <= m ) then
    
            c0 = gm * gc / ( gam * gbm )
            c1 = - gc * ( x - 1.0D+00 )**m / ( ga * gb * rm )
    
            do k = 1, m - 1
              r0 = r0 * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) &
                / ( k * ( k - m ) ) * ( 1.0D+00 - x )
              f0 = f0 + r0
            end do
    
            do k = 1, m
              sp0 = sp0 + 1.0D+00 / ( a + k - 1.0D+00 ) &
                + 1.0D+00 / ( b + k - 1.0D+00 ) - 1.0D+00 / real ( k, kind = 8 )
            end do
    
            f1 = pa + pb + sp0 + 2.0D+00 * el + log ( 1.0D+00 - x )
            hw = f1
    
            do k = 1, 250
    
              sp = sp + ( 1.0D+00 - a ) / ( k * ( a + k - 1.0D+00 ) ) &
                + ( 1.0D+00 - b ) / ( k * ( b + k - 1.0D+00 ) )
    
              sm = 0.0D+00
              do j = 1, m
                sm = sm + ( 1.0D+00 - a ) &
                  / ( ( j + k ) * ( a + j + k - 1.0D+00 ) ) &
                  + 1.0D+00 / ( b + j + k - 1.0D+00 )
              end do
    
              rp = pa + pb + 2.0D+00 * el + sp + sm + log ( 1.0D+00 - x )
    
              r1 = r1 * ( a + m + k - 1.0D+00 ) * ( b + m + k - 1.0D+00 ) &
                / ( k * ( m + k ) ) * ( 1.0D+00 - x )
    
              f1 = f1 + r1 * rp
    
              if ( abs ( f1 - hw ) < abs ( f1 ) * eps ) then
                exit
              end if
    
              hw = f1
    
            end do
    
            hf = f0 * c0 + f1 * c1
    
          else if ( m < 0 ) then
    
            m = - m
            c0 = gm * gc / ( ga * gb * ( 1.0D+00 - x )**m )
            c1 = - ( - 1 )**m * gc / ( gam * gbm * rm )
    
            do k = 1, m - 1
              r0 = r0 * ( a - m + k - 1.0D+00 ) * ( b - m + k - 1.0D+00 ) &
                / ( k * ( k - m ) ) * ( 1.0D+00 - x )
              f0 = f0 + r0
            end do
    
            do k = 1, m
              sp0 = sp0 + 1.0D+00 / real ( k, kind = 8 )
            end do
    
            f1 = pa + pb - sp0 + 2.0D+00 * el + log ( 1.0D+00 - x )
    
            do k = 1, 250
    
              sp = sp + ( 1.0D+00 - a ) &
                / ( k * ( a + k - 1.0D+00 ) ) &
                + ( 1.0D+00 - b ) / ( k * ( b + k - 1.0D+00 ) )
    
              sm = 0.0D+00
              do j = 1, m
                sm = sm + 1.0D+00 / real ( j + k, kind = 8 )
              end do
    
              rp = pa + pb + 2.0D+00 * el + sp - sm + log ( 1.0D+00 - x )
    
              r1 = r1 * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) &
                / ( k * ( m + k ) ) * ( 1.0D+00 - x )
    
              f1 = f1 + r1 * rp
    
              if ( abs ( f1 - hw ) < abs ( f1 ) * eps ) then
                exit
              end if
    
              hw = f1
    
            end do
    
            hf = f0 * c0 + f1 * c1
    
          end if
    
        else
    
          ga = gamma ( a )
          gb = gamma ( b )
          gc = gamma ( c )
          gca = gamma ( c - a )
          gcb = gamma ( c - b )
          gcab = gamma ( c - a - b )
          gabc = gamma ( a + b - c )
          c0 = gc * gcab / ( gca * gcb )
          c1 = gc * gabc / ( ga * gb ) * ( 1.0D+00 - x )**( c - a - b )
          hf = 0.0D+00
          r0 = c0
          r1 = c1
    
          do k = 1, 250
    
            r0 = r0 * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) &
              / ( k * ( a + b - c + k ) ) * ( 1.0D+00 - x )
    
            r1 = r1 * ( c - a + k - 1.0D+00 ) * ( c - b + k - 1.0D+00 ) &
              / ( k * ( c - a - b + k ) ) * ( 1.0D+00 - x )
    
            hf = hf + r0 + r1
    
            if ( abs ( hf - hw ) < abs ( hf ) * eps ) then
              exit
            end if
    
            hw = hf
    
          end do
    
          hf = hf + c0 + c1
    
        end if
    
      else
    
        a0 = 1.0D+00
      
        hf = 1.0D+00
        r = 1.0D+00
    
        do k = 1, 250
    
          r = r * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) &
            / ( k * ( c + k - 1.0D+00 ) ) * x
    
          hf = hf + r
    
          if ( abs ( hf - hw ) <= abs ( hf ) * eps ) then
            exit
          end if
    
          hw = hf
    
        end do
    
        hf = a0 * hf
    
      end if
    
      hygfx = hf
      end function hygfx
    
      
      real(8) function psi ( x )
    
      !*****************************************************************************80
      !
      !! PSI computes the PSI function.
      !
      !  Parameters:
      !
      !    Input, real ( kind = rk ) X, the argument.
      !
        real (8), parameter :: a1 = -0.83333333333333333D-01
        real (8), parameter :: a2 =  0.83333333333333333D-02
        real (8), parameter :: a3 = -0.39682539682539683D-02
        real (8), parameter :: a4 =  0.41666666666666667D-02
        real (8), parameter :: a5 = -0.75757575757575758D-02
        real (8), parameter :: a6 =  0.21092796092796093D-01
        real (8), parameter :: a7 = -0.83333333333333333D-01
        real (8), parameter :: a8 =  0.4432598039215686D+00
        real (8), parameter :: el = 0.5772156649015329D+00
        integer :: k, n
        real (8), parameter :: pi = 3.141592653589793D+00
        real (8) :: ps, s, x, x2, xa
    
        xa = abs ( x )
        s = 0.0D+00
      
        if ( x == aint ( x ) .and. x <= 0.0D+00 ) then
      
          ps = 1.0D+300
      
        else if ( xa == aint ( xa ) ) then
      
          n = int ( xa )
          do k = 1, n - 1
            s = s + 1.0D+00 / real ( k, kind = 8 )
          end do
      
          ps = - el + s
      
        else if ( xa + 0.5D+00 == aint ( xa + 0.5D+00 ) ) then
      
          n = int ( xa - 0.5D+00 )
      
          do k = 1, n
            s = s + 1.0D+00 / real ( 2 * k - 1, kind = 8 )
          end do
      
          ps = - el + 2.0D+00 * s - 1.386294361119891D+00
      
        else
      
          if ( xa < 10.0D+00 ) then
      
            n = 10 - int ( xa )
            do k = 0, n - 1
              s = s + 1.0D+00 / ( xa + real ( k, kind = 8 ) )
            end do
      
            xa = xa + real ( n, kind = 8 )
      
          end if
      
          x2 = 1.0D+00 / ( xa * xa )
      
          ps = log ( xa ) - 0.5D+00 / xa + x2 * ((((((( &
                   a8   &
            * x2 + a7 ) &
            * x2 + a6 ) &
            * x2 + a5 ) &
            * x2 + a4 ) &
            * x2 + a3 ) &
            * x2 + a2 ) &
            * x2 + a1 )
      
          ps = ps - s
      
        end if
      
        if ( x < 0.0D+00 ) then
          ps = ps - pi * cos ( pi * x ) / sin ( pi * x ) - 1.0D+00 / x
        end if
      
        psi = ps
      end function psi
      
end subroutine mergerrate




