      subroutine  debris_cover_main (H, HT, mask, G, G_obs, h_debris, inputlocation_debris, inoutdebris, meltout_debris, h_debrisini, term1_debris_x, term1_debris_y, term2_debris, term3_debris, fdebris, margin, VX, VY, area_debris, actual_mass_debris,expected_mass_debris,mass_ratio_hist,c_debris_out,h_debris_margin,sum_mass_meltout,sir,TMA,actual_mass_hist,dh_debrisdt_yearly,DHDT_yearly,H_old,pert_status)

      use        PARAMETER_LIST

      implicit none

!--   Parameter List in
      real, dimension(NX,NY)               :: HT,H_old
      double precision, dimension(NX,NY,3) :: H
      real,dimension(NX,NY)                :: mask
      real,dimension(NX,NY)                :: margin,h_debris_new_sq,DHDT_yearly
      real,dimension(NX,NY)                :: inoutdebris,c_debris_out,sum_mass_meltout,h_debris_new_margin
      real,dimension(NX,NY)                :: h_debris, h_debris_margin, h_debris_diff,sir,TMA,mass_supra_out_retreat	
      real,dimension(NX,NY)                :: meltout_debris, adv, sr, diff, adv_x, adv_y
      real,dimension(NX,NY)                :: h_debrisini, old_debris_mass2, old_debris_mass3
      real,dimension(NX,NY)                :: term1_debris_x, fl_debris_x, term1a_debris_x, term1b_debris_x, term4_debris_x	
      real,dimension(NX,NY)                :: term1_debris_y, fl_debris_y, term1a_debris_y, term1b_debris_y, term4_debris_y
      real,dimension(NX,NY)                :: term2_debris, advel_x, advel_x_plus, advel_x_min, advel_y_plus, advel_y_min, advel_y, vel_diff, vel_diff_x, vel_diff_y
      real,dimension(NX,NY)                :: term3_debris,flux_adv_out, dh_debrisdt_yearly, hdmean, slopemean_x, slopemean_y
      real,dimension(NX,NY)                :: fdebris, area_debris, actual_mass_debris, new_debris_mass, expected_mass_debris, old_debris_mass
      real,dimension(1)                    :: actual_mass_hist, expected_mass_hist, mass_ratio_hist,actual_mass_debris_before,expected_mass_hist_before, expected_mass_hist_new, mass_ratio_hist_new,F_out_new
      real,dimension(NX,NY)                :: inputlocation_debris,alpha,mass_in_supra, term1_debris_x_plus, term1_debris_x_min, term1_debris_y_plus, term1_debris_y_min
      real,dimension(NX,NY)                :: slope_x_deb, slope_y_deb, T_x_deb, T_y_deb, flux_x_in_deb, flux_x_out_deb, flux_x_deb, flux_y_deb
      real,dimension(NX,NY)                :: flux_y_in_deb, flux_y_out_deb, divergence_flux_deb, flux_grav_out,mass_in_flux
      double precision,dimension(NX,NY,NZ) :: VX,VY
      real                                 :: flux_x_adv, flux_y_adv, flux_y_grav, flux_x_grav,sum_squares_h_debris,sum_h_debris_margin,debris_input_area,depositionvol_debris, length_evac, area_evac
      real                                 :: max_slope, dzdx, dzdy, hsum, slopesum_x, slopesum_y
      real, dimension(4)                   :: slope_dir
      integer, parameter                   :: dp = kind(1.0d0)
      
!--   Parameter List out
      real,dimension(NX,NY) :: G,G_obs

!--   Variables
      integer               :: i,j,it_d,yr_d,iter_mc,iter_smolar,ii, jj,count_h_debris,pert_status, n_up, idir, count_up, k

!-----------------------------------------!
! Run the supraglacial debris cover model !
!-----------------------------------------!

write(*,*),'starting the loop for the supraglacial debris cover...'

yr_d = 1

do while (yr_d.le.numberofyears)

  it_d = 1

  do while (it_d.le.(nyears_d/deltat_d))
                                                                                    
  !------------!
  ! INITIALIZE !
  !------------!
     
  ! Initialize (only for the first time step)

  if((it_d).eq.1)then
        
      do J=1,NY 
         do I=1,NX                                                              
            inoutdebris(I,J) = 0.
            meltout_debris(I,J) = 0.
            inputlocation_debris(I,J) = 0.
            h_debrisini(I,J) = 0.
            adv_x(I,J) = 0.
            adv_y(I,J) = 0.
            term2_debris(I,J) = 0.
            term3_debris(I,J) = 0.
            alpha(I,J) = 0.
            diff(I,J) = 0.
            term4_debris_x(I,J) = 0.
            term4_debris_y(I,J) = 0.
            old_debris_mass(I,J) = 0.
            old_debris_mass2(I,J) = 0.
            old_debris_mass3(I,J) = 0.
            new_debris_mass(I,J) = 0.
            flux_adv_out(I,J) = 0.
            slope_x_deb(I,J) = 0.
            slope_y_deb(I,J) = 0.
            T_x_deb(i, j) = 0.
            T_y_deb(i, j) = 0.
            flux_x_in_deb(i, j) = 0.
            flux_x_out_deb(i, j) = 0.
            flux_y_in_deb(i, j) = 0.
            flux_y_out_deb(i, j) = 0.
            divergence_flux_deb(i,j) = 0.
            flux_grav_out(i, j) = 0.
            h_debris_new_sq(i,j) = 0.
            advel_x_plus(I,J) = 0.
            advel_x_min(I,J) = 0.
            advel_y_plus(I,J) = 0.
            advel_y_min(I,J) = 0.
            vel_diff_x(I,J) = 0.
            vel_diff_y(I,J) = 0.
            vel_diff(I,J) = 0.
            term1_debris_x_plus(I,J) = 0.
            term1_debris_y_plus(I,J) = 0.
            term1_debris_x_min(I,J) = 0.
            term1_debris_y_min(I,J) = 0.
            term1_debris_x(I,J) = 0.
            term1_debris_y(I,J) = 0.
            actual_mass_debris(I,J) = 0.
            h_debris_new_margin(I,J) = 0.
            mass_in_flux(I,J) = 0.
            mass_supra_out_retreat(I,J) = 0.
            if(dh_debrisdt_yearly(I,J).lt.0)then
               dh_debrisdt_yearly(I,J) = 0.
            endif
         end do
      end do

      debris_input_area = 0.
      depositionvol_debris = 0.

   !------------------------------------------------------!
   ! ACTIVE DEBRIS INPUT FROM TOPOGRAPHY IN ABLATION ZONE !
   !------------------------------------------------------!      
      
   ! Debris input (if debris source is situated in ablation zone)

      do J=1,NY 
         do I=1,NX
            if (I.eq.80.and.J.ge.82.and.J.le.96)then
               inputlocation_debris(I,J) = 0     ! No input in ablation zone here...
            else
               inputlocation_debris(I,J) = 0
           end if                            
         end do
       end do

       debris_input_area = sum(inputlocation_debris)*deltax_d*deltax_d                    ! Input area (in m^2)
       depositionvol_debris = deposition_mass_debris_abl / ((1-phi_debris)*rho_debris)    ! Input volume (in m^3/y)

       do J=1,NY
         do I=1,NX       
             if (inputlocation_debris(I,J).gt.0.and.H(I,J,3).gt.0)then
                inoutdebris(I,J) = (depositionvol_debris / debris_input_area)   ! Input rate (in m/y -> m^3/m^2*y)
             else
                inoutdebris(I,J) = 0.
             end if
         end do
      end do

      do J=1,NY
         do I=1,NX
           if (inputlocation_debris(I,J).gt.0.and.H(I,J,3).gt.0)then
              mass_in_flux(I,J) = (1-phi_debris)*rho_debris*(deltax_d*deltax_d*inoutdebris(I,J))
           else
              mass_in_flux(I,J) = 0.
           endif
         enddo
      enddo
      
   !----------------!
   ! DEBRIS MELTOUT !
   !----------------!

   ! Debris input from meltout of englacial debris
       
      do J=1,NY
         do I=1,NX
            if (sum_mass_meltout(I,J).gt.0)then
               meltout_debris(I,J) = -((sum_mass_meltout(I,J)) / ((1-phi_debris)*rho_debris)) / (deltax_d*deltax_d) ! From kg/yr to m/yr
            else
               meltout_debris(I,J) = 0.
            endif
         end do
       end do

       ! Set to zero if no ice
       
       do J=1,NY
          do I=1,NX
             if (H(I,J,3).eq.0)then
                meltout_debris(I,J)=0
             endif
          end do
       end do

   endif  ! endif it_d = 1

   !-----------------------------------!
   ! ADJUSTMENTS BEFORE ADVECTION STEP !
   !-----------------------------------!
   
   ! Reset output matrix for each time step

    flux_x_grav = 0.
    flux_y_grav = 0.

     do J=1,NY
        do I=1,NX
           flux_grav_out(I,J) = 0.
        enddo
     enddo

    ! Supraglacial mass before advection stpes
     
     if(it_d.eq.1)then

      do J=1,NY
         do I=1,NX
            if(h_debris(I,J).gt.0.0.and.H(I,J,3).gt.0)then
               actual_mass_debris(I,J) = (1-phi_debris)*rho_debris*(deltax_d*deltax_d*h_debris(I,J))
            else
               actual_mass_debris(I,J) = 0.
            endif
        end do
     end do

     actual_mass_debris_before = (sum(actual_mass_debris(:,:)))

     ! Put back to zero for now...

      do J=1,NY
         do I=1,NX
            actual_mass_debris(I,J) = 0.
         enddo
      enddo

      do J=1,NY
         do I=1,NX
            if(dh_debrisdt_yearly(I,J).gt.0)then
	       mass_supra_out_retreat(I,J) = ((1-phi_debris)*rho_debris*(deltax_d*deltax_d)*dh_debrisdt_yearly(I,J))
               expected_mass_debris(I,J) = expected_mass_debris(I,J) - ((1-phi_debris)*rho_debris*(deltax_d*deltax_d)*dh_debrisdt_yearly(I,J))
            else
               mass_supra_out_retreat(I,J) = 0.
               expected_mass_debris(I,J) = expected_mass_debris(I,J)
            endif
         enddo
      enddo
     
   endif  ! Endif it_d = 1
     
      !------------------------------------------------!    
      ! STEP 1: DEBRIS THICKNESS EVOLUTION CALCULATION !
      !------------------------------------------------!
 
      ! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ! cccc CALCULATE INITIAL DEBRIS THICKNESS WITHOUT NUMERICAL DIFFUSION CORRECTION cccc
      ! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do J=1,NY 
         do I=1,NX                                                                                                 
    
	   h_debrisini(I,J) = h_debris(I,J)

	   if((I.gt.1).and.(J.gt.1).and.(I.lt.NX).and.(J.lt.NY))then

           ! cccccccccccccccccccccccccccccccccccccccccccccc
           ! cccccccccc START ADVECTION SCHEME cccccccccccc
           ! cccccccccccccccccccccccccccccccccccccccccccccc

           ! cccccccccc UPWIND SCHEME ccccccccc

           ! Advection in x

           if (I.eq.2.or.I.eq.NX-1)then ! First order upstream advection scheme
             if (vx(I,J,1).lt.0)then
                adv_x(I,J) = -vx(I,J,1)*((h_debris(I+1,J)-h_debris(I,J))/(deltax_d)) - h_debris(I,J)*((vx(I+1,J,1)-vx(I,J,1))/(deltax_d))
             elseif (vx(I,J,1).ge.0)then
                adv_x(I,J) = -vx(I,J,1)*((h_debris(I,J)-h_debris(I-1,J))/(deltax_d)) - h_debris(I,J)*((vx(I,J,1)-vx(I-1,J,1))/(deltax_d))
             end if
            else ! Second order upstream advection scheme
             if (vx(I,J,1).lt.0)then
                adv_x(I,J) = -vx(I,J,1)*((-h_debris(I+2,J)+4*h_debris(I+1,J)-3*h_debris(I,J))/(2*deltax_d)) - h_debris(I,J)*((-vx(I+2,J,1)+4*vx(I+1,J,1)-3*vx(I,J,1))/(2*deltax_d))
             elseif (vx(I,J,1).ge.0)then
                adv_x(I,J) = -vx(I,J,1)*((3*h_debris(I,J)-4*h_debris(I-1,J)+h_debris(I-2,J))/(2*deltax_d)) - h_debris(I,J)*((3*vx(I,J,1)-4*vx(I-1,J,1)+vx(I-2,J,1))/(2*deltax_d))
             end if
            end if

           ! Advection in y
         
           if (J.eq.2.or.J.eq.NY-1)then ! First order upstream advection scheme
             if (vy(I,J,1).lt.0)then
                adv_y(I,J) = -vy(I,J,1)*((h_debris(I,J+1)-h_debris(I,J))/(deltay_d)) - h_debris(I,J)*((vy(I,J+1,1)-vy(I,J,1))/(deltay_d))
             elseif (vy(I,J,1).ge.0)then
                adv_y(I,J) = -vy(I,J,1)*((h_debris(I,J)-h_debris(I,J-1))/(deltay_d)) - h_debris(I,J)*((vy(I,J,1)-vy(I,J-1,1))/(deltay_d))
             end if
           else ! Second order upstream advection scheme
              if (vy(I,J,1).lt.0)then
                adv_y(I,J) = -vy(I,J,1)*((-h_debris(I,J+2)+4*h_debris(I,J+1)-3*h_debris(I,J))/(2*deltay_d)) - h_debris(I,J)*((-vy(I,J+2,1)+4*vy(I,J+1,1)-3*vy(I,J,1))/(2*deltay_d))
             elseif (vy(I,J,1).ge.0)then
                adv_y(I,J) = -vy(I,J,1)*((3*h_debris(I,J)-4*h_debris(I,J-1)+h_debris(I,J-2))/(2*deltay_d)) - h_debris(I,J)*((3*vy(I,J,1)-4*vy(I,J-1,1)+vy(I,J-2,1))/(2*deltay_d))
             end if
           end if

           ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
           ! ccccccccccccccccccccccccc ACTIVE INPUT ccccccccccccccccccccc
           ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

           ! Meltout of debris material from the ice (not at marginal pixels)
  
           term2_debris(I,J) = abs(meltout_debris(I,J))

           ! In and output from topography
     
           term3_debris(I,J) = inoutdebris(I,J)

           ! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
           ! ccccccccccccccccccc SLOPE-DRIVEN DISPLACEMENT ccccccccccccccccccccc  
           ! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    

           ! Slopes                                                                                         

           slope_x_deb(i,j) = ((HT(i+1,j)+h_debris(i+1,j)) - (HT(i-1,j)+h_debris(i-1,j))) / (2*deltax_d)
           slope_y_deb(i,j) = ((HT(i,j+1)+h_debris(i,j+1)) - (HT(i,j-1)+h_debris(i,j-1))) / (2*deltay_d)

           ! Transport fluxes (in m^2/y)                                                                

           T_x_deb(i,j) = c_deb_slope * (1-phi_debris)*rho_debris * gav * h_debris(i,j) * slope_x_deb(i,j)
           T_y_deb(i,j) = c_deb_slope * (1-phi_debris)*rho_debris * gav * h_debris(i,j) * slope_y_deb(i,j)

           ! Calculate flux in/out in x-direction                                                        

           if (slope_x_deb(i,j).lt.0)then
              flux_x_in_deb(i,j) = T_x_deb(i-1,j)
              flux_x_out_deb(i,j) = T_x_deb(i,j)
           else
              flux_x_in_deb(i,j) = -T_x_deb(i+1,j)
              flux_x_out_deb(i,j) = -T_x_deb(i,j)
           end if
           if(H(I,J,3).le.0.0)then
              flux_x_in_deb(i,j) = 0.
              flux_x_out_deb(i,j) = 0.
           endif
           if(flux_x_in_deb(I,J).gt.0.0)then
              flux_x_in_deb(I,J) = 0.
           endif
           if(flux_x_out_deb(I,J).gt.0.0)then
              flux_x_out_deb(I,J) = 0.
           endif
           flux_x_deb(i,j) = (flux_x_in_deb(i,j) - flux_x_out_deb(i,j))

           ! Calculate flux in/out in y-direction    

           if (slope_y_deb(i,j).lt.0)then
              flux_y_in_deb(i,j) = T_y_deb(i,j-1)
              flux_y_out_deb(i,j) = T_y_deb(i,j)
           else
              flux_y_in_deb(i,j) = -T_y_deb(i,j+1)
              flux_y_out_deb(i,j) = -T_y_deb(i,j)
           end if
           if(H(I,J,3).le.0.0)then
              flux_y_in_deb(i,j) = 0.
              flux_y_out_deb(i,j) = 0.
           endif
           if(flux_y_in_deb(I,J).gt.0.0)then
              flux_y_in_deb(I,J) = 0.
           endif
           if(flux_y_out_deb(I,J).gt.0.0)then
              flux_y_out_deb(I,J) = 0.
           endif
           flux_y_deb(i,j) = (flux_y_in_deb(i,j) - flux_y_out_deb(i,j))

           ! Flux divergence from slope-driven gravitational transport        

           divergence_flux_deb(i,j) = (flux_x_deb(i,j)/deltax_d + flux_y_deb(i,j)/deltay_d) 
                
           ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
           ! ccccccccccccccccccccc NEW DEBRIS THICKNESS ccccccccccccccccccccc
           ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    

            ! Calculate new debris thickness
  
            h_debris(I,J) = h_debris(I,J) + deltat_d*(adv_x(I,J) + adv_y(I,J) + term2_debris(I,J) + term3_debris(I,J) - divergence_flux_deb(I,J))
           
            ! Adjust false debris thickness

            if (h_debris(I,J).le.1e-5.and.sum_mass_meltout(I,J).eq.0)then
               h_debris(I,J) = 0.
            end if

            if (h_debris(I,J).le.0)then
               h_debris(I,J) = 0.
            end if
            
            if ((sqrt(vx(I,J,1) * vx(I,J,1) + vy(I,J,1) * vy(I,J,1)).eq.0).and.(h_debris(I,J).gt.0))then
               h_debris(I,J) = 0.
            end if

            if (h_debris(I, J) /= h_debris(I, J)) then
               h_debris(I,J) = 0.
            endif

            if (H(I,J,3).eq.0.or.G(I,J).gt.0)then
               h_debris(I,J) = 0.
            endif

           end if
         end do
      end do  
      
     ! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     ! cccccccccccccc REMOVE UNREALISTIC VALUES CHECK cccccccccccccc
     ! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      do J=1,NY
         do I=1,NX
            if (h_debris(I,J).le.1e-5.and.sum_mass_meltout(I,J).eq.0)then
               h_debris(I,J) = 0.
            end if
          end do
      end do

      do J=1,NY
         do I=1,NX      
            if (h_debris(I,J).le.0)then
	       h_debris(I,J) = 0.
            end if
          end do
      end do
            
      do J=1,NY
         do I=1,NX
            if(H(I,J,3).eq.0.or.G(I,J).gt.0)then
               h_debris(I,J) = 0.
             endif
          end do
      end do

      do J=1,NY
         do I=1,NX
            if (h_debris(I, J) /= h_debris(I, J)) then  ! Check for NaN        
               h_debris(I, J) = 0.0  ! Replace NaN with 0
            end if
         enddo
      enddo

      do J=1,NY
         do I=1,NX
            if(margin(I,J).eq.1.and.h_debris(I,J) .gt. 0 .and. &
             h_debris(I-1,J-1) .lt. 1e-5 .and. h_debris(I, J-1) .lt. 1e-5 .and. h_debris(I+1, J-1) .lt. 1e-5 .and. &
             h_debris(I-1,J) .lt. 1e-5 .and. h_debris(I+1, J) .lt. 1e-5 .and. &
             h_debris(I-1,J+1) .lt. 1e-5 .and. h_debris(I, J+1) .lt. 1e-5 .and. h_debris(I+1, J+1) .lt. 1e-5) then
               h_debris(I,J) = 0.
            endif
         enddo
      enddo

  !---------------------------------------------------------------------!
  ! STEP 2: ANTI-DIFFUSION CALCULATION (NUMERICAL DIFFUSION CORRECTION) !
  !---------------------------------------------------------------------!

  do iter_smolar = 1, n_iter_smolar ! Iterations for anti-diffusion correction
    
    ! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ! cccccccccccccc CALCULATE ANTI-DIFFUSION VELOCITIES AFTER SMOLARKIEWICZ cccccccccccc
    ! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    do J=1,NY
        do I=1,NX
	     if((I.gt.1).and.(J.gt.1).and.(I.lt.NX).and.(J.lt.NY).and.H(I,J,3).gt.0.and.h_debris(I,J).gt.0)then

                   ! X-DIRECTION  
                    advel_x_plus(I,J) = (abs(0.5*(vx(I,J,1) + vx(I+1,J,1))) * dx - (0.5*(vx(I,J,1) + vx(I+1,J,1)))**2 * deltat_d) * &
                    ((h_debris(I+1,J) - h_debris(I,J)) / ((h_debris(I+1,J) + h_debris(I,J) + epsilon_d) * dx)) &
                    - ((vx(I,J,1) + vx(I+1,J,1)) * (vy(I,J,1) + vy(I,J+1,1)) * deltat_d) * &
                    ((h_debris(I,J+1) - h_debris(I,J)) / ((h_debris(I,J+1) + h_debris(I,J) + epsilon_d) * dy))
                    advel_x_min(I,J) = (abs(0.5*(vx(I,J,1) + vx(I-1,J,1))) * dx - (0.5*(vx(I,J,1) + vx(I-1,J,1)))**2 * deltat_d) * &
                    ((h_debris(I,J) - h_debris(I-1,J)) / ((h_debris(I,J) + h_debris(I-1,J) + epsilon_d) * dx)) &
                    - ((vx(I,J,1) + vx(I-1,J,1)) * (vy(I,J,1) + vy(I,J-1,1)) * deltat_d) * &
                    ((h_debris(I,J) - h_debris(I,J-1)) / ((h_debris(I,J-1) + h_debris(I,J) + epsilon_d) * dy))

                    ! Y-DIRECTION  
                    advel_y_plus(I,J) = (abs(0.5*(vy(I,J,1) + vy(I,J+1,1))) * dy - (0.5*(vy(I,J,1) + vy(I,J+1,1)))**2 * deltat_d) * &
                    ((h_debris(I,J+1) - h_debris(I,J)) / ((h_debris(I,J+1) + h_debris(I,J) + epsilon_d) * dy)) &
                    - ((vx(I,J,1) + vx(I+1,J,1)) * (vy(I,J,1) + vy(I,J+1,1)) * deltat_d) * &
                    ((h_debris(I+1,J) - h_debris(I,J)) / ((h_debris(I+1,J) + h_debris(I,J) + epsilon_d) * dx))
                    advel_y_min(I,J) = (abs(0.5*(vy(I,J,1) + vy(I,J-1,1))) * dy - (0.5*(vy(I,J,1) + vy(I,J-1,1)))**2 * deltat_d) * &
                    ((h_debris(I,J) - h_debris(I,J-1)) / ((h_debris(I,J) + h_debris(I,J-1) + epsilon_d) * dy)) &
                    - ((vx(I,J,1) + vx(I-1,J,1)) * (vy(I,J,1) + vy(I,J-1,1)) * deltat_d) * &
                    ((h_debris(I,J) - h_debris(I-1,J)) / ((h_debris(I-1,J) + h_debris(I,J) + epsilon_d) * dx))

                    ! ANTI-DIFFUSION VELOCITY  
                    vel_diff_x(I,J) = 0.5 * (advel_x_plus(I,J) + advel_x_min(I,J))
                    vel_diff_y(I,J) = 0.5 * (advel_y_plus(I,J) + advel_y_min(I,J))
                    vel_diff(I,J) = sqrt((vel_diff_x(I,J)**2)+(vel_diff_y(I,J)**2))

                   ! cccccccccccccccccccccccccccccccccccccccccccccccccccccc  
                   ! cccccccccc CALCULATE ANTI-DIFFUSION FLUXES ccccccccccc  
                   ! cccccccccccccccccccccccccccccccccccccccccccccccccccccc  

                   ! X-DIRECTION  
                   term1_debris_x_plus(I,J) = 0.5 * (advel_x_plus(I,J) + abs(advel_x_plus(I,J))) * h_debris(I,J) + &
                          0.5 * (advel_x_plus(I,J) - abs(advel_x_plus(I,J))) * h_debris(I+1,J)
                   term1_debris_x_min(I,J) = 0.5 * (advel_x_min(I,J) + abs(advel_x_min(I,J))) * h_debris(I-1,J) + &
                          0.5 * (advel_x_min(I,J) - abs(advel_x_min(I,J))) * h_debris(I,J)
                   term1_debris_x(I,J) = 0.5 * (term1_debris_x_plus(I,J) - term1_debris_x_min(I,J)) / dx

                   ! Y-DIRECTION  
                   term1_debris_y_plus(I,J) = 0.5 * (advel_y_plus(I,J) + abs(advel_y_plus(I,J))) * h_debris(I,J) + &
                          0.5 * (advel_y_plus(I,J) - abs(advel_y_plus(I,J))) * h_debris(I,J+1)
                   term1_debris_y_min(I,J) = 0.5 * (advel_y_min(I,J) + abs(advel_y_min(I,J))) * h_debris(I,J-1) + &
                          0.5 * (advel_y_min(I,J) - abs(advel_y_min(I,J))) * h_debris(I,J)
                   term1_debris_y(I,J) = 0.5 * (term1_debris_y_plus(I,J) - term1_debris_y_min(I,J)) / dy

                else

                      advel_x_plus(I,J) = 0.
                      advel_x_min(I,J) = 0.
                      advel_y_plus(I,J) = 0.
                      advel_y_min(I,J) = 0.
                      vel_diff_x(I,J) = 0.
                      vel_diff_y(I,J) = 0.
                      vel_diff(I,J) = 0.
                      term1_debris_x_plus(I,J) = 0.
                      term1_debris_y_plus(I,J) = 0.
                      term1_debris_x_min(I,J) = 0.
                      term1_debris_y_min(I,J) = 0.
                      term1_debris_x(I,J) = 0.
                      term1_debris_y(I,J) = 0.

              endif ! endif H > 0 and h_debris > 0
                   
       enddo
   end do

   ! Do the correction step                            

   do J=2,NY-1
       do I=2,NX-1

               IF((H(I,J,3).gt.0.0).and.h_debris(I,J).gt.0)then

               ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc         
               ! cccccccccccccccccccc CORRECTION STEP ccccccccccccccccccccccccccc                  
               ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc                             

               h_debris(I,J) = h_debris(I,J) + deltat_d * (-term1_debris_x(I,J) - term1_debris_y(I,J))

               ! Adjust false debris thickness                                                   
                    if (h_debris(I,J).lt.1e-5.and.sum_mass_meltout(I,J).eq.0)then
                       h_debris(I,J) = 0.
                    end if

                    if(h_debris(I,J).le.0)then
                       h_debris(I,J) = 0.
                    endif
                    
                    if ((sqrt(vx(I,J,1) * vx(I,J,1) + vy(I,J,1) * vy(I,J,1)).eq.0).and.(h_debris(I,J).gt.0))then
                       h_debris(I,J) = 0.
                    end if

                    if (h_debris(I, J) /= h_debris(I, J)) then
                       h_debris(I,J) = 0.
                    endif

                    if (H(I,J,3).eq.0.or.G(I,J).gt.0)then
                       h_debris(I,J) = 0.
                    endif

               endif

      end do
   end do

  end do  ! end do smolar iterations

 ! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 ! cccccccccccccc REMOVE UNREALISTIC VALUES CHECK cccccccccccccc
 ! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 

     do J=1,NY
         do I=1,NX
            if (h_debris(I,J).le.1e-5.and.sum_mass_meltout(I,J).eq.0)then
               h_debris(I,J) = 0.
            end if
          end do
      end do

      do J=1,NY
         do I=1,NX
            if (h_debris(I,J).le.0)then
               h_debris(I,J) = 0.
            end if
          end do
      end do

      do J=1,NY
         do I=1,NX
            if(H(I,J,3).eq.0.or.G(I,J).gt.0)then
               h_debris(I,J) = 0.
             endif
          end do
     end do

     do J=1,NY
         do I=1,NX
            if (h_debris(I, J) /= h_debris(I, J)) then  ! Check for NaN       
               h_debris(I, J) = 0.0  ! Replace NaN with 0       
            end if
         enddo
     enddo

     do J=1,NY
         do I=1,NX
            if(margin(I,J).eq.1.and.h_debris(I,J) .gt. 0 .and. &
             h_debris(I-1,J-1) .lt. 1e-5 .and. h_debris(I, J-1) .lt. 1e-5 .and. h_debris(I+1, J-1) .lt. 1e-5 .and. &
             h_debris(I-1,J) .lt. 1e-5 .and. h_debris(I+1, J) .lt. 1e-5 .and. &
             h_debris(I-1,J+1) .lt. 1e-5 .and. h_debris(I, J+1) .lt. 1e-5 .and. h_debris(I+1, J+1) .lt. 1e-5) then
               h_debris(I,J) = 0.
            endif
         enddo
     enddo

     do J=2,NY-1
        do I=2,NX-1
        ! Check if the cell is at the glacier margin       
          if (margin(I,J) > 0) then
             sum_squares_h_debris = 0.0
             sum_h_debris_margin = 0.0
             count_h_debris = 0

            ! Loop over the 8 neighboring cells        
            do ii = I-1, I+1
                do jj = J-1, J+1
                    if ((ii .ne. I .or. jj .ne. J) .and. h_debris(ii,jj) .gt. 0.0) then
                        sum_squares_h_debris = sum_squares_h_debris + h_debris(ii,jj) ** 2
                        sum_h_debris_margin = sum_h_debris_margin + h_debris(ii,jj)
                        count_h_debris = count_h_debris + 1
                    end if
                end do
            end do

            ! Compute mean of neighboring cells (if neighbors exist)         
            if (count_h_debris .gt. 0) then
               h_debris_new_margin(I,J) = sum_h_debris_margin / count_h_debris     ! Mean value
               h_debris_new_sq(I,J) = sqrt(sum_squares_h_debris)                   ! Sqrt sum of squares
            else
               h_debris_new_margin(I,J) = h_debris(I,J)
               h_debris_new_sq(I,J) = h_debris(I,J)
            end if

          end if
        end do
    end do

    do J=2,NY-1
       do I=2,NX-1
          if (margin(I,J) > 0)then
            ! Check condition: If sum of squares of neighbors is lower than h_debris, update h_debris(I,J) with mean of neighbors
            if (h_debris_new_sq(I,J) .lt. h_debris(I,J)) then
                h_debris(I,J) = h_debris_new_margin(I,J)
             end if
          endif
       end do
    end do
      
 ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 ! cccccccccccc BOUNDARY CONDITION - REMOVAL OF DEBRIS OFF-GLACIER cccccccccccc
 ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 

 do J=1,NY
    do I=1,NX

       if (margin(i, j) > 0 .and. H(i,j,3) > 0.0) then

       ! ---------------------------------------------
       ! ----------- CASE A: lchar_d <= dx -----------
       ! ---------------------------------------------

         if (lchar_d .le. deltax_d)then
          
           ! Slopes       

            if (i > 1 .and. i < NX .and. j > 1 .and. j < NY) then
              slope_x_deb(i,j) = ((HT(i+1,j)+h_debris(i+1,j)) - (HT(i-1,j)+h_debris(i-1,j))) / (2*deltax_d)
              slope_y_deb(i,j) = ((HT(i,j+1)+h_debris(i,j+1)) - (HT(i,j-1)+h_debris(i,j-1))) / (2*deltay_d)
           else
              slope_x_deb(i,j) = 0.
              slope_y_deb(i,j) = 0.
           endif
        
           ! Transport fluxes (in m^2/y)                                    

           T_x_deb(i,j) = c_deb_slope * (1-phi_debris)*rho_debris * gav * h_debris(i,j) * slope_x_deb(i,j)
           T_y_deb(i,j) = c_deb_slope * (1-phi_debris)*rho_debris * gav * h_debris(i,j) * slope_y_deb(i,j)

           ! Calculate flux in/out in x-direction                  

           if (slope_x_deb(i,j).lt.0)then
              flux_x_in_deb(i,j) = T_x_deb(i-1,j)
              flux_x_out_deb(i,j) = T_x_deb(i,j)
           else
              flux_x_in_deb(i,j) = -T_x_deb(i+1,j)
              flux_x_out_deb(i,j) = -T_x_deb(i,j)
           end if
           if(H(I,J,3).le.0.0)then
              flux_x_in_deb(i,j) = 0.
              flux_x_out_deb(i,j) = 0.
           endif
           if(flux_x_in_deb(I,J).gt.0.0)then
              flux_x_in_deb(I,J) = 0.
           endif
           if(flux_x_out_deb(I,J).gt.0.0)then
              flux_x_out_deb(I,J) = 0.
           endif
           flux_x_deb(i,j) = (flux_x_in_deb(i,j) - flux_x_out_deb(i,j))

           ! Calculate flux in/out in y-direction   

           if (slope_y_deb(i,j).lt.0)then
              flux_y_in_deb(i,j) = T_y_deb(i,j-1)
              flux_y_out_deb(i,j) = T_y_deb(i,j)
           else
              flux_y_in_deb(i,j) = -T_y_deb(i,j+1)
              flux_y_out_deb(i,j) = -T_y_deb(i,j)
           end if
           if(H(I,J,3).le.0.0)then
              flux_y_in_deb(i,j) = 0.
              flux_y_out_deb(i,j) = 0.
           endif
           if(flux_y_in_deb(I,J).gt.0.0)then
              flux_y_in_deb(I,J) = 0.
           endif
           if(flux_y_out_deb(I,J).gt.0.0)then
              flux_y_out_deb(I,J) = 0.
           endif
           flux_y_deb(i,j) = (flux_y_in_deb(i,j) - flux_y_out_deb(i,j))

           ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
           ! cccccccccccccccccccccc DEBRIS FLUXES OFF-GLACIER FROM MARGIN (SLOPE-RELATED) ccccccccccccccccc   
           ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     

           ! cccccccccc GRAVITATIONAL DISPLACEMENT OFF-GLACIER INTO THE FORELAND cccccccccccccc
           
                ! Check for debris flux in the x-direction (left or right of the margin)                
                ! Left of the margin (i-1)                                                         
                if (i > 1) then
                   if (H(i-1, j,3) == 0.0 .and. H(i,j,3) > 0.0 .and. slope_x_deb(i,j) > 0.0) then
                      ! Slope is directed into the foreland, calculate debris flux                      
                      flux_x_grav = abs(flux_x_out_deb(i, j)/lchar_d) ! in m/y for each timestep dt      
                      flux_grav_out(i, j) = flux_grav_out(i, j) + flux_x_grav**2  ! Store flux value      
                   end if
                end if
                ! Right of the margin (i+1)                                    
                if (i < NX) then
                   if (H(i+1, j,3) == 0.0 .and. H(i,j,3) > 0.0 .and. slope_x_deb(i,j) < 0.0) then
                      ! Slope is directed into the foreland, calculate debris flux                     
                      flux_x_grav = abs(flux_x_out_deb(i, j)/lchar_d) ! in m/y for each timestep dt      
                      flux_grav_out(i, j) = flux_grav_out(i, j) + flux_x_grav**2  ! Store flux value       
                   end if
                end if
                ! Check for debris flux in the y-direction (below or above the margin)              
                ! Below the margin (j-1)                                                      
                if (j > 1) then
                   if (H(i, j-1,3) == 0.0 .and. H(i,j,3) > 0.0 .and. slope_y_deb(i,j) > 0.0) then
                      ! Slope is directed into the foreland, calculate debris flux                    
                      flux_y_grav = abs(flux_y_out_deb(i, j)/lchar_d) ! in m/y for each timestep dt     
                      flux_grav_out(i, j) = flux_grav_out(i, j) + flux_y_grav**2  ! Store flux value    
                   end if
                end if
                ! Above the margin (j+1)                             
                if (j < NY) then
                   if (H(i, j+1,3) == 0.0 .and. H(i,j,3) > 0.0 .and. slope_y_deb(i,j) < 0.0) then
                      ! Slope is directed into the foreland, calculate debris flux                       
                      flux_y_grav = abs(flux_y_out_deb(i, j)/lchar_d) ! in m/y for each timestep dt   
                      flux_grav_out(i, j) = flux_grav_out(i, j) + flux_y_grav**2  ! Store flux value    
                   end if
                end if
                flux_grav_out(i,j) = sqrt(flux_grav_out(i,j))   ! after summing the squared components, take the square root    
                flux_adv_out(i,j) = 0.

                elseif (lchar_d .gt. deltax_d)then  ! Else lchar_d vs. deltax_d

              ! ---------------------------------------------           
              ! ----------- CASE B: lchar_d > dx ------------   
              ! ---------------------------------------------

                   ! Number of pixels in the upstream band
                   n_up = ceiling(lchar_d / deltax_d)
                   
                   ! Initialize
                   slope_dir(1:4) = 0.0_dp
                   
                   ! Compute directional slopes (negative = uphill)
                   slope_dir(1) = ((HT(I,J)+h_debris(I,J)) - (HT(I,J+1)+h_debris(I,J+1))) / deltay_d   ! North
                   slope_dir(2) = ((HT(I,J)+h_debris(I,J)) - (HT(I,J-1)+h_debris(I,J-1))) / deltay_d   ! South
                   slope_dir(3) = ((HT(I,J)+h_debris(I,J)) - (HT(I+1,J)+h_debris(I+1,J))) / deltax_d   ! East
                   slope_dir(4) = ((HT(I,J)+h_debris(I,J)) - (HT(I-1,J)+h_debris(I-1,J))) / deltax_d   ! West

                   ! Find steepest uphill direction = most negative slope_dir
                   max_slope = 0.0_dp
                   idir = -1

                   do k = 1,4
                      if (slope_dir(k) < 0.0_dp .and. slope_dir(k) < max_slope) then   ! only negative slopes accepted
                         max_slope = slope_dir(k)
                         idir = k
                      end if
                   end do

                   ! idir = -1 means no uphill neighbour exists
                   ! Now initialize

                   hsum = 0.0_dp
                   slopesum_x = 0.0_dp
                   slopesum_y = 0.0_dp
                   count_up = 0
                   dzdx = 0.0_dp
                   dzdy = 0.0_dp

                   ! Slopes    

                   if (i > 1 .and. i < NX .and. j > 1 .and. j < NY) then
                      slope_x_deb(i,j) = ((HT(i+1,j)+h_debris(i+1,j)) - (HT(i-1,j)+h_debris(i-1,j))) / (2*deltax_d)
                      slope_y_deb(i,j) = ((HT(i,j+1)+h_debris(i,j+1)) - (HT(i,j-1)+h_debris(i,j-1))) / (2*deltay_d)
                   else
                      slope_x_deb(i,j) = 0.
                      slope_y_deb(i,j) = 0.
                   endif
                   
                   if (idir.eq.-1)then ! In the case no uphill slopes are found...

                   ! No uphill: use marginal pixel only
                   hdmean(i,j) = h_debris(i,j)
                   slopemean_x(i,j) = slope_x_deb(i,j)
                   slopemean_y(i,j) = slope_y_deb(i,j)
                   count_up = 0
                      
                   elseif (idir.gt.-1)then ! In the case uphill slopes are found...

                   ! Initialize with marginal pixel values
                   hsum = h_debris(i,j)
                   slopesum_x = slope_x_deb(i,j)
                   slopesum_y = slope_y_deb(i,j)
                   count_up = 1
                      
                   ! Walk n_up pixels upstream along idir from marginal pixel

                   ii = I
                   jj = J

                   do k = 1, max(0, n_up - 1)
                      
                      ! Each iteration we move to the upstream pixel along idir
                      select case (idir)
                      case (1)   ! north (increase j)
                         jj = jj + 1
                      case (2)   ! south
                         jj = jj - 1
                      case (3)   ! east (increase i)
                         ii = ii + 1
                      case (4)   ! west
                         ii = ii - 1
                      end select

                     ! Boundary / glacier-mask checks
                      if (ii < 1 .or. ii > NX .or. jj < 1 .or. jj > NY) exit
                      if (H(ii,jj,3) <= 0.0_dp) exit   ! outside glacier → stop including upstream pixels

                      ! Local slope magnitude at this upstream pixel:
                      if (ii > 1 .and. ii < NX .and. jj > 1 .and. jj < NY) then
                         dzdx = ((HT(ii+1,jj) + h_debris(ii+1,jj)) - (HT(ii-1,jj) + h_debris(ii-1,jj))) / (2*deltax_d)
                         dzdy = ((HT(ii,jj+1) + h_debris(ii,jj+1)) - (HT(ii,jj-1) + h_debris(ii,jj-1))) / (2*deltay_d)
                      else
                         exit
                      endif

                      ! Accumulate debris thickness & slope
                      hsum = hsum + h_debris(ii,jj)
                      slopesum_x = slopesum_x + dzdx
                      slopesum_y = slopesum_y + dzdy
                      count_up = count_up + 1

                   end do

                   endif

                   ! If no valid upstream pixels found → use the marginal cell only
                   if (count_up == 0 .or. idir == -1) then
                      hdmean(I,J) = h_debris(I,J)
                      slopemean_x(I,J) = ((HT(i+1,j)+h_debris(i+1,j)) - (HT(i-1,j)+h_debris(i-1,j))) / (2*deltax_d)
                      slopemean_y(I,J) = ((HT(i,j+1)+h_debris(i,j+1)) - (HT(i,j-1)+h_debris(i,j-1))) / (2*deltay_d)
                   else
                      hdmean(I,J) = hsum / real(count_up,dp)
                      slopemean_x(I,J) = slopesum_x / real(count_up,dp)
                      slopemean_y(I,J) = slopesum_y / real(count_up,dp)
                   endif

                   ! RECALCULATE DEBRIS REMOVAL SCHEME WITH DEBRIS THICKNESS AND SLOPE AVERAGED OVER LCHAR_D

                   ! Slopes
                   
                   if (i > 1 .and. i < NX .and. j > 1 .and. j < NY) then
                      slope_x_deb(i,j) = ((HT(i+1,j)+h_debris(i+1,j)) - (HT(i-1,j)+h_debris(i-1,j))) / (2*deltax_d)
                      slope_y_deb(i,j) = ((HT(i,j+1)+h_debris(i,j+1)) - (HT(i,j-1)+h_debris(i,j-1))) / (2*deltay_d)
                   else
                      slope_x_deb(i,j) = 0.
                      slope_y_deb(i,j) = 0.
                   endif
                
                   ! Transport fluxes (in m^2/y)                                                                                                                                                                         

                   T_x_deb(i,j) = c_deb_slope * (1-phi_debris)*rho_debris * gav * hdmean(I,J) * slopemean_x(I,J)
                   T_y_deb(i,j) = c_deb_slope * (1-phi_debris)*rho_debris * gav * hdmean(I,J) * slopemean_y(I,J)

                   ! Calculate flux in/out in x-direction                                                                                                                                                                
                   
                   if (slope_x_deb(i,j).lt.0)then
                      flux_x_in_deb(i,j) = T_x_deb(i-1,j)
                      flux_x_out_deb(i,j) = T_x_deb(i,j)
                   else
                      flux_x_in_deb(i,j) = -T_x_deb(i+1,j)
                      flux_x_out_deb(i,j) = -T_x_deb(i,j)
                   end if
                   if(H(I,J,3).le.0.0)then
                      flux_x_in_deb(i,j) = 0.
                      flux_x_out_deb(i,j) = 0.
                   endif
                   if(flux_x_in_deb(I,J).gt.0.0)then
                      flux_x_in_deb(I,J) = 0.
                   endif
                   if(flux_x_out_deb(I,J).gt.0.0)then
                      flux_x_out_deb(I,J) = 0.
                   endif
                   flux_x_deb(i,j) = (flux_x_in_deb(i,j) - flux_x_out_deb(i,j))

                   ! Calculate flux in/out in y-direction

                   if (slope_y_deb(i,j).lt.0)then
                      flux_y_in_deb(i,j) = T_y_deb(i,j-1)
                      flux_y_out_deb(i,j) = T_y_deb(i,j)
                   else
                      flux_y_in_deb(i,j) = -T_y_deb(i,j+1)
                      flux_y_out_deb(i,j) = -T_y_deb(i,j)
                   end if
                   if(H(I,J,3).le.0.0)then
                      flux_y_in_deb(i,j) = 0.
                      flux_y_out_deb(i,j) = 0.
                   endif
                   if(flux_y_in_deb(I,J).gt.0.0)then
                      flux_y_in_deb(I,J) = 0.
                   endif
                   if(flux_y_out_deb(I,J).gt.0.0)then
                      flux_y_out_deb(I,J) = 0.
                   endif
                   flux_y_deb(i,j) = (flux_y_in_deb(i,j) - flux_y_out_deb(i,j))
                   
                   ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc                                                                                                      
                   ! cccccccccccccccccccccc DEBRIS FLUXES OFF-GLACIER FROM MARGIN (SLOPE-RELATED) ccccccccccccccccc                                                                                                      
                   ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc                                                                                                      

                   ! cccccccccc GRAVITATIONAL DISPLACEMENT OFF-GLACIER INTO THE FORELAND cccccccccccccc                                                                                                                  

                   ! Check for debris flux in the x-direction (left or right of the margin)                                                                                                                         
                   ! Left of the margin (i-1)                                                                                                                                                                       
                   if (i > 1) then
                      if (H(i-1, j,3) == 0.0 .and. H(i,j,3) > 0.0 .and. slope_x_deb(i,j) > 0.0) then
                         ! Slope is directed into the foreland, calculate debris flux                                                                                                                               
                         flux_x_grav = abs(flux_x_out_deb(i, j)/lchar_d) ! in m/y for each timestep dt                                                                                                              
                         flux_grav_out(i, j) = flux_grav_out(i, j) + flux_x_grav**2  ! Store flux value                                                                                                             
                      end if
                   end if
                   ! Right of the margin (i+1)                                                                                                                                                                      
                   if (i < NX) then
                      if (H(i+1, j,3) == 0.0 .and. H(i,j,3) > 0.0 .and. slope_x_deb(i,j) < 0.0) then
                         ! Slope is directed into the foreland, calculate debris flux                                                                                                                               
                         flux_x_grav = abs(flux_x_out_deb(i, j)/lchar_d) ! in m/y for each timestep dt                                                                                                              
                         flux_grav_out(i, j) = flux_grav_out(i, j) + flux_x_grav**2  ! Store flux value                                                                                                             
                      end if
                   end if
                ! Check for debris flux in the y-direction (below or above the margin)                                                                                                                           
                ! Below the margin (j-1)                                                                                                                                                                         
                if (j > 1) then
                   if (H(i, j-1,3) == 0.0 .and. H(i,j,3) > 0.0 .and. slope_y_deb(i,j) > 0.0) then
                      ! Slope is directed into the foreland, calculate debris flux                                                                                                                               
                      flux_y_grav = abs(flux_y_out_deb(i, j)/lchar_d) ! in m/y for each timestep dt                                                                                                              
                      flux_grav_out(i, j) = flux_grav_out(i, j) + flux_y_grav**2  ! Store flux value                                                                                                             
                   end if
                end if
                ! Above the margin (j+1)                                                                                                                                                                         
                if (j < NY) then
                   if (H(i, j+1,3) == 0.0 .and. H(i,j,3) > 0.0 .and. slope_y_deb(i,j) < 0.0) then
                      ! Slope is directed into the foreland, calculate debris flux                                                                                                                               
                      flux_y_grav = abs(flux_y_out_deb(i, j)/lchar_d) ! in m/y for each timestep dt                                                                                                              
                      flux_grav_out(i, j) = flux_grav_out(i, j) + flux_y_grav**2  ! Store flux value                                                                                                             
                   end if
                end if

                flux_grav_out(i,j) = sqrt(flux_grav_out(i,j))   ! after summing the squared components, take the square root                                                                                     
                flux_adv_out(i,j) = 0.

               endif   ! Endif lchar_d vs. deltax_d
                  
             else    ! else (margin(i, j) > 0 .and. H(i,j,3) > 0.0) then (NOT A MARGIN PIXEL)

                flux_grav_out(i,j) = 0.
                flux_adv_out(i,j) = 0.

             endif   ! Endif (margin(i, j) > 0 .and. H(i,j,3) > 0.0) then

     enddo
 enddo
  
 ! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
 ! cccccccccccccc REMOVE UNREALISTIC VALUES CHECK cccccccccccccc            
 ! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     

     do J=1,NY
         do I=1,NX
            if (h_debris(I,J).le.1e-5.and.sum_mass_meltout(I,J).eq.0)then
               h_debris(I,J) = 0.
               flux_grav_out(i, j) = 0.
            end if
          end do
      end do

      do J=1,NY
	 do I=1,NX
            if (h_debris(I,J).le.0)then
               h_debris(I,J) = 0.
               flux_grav_out(i, j) = 0.
            end if
          end do
      end do
      
      do J=1,NY
         do I=1,NX
            if(H(I,J,3).eq.0.or.G(I,J).gt.0)then
               h_debris(I,J) = 0.
               flux_grav_out(i, j) = 0.
             endif
          end do
      end do

      do J=1,NY
         do I=1,NX
            if (h_debris(I, J) /= h_debris(I, J)) then  ! Check for NaN          
               h_debris(I, J) = 0.0  ! Replace NaN with 0                      
               flux_grav_out(i, j) = 0.
            end if
         enddo
      enddo

      do J=1,NY
         do I=1,NX
            if(margin(I,J).eq.1.and.h_debris(I,J) .gt. 0 .and. &
             h_debris(I-1,J-1) .lt. 1e-5 .and. h_debris(I, J-1) .lt. 1e-5 .and. h_debris(I+1, J-1) .lt. 1e-5 .and. &
             h_debris(I-1,J) .lt. 1e-5 .and. h_debris(I+1, J) .lt. 1e-5 .and. &
             h_debris(I-1,J+1) .lt. 1e-5 .and. h_debris(I, J+1) .lt. 1e-5 .and. h_debris(I+1, J+1) .lt. 1e-5) then
               h_debris(I,J) = 0.
               flux_grav_out(i, j) = 0.
            endif
         enddo
      enddo

      do J=1,NY
         do I=1,NX
            if(flux_grav_out(i, j).lt.1e-5)then
               flux_grav_out(i, j) = 0.
            endif      
         enddo
      enddo
 
  ! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
  ! cccccccccccccccccccccccc SAVE TERMS FOR MASS CONSERVATION ccccccccccccccccccccc                                  
  ! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc                  

  do J=1,NY
      do I=1,NX
         ! New mass is from input from topography and meltout of debris-loaded ice (if exists)                            
         new_debris_mass(I,J) = abs(mass_in_flux(I,J)) + abs(sum_mass_meltout(I,J))
         ! Old mass is gravitational displacement off-glacier from margin                         
         old_debris_mass(I,J) = ((1-phi_debris)*rho_debris*(deltax_d*lchar_d)*(flux_grav_out(I,J)))
         ! Expected debris mass is new mass minus old minus for every time step                                                    
         expected_mass_debris(I,J) = expected_mass_debris(I,J) + new_debris_mass(I,J)*deltat_d - old_debris_mass(I,J)*deltat_d
      end do
  end do

  ! Adjust false values
   
  do J=1,NY
     do I=1,NX
         if(H(I,J,3).eq.0.or.G(I,J).gt.0.or.h_debris(I,J).le.0.0)then
            new_debris_mass(I,J) = 0.
            old_debris_mass(I,J) = 0.
            expected_mass_debris(I,J) = 0.
          endif
      end do
  end do

  do J=1,NY
       do I=1,NX
          if (expected_mass_debris(I,J) /= expected_mass_debris(I, J)) then  ! Check for NaN            
               expected_mass_debris(I,J) = 0.0  ! Replace NaN with 0      
          end if
      enddo
  enddo
  
  ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc       
  ! ccccccccccccccccc UPDATE TIME STEP (1/DT) cccccccccccccccccccccc            
  ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc          

  it_d = it_d + 1

  end do   ! end do main loop iterations
  
  !----------------------------------!    
  ! STEP 3: ENSURE MASS CONSERVATION ! 
  !----------------------------------!

  ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! ccccccccc GET TERMS FOR MASS CONSERVATION CORRECTION ccccccccccc
  ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  ! Actual debris mass on top of glacier

   do J=1,NY
       do I=1,NX
          if(h_debris(I,J).gt.0.and.H(I,J,3).gt.0)then
             actual_mass_debris(I,J) = (1-phi_debris)*rho_debris*(deltax_d*deltax_d*h_debris(I,J))
          else
             actual_mass_debris(I,J) = 0.
          endif
      end do
   end do

   ! Get terms
  
   actual_mass_hist = (sum(actual_mass_debris(:,:)))         ! Actual debris mass after 1/dt iterations 
   expected_mass_hist = (sum(expected_mass_debris(:,:)))     ! Expected debris mass after 1/dt iterations
   expected_mass_hist_before = expected_mass_hist            ! Save this expected debris mass before correction
   mass_ratio_hist = (expected_mass_hist/actual_mass_hist)   ! Mass ratio before correction

!   write(*,*),'------ before correction --------'
!   write(*,*),'sum actual supra mass before adv= ',actual_mass_debris_before
!   write(*,*),'sum actual supra mass before corr= ',sum(actual_mass_debris)                                                                        
!   write(*,*),'sum expected supra mass before corr= ',sum(expected_mass_debris)                        
!   write(*,*),'mass ratio hist supra before= ',mass_ratio_hist(1)

   ! Mass conservation iterations
   
   do iter_mc = 1, max_iters_mc

    ! Initial guess of off-glacier deposition mass                              
    F_out_new = 0.0
    do i = 1, NX
        do j = 1, NY
            if(old_debris_mass(I,J).ne.0)then
              F_out_new = F_out_new + (mass_ratio_hist(1) * old_debris_mass(I,J))
            endif
        end do
     end do

    ! Update mass ratio with initial guess of new melt-out mass                                                
    expected_mass_hist_new = sum(actual_mass_debris_before) + sum(new_debris_mass) - F_out_new
    mass_ratio_hist_new = expected_mass_hist_new / actual_mass_hist
    
    ! Check for convergence                                                                         
    mass_ratio_hist(1) = mass_ratio_hist_new(1)                ! Mass ratio after correction                            
    expected_mass_hist(1) = expected_mass_hist_new(1)          ! Total expected mass after correction               
    if (abs(mass_ratio_hist_new(1) - mass_ratio_hist(1)) < 1.0E-4) exit

   end do ! end do iterations mass conservation            
   
  ! Ensure mass conservation

  if(mass_ratio_hist(1).gt.0.0)then
    do J=1,NY
      do I=1,NX
           h_debris(I,J) = h_debris(I,J)*mass_ratio_hist(1)
       end do
    end do
  endif
 
  ! Check if it worked

  do J=1,NY
       do I=1,NX
          actual_mass_debris(I,J) = (1-phi_debris)*rho_debris*(deltax_d*deltax_d*h_debris(I,J))
       end do
  end do

  ! Calculate sum of masses
  
   actual_mass_hist = (sum(actual_mass_debris(:,:)))
   mass_ratio_hist = (expected_mass_hist/actual_mass_hist)   ! Should be 1 now...

  ! Update the expected debris mass

  do J=1,NY
     do I=1,NX
        expected_mass_debris(I,J) = expected_mass_debris(I,J) * (expected_mass_hist(1)/expected_mass_hist_before(1))
     enddo
  enddo
   
   write(*,*),'------ after correction --------'
   write(*,*),'sum actual supra mass before= ',sum(actual_mass_debris_before)
   write(*,*),'sum actual supra mass after= ',sum(actual_mass_debris)
   write(*,*),'sum expected supra mass after= ',sum(expected_mass_debris)
   write(*,*),'mass ratio hist supra after= ',mass_ratio_hist(1)
   write(*,*),'mass out supra =',F_out_new
   write(*,*),'mass retreat supra =',sum(mass_supra_out_retreat)
   write(*,*),'mass in supra = ',sum(sum_mass_meltout) + sum(mass_in_flux)

   ! Area and length of debris evacuation

   length_evac = 0.
   area_evac = 0.
   
   do J = 1, NY
      do I = 1, NX
         ! Count only margin cells in ablation area (G < 0) with debris
         if (margin(I,J) .gt. 0 .and. G(I,J) .le. 0 .and. h_debris(I,J) .gt. 0) then
            length_evac = length_evac + deltax_d    ! each margin cell ≈ dx of front
         end if
      end do
   end do

   ! Evacuation area: length * L_d*
   area_evac = (length_evac * lchar_d) / 1e6          ! [km^2]

   write(*,*),'length evac = ', length_evac
   write(*,*),'area evac = ', area_evac
   
  !-----------------------------------------------------------!        
  ! STEP 4: DETERMINE DEBRIS-RELATED MELT-MODIFICATION FACTOR !      
  !-----------------------------------------------------------!

  ! Calculate the debris-related melt-modification factor fdebris                                                                                 

   do J=1,NY
      do I=1,NX
         if(h_debris(I,J).le.0.015)then
            fdebris(I,J) = (26.667*h_debris(I,J)) + 1
         else if(h_debris(I,J).gt.0.015.and.h_debris(I,J).le.0.04)then
            fdebris(I,J) = (-16*h_debris(I,J)) + 1.64
         else if(h_debris(I,J).gt.0.04)then
            fdebris(I,J) = 0.1061*(h_debris(I,J)**(-0.7205))-0.07925
         end if
      end do
   end do

!   do J=1,NY
!      do I=1,NX
!         if(h_debris(I,J).le.0.00375)then
!            fdebris(I,J) = (26.667*h_debris(I,J)) + 1
!         else if(h_debris(I,J).gt.0.00375.and.h_debris(I,J).le.0.01025)then
!            fdebris(I,J) = ((-16*h_debris(I,J)) + 1.16)
!         else if(h_debris(I,J).gt.0.01025)then
!            fdebris(I,J) = (0.1061*(h_debris(I,J)**(-0.7205))-0.225)*0.375
!         end if
!      end do
!   end do
   
   do J=1,NY
      do I=1,NX
         if(fdebris(I,J).le.1e-3)then
            fdebris(I,J)=1e-3
         end if
      end do
   end do

  !--------------------------------------------------!  
  ! STEP 5: DETERMINE FRACTIONAL DEBRIS-COVERED AREA !
  !--------------------------------------------------!

  ! Calculate fractional debris-covered area

   do J=1,NY
      do I=1,NX
         area_debris(I,J)=(1-exp(-20*h_debris(I,J)))
      end do
  end do

!fdebris = 1.0 !!!!!!!!!!
!sir = flux_grav_out !!!!!!!!!
!TMA = vel_diff !!!!!!!!!
TMA = dh_debrisdt_yearly
sir = mass_supra_out_retreat

! Update time step for next year

yr_d = yr_d + 1

end do

      end subroutine debris_cover_main
