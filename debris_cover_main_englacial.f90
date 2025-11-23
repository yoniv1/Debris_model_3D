      subroutine  debris_cover_main_englacial(H, HT, G, G_obs, h_debris, margin, VX, VY, VZ, VZAGE, meltout_debris, inputlocation_debris_acc, ui, vi, wi, za, ti, ti_new,actual_mass_hist_eng, expected_mass_hist_eng, expected_debris_mass_eng,c_debris_out,DHDT_yearly,ti_add,t_s,pert_status,H_old,TIME,sum_mass_meltout,mass_in_old,margin_engl,max_ti_hist, inoutdebris_acc,DHTDT_yearly,dc_ddt_yearly,mask)

      use        PARAMETER_LIST

      implicit none

!--   Parameter List in                                                                                                                                                                  
      real, dimension(NX,NY)               :: HT,DHDT_yearly, H_old,mask
      double precision                     :: TIME
      double precision, dimension(NX,NY,3) :: H
      real,dimension(NX,NY)                :: margin,margin_engl
      real,dimension(NX,NY)                :: h_debris, c_debris_out,count_debris,cells_debris,mass_engl_out_retreat,dc_ddt_yearly
      real,dimension(NX,NY)                :: meltout_debris, inoutdebris_acc,mass_in,mass_out,loop_number,mass_in_old,mass_melt_out,DHTDT_yearly       
      real,dimension(NX,NY)                :: inputlocation_debris_acc, t_s,sum_C,ti_corr,th_per_piece,meltout_pieces,sum_c_meltout,sum_mass_meltout
      double precision,dimension(NX,NY,NZ) :: VX,VY,VZ,VZAGE, ui, vi, wi, ti, ti_new, tinit, advx_eng, advy_eng, advz_eng, advel_x_plus, advel_x_min, advel_y_plus, advel_y_min
      double precision,dimension(NX,NY,NZ) :: advel_z_plus, advel_z_min, ti_diff, vel_diff, actual_mass_debris_eng, sink_debris, new_debris_mass_eng, expected_debris_mass_eng
      double precision,dimension(NX,NY,NZ) :: old_debris_mass_eng, term1_debris_eng_x, term1_debris_eng_y, term1_debris_eng_z, new_C,ti_add,ti_add2,ti_add3,c_meltout,ti_hist
      double precision,dimension(NX,NY,NZ) :: latx_eng,laty_eng,temp_eng,advel_x,advel_y,advel_z, vel_diff_x, vel_diff_y, vel_diff_z, term1_debris_eng_x_plus,term1_debris_eng_x_min
      double precision,dimension(NX,NY,NZ) :: term1_debris_eng_y_plus,term1_debris_eng_y_min,term1_debris_eng_z_plus, term1_debris_eng_z_min,actual_mass_debris_eng_before
      real,dimension(NZ)                   :: za, za_eng, dzeta, dza_eng
      real,dimension(1)                    :: actual_mass_hist_eng, expected_mass_hist_eng, mass_ratio_hist_eng, expected_mass_before_eng,mass_ratio_hist_eng_new,expected_mass_hist_eng_new
      
!--   Parameter List out                                                                                                                                                                 
      real,dimension(NX,NY) :: G,G_obs

!--   Variables                                                                                                                                                                          
      integer               :: i,j,k,yr_deng, it_deng,it_smolar_eng,c,pert_status,W,iter_smolar,iter_mc
      real                  :: startTime, stopTime,max_ti_hist,debris_input_area, deposition_vol_debris, F_out_new
      real, allocatable     :: ti_corr2(:,:)

!---------------------------!                              
! FLIP MATRICES UPSIDE DOWN !                                                                          
!---------------------------!

! Flip velocity matrices    
                                                   
do K = 1, NZ
   ui(:,:,K) = vx(:,:,NZ+1 - K)
   vi(:,:,K) = vy(:,:,NZ+1 - K)
   wi(:,:,K) = vzage(:,:,NZ+1 - K)
   za_eng(K) = za(NZ+1 - K)
end do

! Adjust zeros for vertical velocity                                       

do J=1,NY
   do I=1,NX
      do K=1,NZ
          if(H(I,J,3).gt.0.and.G(I,J).le.0.and.wi(I,J,K).eq.0)then
              wi(I,J,K) = 1e-5
           endif
       end do
   end do
end do

!----------------------------------------------!        
! RUN THE ENGLACIAL DEBRIS CONCENTRATION MODEL !
!----------------------------------------------!

write(*,*),'starting the loop for the englacial debris concentration...'

yr_deng = 1

do while (yr_deng.le.numberofyears_eng)

  it_deng = 1

  do while (it_deng.le.(nyears_deng/dt_eng))

  !---------------------------------!                        
  ! DEBRIS INPUT AND INITIALIZATION !                      
  !---------------------------------!

  if(it_deng.eq.1)then  ! Only for the first time step
     
  ! cccccccccccccccccccccccccccccccccccccccccccccc 
  ! ccccccccccc INITIALIZE MATRICES cccccccccccccc  
  ! cccccccccccccccccccccccccccccccccccccccccccccc         

  ! Initialize 1D
     
  expected_mass_before_eng = 0.   
  actual_mass_hist_eng = 0.
  mass_ratio_hist_eng = 0.
  
  ! Initialize 2D

      do J=1,NY 
         do I=1,NX      
            inoutdebris_acc(I,J) = 0.
            mass_out(I,J) = 0.
            mass_in(I,J) = 0.
            sink_debris(I,J,1) = 0.
            mass_melt_out(I,J) = 0.
            mass_engl_out_retreat(I,J) = 0.
            if(dc_ddt_yearly(I,J).lt.0)then
               dc_ddt_yearly(I,J) = 0.
            endif
         end do
       end do

  ! Initialize 3D                                                                                                

      do J=1,NY
         do I=1,NX
            do K=1,NZ
              advx_eng(I,J,K) = 0.
              advy_eng(I,J,K) = 0.
              advz_eng(I,J,K) = 0.
              latx_eng(I,J,K) = 0.
              laty_eng(I,J,K) = 0.
              ti_add2(I,J,K) = 0.
              ti_add(I,J,K) = 0.
              advel_x(I,J,K) = 0.
              advel_y(I,J,K) = 0.
              advel_z(I,J,K) = 0.
              advel_x_plus(I,J,K) = 0.
              advel_x_min(I,J,K) = 0.
              advel_y_plus(I,J,K) = 0.
              advel_y_min(I,J,K) = 0.
              advel_z_plus(I,J,K) = 0.
              advel_z_min(I,J,K) = 0.
              vel_diff(I,J,K) = 0.
              vel_diff_x(I,J,K) = 0.
              vel_diff_y(I,J,K) = 0.
              vel_diff_z(I,J,K) = 0.
              term1_debris_eng_x_plus(I,J,K) = 0.
              term1_debris_eng_x_min(I,J,K) = 0.
              term1_debris_eng_y_plus(I,J,K) = 0.
              term1_debris_eng_y_min(I,J,K) = 0.
              term1_debris_eng_z_plus(I,J,K) = 0.
              term1_debris_eng_z_min(I,J,K) = 0.
              actual_mass_debris_eng(I,J,K) = 0.
              new_debris_mass_eng(I,J,K) = 0.
           end do
         end do
      end do

  ! Zeta coordinate   
  
      do K=1,NZ-1
          dzeta(K)=za_eng(K+1)-za_eng(K)
      enddo

      dzeta(NZ) = dzeta(NZ-1)
      dza_eng = dzeta

  endif ! endif it_deng = 1
      
  ! cccccccccccccccccccccccccccccccccccccccccccccc
  ! ccccccccccccc DEFINE THE INPUT ccccccccccccccc                       
  ! cccccccccccccccccccccccccccccccccccccccccccccc               

  ! Debris input                  

      do J=1,NY
         do I=1,NX

             if (I.eq.6.and.J.ge.83.and.J.le.95)then
                inputlocation_debris_acc(I,J) = 1
             end if
             if (I.eq.7.and.J.ge.83.and.J.le.95)then
                inputlocation_debris_acc(I,J) = 1
             end if

         end do
       end do

       debris_input_area = sum(inputlocation_debris_acc)*deltax_d*deltax_d                ! Input area (in m^2)                                          
       deposition_vol_debris = deposition_mass_debris_acc  / ((1-phi_debris)*rho_debris)  ! Input volume from input mass (in m^3/y)   

       do J=1,NY
         do I=1,NX

             if (inputlocation_debris_acc(I,J).gt.0)then
                inoutdebris_acc(I,J) = (deposition_vol_debris / debris_input_area)    ! Input rate (in m/y -> m^3/(y*m^2))     
             else
                inoutdebris_acc(I,J) = 0.
             end if

         end do
      end do

   ! Recalculate to concentrations      

      do J=1,NY
         do I=1,NX
           if(H(I,J,3).gt.0.and.G(i,j).gt.0)then
              tinit(I,J,NZ) = (inoutdebris_acc(I,J)*(1-phi_debris)*rho_debris*deltax_d*deltax_d) / (G(I,J)*deltax_d*deltax_d)
           else
              tinit(I,J,NZ) = 0.
           endif
         enddo
      enddo

  ! Set top boundary (input) condition

     do J=1,NY
        do I=1,NX
              if (inputlocation_debris_acc(I,J).gt.0)then
                 t_s(I,J) = tinit(I,J,NZ)
              else
                 t_s(I,J) = 0.
              endif
        end do
     end do

  ! cccccccccccccccccccccccccccccccccccccccccccccc 
  ! cccccc APPLY INPUT BOUNDARY CONDITION cccccccc  
  ! cccccccccccccccccccccccccccccccccccccccccccccc    

  ! Initialization of ti and ti_new

      do J=1,NY
         do I=1,NX
                ti(I,J,NZ) = t_s(I,J)
	        ti_new(I,J,NZ) = t_s(I,J)
         end do
      end do

   ! Mass input   

      do J=1,NY
         do I=1,NX
              mass_in(I,J) = (inoutdebris_acc(I,J)*(1-phi_debris)*rho_debris*deltax_d*deltax_d)
         end do
      end do

  ! cccccccccccccccccccccccccccccccccccccccccccccc 
  ! ccccccccccc CORRECTION FOR DHDT cccccccccccccc     
  ! cccccccccccccccccccccccccccccccccccccccccccccc 

  if(it_deng.eq.1.and.TIME.gt.tdebris+1)then

  ! Initialize

      do J=1,NY
         do I=1,NX
           do K=1,NZ
              ti_add2(I,J,K) = 0
           end do
         end do
      end do

   ! Take into account dHdt

      do J=1,NY
         do I=1,NX
           do K=1,NZ
              if(ti(I,J,K).gt.0)then
                ti_add2(I,J,K) = -ti(I,J,K) * (NZ/H(I,J,3)) * (-DHDT_yearly(I,J)/NZ)
              else
                ti_add2(I,J,K) = 0
           end if
         end do
       end do
     end do

  ! Delete false pixels

      do J=1,NY
         do I=1,NX
           do K=1,NZ
              if(H(I,J,3).eq.0)then
                ti_add2(I,J,K) = 0
              else
                ti_add2(I,J,K) = ti_add2(I,J,K)
           end if
         end do
       end do
     end do

   ! Adjust debris concentration accordingly

      ti = ti + ti_add2
      ti_new = ti_new + ti_add2

       do J=1,NY
         do I=1,NX
           if(G(I,J).gt.0.and.inputlocation_debris_acc(I,J).gt.0)then
              t_s(I,J) = ti(I,J,NZ)
           else
              t_s(I,J) = 0.
           endif
         end do
       end do

   ! Mass input     

    do J=1,NY
        do I=1,NX
              mass_in(I,J) = (inoutdebris_acc(I,J)*(1-phi_debris)*rho_debris*deltax_d*deltax_d)
        end do
    end do

    do J=1,NY
       do I=1,NX
          do K=1,NZ
             actual_mass_debris_eng_before(I,J,K) = 0.
          end do
      end do
   end do

   endif ! endif it_deng = 1

   ! ccccccccccccccccccccccccccccccccccccccccccccccc
   ! cccccc ADJUSTMENTS BEFORE ADVECTION STEP cccccc
   ! ccccccccccccccccccccccccccccccccccccccccccccccc

   if(it_deng.eq.1)then
      
     do J=2,NY-1
       do I=2,NX-1
          do K=2,NZ-1
             if(ti(I,J,K).gt.0.0.and.H(I,J,3).gt.0)then
                actual_mass_debris_eng_before(I,J,K) = (H(I,J,3)/nz)*(deltax_d*deltax_d*ti(I,J,K))
             else
                actual_mass_debris_eng_before(I,J,K) = 0.
             endif
           end do
       end do
    end do

    write(*,*),'actual mass before loop = ',sum(actual_mass_debris_eng_before)
    write(*,*),'expected mass before loop = ',sum(expected_debris_mass_eng(:,:,NZ-1))
    
    do J=1,NY
       do I=1,NX
          if(dc_ddt_yearly(I,J).gt.0)then
            mass_engl_out_retreat(I,J) = (H_old(I,J))*(deltax_d*deltax_d*dc_ddt_yearly(I,J))
            expected_debris_mass_eng(I,J,NZ-1) = expected_debris_mass_eng(I,J,NZ-1) - (H_old(I,J)*(deltax_d*deltax_d*dc_ddt_yearly(I,J)))
         else
            mass_engl_out_retreat(I,J) = 0.
            expected_debris_mass_eng(I,J,NZ-1) = expected_debris_mass_eng(I,J,NZ-1)
          endif
      end do
    end do

    write(*,*),'actual mass engl before loop = ',sum(actual_mass_debris_eng_before)
    write(*,*),'expected mass engl before loop = ',sum(expected_debris_mass_eng(:,:,NZ-1))
    write(*,*),'retreat mass engl before loop = ',sum(mass_engl_out_retreat)
    
   endif ! endif it_deng = 1
   
   !------------------------------------------!    
   ! STEP 1: DEBRIS CONCENTRATION CALCULATION !
   !------------------------------------------!

      do J=1,NY 
         do I=1,NX
	   do K=1,NZ                                                                                                 
    
	   ti_new(I,J,K) = ti(I,J,K)

	   IF((H(I,J,3).ge.0).and.(I.gt.1).and.(J.gt.1).and.(K.gt.1).and.(I.lt.NX).and.(J.lt.NY).and.(K.lt.NZ))then

            ! cccccccccccccccccccccccccccccccccccccccccccccc
            ! cccccccccc START ADVECTION SCHEME cccccccccccc
            ! cccccccccccccccccccccccccccccccccccccccccccccc

            ! cccccccccccccccccccccccccccccccccccccccccccccc                                                                                    
            ! cccccccccccccccc Advection in X cccccccccccccc                   
            ! cccccccccccccccccccccccccccccccccccccccccccccc                

             if ((I.eq.2).or.(K.eq.2).or.(I.eq.NX-1).or.(K.eq.NZ-1)) then
               if (ui(I,J,K).lt.0) then
                  advx_eng(I,J,K) = ui(I,J,K) * (ti(I+1,J,K)-ti(I,J,K))/deltax_d + ti(I,J,K) * (ui(I+1,J,K)-ui(I,J,K))/deltax_d
               else if (ui(I,J,K).ge.0) then
                  advx_eng(I,J,K) = ui(I,J,K) * (ti(I,J,K)-ti(I-1,J,K))/deltax_d + ti(I,J,K) * (ui(I,J,K)-ui(I-1,J,K))/deltax_d
               endif
             else
               if (ui(I,J,K).ge.0) then
                  advx_eng(I,J,K) = ui(I,J,K) * (((3*ti(I,J,K))-(4*ti(I-1,J,K))+ti(I-2,J,K)))/(2*deltax_d) + ti(I,J,K) * (((3*ui(I,J,K))-(4*ui(I-1,J,K))+ui(I-2,J,K)))/(2*deltax_d)
               else if (ui(I,J,K).lt.0) then
                  advx_eng(I,J,K) = ui(I,J,K) * ((-ti(I+2,J,K)+(4*ti(I+1,J,K))-(3*ti(I,J,K))))/(2*deltax_d) + ti(I,J,K) * ((-ui(I+2,J,K)+(4*ui(I+1,J,K))-(3*ui(I,J,K))))/(2*deltax_d)
               endif
             endif

            ! cccccccccccccccccccccccccccccccccccccccccccccc                                        
            ! cccccccccccccccc Advection in Y cccccccccccccc                                             
            ! cccccccccccccccccccccccccccccccccccccccccccccc
            
             if ((J.eq.2).or.(K.eq.2).or.(J.eq.NY-1).or.(K.eq.NZ-1)) then
               if (vi(I,J,K).lt.0) then
                  advy_eng(I,J,K) = vi(I,J,K) * (ti(I,J+1,K)-ti(I,J,K))/deltay_d + ti(I,J,K) * (vi(I,J+1,K)-vi(I,J,K))/deltay_d 
               else if (vi(I,J,K).ge.0) then
                  advy_eng(I,J,K) = vi(I,J,K) * (ti(I,J,K)-ti(I,J-1,K))/deltay_d + ti(I,J,K) * (vi(I,J,K)-vi(I,J-1,K))/deltay_d
               endif
             else
               if (vi(I,J,K).ge.0) then
                  advy_eng(I,J,K) = vi(I,J,K) * (((3*ti(I,J,K))-(4*ti(I,J-1,K))+ti(I,J-2,K)))/(2*deltay_d) + ti(I,J,K) * (((3*vi(I,J,K))-(4*vi(I,J-1,K))+vi(I,J-2,K)))/(2*deltay_d)
               else if (vi(I,J,K).lt.0) then
                  advy_eng(I,J,K) = vi(I,J,K) * ((-ti(I,J+2,K)+(4*ti(I,J+1,K))-(3*ti(I,J,K))))/(2*deltay_d) + ti(I,J,K) * ((-vi(I,J+2,K)+(4*vi(I,J+1,K))-(3*vi(I,J,K))))/(2*deltay_d)
               endif
             endif

            ! cccccccccccccccccccccccccccccccccccccccccccccc                                                         
            ! cccccccccccccccc Advection in Z cccccccccccccc                                             
            ! cccccccccccccccccccccccccccccccccccccccccccccc                    

            if ((K.eq.2).or.(K.eq.NZ-1)) then
               if (wi(I,J,K).lt.0) then
                  advz_eng(I,J,K) = (wi(I,J,K)/H(I,J,3)) * (ti(I,J,K+1)-ti(I,J,K))/(za_eng(K+1) - za_eng(K)) + (ti(I,J,K)/H(I,J,3)) * (wi(I,J,K+1)-wi(I,J,K))/(za_eng(K+1) - za_eng(K))
               else if (wi(I,J,K).ge.0) then
                  advz_eng(I,J,K) = (wi(I,J,K)/H(I,J,3)) * (ti(I,J,K)-ti(I,J,K-1))/(za_eng(K) - za_eng(K-1)) + (ti(I,J,K)/H(I,J,3)) * (wi(I,J,K)-wi(I,J,K-1))/(za_eng(K) - za_eng(K-1))
               endif
            else
               if (wi(I,J,K).ge.0) then
                  advz_eng(I,J,K) = (wi(I,J,K)/H(I,J,3)) * (((3*ti(I,J,K))-(4*ti(I,J,K-1))+ti(I,J,K-2)))/(za_eng(K) - za_eng(K-2)) + (ti(I,J,K)/H(I,J,3)) * (((3*wi(I,J,K))-(4*wi(I,J,K-1))+wi(I,J,K-2)))/(za_eng(K) - za_eng(K-2))
               else if (wi(I,J,K).lt.0) then
                  advz_eng(I,J,K) = (wi(I,J,K)/H(I,J,3)) * (-ti(I,J,K+2)+(4*ti(I,J,K+1))-(3*ti(I,J,K)))/(za_eng(K+2) - za_eng(K)) + (ti(I,J,K)/H(I,J,3)) * ((-wi(I,J,K+2)+(4*wi(I,J,K+1))-(3*wi(I,J,K))))/(za_eng(K+2) - za_eng(K))
               endif
            endif

          ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc             
          ! ccccccccc Longitudonal/lateral ice thickness changes ccccccccc    
          ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    

          latx_eng(I,J,K) = (ti(I,J,K)/H(I,J,3))*ui(I,J,K)*((H(I+1,J,3)-H(I-1,J,3)) /(2*deltax_d))
          laty_eng(I,J,K) = (ti(I,J,K)/H(I,J,3))*vi(I,J,K)*((H(I,J+1,3)-H(I,J-1,3)) /(2*deltay_d))

          ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          ! cccccccccccccccccccc END OF ADVECTION SCHEME ccccccccccccccccccc
          ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  
          ! Calculate new debris concentration
  
           ti(I,J,K) = ti(I,J,K) + dt_eng * (-advx_eng(I,J,K) - advy_eng(I,J,K) + advz_eng(I,J,K) - latx_eng(I,J,K) - laty_eng(I,J,K))

          ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
          ! cccccccccccccccc UPDATE BOUNDARY CONDITIONS cccccccccccccccccccc 
          ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
                      
              ! Boundary conditions on top                                                                                          
               if (G(I,J).gt.0)then
                ti(I,J,NZ) = t_s(I,J)
                ti_new(I,J,NZ) = t_s(I,J)
               endif
               if (G(I,J).lt.0)then
                 ti(I,J,NZ) = ti(I,J,NZ-1)
                 ti_new(I,J,NZ) = ti_new(I,J,NZ-1)
               end if
              ! Boundary conditions on sides                                                                                                       
                ti_new(NX,J,K) = 0.0
                ti(NX,J,K) = 0.0
                ti_new(I,NY,K) = 0.0
                ti(I,NY,K) = 0.0
                ti_new(1,J,K) = 0.0
                ti(1,J,K) = 0.0
                ti_new(I,1,K) = 0.0
                ti(I,1,K) = 0.0
              ! Boundary conditions at bottom                                                             
                ti(I,J,1) = 0.0
                ti_new(I,J,1) = 0.0

            ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            ! cccccccccccccccc ADJUST FALSE CONCENTRATION cccccccccccccccccccc  
            ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
                     
               if ((sqrt(vx(I,J,1) * vx(I,J,1) + vy(I,J,1) * vy(I,J,1) + vz(I,J,1) * vz(I,J,1)).eq.0).and.(ti(I,J,K).gt.0))then
                    ti(I,J,K) = 0.
		    ti_new(I,J,K) = 0.
               end if

               if (H(I,J,3).eq.0)then
                    ti(I,J,K) = 0.
		    ti_new(I,J,K) = 0.
               end if

               if (ti(I,J,K) /= ti(I,J,K))then
                    ti(I,J,K) = 0.
                    ti_new(I,J,K) = 0.
               end if
 
               if(ti(I,J,K).lt.1e-3.or.ti(I,J,K).gt.1000)then
                  ti(I,J,K) = 0.
                  ti_new(I,J,K) = 0.
               end if

           end if
         end do
      end do
    end do

    ! Update boundary condition at the top             

     do J=1,NY
        do I=1,NX
           do K=1,NZ
               if (G(I,J).le.0)then
                ti(I,J,NZ) = ti(I,J,NZ-1)
                ti_new(I,J,NZ) = ti_new(I,J,NZ-1)
                ti(I,J,1) = 0.
                ti_new(I,J,1) = 0.
               else
                ti(I,J,NZ) = ti(I,J,NZ)
                ti_new(I,J,NZ) = ti_new(I,J,NZ)
                ti(I,J,1) = 0.
                ti_new(I,J,1) = 0.
               endif
            end do
        end do
     end do

     ! Set to zero if no ice

     do J=1,NY
        do I=1,NX
           do K=1,NZ
             if(H(I,J,3).eq.0.or.ti(I,J,K).lt.1e-3)then
                ti(I,J,K) = 0.
                ti_new(I,J,K) = 0.
             endif
          end do
       end do
    end do

  !-----------------------------------!    
  ! STEP 2: ANTI-DIFFUSION CORRECTION !
  !-----------------------------------!

  do iter_smolar = 1, n_iter_smolar ! Iterations for anti-diffusion correction
       
    do J=2,NY-1
       do I=2,NX-1
          do K=2,NZ-1

                ! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc            
                ! cccccccccccccc CALCULATE ANTI-DIFFUSION VELOCITIES AFTER SMOLARKIEWICZ cccccccccccc                                            
                ! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc                 
               
	        IF((H(I,J,3).gt.0.0).and.ti(I,J,K).gt.1e-3)then

                    ! X-DIRECTION 
                    advel_x_plus(I,J,K) = (abs(0.5*(ui(I,J,K) + ui(I+1,J,K))) * dx - (0.5*(ui(I,J,K) + ui(I+1,J,K)))**2 * dt_eng) * &
                    ((ti(I+1,J,K) - ti(I,J,K)) / ((ti(I+1,J,K) + ti(I,J,K) + epsilon_d) * dx)) &
                    - ((ui(I,J,K) + ui(I+1,J,K)) * (vi(I,J,K) + vi(I,J+1,K)) * dt_eng) * &
                    ((ti(I,J+1,K) - ti(I,J,K)) / ((ti(I,J+1,K) + ti(I,J,K) + epsilon_d) * dy)) &
                    - ((ui(I,J,K) + ui(I+1,J,K)) * (wi(I,J,K) + wi(I,J,K+1)) * dt_eng) * &
                    ((ti(I,J,K+1) - ti(I,J,K)) / ((ti(I,J,K+1) + ti(I,J,K) + epsilon_d) * (-H(I,J,3) * dzeta(K))))
                     advel_x_min(I,J,K) = (abs(0.5*(ui(I,J,K) + ui(I-1,J,K))) * dx - (0.5*(ui(I,J,K) + ui(I-1,J,K)))**2 * dt_eng) * &
                    ((ti(I,J,K) - ti(I-1,J,K)) / ((ti(I,J,K) + ti(I-1,J,K) + epsilon_d) * dx)) &
                    - ((ui(I,J,K) + ui(I-1,J,K)) * (vi(I,J,K) + vi(I,J-1,K)) * dt_eng) * &
                    ((ti(I,J,K) - ti(I,J-1,K)) / ((ti(I,J-1,K) + ti(I,J,K) + epsilon_d) * dy)) &
                    - ((ui(I,J,K) + ui(I-1,J,K)) * (wi(I,J,K) + wi(I,J,K-1)) * dt_eng) * &
                    ((ti(I,J,K) - ti(I,J,K-1)) / ((ti(I,J,K-1) + ti(I,J,K) + epsilon_d) * (-H(I,J,3) * dzeta(K))))
                    
                    ! Y-DIRECTION 
                    advel_y_plus(I,J,K) = (abs(0.5*(vi(I,J,K) + vi(I,J+1,K))) * dy - (0.5*(vi(I,J,K) + vi(I,J+1,K)))**2 * dt_eng) * &
                    ((ti(I,J+1,K) - ti(I,J,K)) / ((ti(I,J+1,K) + ti(I,J,K) + epsilon_d) * dy)) &
                    - ((ui(I,J,K) + ui(I+1,J,K)) * (vi(I,J,K) + vi(I,J+1,K)) * dt_eng) * &
                    ((ti(I+1,J,K) - ti(I,J,K)) / ((ti(I+1,J,K) + ti(I,J,K) + epsilon_d) * dx)) &
                    - ((vi(I,J,K) + vi(I,J+1,K)) * (wi(I,J,K) + wi(I,J,K+1)) * dt_eng) * &
                    ((ti(I,J,K+1) - ti(I,J,K)) / ((ti(I,J,K+1) + ti(I,J,K) + epsilon_d) * (-H(I,J,3) * dzeta(K))))
                     advel_y_min(I,J,K) = (abs(0.5*(vi(I,J,K) + vi(I,J-1,K))) * dy - (0.5*(vi(I,J,K) + vi(I,J-1,K)))**2 * dt_eng) * &   
                    ((ti(I,J,K) - ti(I,J-1,K)) / ((ti(I,J,K) + ti(I,J-1,K) + epsilon_d) * dy)) &
                    - ((ui(I,J,K) + ui(I-1,J,K)) * (vi(I,J,K) + vi(I,J-1,K)) * dt_eng) * &
                    ((ti(I,J,K) - ti(I-1,J,K)) / ((ti(I-1,J,K) + ti(I,J,K) + epsilon_d) * dx)) &
                    - ((vi(I,J,K) + vi(I,J-1,K)) * (wi(I,J,K) + wi(I,J,K-1)) * dt_eng) * &
                    ((ti(I,J,K) - ti(I,J,K-1)) / ((ti(I,J,K-1) + ti(I,J,K) + epsilon_d) * (-H(I,J,3) * dzeta(K))))
                    
                    ! Z-DIRECTION
                    advel_z_plus(I,J,K) = (abs(0.5*(wi(I,J,K) + wi(I,J,K+1))) * (-H(I,J,3) * dzeta(K)) - (0.5*(wi(I,J,K) + wi(I,J,K+1)))**2 * dt_eng) * &
                    ((ti(I,J,K+1) - ti(I,J,K)) / ((ti(I,J,K+1) + ti(I,J,K) + epsilon_d) * (-H(I,J,3) * dzeta(K)))) &
                    - ((ui(I,J,K) + ui(I+1,J,K)) * (wi(I,J,K) + wi(I,J,K+1)) * dt_eng) * &
                    ((ti(I+1,J,K) - ti(I,J,K)) / ((ti(I+1,J,K) + ti(I,J,K) + epsilon_d) * dx)) &
                    - ((vi(I,J,K) + vi(I,J+1,K)) * (wi(I,J,K) + wi(I,J,K+1)) * dt_eng) * &
                    ((ti(I,J+1,K) - ti(I,J,K)) / ((ti(I,J+1,K) + ti(I,J,K) + epsilon_d) * dy))
                     advel_z_min(I,J,K) = (abs(0.5*(wi(I,J,K) + wi(I,J,K-1))) * (-H(I,J,3) * dzeta(K)) - (0.5*(wi(I,J,K) + wi(I,J,K-1)))**2 * dt_eng) * &
                    ((ti(I,J,K) - ti(I,J,K-1)) / ((ti(I,J,K) + ti(I,J,K-1) + epsilon_d) * (-H(I,J,3) * dzeta(K)))) &
                    - ((ui(I,J,K) + ui(I-1,J,K)) * (wi(I,J,K) + wi(I,J,K-1)) * dt_eng) * &
                    ((ti(I,J,K) - ti(I-1,J,K)) / ((ti(I-1,J,K) + ti(I,J,K) + epsilon_d) * dx)) &
                    - ((vi(I,J,K) + vi(I,J-1,K)) * (wi(I,J,K) + wi(I,J,K-1)) * dt_eng) * &
                    ((ti(I,J,K) - ti(I,J-1,K)) / ((ti(I,J-1,K) + ti(I,J,K) + epsilon_d) * dy))
                    
                    ! ANTI-DIFFUSION VELOCITY
                    vel_diff_x(I,J,K) = 0.5 * (advel_x_plus(I,J,K) + advel_x_min(I,J,K))
                    vel_diff_y(I,J,K) = 0.5 * (advel_y_plus(I,J,K) + advel_y_min(I,J,K))
                    vel_diff_z(I,J,K) = 0.5 * (advel_z_plus(I,J,K) + advel_z_min(I,J,K))
                    vel_diff(I,J,K) = sqrt((vel_diff_x(I,J,K)**2)+(vel_diff_y(I,J,K)**2)+(vel_diff_z(I,J,K)**2))

                   ! cccccccccccccccccccccccccccccccccccccccccccccccccccccc   
                   ! cccccccccc CALCULATE ANTI-DIFFUSION FLUXES ccccccccccc   
                   ! cccccccccccccccccccccccccccccccccccccccccccccccccccccc     
               
                   ! X-DIRECTION
                   term1_debris_eng_x_plus(I,J,K) = 0.5 * (advel_x_plus(I,J,K) + abs(advel_x_plus(I,J,K))) * ti(I,J,K) + &
                          0.5 * (advel_x_plus(I,J,K) - abs(advel_x_plus(I,J,K))) * ti(I+1,J,K)
                   term1_debris_eng_x_min(I,J,K) = 0.5 * (advel_x_min(I,J,K) + abs(advel_x_min(I,J,K))) * ti(I-1,J,K) + &
                          0.5 * (advel_x_min(I,J,K) - abs(advel_x_min(I,J,K))) * ti(I,J,K)
                   term1_debris_eng_x(I,J,K) = 0.5 * (term1_debris_eng_x_plus(I,J,K) - term1_debris_eng_x_min(I,J,K)) / dx
 
                   ! Y-DIRECTION
                   term1_debris_eng_y_plus(I,J,K) = 0.5 * (advel_y_plus(I,J,K) + abs(advel_y_plus(I,J,K))) * ti(I,J,K) + &
                          0.5 * (advel_y_plus(I,J,K) - abs(advel_y_plus(I,J,K))) * ti(I,J+1,K)
                   term1_debris_eng_y_min(I,J,K) = 0.5 * (advel_y_min(I,J,K) + abs(advel_y_min(I,J,K))) * ti(I,J-1,K) + &
                          0.5 * (advel_y_min(I,J,K) - abs(advel_y_min(I,J,K))) * ti(I,J,K)
                   term1_debris_eng_y(I,J,K) = 0.5 * (term1_debris_eng_y_plus(I,J,K) - term1_debris_eng_y_min(I,J,K)) / dy

                   ! Z-DIRECTION
                   term1_debris_eng_z_plus(I,J,K) = 0.5 * (advel_z_plus(I,J,K) + abs(advel_z_plus(I,J,K))) * ti(I,J,K) + &
                          0.5 * (advel_z_plus(I,J,K) - abs(advel_z_plus(I,J,K))) * ti(I,J,K+1)
                   term1_debris_eng_z_min(I,J,K) = 0.5 * (advel_z_min(I,J,K) + abs(advel_z_min(I,J,K))) * ti(I,J,K-1) + &
                          0.5 * (advel_z_min(I,J,K) - abs(advel_z_min(I,J,K))) * ti(I,J,K)
                   term1_debris_eng_z(I,J,K) = 0.5 * (term1_debris_eng_z_plus(I,J,K) - term1_debris_eng_z_min(I,J,K)) / (H(I,J,3)*dzeta(K))

                else
                   
                      advel_x_plus(I,J,K) = 0.
                      advel_x_min(I,J,K) = 0.
                      advel_y_plus(I,J,K) = 0.
                      advel_y_min(I,J,K) = 0.
                      advel_z_plus(I,J,K) = 0.
                      advel_z_min(I,J,K) = 0.
                      vel_diff_x(I,J,K) = 0.
                      vel_diff_y(I,J,K) = 0.
                      vel_diff_z(I,J,K) = 0.
                      vel_diff(I,J,K) = 0.
                      term1_debris_eng_x_plus(I,J,K) = 0.
                      term1_debris_eng_y_plus(I,J,K) = 0.
                      term1_debris_eng_z_plus(I,J,K) = 0.
                      term1_debris_eng_x_min(I,J,K) = 0.
                      term1_debris_eng_y_min(I,J,K) = 0.
                      term1_debris_eng_z_min(I,J,K) = 0.
                      term1_debris_eng_x(I,J,K) = 0.
                      term1_debris_eng_y(I,J,K) = 0.
                      term1_debris_eng_z(I,J,K) = 0.

               endif ! endif H > 0 and ti > 0

         end do
      end do
   end do

  ! Do the correction step

   do J=2,NY-1
       do I=2,NX-1
          do K=2,NZ-1

               IF((H(I,J,3).gt.0.0).and.ti(I,J,K).gt.1e-3)then
             
               ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   
               ! cccccccccccccccccccc CORRECTION STEP ccccccccccccccccccccccccccc     
               ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc                           

                    ti(I,J,K) = ti(I,J,K) + dt_eng * (-term1_debris_eng_x(I,J,K) - term1_debris_eng_y(I,J,K) + term1_debris_eng_z(I,J,K))
                    ! Boundary conditions on top                                                        
                    if (G(I,J).ge.0)then
                       ti(I,J,NZ) = t_s(I,J)
                    endif
                    if (G(I,J).lt.0)then
                       ti(I,J,NZ) = ti(I,J,NZ-1)
                    end if
                    ! Boundary conditions on sides                                                                                     
                    ti(NX,J,K) = 0.0
                    ti(I,NY,K) = 0.0
                    ti(1,J,K) = 0.0
                    ti(I,1,K) = 0.0
                    ! Boundary conditions at bottom                                            
                    ti(I,J,1) = 0.0
                    ! Adjust false concentrations                
                    if ((sqrt(vx(I,J,1) * vx(I,J,1) + vy(I,J,1) * vy(I,J,1) + vz(I,J,1) * vz(I,J,1)).eq.0).and.(ti(I,J,K).gt.0))then
                       ti(I,J,K) = 0.
                    end if
                    if (H(I,J,3).eq.0)then
                       ti(I,J,K) = 0.
                    end if
                    if (ti(I,J,K) /= ti(I,J,K))then
                       ti(I,J,K) = 0.
                    end if
                    if(ti(I,J,K).lt.1e-3.or.ti(I,J,K).gt.1000)then
                       ti(I,J,K) = 0.
                    end if
                    ti_new(I,J,K) = ti(I,J,K)

                 endif
               
         end do
      end do
   end do

 enddo ! enddo smolar iterations 
   
   ! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
   ! ccccccccccccccccccccccc SAVE FOR MASS CONSERVATION cccccccccccccccccccccc    
   ! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    

   ! Update boundary conditions
   
    do J=1,NY
      do I=1,NX
          do K=1,NZ
               if (G(I,J).ge.0)then
                ti(I,J,NZ) = t_s(I,J)
                ti_new(I,J,NZ) = t_s(I,J)
                ti(I,J,1) = 0.
                ti_new(I,J,1) = 0.
               endif
               if (G(I,J).lt.0)then
                 ti(I,J,NZ) = ti(I,J,NZ-1)
                 ti_new(I,J,NZ) = ti_new(I,J,NZ-1)
                 ti(I,J,1) = 0.
                 ti_new(I,J,1) = 0.
               end if
         end do
      end do
    end do

   ! Adjust false concentrations
    
   do J=1,NY
      do I=1,NX
         do K=1,NZ
               if ((sqrt(vx(I,J,1) * vx(I,J,1) + vy(I,J,1) * vy(I,J,1) + vz(I,J,1) * vz(I,J,1)).eq.0).and.(ti(I,J,K).gt.0))then
                    ti(I,J,K) = 0.
                    ti_new(I,J,K) = 0.
               end if
               if (H(I,J,3).eq.0)then
                    ti(I,J,K) = 0.
                    ti_new(I,J,K) = 0.
               end if
               if (ti(I,J,K) /= ti(I,J,K))then
                    ti(I,J,K) = 0.
                    ti_new(I,J,K) = 0.
               end if
               if(ti(I,J,K).lt.1e-3.or.ti(I,J,K).gt.1000)then
                    ti(I,J,K) = 0.
                    ti_new(I,J,K) = 0.
               end if
          enddo
       enddo
   enddo

   ! Update sinks (debris meltout)
   
   do J=1,NY
     do I=1,NX
        if (G(I,J).lt.0.and.H(I,J,3).gt.0.and.ti(I,J,NZ-1).gt.1e-3)then
             sink_debris(I,J,1) = (-(ti(I,J,NZ-1)*G(I,J)*dx*dx)) / ((1-phi_debris)*rho_debris*dx*dx) ! in m/yr
             mass_out(I,J) = -(ti(I,J,NZ-1)*G(I,J)*dx*dx) ! in kg/yr
        else
             sink_debris(I,J,1) = 0.
             mass_out(I,J) = 0.
        endif
        if(I.eq.1.or.I.eq.NX.or.J.eq.1.or.J.eq.NY.or.H(I,J,3).eq.0.or.sink_debris(I,J,1).lt.1e-5)then
             sink_debris(I,J,1) = 0
             mass_out(I,J) = 0.
         endif
     end do
   end do

   ! Save terms for mass conservation
   
   do J=2,NY-1
      do I=2,NX-1
         do K=2,NZ-1
             actual_mass_debris_eng(I,J,K) = (H(I,J,3)/nz)*(deltax_d*deltax_d*ti(I,J,K))
             new_debris_mass_eng(I,J,NZ-1) = (mass_in(I,J))*dt_eng - (mass_out(I,J))*dt_eng
             expected_debris_mass_eng(I,J,K) = expected_debris_mass_eng(I,J,K) + new_debris_mass_eng(I,J,K)
          end do
     end do
   end do

   do J=1,NY
      do I=1,NX
        do K=1,NZ
            if(I.eq.1.or.I.eq.NX.or.J.eq.1.or.J.eq.NY.or.H(I,J,3).eq.0.or.ti(I,J,K).lt.1e-3)then
              actual_mass_debris_eng(I,J,K) = 0.
              new_debris_mass_eng(I,J,K) = 0.
              expected_debris_mass_eng(I,J,K) = 0.
            end if
        end do
     end do
   end do

  ! cccccccccccccccccccccccccccccccccccccccccc     
  ! cccccccccc END OF THE LOOP ccccccccccccccc                   
  ! cccccccccccccccccccccccccccccccccccccccccc             

  ! Upload debris time step                     

  it_deng = it_deng + 1

  end do ! end do main loop for one iteration                                                      
  
  ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! ccccccccccccccccccc END OF MAIN ADVECTION MODULE ccccccccccccccc
  ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  
  !----------------------------------!
  ! STEP 3: ENSURE MASS CONSERVATION !
  !----------------------------------!

  ! Ensure mass conservation

  actual_mass_hist_eng = sum(actual_mass_debris_eng)                      ! Actual debris mass after 1/dt iterations
  expected_mass_hist_eng = sum(expected_debris_mass_eng(:,:,NZ-1))        ! Expected debris mass after 1/dt iterations
  expected_mass_before_eng = expected_mass_hist_eng                       ! Save this expected debris mass before correction
  mass_ratio_hist_eng = (expected_mass_hist_eng/actual_mass_hist_eng)     ! Mass ratio before correction

  if(TIME.gt.tdebris+1)then
 
   do iter_mc = 1, max_iters_mc

    ! Initial guess of melt-out mass
    F_out_new = 0.0
    do i = 1, NX
        do j = 1, NY
            if(mass_out(I,J).ne.0)then
              F_out_new = F_out_new + (mass_ratio_hist_eng(1) * mass_out(I,J))
            endif
        end do
    end do

    ! Update mass ratio with initial guess of new melt-out mass
    expected_mass_hist_eng_new = sum(actual_mass_debris_eng_before) + sum(mass_in) - F_out_new
    mass_ratio_hist_eng_new = expected_mass_hist_eng_new / actual_mass_hist_eng

    ! Check for convergence
    mass_ratio_hist_eng(1) = mass_ratio_hist_eng_new(1)                ! Mass ratio after correction
    expected_mass_hist_eng(1) = expected_mass_hist_eng_new(1)          ! Total expected mass after correction
    if (abs(mass_ratio_hist_eng_new(1) - mass_ratio_hist_eng(1)) < 1.0E-4) exit
    
   end do ! end do iterations mass conservation

  endif  ! endif tdebris+1
  
  ! Multiply with mass conservation factor

  if (mass_ratio_hist_eng(1).gt.0.0)then
    do J=2,NY-1
       do I=2,NX-1
          do K=2,NZ-1
             ti(I,J,K) = ti(I,J,K)*mass_ratio_hist_eng(1)
          end do
       end do
     end do
  endif

  ! Check if it worked

  do J=2,NY-1
      do I=2,NX-1
          do K=2,NZ-1
           actual_mass_debris_eng(I,J,K) = (H(I,J,3)/nz)*(deltax_d*deltax_d*ti(I,J,K))
          end do
      end do
  end do

  actual_mass_hist_eng = sum(actual_mass_debris_eng)                   ! Mass after correction
  mass_ratio_hist_eng = (expected_mass_hist_eng/actual_mass_hist_eng)  ! New mass ratio... should be 1 now...

  ! Update the expected debris mass

    do J=1,NY
       do I=1,NX
          do K=1,NZ
             expected_debris_mass_eng(I,J,K) = expected_debris_mass_eng(I,J,K) * (expected_mass_hist_eng(1)/expected_mass_before_eng(1))
          end do
       end do
    end do
  
  ! Update boundary conditions 
  
  do J=1,NY
      do I=1,NX
          do K=1,NZ
               if (G(I,J).ge.0)then
                ti(I,J,NZ) = t_s(I,J)
                ti_new(I,J,NZ) = t_s(I,J)
                ti(I,J,1) = 0.
                ti_new(I,J,1) = 0.
               endif
               if (G(I,J).lt.0)then
                ti(I,J,NZ) = ti(I,J,NZ-1)
                ti_new(I,J,NZ) = ti_new(I,J,NZ-1)
                ti(I,J,1) = 0.
                ti_new(I,J,1) = 0.
               end if
           end do
      end do
  end do

  ! Update sources and sinks

   do J=1,NY
        do I=1,NX
              mass_in(I,J) = (inoutdebris_acc(I,J)*(1-phi_debris)*rho_debris*deltax_d*deltax_d)
        end do
   end do
     
   do J=1,NY
       do I=1,NX
            if (G(I,J).lt.0.and.H(I,J,3).gt.0.and.ti(I,J,NZ-1).gt.1e-3)then
               sink_debris(I,J,1) = (-(ti(I,J,NZ-1)*G(I,J)*dx*dx)) / ((1-phi_debris)*rho_debris*dx*dx)
               mass_out(I,J) = -(ti(I,J,NZ-1)*G(I,J)*dx*dx)
            else
               sink_debris(I,J,1) = 0.
               mass_out(I,J) = 0.
            endif
            if(I.eq.1.or.I.eq.NX.or.J.eq.1.or.J.eq.NY.or.H(I,J,3).eq.0.or.sink_debris(I,J,1).lt.1e-5)then
               sink_debris(I,J,1) = 0
               mass_out(I,J) = 0.
            endif
       end do
   end do

   write(*,*),'----- after MC correction -----'
   write(*,*),'mass in engl = ',sum(mass_in)
   write(*,*),'mass out engl = ',sum(mass_out)
   write(*,*),'mass retreat engl = ',sum(mass_engl_out_retreat)
   write(*,*),'mass_ratio_hist_eng after = ',mass_ratio_hist_eng(1)
   write(*,*),'actual_mass_hist_eng = ',actual_mass_hist_eng
   write(*,*),'actual mass hist eng before = ',sum(actual_mass_debris_eng_before)
   write(*,*),'expected_mass_hist_eng = ',expected_mass_hist_eng
   
yr_deng = yr_deng + 1

end do ! end do loop full year

ti_new = ti

! cccccccccccccccccccccccccccccccccccccccccc   
! cccccccc PREPARE FOR OUTPUT cccccccccccccc    
! cccccccccccccccccccccccccccccccccccccccccc    

! Total mass that melts out near the surface   

do J=1,NY
   do I=1,NX
      if(G(I,J).le.0)then
         mass_melt_out(I,J) = mass_out(I,J)
      else
         mass_melt_out(I,J) = 0.
      endif
      if(I.eq.1.or.I.eq.NX.or.J.eq.1.or.J.eq.NY.or.H(I,J,3).eq.0.or.sink_debris(I,J,1).lt.1e-5)then
         mass_melt_out(I,J) = 0.
      endif
   enddo
enddo

sum_mass_meltout = mass_melt_out
ti_hist = ti
ti_add = vel_diff
ti_add(:,:,1)=dc_ddt_yearly
ti_add(:,:,2)=mass_engl_out_retreat

      end subroutine debris_cover_main_englacial
