# 3D englacial and supraglacial debris advection model

## General information

The model code is a 3D englacial and supraglacial debris transport model. It is written in Fortran90. The model is used in Verhaegen and Huybrechts (2025), where it is applied to a synthetic, idealized glacier. It builds further upon earlier work in Anderson
& Anderson (2016) and Verhaegen et al. (2020).

## Model structure

The model code consists of the following elemennts:

a) debris_cover_main_supraglacial.f90: supraglacial debris advection model

This model code solves for the debris thickness, including advection by surface ice flow, melt-out of englacial debris, gravitational transport, and off-glacier debris removal at the margin.
The advection submodule uses a second-order upwind transport scheme, a Smolarkiewicz-based anti-diffusion correction and also uses a mass conservation correction.

b) debris_cover_main_englacial.f90: englacial debris advection model

This model advects debris mass concentration within the ice column, using a second-order upwind scheme + Smolarkiewicz anti-diffusion correction and also uses a mass conservation correction. It provides the melt-out mass that is used as input into the supraglacial debris model.

## Model input and requirements

The model requieres certain input variables and parameters.

### Input variables

The following input variables are needed to run the model:

- H (ice thickness - m)
- HT (surface elevation - m)
- G (surface mass balance - m i.e. yr-1)
- mask (ice mask - 0 or 1)
- VX (horizontal velocity x - m yr-1)
- VY (horizontal velocity y - m yr-1)
- VZ (vertical velocity z - m yr-1)
- DHTDT_yearly (surface elevation change rate - m yr-1)
- DHDT_yearly (ice thickness change rate - m yr-1)
- margin (marginal pixels of glacier - 0 or 1)

### Debris input locations

Debris input locations are X,Y positions of the pixels within the domain and should be set in debris_cover_main_englacial.f90 (if in the accumulation zone) or debris_cover_main_supraglacial.f90 (if in the ablation zone).

### Example input parameter values

The following parameter values were inserted to obtain the debris-covered reference glacier in Verhaegen & Huybrechts (2025):

      INTEGER, PARAMETER   :: NY=178                               ! Number of pixels in y (-)
      INTEGER, PARAMETER   :: NX=300                               ! Number of pixels in x (-)
      INTEGER, PARAMETER   :: NZ=21                                ! Number of layers in z (-)
      INTEGER, PARAMETER   :: DX=25                                ! Spatial resolution ice flow model x (m)
      INTEGER, PARAMETER   :: DY=DX                                ! Spatial resolution ice flow model y (m)
      integer,parameter    :: tdebris=5                            ! Time of debris source release (yr)                                                                                                                                  
      integer,parameter    :: deltax_d=DX                          ! Spatial resolution debris model x (m)                                                                                                                                       
      integer,parameter    :: deltay_d=deltax_d                    ! Spatial resolution debris model y (m)                                                                                                                                       
      real,parameter       :: deltat_d=0.005                       ! Time step supraglacial debris model (yr)                                                                                                                                    
      real,parameter       :: dt_eng=0.005                         ! Time step englacial debris model (yr)                                                                                                                               
      integer,parameter    :: dtdiag_d=1                           ! Diagnostic time step debris models (yr)                                                                                                                         
      integer,parameter    :: nyears_d=1                           ! Number of years supragglacial debris (yr)                                                                                                                              
      integer,parameter    :: nyears_deng=1                        ! Number of years englacial debris (yr)                                                                                                                                   
      integer,parameter    :: numberofyears=1                      ! Number of years supraglacial debris (yr)                                                                                                                               
      integer,parameter    :: numberofyears_eng=1                  ! Number of years englacial debris (yr)                                                                                                                                  
      real,parameter       :: deposition_mass_debris_abl=0.0       ! Deposition mass in kg/yr (ABLATION AREA)                                                                                                                        
      real,parameter       :: deposition_mass_debris_acc=2000000.  ! Deposition mass in kg/yr (ACCUMULATION AREA)                                                                                                                    
      real,parameter       :: phi_debris=0.3                       ! Debris cover porosity (-)                                                                                                                                        
      real,parameter       :: rho_debris=2650                      ! Debris cover density (kg m-3)                                                                                                                                           
      real,parameter       :: epsilon_d = 1e-15                    ! Ensuring finite values (-)                                                                                                                                        
      real,parameter       :: n_iter_smolar = 1                    ! Number of iterations for numerical diffusion correction (-)                                                                                                        
      real,parameter       :: c_deb_slope = 0.0001                 ! Debris removal efficiency parameter (Pa-1 m2 yr-1)                                                                                                                            
      real,parameter       :: max_iters_mc = 10                    ! For mass conservation check (-)                                                                                                                                    
      real,parameter       :: lchar_d = 25                         ! Characteristic marginal length scale for debris removal (m)                                                                                                              


## References

Anderson, L. S. and Anderson, R. S.: Modeling debris-covered glaciers: response to steady debris deposition, The Cryosphere, 10, 1105–1124, https://doi.org/10.5194/tc-10-1105-2016, 2016.

Verhaegen, Y., Huybrechts, P., Rybak, O., and Popovnin, V. V.: Modelling the evolution of Djankuat Glacier, North Caucasus, from 1752 until 2100 CE, The Cryosphere, 14, 4039–4061, https://doi.org/10.5194/tc-14-4039-2020, 2020.

## Citation

If you use this model, please cite:

Verhaegen, Y. & Huybrechts, P. (2025). Coupling debris transport to 3D higher-order ice flow dynamics to model the behavior and climate change response of debris-covered glaciers. Journal of Geophysical Research – Earth Surface.

The corresponding Zenodo folder can be found here: http://doi.org/10.5281/zenodo.17693004.

