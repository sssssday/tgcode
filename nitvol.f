      subroutine nitvol

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine estimates daily mineralization (NH3 to NO3)
!!    and volatilization of NH3

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name          |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    curyr         |none          |current year of simulation
!!    hru_dafr(:)   |km**2/km**2   |fraction of watershed area in HRU
!!    ihru          |none          |HRU number
!!    nyskip        |none          |number of years to skip output
!!                                 |summarization and printing
!!    sol_fc(:,:)   |mm H2O        |amount of water available to plants in soil
!!                                 |layer at field capacity (fc - wp)
!!    sol_nh3(:,:)  |kg N/ha       |amount of nitrogen stored in the ammonium
!!                                 |pool in soil layer
!!    sol_nly(:)    |none          |number of layers in soil profile
!!    sol_no3(:,:)  |kg N/ha       |amount of nitrogen stored in the
!!                                 |nitrate pool in soil layer
!!    sol_st(:,:)   |mm H2O        |amount of water stored in the soil layer
!!                                 |on any given day (less wp water)
!!    sol_tmp(:,:)  |deg C         |daily average temperature of soil layer
!!    sol_wpmm(:,:) |mm H20        |water content of soil at -1.5 MPa (wilting
!!                                 |point)
!!    sol_z(:,:)    |mm            |depth to bottom of soil layer
!!    wshd_nitn     |kg N/ha       |average annual amount of nitrogen moving
!!                                 |from the NH3 to the NO3 pool by
!!                                 |nitrification in the watershed
!!    wshd_voln     |kg N/ha       |average annual amount if nitrogen lost by
!!                                 |ammonia volatilization in watershed
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name          |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    sol_nh3(:,:)  |kg N/ha       |amount of nitrogen stored in the ammonium
!!                                 |pool in soil layer
!!    sol_no3(:,:)  |kg N/ha       |amount of nitrogen stored in the
!!                                 |nitrate pool in soil layer
!!    wshd_nitn     |kg N/ha       |average annual amount of nitrogen moving
!!                                 |from the NH3 to the NO3 pool by
!!                                 |nitrification in the watershed
!!    wshd_voln     |kg N/ha       |average annual amount if nitrogen lost by
!!                                 |ammonia volatilization in watershed
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    akn         |
!!    akv         |
!!    cecf        |none          |volatilization CEC factor
!!    dmidl       |
!!    dpf         |
!!    j           |none          |HRU number
!!    k           |none          |counter (soil layer)
!!    rnit        |kg N/ha       |amount of nitrogen moving from the NH3 to the
!!                               |NO3 pool (nitrification) in the layer
!!    rnv         |
!!    rvol        |kg N/ha       |amount of nitrogen lost from the NH3 pool due
!!                               |to volatilization
!!    sw25        |
!!    swf         |
!!    swwp        |
!!    tf          |
!!    xx          |
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Exp, Max

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~
      use parm

      integer :: j, k
      real :: sw25, swwp, swf, xx, dmidl, dpf, akn, akv, rnv, rnit, rvol
      real :: tf 
      real :: cecf = 0.15
      real :: turnovfrac, N2Oadjust_wp, N2Oadjust_fc
      real :: dDO_wp, dDO_fc
      real :: swcfrac_fc
      real :: swcfrac_wp
      real :: swcfrac_dDO
      real :: porespace   
      real :: porespace_dDO        
      real :: wfpslyr_fc
      real :: wfpslyr_wp
      real :: wfpslyr_dDO
               

      real :: dDO
      
      real :: rnit_to_n2o
        
      

      j = 0
      j = ihru 
      
      N2Oadjust_wp = 0.002
      N2Oadjust_fc = 0.015
      dDO_wp = 0.
      dDO_fc = 0.
      
      
      swcfrac_fc = 0.
      swcfrac_wp = 0.
      swcfrac_dDO = 0.
      porespace = 0.   
      porespace_dDO = 0.        
      wfpslyr_fc = 0.
      wfpslyr_wp = 0.
      wfpslyr_dDO = 0.
               
      dDO_fc = 0.
      dDO_wp = 0.
      dDO = 0.
      rnit_to_n2o = 0.
      

      do k = 1, sol_nly(j)
        tf = 0.
        tf = .41 * (sol_tmp(k,j) - 5.) / 10.
        
        if (k == 1)cal_temp(11) = 0

        if (sol_nh3(k,j) > 0. .and. tf >= 0.001) then
          sw25 = 0.
          swwp = 0.
          sw25 = sol_wpmm(k,j) + 0.25 * sol_fc(k,j)
          swwp = sol_wpmm(k,j) + sol_st(k,j)
          if (swwp < sw25) then
            swf = 0.
            swf = (swwp - sol_wpmm(k,j)) /(sw25 - sol_wpmm(k,j))
          else
            swf = 1.
          endif

          if (k == 1) then
            xx = 0.
          else
            xx = 0.
            xx = sol_z(k-1,j)
          endif


        	if (k == 1) then
	       sol_thick = sol_z(k,j)
	    
	      else	
	       sol_thick = sol_z(k,j) - sol_z(k-1,j)
	    
	     end if
      

          dmidl = 0.
          dpf = 0.
          akn = 0.
          akv = 0.
          rnv = 0.
          rnit = 0.
          rvol = 0.
          dmidl = (sol_z(k,j) + xx) / 2.
          dpf = 1. - dmidl / (dmidl + Exp(4.706 - .0305 * dmidl))
          akn = tf * swf
          akv = tf * dpf * cecf
          rnv = sol_nh3(k,j) * (1. - Exp(-akn - akv))
          rnit = 1. - Exp(-akn)
          rvol = 1. - Exp(-akv)

          !! calculate nitrification (NH3 => NO3)
	    !! apply septic algorithm only to active septic systems
          if(k/=i_sep(j).or.isep_opt(j)/= 1) then  ! J.Jeong for septic, biozone layer
             if (rvol + rnit > 1.e-6) then
               rvol = rnv * rvol / (rvol + rnit)
               rnit = rnv - rvol
               if (rnit < 0.) rnit = 0.
               
               
              !! sol_nh3(k,j) = Max(1.e-6, sol_nh3(k,j) - rnit)
               
             rnit = Min(rnit, sol_nh3(k,j))     !! qichun
             sol_nh3(k,j) = sol_nh3(k,j)- rnit  !! qicun
             
 !! n2o emission from nitrified N
 

            swcfrac_fc = (sol_fc(k, j)+sol_wpmm(k,j)) / sol_thick
            swcfrac_wp = sol_wpmm(k,j)/ sol_thick
            swcfrac_dDO = (sol_fc(1, j)+sol_wpmm(1,j)) / sol_z(1,j)
            porespace = 1.0 - sol_bd(k,j)/2.56      
            porespace_dDO = 1.0 - sol_bd(1,j)/2.56           
            wfpslyr_fc = swcfrac_fc/porespace
            wfpslyr_wp = swcfrac_wp/porespace
            wfpslyr_dDO = swcfrac_dDO/porespace_dDO
               
            dDO_fc = diffusiv(swcfrac_fc, sol_bd(k,j), wfpslyr_fc)
            dDO_wp =  diffusiv(swcfrac_wp, sol_bd(k,j), wfpslyr_wp)
            dDO =  diffusiv(swcfrac_dDO, sol_bd(1,j), wfpslyr_dDO)
        
         
           
             
           if (rnit > 4.) rnit = 4.
           if (rnit > 1.0E-30) then
            turnovfrac = (N2Oadjust_wp - N2Oadjust_fc) /
     &                 (dDO_wp - dDO_fc) *
     &                  (dDO - dDO_wp) + N2Oadjust_wp
            turnovfrac = max(turnovfrac, N2Oadjust_wp)
            turnovfrac = min(turnovfrac, N2Oadjust_fc)
            rnit_to_n2o = rnit * turnovfrac
            rnit = rnit - rnit_to_n2o
            N2O(j)  = N2O(j) + rnit_to_n2o*0.1  !! unit conversion
            
            if (k==1)cal_temp(2) = 0. 
            cal_temp(2) = cal_temp(2)+ rnit_to_n2o*0.1
            
            end if
             
             
             endif
             if (sol_nh3(k,j) < 0.) then
               rnit = rnit + sol_nh3(k,j)
               sol_nh3(k,j) = 0.
             endif
             
             sol_no3(k,j) = sol_no3(k,j) + rnit
             

             !! calculate ammonia volatilization
             
             rvol = Min(rvol, sol_nh3(k,j))   !! qichun
             !!rvol = 0.
             sol_nh3(k,j) = sol_nh3(k,j)- rvol !! qicun
             
             if (sol_nh3(k,j) < 0.) then
               rvol = rvol + sol_nh3(k,j)
               sol_nh3(k,j) = 0.
             endif

             !! summary calculations
             if (curyr > nyskip) then
               wshd_voln = wshd_voln + rvol * hru_dafr(j)
               wshd_nitn = wshd_nitn + rnit * hru_dafr(j)
             end if
          end if
        end if

      cal_temp(11) = cal_temp(11) +  rvol

      end do
      
      


      return
      end