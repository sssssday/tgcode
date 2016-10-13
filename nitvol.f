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
      !!implicit none
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
      real :: no_n2o_ratio
      real :: tem_no
      real :: krain_no  !! rain impacts on n2o emission
      real :: canopy_red
      real :: no_adsorb  
      real :: rel_wc
      real :: avgstemp
      real :: fNwfps
      real :: fNsoilt
      real :: fNph 
      real :: fNnh4
      real :: abiotic
      real :: absoluteMaxRate
      real :: grams_soil
      real :: nh4_conc 
      real :: swclimit
      real :: deltamin
      real :: maxt
      real :: base_flux 
      real :: MaxRate !! faction of ammonia that is nitrified = 0.15;
      real :: sol_thick
      integer :: mon
      
      real :: totnh4
      
      !! define functions
      
      real :: f_gen_poisson_density
      
      
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
      no_n2o_ratio = 0.
      tem_no = 0.
      krain_no  = 0.
      canopy_red = 0.
      no_adsorb = 0.
      
      
       rel_wc = 0.
       avgstemp = 0.
       fNwfps = 0.
       fNsoilt = 0.
       fNph  = 0.
       fNnh4 = 0.
       abiotic = 0.
       absoluteMaxRate = 0.
       grams_soil = 0.
       nh4_conc = 0. 
       swclimit = 0.
       deltamin = 0.
       maxt = 0.
       base_flux = 0.
       MaxRate  = 0. !! faction of ammonia that is nitrified = 0.15;
       sol_thick = 0.
       mon  = 0
       totnh4 = 0.
      

      do k = 1, sol_nly(j)
        tf = 0.
        tf = 0.41 * (sol_tmp(k,j) - 5.) / 10.
        totnh4 = totnh4 + sol_nh3(k,j)
        

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
          cal_temp(2) = 0
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
          dpf = 1. - dmidl / (dmidl + Exp(4.706 - 0.0305 * dmidl))
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
               

             
             
              rnit = Min(rnit, 4.* (sol_nh3(k,j)/totnh4))   !! need to double check
              
              
              rnit = Min(rnit, sol_nh3(k,j))     !! qichun
             
             
              
             !!sol_nh3(k,j) = sol_nh3(k,j)- rnit  !! qicun
             
             
             
             
 !! nitrificaiton algorithm from Daycent~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
      rnit = 0.
 
      grams_soil = sol_bd(k,j)* sol_thick/1000. * 1000000. !! g/soil /m2
       
      nh4_conc =  sol_nh3(k,j) / 10. / grams_soil* 1000000. !! ppm nh4
      MaxRate  = 0.15  !! original is 0.15
      
      !! calculate water factor
      deltamin = 0.042  !! 0.042   in daycent this value ranges from 0.012 to 0.062
      
      swclimit = sol_wpmm(k,j)/sol_thick - deltamin   
      
      if (swclimit < 0.) swclimit = 0.
      
      
      rel_wc = ((sol_st(k, j)+sol_wpmm(k,j)) / sol_thick - swclimit) / 
     & ((sol_fc(k, j)+sol_wpmm(k,j)) / sol_thick -swclimit)  
      
      
       if (rel_wc < 0.) rel_wc = 0.
       if (rel_wc > 1.) rel_wc = 1.
       
       fNwfps = 1.0/(1.0 + 30.0 * exp(-9.0 * rel_wc))
       
      
     
      
      !! calculate temperature factor  
      
      
      !! change later
      !! i = 1
      do mon = 1,12
       if (tmpmx(mon,1) > maxt) maxt = tmpmx(mon,1)
      end do 
      
      
            
 
      if (maxt .ge. 35.0) then
      
      
         fNsoilt = f_gen_poisson_density(sol_tmp(k,j), maxt)
            
      else
         
         fNsoilt = f_gen_poisson_density(sol_tmp(k,j) + (35. - maxt), 
     &     maxt)
      
             
      endif
      
   !! calculate soil pH factor      
      
   
      fNph = 0.56 + (1./ 3.14159) * 
     &  atan(3.1415926 * 0.45 * (sol_ph(k,j) - 5))
 
   !! calculate ammonia factor     
   
       fNnh4 = 1.0 - exp(-0.0105 * nh4_conc)
      
   !! calculate nitrification rate
       base_flux = 0.1/10000.0 !! origianl is 0.1
       abiotic = max(fNwfps * fNsoilt, 0.03)
   
       absoluteMaxRate = 
     & min(0.4*(sol_nh3(k,j)/totnh4), sol_nh3(k,j)*0.1 * MaxRate) !! here is g N /m2   
     
       rnit = 10. * (absoluteMaxRate * fNph * abiotic  +  base_flux ) !!change unit to kg/ha
       
       
        rnit = Min(rnit, sol_nh3(k,j))     !! qichun
       cal_temp(4) = cal_temp(4) + rnit
       
       sol_nh3(k,j) = sol_nh3(k,j) - rnit

    !!  if (sol_nh3(k,j)> rnit)then
    !!    sol_nh3(k,j) = sol_nh3(k,j) - rnit
    !!   else 
    !!    rnit = sol_nh3(k,j)
   !!     sol_nh3(k,j) = 0.
   !!   end if
 
 !! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
             
             
 !! n2o emission from nitrified N
 

            swcfrac_fc = (sol_fc(k, j)+sol_wpmm(k,j)) / sol_thick
            swcfrac_wp = sol_wpmm(k,j)/ sol_thick
            swcfrac_dDO = (sol_fc(1, j)+sol_wpmm(1,j)) / sol_z(1,j)
            porespace = 1.0 - sol_bd(k,j)/2.56      
            porespace_dDO = 1.0 - sol_bd(1,j)/2.56           
            wfpslyr_fc = swcfrac_fc/porespace
            wfpslyr_wp = swcfrac_wp/porespace
            wfpslyr_dDO = (sol_st(2, j)+sol_wpmm(2,j)) / sol_z(2,j)
     &        /porespace_dDO
               
            dDO_fc = diffusiv(swcfrac_fc, sol_bd(k,j), wfpslyr_fc)
            dDO_wp =  diffusiv(swcfrac_wp, sol_bd(k,j), wfpslyr_wp)
            dDO =  diffusiv(swcfrac_dDO, sol_bd(1,j), wfpslyr_dDO)
        
         
           
             
         
           if (rnit > 1.0E-30) then
               turnovfrac = (N2Oadjust_wp - N2Oadjust_fc) /
     &                 (dDO_wp - dDO_fc) *
     &                  (dDO - dDO_wp) + N2Oadjust_wp
                turnovfrac = max(turnovfrac, N2Oadjust_wp)
                turnovfrac = min(turnovfrac, N2Oadjust_fc)
                rnit_to_n2o = rnit * turnovfrac
            
                 
                if (rnit > rnit_to_n2o)then 
                    rnit = rnit - rnit_to_n2o
                else
                rnit_to_n2o = rnit
                rnit = 0.
                end if    
            
  !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !! calculated NO from N2O
               krain_no = 1.0
               no_n2o_ratio = 8.0 + (18.0*atan(0.75*3.1415926*(
     &         10.*dDO-1.86))) /3.1415926   
         
          !!print*, 'no_n2o_ratio', no_n2o_ratio
         
             if (idplt(j) == 4) no_n2o_ratio = no_n2o_ratio* 0.5
            tem_no = no_n2o_ratio * rnit_to_n2o * krain_no 
            
            !!tem_no = 0.
            
            if (tem_no .le. rnit) then
            
                rnit = rnit - tem_no 
                    
            else
                
                tem_no = rnit + min ( sol_nh3(k,j), 
     &          (tem_no - rnit))
         
                sol_nh3(k,j) = sol_nh3(k,j) - min ( sol_nh3(k,j), 
     &          (tem_no - rnit))
         
                rnit = 0.
                
            endif    
              
  !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~       
                  
         
            
            if (rnit_to_n2o < 1.0E-30) rnit_to_n2o = 0.
            cal_temp(2) = cal_temp(2) + rnit_to_n2o*0.1
            
            N2O(j)  = N2O(j) + rnit_to_n2o*0.1  !! unit conversion
            
            if (tem_no < 1.0E-30) tem_no = 0.
            NO(j)  = NO(j) + tem_no * 0.1
            
            
            
            
            endif
             
             
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

     

      end do
      
      

      !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
        !! reduce no emission by plant obsorption
         
         if (laiday(j) > 0.) then
            if (laiday(j) > 8.) then
                canopy_red = 0.4428
            else
                canopy_red = (0.0077 * laiday(j)**2.0  -0.13 * laiday(j)
     &         + 0.99)
            
            endif   
         
           no_adsorb = NO(j) * (1.-canopy_red)     
           plantn(j) = plantn(j) + no_adsorb 
           NO(j) = NO(j) - no_adsorb 
         
         end if
         
      
      
      !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   





      return
      end
      
      
!!nitrification fuctions

      function f_gen_poisson_density (x1, x2)
      
      use parm
      
      implicit none
      
      real :: x1, x2
      real :: aa1, aa2, aa3, aa4 
      real :: f_gen_poisson_density
      real :: temp1, temp2, temp3
      

      
      aa1 =  35.
      aa2 = -5.
      aa3 = 4.5
      aa4 = 7.
      
      if ( x2 .ge. 35.) aa1 = x2
      
      if (aa2 == aa1) then
         f_gen_poisson_density = 0.
         
      else
      
          temp1 = (aa2-x1)/(aa2 -aa1)
            
          if (temp1 .le. 0.) then
             f_gen_poisson_density = 0.
             
          else
           temp2 = 1.0 - (temp1**aa4)
           temp3 = temp1** aa3
           f_gen_poisson_density = exp(aa3 * temp2 / aa4 * temp3)
      
          end if
       end if    
      return 
      end function
     