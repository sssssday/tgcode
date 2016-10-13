      subroutine methane ()
      
      use parm   
      implicit none
      integer :: lyr, j
      integer :: idf
      real ::  soildthmax, soildthmin
      real :: CH4DEPTH 
      real :: bulkdensity  !! g/cm3
      real :: fieldcapacity !! cm h2o/cm soil
      real :: soiltemp !! degree 
      real :: soilwater !! cm water
      real :: wfps
      real :: percentlayer
      real :: PARTDENS
      real :: swcfrac
      real :: porespace
      real :: wfpslyr
      real :: CH4max
      real :: temp_adjust
      real :: wfps_adjust
      real :: Wmin
      real :: Wopt
      real :: Wmax
      real :: temp
      real :: Dopt
      real :: agri_adjust
      real :: watr_adjust
      real :: diffusiv
      real :: sol_thick
      real :: tsoilwater
 !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ variables from DLEM

      real :: atmCH4           !!(ppm) atmosphere CH4 content
      real :: production       !!rate of methane production (gC/m3/day)
      real :: ftemp,fpH,fdoc	        !!W,fEh,SI,VI unused variable deleted by ZC 3/23/07
      real :: fmoist_prod
      real :: fmoist_oxid
      real :: oxidize_air              !!oxidized ch4 from Atmosphere         gC/m3
      real :: oxidize_soil            !!oxidized ch4 from production in soil gC/m3
      real :: avDOC          !!(gC/m3)(concentration of DOC in top 50 cm soil layer that can be used for CH4 microbes)
      real :: Sita_ch4   !!(gC/m3)
      real :: pHmax !! changed by xiaofeng and ZC 4/16/07
      real :: pHmin !! changed by xiaofeng and ZC 4/16/07
      real :: maxsoilch4conc        !!//(gC/m3 h20)  // Kettunen 2003
      real :: planttrans                !!ratio of transport through plant
      real :: exchange_soilair          !!ratio of exchange between soil and air
      real :: emission
      real :: ch4depthd   !!assume ch4 only locates at the upper soil of 0.5 meter
      
      real :: pH
     
      real :: Q10
      real :: tempt 
      
      real :: fairt

      real :: total_doc
      real :: total_doc_top
      real :: temp2
      
      real :: wvwc
          
      real :: wpot_sati
      real :: vsat 
      real :: vfc 
      
      real :: PRMT_21, DK, X1, XX, V, X3
      integer :: kk
      real :: sut
      
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine calculate methane production 
!! DESCRIPTION:
!!    Methane oxidation is a function of soil water content, temperature,
!!    porosity, and field capacity.
!!
!!  REFERENCE:
!!    General CH4 Oxidation Model and Comparison of CH4 Oxidation in Natural
!!    and Managed Systems.  S.J. Del Grosso, W.J. Parton, A.R. Mosier, D.S.
!!    Ojima, C.S. Potter, W. Borken, R. Brumme, K. Butterbach-Bahi, P.M.
!!    Crill, K. Dobbie, and K.A. Smith.  2000.  Global Biogeochemical Cycles
!!    14:999-1019.

    !!   1,  Compute a weighted average for soil temperature, field capacity, */
      !! bulk density, water filled pore space, and volumetric soil water */
      !! content in top 15 cm of soil profile */
      j = ihru
      soildthmax = 0.
      soildthmin = 0.
      CH4DEPTH = 150. 
      PARTDENS = 2.65 !! g/cm3
      bulkdensity = 0.
      fieldcapacity = 0.
      soiltemp = 0.
      soilwater = 0.
      wfps = 0.
      percentlayer = 0.
      swcfrac = 0.
      porespace = 0.
      wfpslyr = 0.
      CH4max = 0.
      temp_adjust = 0.
      wfps_adjust = 0.
      sol_thick = 0.
      
      Wmin = 0.
      Wopt = 0.
      Wmax = 0.
      temp = 0.
      Dopt = 0.
      agri_adjust =0.
      watr_adjust = 0.0
      tsoilwater = 0.0
      
      idf = idplt(j)
      
  !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  


       atmCH4 = 1.8             
       production = 0.0         
       oxidize_air = 0.0              
       oxidize_soil = 0.0            
       pHmax = 10.0 
       pHmin = 4.0 
       maxsoilch4conc = atmCH4 * 12. / 22400. * 0.035       
       planttrans = 0.68                 
       exchange_soilair = 0.3           

       ch4depthd = 0.5  
       !! Sita_ch4 = soil_ch4 / ch4depth      
       
       Q10 = 2.5
       
       
             
      do lyr = 1, 2 !sol_nly(j)
      
      
       if (lyr == 1) then
	    sol_thick = sol_z(lyr,j)
	else	
	    sol_thick = sol_z(lyr,j) - sol_z(lyr-1,j)
	end if
      
       tempt = min(30.0, sol_tmp(lyr,j))
       pH = sol_ph(lyr,j)
        if(pH<pHmin .or. pH>pHmax) then 
            fpH = 0.
        else if (pH <= 7.0) then
            fpH = 1.02/(1+1000000 * exp(-2.5 * pH))
        else  
            fpH = 1.02/(1+1000000 * exp(-2.5 * (14.0-pH)))
        endif
       fpH = 10.**(-0.2235*pH*pH + 2.7727*pH - 8.6)
       
        if(tempt<-5.0) then
            ftemp = 0.0
        else 
            ftemp = Q10**((tempt-30.0) / 10.0)
        end if    
        
        fairt = Q10**((tmpav(j)-30.0) / 10.0)
        !! XXXXXXXXXXXXXXXXXXXXXX
          PRMT_21 = 0.  !KOC FOR CARBON LOSS IN WATER AND SEDIMENT(500._1500.) KD = KOC * C
          PRMT_21 = 1500.
          sol_WOC(lyr,j) = sol_LSC(lyr,j)+sol_LMC(lyr,j)+sol_HPC(lyr,j)
     &                      +sol_HSC(lyr,j)+sol_BMC(lyr,j) 
          DK=.0001*PRMT_21*sol_WOC(lyr,j)
          !X1=PO(LD1)-S15(LD1)
          X1 = sol_por(lyr,j)*sol_thick-sol_wpmm(lyr,j) !mm
          IF (X1 <= 0.) THEN
            X1 = 0.01
          END IF
          XX=X1+DK
          !V=QD+Y4
          V = sol_st(lyr,j) !+ sol_wpmm(lyr,j)
          X3=0.
          IF(V>1.E-10)THEN
              X3=sol_BMC(lyr,j)*(1.-EXP(-V/XX)) !loss of biomass C
          END IF
        
        total_doc = X3/10. !From Kg/ha C ot g/m2 C !DOM.C()

        total_doc_top = total_doc * 0.8  
        avDOC = total_doc_top
        fdoc = avDOC/(avDOC+ 1.)!! pftp[ecotype].Kmch4pro)
        
        wvwc = (sol_st(lyr,j)+ sol_wpmm(lyr,j))
     &  /sol_thick
        
        wpot_sati = -10.*10.**(1.88-0.0131*sol_sand(lyr,j))
        vsat = 0.489-0.00126*sol_sand(lyr,j)
        vfc = vsat*(-3.3e3/wpot_sati)**
     &   (-1./(2.91+0.159*sol_clay(lyr,j)))
        temp2 = (wvwc-vfc)/(vsat-vfc)

        if(wvwc < vfc) then
                fmoist_prod = 0.0
        else
       
                if(wvwc >= vsat) then
                    fmoist_prod = 1.0
                else 
                    fmoist_prod = temp2 * temp2 * 0.367879 * exp(temp2)
                end if    
        end if
       !!moist_prod = 1.
       !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         	    kk =0 
	          if (lyr == 1) then
	            kk = 2
	          else
	            kk = lyr
	          end if
       
       
                     X1=sol_st(lyr,j)
	          IF(X1<0.)THEN
	              SUT=.1*(sol_st(kk,j) /sol_wpmm(lyr,j))**2
	          ELSE
	              SUT = .1 + .9 * Sqrt(sol_st(lyr,j) / sol_fc(lyr,j))
              END IF             
              sut = min(1., sut)
              sut = max(.05, sut)
              
             !  fmoist_prod = sut
       !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
         if (lyr == 1)production = 0.
         
        production = production + min(0.4            !!pftp[ecotype].CH4ProMaxRate 
     &   * fdoc * ftemp* fpH * fmoist_prod, avDOC)
     & * sol_thick/1000.    !!(gC/m3 soil)  

        
        
        
        sol_BMC(lyr,j) = sol_BMC(lyr,j) - production*10. !remove produced CH4 from BMC pool
     
      end do
      
      
     
   !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
      
      
      do lyr = 1, sol_nly(j)
      
       if (lyr == 1) then
	    sol_thick = sol_z(lyr,j)
	else	
	    sol_thick = sol_z(lyr,j) - sol_z(lyr-1,j)
	end if
      
         soildthmin = soildthmax
         soildthmax = soildthmax + sol_thick 
        swcfrac = (sol_st(lyr, j)+ sol_wpmm(lyr,j)) / sol_thick 
        porespace = 1.0 - sol_bd(lyr,j)/PARTDENS !! orginal unit for bulk density is g/cm3 in daycent. check if SWAT variable match this value
        wfpslyr = swcfrac/porespace
      
        if (soildthmin < CH4DEPTH) then
          if (soildthmax <= CH4DEPTH)then
            bulkdensity  = bulkdensity + sol_bd(lyr,j) *  sol_thick / 
     &                         CH4DEPTH
            fieldcapacity = fieldcapacity + (sol_fc (lyr,j)+
     &       sol_wpmm(lyr,j))*0.1 * sol_thick / CH4DEPTH
            soiltemp = soiltemp + sol_tmp(lyr,j)* sol_thick / 
     &                         CH4DEPTH
            soilwater = soilwater + (sol_st(lyr, j)+sol_wpmm(lyr,j))
     &        *0.1* sol_thick /CH4DEPTH
            wfps = wfps + wfpslyr*sol_thick / CH4DEPTH
            else if ((soildthmax - soildthmin) > 0.0) then
            percentlayer = (CH4DEPTH - soildthmin) / 
     &                      (soildthmax - soildthmin)
     
            bulkdensity  = bulkdensity + sol_bd(lyr,j) *  sol_thick/ 
     &                         CH4DEPTH*percentlayer
            fieldcapacity = fieldcapacity + (sol_fc (lyr,j) +  
     &      sol_wpmm(lyr,j))*0.1 * sol_thick / CH4DEPTH*percentlayer
            soiltemp = soiltemp + sol_tmp(lyr,j)* sol_thick / 
     &                         CH4DEPTH*percentlayer
            soilwater = soilwater + (sol_st(lyr, j)+ sol_wpmm(lyr,j)) 
     &      *0.1* sol_thick /CH4DEPTH*percentlayer
            wfps = wfps + wfpslyr*sol_thick /CH4DEPTH*percentlayer
       
          endif
        endif
      
      end do
      tsoilwater = soilwater

     !! Convert from water filled pore space to volumetric water */
      soilwater = wfps  * (1.0 - (bulkdensity / PARTDENS))
      soilwater = wfps *100.0
      
 !! CH4 oxidation for a deciduous system */
      if (idf.eq.7) then  !! XXXXXXXXXXXXXXXXXXXXXXXXXXXneed to check which category is deciduous forest
        CH4max = 40.0 - 18.3 * bulkdensity
        temp_adjust = 0.0209 * soiltemp + 0.845
        !!Use bounded value for wfps_adjust if wfps falls below a critical */
        !! value, cak - 11/12/02 */
        if (wfps <= 0.05) then
          wfps_adjust = 0.1;
         else 
          wfps_adjust = ((10.0 * wfps - 0.5) / (1.84 - 0.5))** (0.13)
          wfps_adjust = wfps_adjust * ((10.0 * wfps - 55) / 
     &       (1.84 - 55))**(0.13 * (55 - 1.84) / (1.84 - 0.5))
          wfps_adjust = max(0.1, wfps_adjust);
        endif
         CH4(j) = CH4max * wfps_adjust * temp_adjust;

       else 
        !! CH4 oxidation for a grassland/coniferous/tropical system */
        Wmin = 3.0 * fieldcapacity - 0.28
        Wopt = 6.3 * fieldcapacity - 0.58
        Wmax = 10.6 * fieldcapacity + 1.9
        temp = Wopt * 0.1 / (1.0 - (bulkdensity / PARTDENS))
        temp = tsoilwater*10./CH4DEPTH/(1.0 - (bulkdensity / PARTDENS))
        Dopt = diffusiv(fieldcapacity*10./CH4DEPTH, bulkdensity, temp)
        
        CH4max = 53.8 * Dopt + 0.58
        if (((soilwater) < Wmin) .or. 
     &    ((soilwater) > Wmax)) then 
          watr_adjust = 0.1
         else 
          watr_adjust = (((soilwater - Wmin) / 
     &     (Wopt - Wmin))** 0.4) *
     &                   (((soilwater - Wmax) / (Wopt - Wmax))**
     &                        ((0.4 * (Wmax - Wopt)) / (Wopt - Wmin)))
     
     
          watr_adjust = max(0.1, watr_adjust)
        endif
        if (idf <7 ) then
          if (Dopt < 0.15) then
            agri_adjust = 0.9
           else if (Dopt > 0.28) then
            agri_adjust = 0.28
           else 
            agri_adjust = -4.6 * Dopt + 1.6
          endif
         else 
          agri_adjust = 1.0
        endif
        temp_adjust = (soiltemp * max(0.11, Dopt) * 0.095) + 0.9
        CH4(j) = production*10000. - CH4max * watr_adjust * temp_adjust 
     &  * agri_adjust
      endif
      return
      end
