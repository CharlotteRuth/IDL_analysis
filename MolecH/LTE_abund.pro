;base ='/astro/net/scratch2/christensen/MolecH/12M/Disk_Iso_1e6g2/MW_disk'
;step = '00004'
;msol_per_sysmass = 1.36e17
;kpc_per_syslength = 1e5
;LTE_abund_master,base = base,step = step,msol_per_sysmass = msol_per_sysmass, kpc_per_syslength = kpc_per_syslength,/calc_CO

;filename ='/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g3HBWK/steps/h516.cosmo25cmb.1536g3HBWK.00512.dir/h516.cosmo25cmb.1536g3HBWK.00512.halo.1'
;dir ='/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g3HBWK/steps/h516.cosmo25cmb.1536g3HBWK.00512.dir/'
;filename = 'h516.cosmo25cmb.1536g3HBWK.00512.halo.1'
;msol_per_sysmass = 2.310e15
;kpc_per_syslength        = 25000.
;LTE_abund_master,filename = filename,msol_per_sysmass = msol_per_sysmass, kpc_per_syslength = kpc_per_syslength,/molecularH

;dir = '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g8HBWK/steps/h516.cosmo25cmb.1536g9HBWK.00408.dir/'
;filename = 'h516.cosmo25cmb.1536g9HBWK.00408.halo.1'
;msol_per_sysmass = 2.310e15
;kpc_per_syslength        = 25000.
;LTE_abund_master,filename = filename,msol_per_sysmass = msol_per_sysmass, kpc_per_syslength = kpc_per_syslength,/molecularH

;dir ='/astro/net/scratch2/christensen/MolecH/12M/Disk_Collapse_1e6_H2_00100/'
;filename = 'Disk_Collapse.H2.ch0.c100.00100.00010'
;msol_per_sysmass = '2.362e5'
;kpc_per_syslength        = 1
;LTE_abund_master,filename = filename,msol_per_sysmass =msol_per_sysmass, kpc_per_syslength = kpc_per_syslength,/IsoDisk,/molecularH

;dir = '/astro/net/scratch2/christensen/MolecH/11M/Disk_Iso_1e5_zsol/uvtest/'
;filename = 'MW_disk.00001';'MW_disk.00010.db'
;dir = '/astro/net/scratch2/christensen/MolecH/11M/Disk_Iso_1e5_zsol/'
;filename = 'MW_disk.00010.db.LTE'
;msol_per_sysmass = 1.36e17
;kpc_per_syslength = 1e5
;LTE_abund_master,filename = filename,msol_per_sysmass =msol_per_sysmass, kpc_per_syslength = kpc_per_syslength,/IsoDisk,/molecularH

;dir ='/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g3HBWK_00504/steps/h516.cosmo25cmb.3072g3HBWK.00504.00002.dir/'
;filename = 'h516.cosmo25cmb.3072g3HBWK.00504.00002.halo.1'

;dir='/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g3HBWK_00504/steps/h516.cosmo25cmb.3072g6HBWK.00504.00008.dir/'
;filename = 'h516.cosmo25cmb.3072g6HBWK.00504.00008.halo.1'
;msol_per_sysmass = 2.310e15
;kpc_per_syslength = 25000.0
;LTE_abund_master,filename = filename,msol_per_sysmass =msol_per_sysmass, kpc_per_syslength = kpc_per_syslength

PRO LTE_abund_master,filename = filename,base = base, step = step, msol_per_sysmass = msol_per_sysmass, kpc_per_syslength = kpc_per_syslength, calc_CO = calc_CO, molecularH = molecularH, IsoDisk = IsoDisk
cm_per_kpc = 3.0857d21
gm_per_msol = 1.989d33
amu_per_gm = 6.022142d23
H_per_gm = 5.9753790e+23
molec_weight = (0.76*1 + 0.24*4.0)
ZSOLAR = 0.0130215
nbins = 100
loadct,39

IF KEYWORD_SET(base) AND KEYWORD_SET(step) THEN filename = base+'.'+step

rtipsy,filename,h,g,d,s
IF NOT KEYWORD_SET(ISODISK) THEN kpc_per_syslength = kpc_per_syslength*h.time
lengthunit =  kpc_per_syslength*cm_per_kpc ;system length unit in cm (=1kpc)
massunit   = msol_per_sysmass*gm_per_msol ;system mass unit in gm
dens_convert =  msol_per_sysmass * gm_per_msol/kpc_per_syslength^3/cm_per_kpc^3 ;massunit*amu_per_gm/molec_weight/lengthunit^3 ;Converts to grams/cm^3
t = g.tempg
rho = g.dens*dens_convert
zmetal = g.zmetal
lw_const = 0.5/kpc_per_syslength/kpc_per_syslength*1.2363437d-29;/h.time/h.time
readcol,filename+".HI",HI_prime,/silent
HI = HI_prime[1:h.ngas]
H2_prime = HI_prime
if keyword_set(molecularH) then readarr,filename+".H2",h,H2,/ascii,part = 'gas'

IF (FILE_TEST(filename+'.smoothlength')) THEN BEGIN
    readarr,filename+'.smoothlength',h,smoothlengths_L,/ascii,part = 'gas'
    smoothlengths = (smoothlengths_L*kpc_per_syslength*cm_per_kpc)
;    length = smoothlengths
ENDIF ELSE smoothlengths = g.tempg*0 + 1e20
IF (FILE_TEST(filename+'.shear')) THEN BEGIN
    readarr,filename+'.shear',h,mach,/ascii,part = 'gas'
    column = (1.0/mach*kpc_per_syslength*cm_per_kpc)
    length = column
ENDIF
IF (FILE_TEST(filename+'.correL')) then BEGIN
    readarr,filename+'.correL',h,column_L,/ascii,part = 'gas'
    column = column_L*cm_per_kpc*kpc_per_syslength
;    length = column
;;    mach = smoothlengths/column
ENDIF
;length = length*0 + 1e20
IF (FILE_TEST(filename+'.lw')) THEN BEGIN
    readarr,filename+".lw",h,lw_prime,/ascii,part = 'gas'
;    lw = lw_prime*1d30/cm_per_kpc^2/kpc_per_syslength^2*1.1d8*12.87*1.6021746d-12*4.0*!PI*4.13e-15
    lw = lw_prime[0:h.ngas - 1]*lw_const
;*1d30/cm_per_kpc^2/kpc_per_syslength^2/1e6
ENDIF

window,0
!p.multi = [0,2,2]
minlw = alog10(1d20*lw_const)
maxlw = alog10(2d34*lw_const)
lw_cut = lw
IF ((where(lw gt 2d31*lw_const))[0] ne -1) THEN lw_cut[where(lw gt 2d31*lw_const)] = 2d31*lw_const
colors_lw = (alog10(lw) - minlw)*240/(maxlw - minlw)
ind = where(alog10(lw) lt minlw)
if (ind[0] ne -1) THEN colors_lw[ind] = 0
ind = where(alog10(lw) gt maxlw)
if (ind[0] ne -1) THEN colors_lw[ind] = 240
y_lwhist = histogram(alog10(lw), min = minlw, max = maxlw, nbins = nbins,locations = x_lwhist)
plot,x_lwhist - alog10(lw_const),y_lwhist/FLOAT(N_ELEMENTS(lw)),psym = 10,xtitle = 'Log Lyman Werner Radiation',charsize = 2

minz = MIN(g[where(g.zmetal gt 1e-10)].zmetal/zsolar) ;1e-2
maxz = MAX(g.zmetal/zsolar)
colors_z = (g.zmetal/zsolar - minz)*240/(maxz - minz)
ind = where(g.zmetal/zsolar lt minz)
if (ind[0] ne -1) THEN colors_z[ind] = 0
ind = where(g.zmetal/zsolar gt maxz)
if (ind[0] ne -1) THEN colors_z[ind] = 240
y_zhist = histogram(g.zmetal/zsolar, min = minz, max = maxz, nbins = nbins,locations = x_zhist)
plot,x_zhist,y_zhist,psym = 10,xtitle = 'Metallicity [z/z_solar]',charsize = 2

minL = 18.0
maxL = 22.0
colors_L = (alog10(smoothlengths) - minL)/(maxL - minL)*240
ind = where(alog10(smoothlengths) lt minL)
if (ind[0] ne -1) THEN colors_L[ind] = 0
ind = where(alog10(smoothlengths) gt maxL)
if (ind[0] ne -1) THEN colors_L[ind] = 240

minL = 18.0
maxL = 22.0
colors_c = (alog10(column) - minL)/(maxL - minL)*240
ind = where(alog10(column) lt minL)
if (ind[0] ne -1) THEN colors_c[ind] = 0
ind = where(alog10(column) gt maxL)
if (ind[0] ne -1) THEN colors_c[ind] = 240
y_Lhist = histogram(alog10(column), min = minL, max = maxL, nbins = nbins,locations = x_Lhist)
plot,x_Lhist,y_Lhist,psym = 10,xtitle = 'Log Column Length [cm]',charsize = 2

y_Lhist = histogram(alog10(smoothlengths), min = minL, max = maxL, nbins = nbins,locations = x_Lhist)
oplot,x_Lhist,y_Lhist,psym = 10,linestyle = 2
legend,['Column Length','Smooth Length'],linestyle = [0,2]

minD = -2
maxD = 3
colors_D = (alog10(H_per_gm*rho) - minD)/(maxD - minD)*240
ind = where(alog10(H_per_gm*rho) lt minD)
if (ind[0] ne -1) THEN colors_D[ind] = 0
ind = where(alog10(H_per_gm*rho) gt maxD)
if (ind[0] ne -1) THEN colors_D[ind] = 240
y_Dhist = histogram(alog10(H_per_gm*rho), min = minD, max = maxD, nbins = nbins,locations = x_Dhist)
plot,x_Dhist,y_Dhist,psym = 10,xtitle = 'Log Density [amu/cc]',charsize = 2

minT = 1
maxT = 4.5
colors_T = (alog10(g.tempg) - minT)/(maxT - minT)*240
ind = where(alog10(g.tempg) lt minT)
if (ind[0] ne -1) THEN colors_T[ind] = 0
ind = where(alog10(g.tempg) gt maxT)
if (ind[0] ne -1) THEN colors_T[ind] = 240
;y_Dhist = histogram(alog10(H_per_gm*rho), min = minD, max = maxD, nbins = nbins,locations = x_Dhist)
;plot,x_Dhist,y_Dhist,psym = 10,xtitle = 'Log Density [amu/cc]',charsize = 2
stop

;mach = fltarr(N_ELEMENTS(t)) + 1.0      ;*******************THIS ASSUMES A MACH SPEED OF 1.0 TO CALCULATE A TURBULENT LENGTH SCALE (h) -- PROBABLY WRONG***************
all = fltarr(9,N_ELEMENTS(t))
all_smooth = fltarr(9,N_ELEMENTS(t))
;all_lw = fltarr(9,N_ELEMENTS(t))
All_smooth_lw = fltarr(9,N_ELEMENTS(t))
All_smooth_clump = fltarr(9,N_ELEMENTS(t))
All_clump = fltarr(9,N_ELEMENTS(t))
IF KEYWORD_SET(calc_CO) THEN  co = fltarr(2,N_ELEMENTS(t))

length = column ;smoothlengths

for i=0LL,LONG64(N_ELEMENTS(t)) - 1 do  BEGIN
    if (g[i].dens gt 0.01) THEN BEGIN
;    all[*,i] = LTE_abund(t[i],g[i].dens*dens_convert,zmetal[i],length[i],lw[i],10)
        if i eq 57507 then all[*,i] = LTE_abund(t[i],g[i].dens*dens_convert,zmetal[i],column[i],lw[i],10,/verbose) else all[*,i] = LTE_abund(t[i],g[i].dens*dens_convert,zmetal[i],column[i],lw[i],10)
        all_smooth[*,i] = LTE_abund(t[i],g[i].dens*dens_convert,zmetal[i],smoothlengths[i],lw[i],10)
;    all_lw[*,i] = LTE_abund(t[i],g[i].dens*dens_convert,zmetal[i],column[i],lw_cut[i],10)
        all_smooth_lw[*,i] = LTE_abund(t[i],g[i].dens*dens_convert,zmetal[i],smoothlengths[i],lw_cut[i],10)
        all_smooth_clump[*,i] = LTE_abund(t[i],g[i].dens*dens_convert,zmetal[i],smoothlengths[i],lw[i],1000)
        all_clump[*,i] = LTE_abund(t[i],g[i].dens*dens_convert,zmetal[i],smoothlengths[i],lw[i],100)
        IF KEYWORD_SET(calc_co) THEN co[*,i] = calcCO(t[i],g[i].dens*dens_convert,zmetal[i],length[i])
        if (i mod 100000) eq 0 THEN print,i
        ENDIF
ENDFOR

Y_HI = REFORM(all_clump[0,*])
Y_H2 = REFORM(all_clump[1,*]) 
Y_HI_smooth = REFORM(all_smooth[0,*])
Y_H2_smooth = REFORM(all_smooth[1,*])
;Y_HI = REFORM(all_lw[0,*])
;Y_H2 = REFORM(all_lw[1,*])
;Y_HI = REFORM(all_smooth_lw[0,*])
;Y_H2 = REFORM(all_smooth_lw[1,*])   
print,''
;print,'LTE Range'
;print,'HI: ',MINMAX(Y_HI/total_H)
;print,'H2: ',MINMAX(Y_H2*2/total_H)
;print,'LOG(HI): ',aLOG10(MINMAX(Y_HI/total_H))
;print,'LOG(H2): ',ALOG10(MINMAX(Y_H2*2/total_H))

window,4,xsize = 800,ysize = 800
colors = colors_T
IF KEYWORD_SET(molecularH) THEN !p.multi = [0,2,2] ELSE !p.multi = [0,2,1]
plot,rho*5.9753790d23,2.0*Y_H2/(2.0*Y_H2 + Y_HI),psym = 3,/xlog,yrange = [1e-6,1],xrange = [1e-2,1e4],/ylog
FOR i = 0L, N_ELEMENTS(g.dens) - 1 DO $
  oplot,[rho[i],rho[i]]*5.9753790d23,2.0*[Y_H2[i]/(2.0*Y_H2[i] + Y_HI[i]),Y_H2[i]/(2.0*Y_H2[i] + Y_HI[i])],psym = 3,color = colors[i]

plot,rho*5.9753790d23*length,2.0*Y_H2/(2.0*Y_H2 + Y_HI),psym = 3,/xlog,/ylog,yrange = [1e-6,1],xrange = [1e17,1e24]
FOR i = 0L, N_ELEMENTS(g.dens) - 1 DO $
  oplot,[rho[i]*length[i],rho[i]*length[i]]*5.9753790d23,2.0*[Y_H2[i]/(2.0*Y_H2[i] + Y_HI[i]),Y_H2[i]/(2.0*Y_H2[i] + Y_HI[i])],psym = 3,color = colors[i]

IF KEYWORD_SET(molecularH) THEN BEGIN
    plot,rho*5.9753790d23,2.0*H2/(2.0*H2 + HI),psym = 3,/xlog,yrange = [1e-6,1],xrange = [1e-2,1e4],/ylog
    FOR i = 0L, N_ELEMENTS(g.dens) - 1 DO $
      oplot,[rho[i],rho[i]]*5.9753790d23,2.0*[H2[i]/(2.0*H2[i] + HI[i]),H2[i]/(2.0*H2[i] + HI[i])],psym = 3,color = colors[i]

    plot,rho*5.9753790d23*length,2.0*H2/(2.0*H2 + HI),psym = 3,/xlog,/ylog,yrange = [1e-6,1],xrange = [1e17,1e24]
    FOR i = 0L, N_ELEMENTS(g.dens) - 1 DO $
      oplot,[rho[i]*column[i],rho[i]*length[i]]*5.9753790d23,2.0*[H2[i]/(2.0*H2[i] + HI[i]),H2[i]/(2.0*H2[i] + HI[i])],psym = 3,color = colors[i]
ENDIF
stop

HI_prime[1:h.ngas] = Y_HI
H2_prime[1:h.ngas] = Y_H2
openw,1,filename+".HI_LTE"
openw,2,filename+".H2_LTE"
IF KEYWORD_SET(calc_CO) THEN BEGIN
    co1 = HI_prime
    co2 = HI_prime
    co1[1:h.ngas] = reform(co[0,*])
    co2[1:h.ngas] = reform(co[1,*])
    openw,3,filename+".CO"
    openw,4,filename+".I_CO"
ENDIF
printf,1,LONG(HI_prime[0]),FORMAT = '(I)'
printf,2,LONG(H2_prime[0]),FORMAT = '(I)'
IF KEYWORD_SET(calc_CO) THEN printf,3,FIX(co1[0]),FORMAT = '(I)'
IF KEYWORD_SET(calc_CO) THEN printf,4,FIX(co2[0]),FORMAT = '(I)'
for i=1LL,LONG64(N_ELEMENTS(HI_prime)) - 1 do  BEGIN
    printf,1,HI_prime[i]
    printf,2,H2_prime[i]
    IF KEYWORD_SET(calc_CO) THEN    printf,3,co1[i]
    IF KEYWORD_SET(calc_CO) THEN    printf,4,co2[i]
ENDFOR
close,1
close,2
IF KEYWORD_SET(calc_CO) THEN close,3
IF KEYWORD_SET(calc_CO) THEN close,4
sortrho = SORT(rho)
sortz = SORT(zmetal)
IF KEYWORD_SET(calc_CO) THEN plot,rho[sortrho]*5.9753790d23,co[0,sortrho],psym = 3,/xlog
IF KEYWORD_SET(calc_CO) THEN lot,zmetal[sortz]/0.0177,co[0,sortz]


END

FUNCTION LTE_abund,T,rho,zmetal,length,lw,clump,verbose = verbose
;LTE_abund,temperature,density,metallicity,length radiation travels over in cm ( = c/gas_shear)
CL_B_gm  = 6.022e23*(938.7830/931.494)
en_B = CL_B_gm*rho     
EPS = 1e-5
mu_metal = 17.600299999999997

if (ZMetal le 0.1) then  yHe = (0.236 + 2.1*ZMetal)/4.0 else yHe = (-0.446*(ZMetal - 0.1)/0.9 + 0.446)/4.0
yH = 1.0 - yHe*4.0 - ZMetal;
yH2 = 1d-6;0                      
yHI = 1d0;0                      
yHII = 0                        
;yHe
yHeI = 0                        
yHeII = 0                   
yHeIII = 0       
YeMax = yH + yHe*2               

R_Radr_HII = clRateRadrHII(T)  ;radiative recombination
R_Phot_HI = 9.300920e-13; photon dissociation from UV background radation at z = 0
R_Coll_HI = clRateCollHI(T) ; collisional dissociation
R_TotrForm_H2 = clRateTotrFormH2(ZMetal, clump); Total rate of formation of H2
R_Coll_e_H2 = clRateColl_e_H2(T); collisional dissociation by electrons
R_Coll_H_H2 = clRateColl_H_H2(T); collisional dissociation by H atoms
R_Coll_H2_H2 = clRateColl_H2_H2(T); collisional dissociation by H2
R_Phot_H2 = lw + 10^(-12.518825);9.38052e-11;(Corresponds to 1e6 photons/str) -1.591700e-11; rate of photon dissociation from UV background radation at z = 0 ;0.5*lw*1.175577d-16

R_Phot_HeI = 2.850000e-14 ; rate of photon dissociation from UV background radation at z = 0
R_Totr_HeII =clRateRadrHeII(T) + clRateDielHeII(T) ; Total rate of recombination of HeII
R_Phot_HeII = 2.940000e-16 ; rate of photon dissociation from UV background radation at z = 0
R_Radr_HeIII =clRateRadrHeII(T); Rate of radative recombination
R_Coll_HeI = clRateCollHeI(T) ; Rate of collisional dissociation
R_Coll_HeII = clRateCollHeII(T) ; Rate of collisional dissociation

rcirrHeI  = (R_Coll_HeI)/(R_Totr_HeII) ;
rcirrHeII = (R_Coll_HeII)/(R_Radr_HeIII) ;
rpirrHeI  = (R_Phot_HeI)/(R_Totr_HeII * en_B) ;
rpirrHeII = (R_Phot_HeII)/(R_Radr_HeIII * en_B) ;


S_H2 = 1.0 ;S_H2(,h)         ;self shielding
S_d  = S_d(yH*en_B,0,length,zmetal) ;shielding from dust

for  i=0, 20 do begin
    yHI_old   = yHI             
    yHeII_old = yHeII           
    yH2_old   = yH2             
    yHII_old = yH - yHI - 2.0*yH2;
   
    ye = (yeMax-(yHI + 2 * yH2 + 2 * yHeI + yHeII)) ; //Free electrons

    if (ye le 0) then begin
        y_e = 0                 
        yHII = 0    
        yHeI = yHe             
        yHeII = 0              
        yHeIII = 0       
        if (s_d eq 0 OR s_H2 eq 0) THEN BEGIN
            yHI = 0             ;
            yHII = 0            ;
            yH2 = yH/2.0        ;
        ENDIF ELSE BEGIN
        fHI = 2.0*(R_TotrForm_H2*en_B*(yHI_old + yH2_old))/ $
                 (R_Coll_H_H2*en_B*yHI_old  +  $
                   R_Coll_H2_H2*en_B*yH2_old +  $ 
                   (R_Phot_H2)*S_d*S_H2)
        yHI = yH/(1 + fHI)
        yH2 = (yH - yHI)/2.
        ENDELSE
    endif else begin
        rye = 1/ye              
    
        IF (s_d le 1d-150) THEN BEGIN
            yHI = 0
            yH2 = yH/2.0
            yHII = 0
        ENDIF ELSE BEGIN
            fHII = (R_Radr_HII*en_B*ye)/(R_Phot_HI*S_d + R_Coll_HI*en_B*ye) 
            fHI = 2.0*(R_TotrForm_H2*en_B*(yHI_old + yH2_old))/ $
                     (R_Coll_e_H2*en_B*ye       +  $
                       R_Coll_H_H2*en_B*yHI_old  +  $ 
                       R_Coll_H2_H2*en_B*yH2_old +  $
                       (R_Phot_H2)*S_d*S_H2)
            rfH  =  1 / ( 1 + fHII * (1 + fHI) )  
            yHII =  yH * rfH        
            yHI  =  yH * rfH * fHII 
            yH2  = (yH - yHI - yHII)/2.0 
        ENDELSE

        fHeI  = rcirrHeI + rpirrHeI * rye 
        fHeII = rcirrHeII+rpirrHeII * rye 
        rfHe  = 1 / ( 1 + fHeI * (1 + fHeII) ) 
        yHeI  = yHe * rfHe      
        yHeII = yHe * fHeI * rfHe 
        yHeIII = yHe / ((1.0/fHeI+1.0)/fHeII+1.0) 
    endelse
    IF(~FINITE(YHI) OR ~FINITE(YH2)) THEN BEGIN
        print,rho*5.9753790e+23,T
        print,yHI, yH2
        print,yHI_old,yH2_old
        print,s_d,s_H2
        yHI = yHI_old
        yH2 = yH2_old
        stop
        break
    ENDIF
    if ( ABS(yHeII_old-yHeII) lt EPS * yHeII AND abs(yHI_old-yHI) lt EPS * yHI ) then  break

    if (YH2 lt 0) then begin
        YH2 = 0                 ;
        YHII = yH - yHI
    ENDIF
    if (YHI lt 0) then yHI = 0;
    if (YHII lt 0) then begin
        YHII = 0 
        YH2 = (YH - yHI)/2;
    ENDIF
    if (YH2 gt yH/2.0) then YH2 = yH/2.0 
    if (YHI gt yH) then YHI = yH 
    if (YHII gt yH) then YHII = yH
    min = 1e-12
    if (YH2 lt min) then YH2 = min
    if (YHI lt min) then YHI = min
    if (YHII lt min) then YHII = min
;    if (en_B gt 8 AND T lt 5000) then BEGIN
;        print,strtrim(i,2)," enB: ",strtrim(en_B,2)," T: ",strtrim(T,2)," HI: ",strtrim(YHI,2),";  H2: ",strtrim(YH2,2),"; HII: ",strtrim(YHII,2),";"
;        stop
;    ENDIF
    S_H2 = S_H2(yH2*en_B,length) ;self shielding
    S_d  = S_d(yHI*en_B,yH2*en_B,length,zmetal) ;shielding from dust

endfor
Y_e = ye                        

Y_HI = yHI                      
Y_HII = yHII                   
Y_H2 = yH2                      
Y_HeI = yHeI                    
Y_HeII = yHeII                  
Y_HeIII = yHeIII
IF(~FINITE(Y_HI) OR ~FINITE(Y_H2)) THEN stop
;print,"HI: ",Y_HI,"; H2: ",Y_H2,"; HII: ",Y_HII,";"
if KEYWORD_SET(verbose) THEN stop
RETURN,[Y_HI,Y_H2,Y_HII,Y_HeI,Y_HeII,Y_HeIII,Y_e,s_d,s_H2]
;return,number fraction of HI, number fraction of H2, number fraction of HeI,  number fraction of HeII,  number fraction of HeIII,  number fraction of electrons, dust shielding parameter, self-shielding parameter
;In general, by number fraction, I mean (the number of HI)/(total number of atoms)  
;This is the same format as the XXXXX.HI files outputted from gasoline are in
;For electrons, Y_e is the average number of free electrons per nuclei
end


;CL_Rgascode      =   8.2494e7
;CL_Eerg_gm_degK    =    CL_Rgascode
;CL_ev_degK         =    1.0/1.1604e4
;CL_Eerg_gm_ev      =    CL_Eerg_gm_degK/CL_ev_degK
;CL_Eerg_gm_degK3_2   =  1.5*CL_Eerg_gm_degK
;CL_MAX_NEG_EXP_ARG  =   -500.
;CL_B_gm        =    (6.022e23*(938.7830/931.494))  ;/*Avegadro's Number * Mass_Hydrogen/Energy_AMU */
;CL_k_Boltzmann   =  1.38066e-16
;CL_eV_erg       =   1.60219e-12
;CL_eV_per_K     =   (CL_k_Boltzmann/CL_eV_erg)

; /*-----------------------------------------------------------------
; *     Collisional Ionization rates
; *-----------------------------------------------------------------*/
;/*     H + e- -> H+ + 2e-  Janev et al. 1987 (Abel 1996) */
function clRateCollHI,T
CL_MAX_NEG_EXP_ARG = -500.
CL_k_Boltzmann   =  1.38066d-16
CL_eV_erg       =   1.60219d-12
CL_eV_per_K     =   (CL_k_Boltzmann/CL_eV_erg)

TL = alog(T*CL_eV_per_K)        
arg = -32.713967867 + TL*(13.536556     + TL*(-5.73932875 + $
                                              TL*(1.56315498 + $
                                                  TL*(-0.2877056     + TL*(3.48255977e-2 + TL*(-2.63197617e-3 + $
                                                                                               TL*(1.11954395e-4  + TL*(-2.03914985e-6)))))))) 
if (arg lt CL_MAX_NEG_EXP_ARG) then return, 0 
return, DOUBLE(exp ( arg ))     
end
  
;/*     He + e- -> He+ + 2e-  Janev et al. 1987 (Abel 1996) */
function clRateCollHeI,T 
CL_MAX_NEG_EXP_ARG = -500.
CL_k_Boltzmann   =  1.38066d-16
CL_eV_erg       =   1.60219d-12
CL_eV_per_K     =   (CL_k_Boltzmann/CL_eV_erg)

TL = alog(T*CL_eV_per_K)        
arg = -44.09864886  + TL*(23.91596563   + TL*(-10.7532302 +  $
                      TL*(3.05803875 +  $
                      TL*(-0.56851189    + TL*(6.79539123e-2 + TL*(-5.00905610e-3 + $
                      TL*(2.06723616e-4  + TL*(-3.64916141e-6)))))))) 
if (arg lt CL_MAX_NEG_EXP_ARG) then return, 0 
return, DOUBLE(exp( arg ))
end

;/*     He+ + e- -> He++ + 2e- Aladdin Database 1989 (Abel 1996) */
function clRateCollHeII,T  
CL_MAX_NEG_EXP_ARG = -500.
CL_k_Boltzmann   =  1.38066e-16
CL_eV_erg       =   1.60219e-12
CL_eV_per_K     =   (CL_k_Boltzmann/CL_eV_erg)

TL = alog(T*CL_eV_per_K)        
arg = -68.71040990  + TL*(43.93347633 + TL*(-18.4806699 + $
                      TL*(4.70162649  + TL*(-0.76924663 + $
                      TL*(8.113042e-2 + TL*(-5.32402063e-3 +  $
                      TL*(1.97570531e-4  + TL*(-3.16558106e-6))))))))
if (arg lt CL_MAX_NEG_EXP_ARG) then return, 0
return, DOUBLE(exp( arg ))
end

;//Proton Collision
;//Lepp & Shull, 1983
function clRateColl_H2_H2,T
CL_MAX_NEG_EXP_ARG = -500.
CL_k_Boltzmann   =  1.38066e-16
CL_eV_erg       =   1.60219e-12
CL_eV_per_K     =   (CL_k_Boltzmann/CL_eV_erg)

if (T lt 7291) then ratecollH2 = 5.22e-14*exp(-3.22e4/T) else ratecollH2 = 3.17e-15*exp(-4060./T - (7500./T)*(7500./T))
return, DOUBLE(ratecollH2)
end

;//Electron Collision
;//Donahue & Shull, Abel 1997
function clRateColl_e_H2,T
CL_MAX_NEG_EXP_ARG = -500.
CL_k_Boltzmann   =  1.38066e-16
CL_eV_erg       =   1.60219e-12
CL_eV_per_K     =   (CL_k_Boltzmann/CL_eV_erg)

ratecollH2 = 5.6e-11*SQRT(CL_eV_per_K*T)*exp(-8.8/(CL_eV_per_K*T)) ; //Is 8.8eV the right number?  Abel has 102124
return, DOUBLE(ratecollH2)
end

;//Neutral H
;//Dove and Mandy 1986, Abel 1997
function clRateColl_H_H2,T
;  //Donahue & Shull  ratecollH2 = 6.11e-14*exp(-4.48*CL_eV_per_K/T);
CL_MAX_NEG_EXP_ARG  =-500.
CL_k_Boltzmann   =  1.38066e-16
CL_eV_erg       =   1.60219e-12
CL_eV_per_K     =   (CL_k_Boltzmann/CL_eV_erg)

ratecollH2 = 1.067e-10*(CL_eV_per_K*T)^2.012*exp(-1.0*(4.463/T/CL_eV_per_K)*(1+0.2472*CL_eV_per_K*T)^3.512) 
return, DOUBLE(ratecollH2)      
end

;'/*-----------------------------------------------------------------
; *     Radiative Recombination rates
; *-----------------------------------------------------------------*/
;/*     H+ + e- -> H + gam  Verner & Ferland 1996 */
function clRateRadrHII,T 
Tsq = sqrt(T)                   
return, DOUBLE(7.982e-11/( Tsq*0.563615 * (1+Tsq*0.563615)^0.252 * (1+Tsq*1.192167e-3)^1.748)) 
end

;/*     He+ + e- -> He + gam  radiative  Verner & Ferland 1996 */
function clRateRadrHeII,T 
lgT = alog10(T)                  
lg_rec = -10.3217 + 0.749031* lgT - 1.66439*lgT^2 + 1.13114 * lgT^3 -0.459420 * lgT^4.0 + 0.114337 * lgT^5.0 -0.0174612 * lgT^6.0 + 0.00158680 * lgT^7.0 -7.86098e-05 * lgT^8 +  1.63398e-06 * lgT^9.0 
return, DOUBLE(10.0^lg_rec)      
end

function clRateTotrFormH2,z,clump
ZSOLAR = 0.0130215
;clump = 10.0   ; /*Ranges from 2-10 to 30-100, Gendin et al 2008 CC*/ 
Rate_dust = 3.5e-17*z/ZSOLAR*clump ; /*Formation rate coefficient of molecular hydrogen on dust, Gnedin et al 2008, Wolfire 2008, unit of cc per s CC*/  
return, DOUBLE(Rate_dust)       
end

function clRateGasFormH2,T 
  return, 0;
end

function clRateDielHeII,T
CL_MAX_NEG_EXP_ARG  =   -500.
T_inv = 1.0/T
arg = -4.7e5*T_inv              
if (arg lt CL_MAX_NEG_EXP_ARG)THEN return, 0 
return, DOUBLE(1.9e-3*T^(-1.5)*exp(arg)*(1+0.3*exp(-9.4e4*T_inv))) ;
end

function clRateChtrHeII,T
T_4 = T/1e4                      
return, 0.                       
end

;/*     He++ + e- -> He+ + gam  Verner & Ferland 1996 */
function clRateRadrHeIII,T
Tsq = sqrt(T)                   ;
return, DOUBLE(1.891e-10/( Tsq*0.326686 * (1+Tsq*0.326686)^0.2476 * (1+Tsq*6.004084e-4)^1.7524)) ;
end

function S_H2,yH2,h
omega_H2 = 0.2
x = yH2*h/5d14
return, (1 - omega_H2)/(1 + x)^2 + omega_H2/SQRT(1 + x)*exp(-0.00085*SQRT(1 + x))
end

function S_d,yHI,yH2,z,h
ZSOLAR = 0.0130215
sigma_d = 2d-21
return, exp(-1.0*sigma_d*z/ZSOLAR*(yHI*h + 2.0*yH2*h))
end

