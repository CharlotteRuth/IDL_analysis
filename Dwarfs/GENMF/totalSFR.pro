; simple example idl script to generate and plot mass functions 

; generate Reed 06 mass functions
;$./genmf 0.23831863389003569 0.76168136610996429 .74 0. reed06-z=0.mf 0
;$./genmf 0.23831863389003569 0.76168136610996429 .74 30. reed06-z=30.mf 0

; ; generate corresponding Sheth and Tormen fit
;$./genmf 0.23831863389003569 0.76168136610996429 .74 0. ST-z=0.mf 2
;$./genmf 0.23831863389003569 0.76168136610996429 .74 30. ST-z=30.mf 2

pro totalSFR

; read them
readcol,'PS_0.mf',x0,y0,format='(d,d)'
loadct,39

tengalaxies = where(x0 lt 10.25 AND x0 gt 9.75)
N_galaxies = 0
FOR i = 0, N_ELEMENTS(tengalaxies)-2 DO BEGIN
    N_galaxies = N_galaxies + (10^y0[tengalaxies[i]]+10^y0[tengalaxies[i+1]])/2*(x0[tengalaxies[i]] - x0[tengalaxies[i+1]])
    print,x0[tengalaxies[i]],x0[tengalaxies[i+1]],y0[tengalaxies[i]],y0[tengalaxies[i+1]],(10^y0[tengalaxies[i]]+10^y0[tengalaxies[i+1]])/2*(x0[tengalaxies[i]] - x0[tengalaxies[i+1]]),N_galaxies
ENDFOR
print,N_galaxies

bins = [6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,15.5,16.5]
N_galaxies_array = fltarr(N_ELEMENTS(bins)-1)
FOR i1 = 0,N_ELEMENTS(bins) -2 DO BEGIN
    galaxies = where(x0 lt bins[i1+1] AND x0 gt bins[i1])
    FOR i = 0, N_ELEMENTS(tengalaxies)-2 DO BEGIN
        N_galaxies_array[i1] = N_galaxies_array[i1] + (10^y0[galaxies[i]]+10^y0[galaxies[i+1]])/2*(x0[galaxies[i]] - x0[galaxies[i+1]])
        print,x0[galaxies[i]],x0[galaxies[i+1]],y0[galaxies[i]],y0[galaxies[i+1]],(10^y0[galaxies[i]]+10^y0[galaxies[i+1]])/2*(x0[galaxies[i]] - x0[galaxies[i+1]]),N_galaxies_array[i1]
    ENDFOR
ENDFOR
print,N_galaxies_array
base = '~/Scratch1/DwarfResearch/ResMassTests/'

rtipsy,base+'5E1R/9M_k/o9M_1.00300',h9,g9,d9,s9
rtipsy,base+'1E2R/10M_k/o10M_1.00300',h10,g10,d10,s10
rtipsy,base+'1E3R/11M_k/o11M_1.00300',h11,g11,d11,s11
rtipsy,base+'1E4R/12M_k/o12M_1.00300',h12,g12,d12,s12
rtipsy,base+'1E5R/13M_k/o13M_1.00300',h13,g13,d13,s13

IF(N_ELEMENTS(s9) gt 0 ) THEN s9.mass = s9.mass*N_galaxies_array[2]
s10.mass = s10.mass*N_galaxies_array[3]
s11.mass = s11.mass*N_galaxies_array[4]
s12.mass = s12.mass*N_galaxies_array[5]
s13.mass = s13.mass*N_galaxies_array[6]

IF(N_ELEMENTS(s9) gt 0 ) THEN s_poor = [s9,s10,s11,s12,s13] ELSE s_poor = [s10,s11,s12,s13]

rtipsy,base+'1E5R/9M_k/o9M_1.00300',h9,g9,d9,s9
rtipsy,base+'1E6R/10M_k/o10M_1.00300',h10,g10,d10,s10
rtipsy,base+'1E5R/11M_k/o11M_1.00300',h11,g11,d11,s11
rtipsy,base+'1E6R/12M_k/o12M_1.00300',h12,g12,d12,s12
rtipsy,base+'1E5R/13M_k/o13M_1.00300',h13,g13,d13,s13

IF(N_ELEMENTS(s9) gt 0 ) THEN s9.mass = s9.mass*N_galaxies_array[2]
s10.mass = s10.mass*N_galaxies_array[3]
s11.mass = s11.mass*N_galaxies_array[4]
s12.mass = s12.mass*N_galaxies_array[5]
s13.mass = s13.mass*N_galaxies_array[6]

IF(N_ELEMENTS(s9) gt 0 ) THEN s_best = [s9,s10,s11,s12,s13] ELSE s_bset = [s10,s11,s12,s13]

set_plot,'x';,'ps'
;device,filename = "universalSFR_1e6.eps"
binsize=1.5e9
timeunit=1e9
sfr,s_best,massu=2.325e5,time=timeunit,OVERPLOT=0,xtitle = "t [Gyr]",ytitle = "SFR [solar masses year^-1 Mpc^-3]",binsize = binsize
sfr,s_poor,massu=2.325e5,time=timeunit,OVERPLOT=1,linestyle = 1,binsize = binsize
;;device,/close
stop
rtipsy,base+'1E5R/10M_k/o10M_1.00300',h10,g10,d10,s10
rtipsy,base+'1E5R/12M_k/o12M_1.00300',h12,g12,d12,s12

s10.mass = s10.mass*N_galaxies_array[3]
s12.mass = s12.mass*N_galaxies_array[5]
s_best = [s9,s10,s11,s12,s13]

;device,filename = "universalSFR_1e5.eps"
binsize=1.5e9
timeunit=1e9
sfr,s_best,massu=2.325e5,time=timeunit,OVERPLOT=0,xtitle = "t [Gyr]",ytitle = "SFR [solar masses year^-1 Mpc^-3]",binsize = binsize
sfr,s_poor,massu=2.325e5,time=timeunit,OVERPLOT=1,linestyle = 1,binsize = binsize
;device,/close




end
