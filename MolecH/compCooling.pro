;also see coolingCurveDetail

;Charlotte Christensen: 4/15/10
;This program is written to analyze the cooling debugging information
;from gasoline contained in coolingdebug2.txt 

;To create the input for this program, run gasoline with 
;CODE_DEF = -DGASOLINE -DCHANGESOFT -DCOOLING_METAL -DCOOLDEBUGOUT -DCOOLDEBUG
;To avoid extraneous heating, I recomend running it on glass_cube.std
;with glass_cube.param with no star formation, no gravity and no
;external heating ( -DNOEXTHEAT )

;Main Program
;compCooling,dirs,key = key,outplot = outplot,color = color,sym = sym,dosteps = dosteps

;To run:
;compCooling,dirs
;dirs: array of directories containing the cooldebug2.txt files you
;would like to analyze

;Optional inputs
;key: the labels for the legend,
;outplot: if set to the name of the file you would like a color = color,sym = sym
;color: array of colors used to plot the different cooling curves
;sym: array of symbols used to plot the different cooling curves

;compCooling_master can be modified to set up the input for compCooling.pro

FUNCTION readCoolDebug,filename
filecool = {cooldebug, $
          T: DOUBLE(0.0), $
          energy: 0.0, $
          en_B: 0.0, $
          n_e: 0.0, $
          H2: 0.0, $
          HI: 0.0, $
          HII: 0.0, $
          HeI: 0.0, $
          HeII: 0.0, $
          HeIII: 0.0, $
          Sheild: 0.0, $
          Edot: 0.0, $
          Comp: 0.0, $
          Brems: 0.0, $
          Diel: 0.0, $
          Radr: 0.0, $
          LineHI: 0.0, $
          LineHeI: 0.0, $
          LineHeII: 0.0, $
          LineH2_H: 0.0, $
          LineH2_H2: 0.0, $
          LineH2_He: 0.0, $
          LineH2_e: 0.0, $
          LineH2_HII: 0.0, $
          Coll: 0.0, $
          CollH2: 0.0, $
;          CollH2_e: 0.0, $
;          CollH2_H: 0.0, $
;          CollH2_H2: 0.0, $
          LowT: 0.0, $
          MCool: 0.0, $
          Cool: 0.0, $
          Phot: 0.0, $
          PhotH2: 0.0, $
          MHeat: 0.0 $
}
print,filename
readcol,filename, T, energy, en_B, n_e, H2, HI, HII, HeI, HeII, HeIII, Sheild, Edot, $
                  Comp, Brems, Diel, Radr, LineHI, LineHeI, LineHeII, LineH2_H, LineH2_H2, LineH2_He, LineH2_e, LineH2_HII, Coll, CollH2_e, CollH2_H, CollH2_H2, LowT, MCool, Cool, Phot, PhotH2, MHeat,/silent
; readcol,filename, T, energy, en_B, n_e, H2, HI, HII, HeI, HeII, HeIII, Sheild, Edot, $
;                  Comp, Brems, Diel, Radr, LineHI, LineHeI, LineHeII, LineH2_H, LineH2_H2, LineH2_He, LineH2_e, LineH2_HII, Coll, CollH2, LowT, MCool, Cool, Phot, PhotH2, MHeat,/silent
;readcol,filename, T,         en_B, n_e, H2, HI, HII, HeI, HeII, HeIII, Sheild, Edot, $
;                  Comp, Brems, Diel, Radr, LineHI,                    LineH2_H, LineH2_H2, LineH2_He, LineH2_e, LineH2_HII, Coll, CollH2, LowT, MCool, Cool, Phot, PhotH2, MHeat

filecool_arr = REPLICATE(filecool,N_ELEMENTS(en_B))

FOR i = 0L, LONG(N_ELEMENTS(en_B) - 1) DO BEGIN
    filecool_arr[i].LowT = LowT[i]
    filecool_arr[i].energy = energy[i]
    filecool_arr[i].MCool = MCool[i]
    filecool_arr[i].Cool = Cool[i]
    filecool_arr[i].Phot = Phot[i]
    filecool_arr[i].PhotH2 = PhotH2[i]
    filecool_arr[i].MHeat = MHeat[i]
    filecool_arr[i].T = T[i]
    filecool_arr[i].en_B = en_B[i]
    filecool_arr[i].n_e = n_e[i]
    filecool_arr[i].H2 = H2[i]
    filecool_arr[i].HI = HI[i]
    filecool_arr[i].HII = HII[i]
    filecool_arr[i].HeI = HeI[i]
    filecool_arr[i].HeII = HeII[i]
    filecool_arr[i].HeIII = HeIII[i]
    filecool_arr[i].Sheild = Sheild[i]
    filecool_arr[i].Edot = Edot[i]
    filecool_arr[i].Comp = Comp[i]
    filecool_arr[i].Brems = Brems[i]
    filecool_arr[i].Diel = Diel[i]
    filecool_arr[i].Radr = Radr[i]
    filecool_arr[i].LineHI = LineHI[i]
    filecool_arr[i].LineHeI = LineHeI[i]   
    filecool_arr[i].LineHeII = LineHeII[i]
    filecool_arr[i].LineH2_H = LineH2_H[i] 
    filecool_arr[i].LineH2_H2 = LineH2_H2[i] 
    filecool_arr[i].LineH2_He = LineH2_He[i] 
    filecool_arr[i].LineH2_e = LineH2_e[i] 
    filecool_arr[i].LineH2_HII = LineH2_HII[i] 
    filecool_arr[i].Coll = Coll[i]
;    filecool_arr[i].CollH2 = CollH2[i]
    filecool_arr[i].CollH2 = CollH2_e[i] + CollH2_H[i] + CollH2_H2[i]
;    filecool_arr[i].CollH2_e = CollH2_e[i]
;    filecool_arr[i].CollH2_H = CollH2_H[i]
;    filecool_arr[i].CollH2_H2 = CollH2_H2[i]

ENDFOR

RETURN,filecool_arr
END

PRO compCooling,dirs,key = key,outplot = outplot,color = color,sym = sym,dosteps = dosteps,thick = thick, linestyle = linestyle

IF KEYWORD_SET (outplot) THEN BEGIN
    loadct,39
    set_plot,'ps'
    !P.CHARTHICK=4
    !X.THICK=4
    !Y.THICK=4
    !p.charsize=2.25
    legendsize = 1.2
    device,filename = outplot,/color,bits_per_pixel= 8,/times
ENDIF ELSE BEGIN
    loadct,39
    set_plot,'x'
    !P.CHARTHICK=1.5
    !X.THICK=1.5
    !Y.THICK=1.5
    !p.charsize=1.0
    !x.charsize=1.5
    !y.charsize=1.5 
    legendsize = 1
    window,0
ENDELSE

IF NOT KEYWORD_SET (color) THEN color = (findgen(N_ELEMENTS(dirs))+1)/N_ELEMENTS(dirs)*240
IF NOT KEYWORD_SET (sym)   THEN sym =  findgen(N_ELEMENTS(dirs))
IF NOT KEYWORD_SET (thick) THEN thick = (findgen(N_ELEMENTS(dirs))+1)
IF NOT KEYWORD_SET (linestyle)   THEN linestyle = (findgen(N_ELEMENTS(dirs))+1)

;IF KEYWORD_SET(dosteps) THEN BEGIN
;    FOR i = 0, N_ELEMENTS(dirs) -1 DO BEGIN
;        spawn,'ls '+dirs[i]+'glass_highd.0*0',steps
;        temps = findgen(N_ELEMENTS(steps))
;        FOR j = 0, N_ELEMENTS(steps) - 1 DO BEGIN
;            rtipsy,steps[j],h,g,d,s
;            print,steps[j]
;            temps[j] = MEAN(g.tempg)
;        ENDFOR  
;        IF i eq 0 THEN plot,temps[*], /ylog, xtitle = 'Output', ytitle = 'Temperature [K]' ;yrange = [9900,15000],psym = sym[0]
;        oplot,temps[*],color = color[i],psym = sym[i]
;    ENDFOR
;    IF KEYWORD_SET (key) THEN legend,key,color = color,psym = sym ELSE legend,STRTRIM(findgen(N_ELEMENTS(dirs)),2),color = color,psym = sym
;ENDIF 

CL_B_gm = (6.022d23*(938.7830/931.494))
CL_Rgascode    =     8.2494e7
CL_Eerg_gm_degK =    CL_Rgascode
CL_ev_degK       =   1.0/1.1604e4
CL_Eerg_gm_ev    =   CL_Eerg_gm_degK/CL_ev_degK
CL_Eerg_gm_degK3_2 = 1.5*CL_Eerg_gm_degK

cool_arr0 = readCoolDebug(dirs[0]+'cooldebug_table.txt')
plot,cool_arr0.T,-1.0/cool_arr0.en_B/CL_B_gm*cool_arr0.edot,/xlog,/ylog,xtitle = 'Temperature [K]',ytitle = 'Edot',xrange = [1e1,1e8],yrange = [1e-30,1e-20],xstyle = 1, ystyle = 1,psym = sym[0],linestyle = linestyle[0], thick = thick[0]
oplot,cool_arr0.T,-1.0/cool_arr0.en_B/CL_B_gm*cool_arr0.edot,psym = sym[0], color = color[0],linestyle = linestyle[0], thick = thick[0]

FOR i = 1, N_ELEMENTS(dirs) - 1 DO BEGIN
    CASE i OF
        1: BEGIN
            cool_arr1 = readCoolDebug(dirs[i]+'cooldebug_table.txt')
            oplot,cool_arr1.T,-1.0/cool_arr1.en_B/CL_B_gm*cool_arr1.edot,psym = sym[i], color = color[i],linestyle = linestyle[i], thick = thick[i]
        END
        2: BEGIN
            cool_arr2 = readCoolDebug(dirs[i]+'cooldebug_table.txt')
            oplot,cool_arr2.T,-1.0/cool_arr2.en_B/CL_B_gm*cool_arr2.edot,psym = sym[i], color = color[i],linestyle = linestyle[i], thick = thick[i]
        END
        3: BEGIN
            cool_arr3 = readCoolDebug(dirs[i]+'cooldebug_table.txt')
            oplot,cool_arr3.T,-1.0/cool_arr3.en_B/CL_B_gm*cool_arr3.edot,psym = sym[i], color = color[i],linestyle = linestyle[i], thick = thick[i]
        END
        4: BEGIN
            cool_arr4 = readCoolDebug(dirs[i]+'cooldebug_table.txt')
            oplot,cool_arr4.T,-1.0/cool_arr4.en_B/CL_B_gm*cool_arr4.edot,psym = sym[i], color = color[i],linestyle = linestyle[i], thick = thick[i]
        END
        5: BEGIN
            cool_arr5 = readCoolDebug(dirs[i]+'cooldebug_table.txt')
            oplot,cool_arr5.T,-1.0/cool_arr5.en_B/CL_B_gm*cool_arr5.edot,psym = sym[i], color = color[i],linestyle = linestyle[i], thick = thick[i]
        END
        6: BEGIN
            cool_arr6 = readCoolDebug(dirs[i]+'cooldebug_table.txt')
            oplot,cool_arr6.T,-1.0/cool_arr6.en_B/CL_B_gm*cool_arr6.edot,psym = sym[i], color = color[i],linestyle = linestyle[i], thick = thick[i]
        END
        7: BEGIN
            cool_arr7 = readCoolDebug(dirs[i]+'cooldebug_table.txt')
            oplot,cool_arr7.T,-1.0/cool_arr7.en_B/CL_B_gm*cool_arr7.edot,psym = sym[i], color = color[i],linestyle = linestyle[i], thick = thick[i]
        END
        8: BEGIN
            cool_arr8 = readCoolDebug(dirs[i]+'cooldebug_table.txt')
            oplot,cool_arr8.T,-1.0/cool_arr8.en_B/CL_B_gm*cool_arr8.edot,psym = sym[i], color = color[i],linestyle = linestyle[i], thick = thick[i]
        END
        9: BEGIN
            cool_arr9 = readCoolDebug(dirs[i]+'cooldebug_table.txt')
            oplot,cool_arr9.T,-1.0/cool_arr9.en_B/CL_B_gm*cool_arr9.edot,psym = sym[i], color = color[i],linestyle = linestyle[i], thick = thick[i]
        END
        10: BEGIN
            cool_arr10 = readCoolDebug(dirs[i]+'cooldebug_table.txt')
            oplot,cool_arr10.T,-1.0/cool_arr10.en_B/CL_B_gm*cool_arr10.edot,psym = sym[i], color = color[i],linestyle = linestyle[i], thick = thick[i]
        END
        11: BEGIN
            cool_arr11 = readCoolDebug(dirs[i]+'cooldebug_table.txt')
            oplot,cool_arr11.T,-1.0/cool_arr11.en_B/CL_B_gm*cool_arr11.edot,psym = sym[i], color = color[i],linestyle = linestyle[i], thick = thick[i]
        END
        12: BEGIN
            cool_arr12 = readCoolDebug(dirs[i]+'cooldebug_table.txt')
            oplot,cool_arr12.T,-1.0/cool_arr12.en_B/CL_B_gm*cool_arr12.edot,psym = sym[i], color = color[i],linestyle = linestyle[i], thick = thick[i]
        END
        13: BEGIN
            cool_arr13 = readCoolDebug(dirs[i]+'cooldebug_table.txt')
            oplot,cool_arr13.T,-1.0/cool_arr13.en_B/CL_B_gm*cool_arr13.edot,psym = sym[i], color = color[i],linestyle = linestyle[i], thick = thick[i]
        END
        14: BEGIN
            cool_arr14 = readCoolDebug(dirs[i]+'cooldebug_table.txt')
            oplot,cool_arr14.T,-1.0/cool_arr14.en_B/CL_B_gm*cool_arr14.edot,psym = sym[i], color = color[i],linestyle = linestyle[i], thick = thick[i]
        END
        15: BEGIN
            cool_arr15 = readCoolDebug(dirs[i]+'cooldebug_table.txt')
            oplot,cool_arr15.T,-1.0/cool_arr15.en_B/CL_B_gm*cool_arr15.edot,psym = sym[i], color = color[i],linestyle = linestyle[i], thick = thick[i]
        END
        16: BEGIN
            cool_arr16 = readCoolDebug(dirs[i]+'cooldebug_table.txt')
            oplot,cool_arr16.T,-1.0/cool_arr16.en_B/CL_B_gm*cool_arr16.edot,psym = sym[i], color = color[i],linestyle = linestyle[i], thick = thick[i]
        END
        17: BEGIN
            cool_arr17 = readCoolDebug(dirs[i]+'cooldebug_table.txt')
            oplot,cool_arr17.T,-1.0/cool_arr17.en_B/CL_B_gm*cool_arr17.edot,psym = sym[i], color = color[i],linestyle = linestyle[i], thick = thick[i]
        END
        18: BEGIN
            cool_arr18 = readCoolDebug(dirs[i]+'cooldebug_table.txt')
            oplot,cool_arr18.T,-1.0/cool_arr18.en_B/CL_B_gm*cool_arr18.edot,psym = sym[i], color = color[i],linestyle = linestyle[i], thick = thick[i]
        END
        19: BEGIN
            cool_arr19 = readCoolDebug(dirs[i]+'cooldebug_table.txt')
            oplot,cool_arr19.T,-1.0/cool_arr19.en_B/CL_B_gm*cool_arr19.edot,psym = sym[i], color = color[i],linestyle = linestyle[i], thick = thick[i]
        END
        ENDCASE
ENDFOR
IF KEYWORD_SET (key) THEN legend,key,color = color,thick = thick,linestyle = linestyle,psym = sym,/bottom,/right,charsize = legendsize ELSE legend,STRTRIM(findgen(N_ELEMENTS(dirs)),2),color = color,psym = sym,thick = thick,linestyle = linestyle,/bottom,/left,charsize = legendsize
IF KEYWORD_SET (outplot) THEN device,/close

;window,1
;plot,cool_arr0.energy,-1.0*cool_arr0.MCool,/ylog,/xlog,xrange = [1e12,1e15]
;oplot,cool_arr1.energy,-1.0*cool_arr1.MCool,color = color[1],linestyle = 2
stop
END



PRO compCooling_master
root = "/astro/net/scratch1/christensen/MolecH/ShockTube/st_glass/"
root = "/astro/net/scratch1/jillian/runforcharlotte/"
root = "/astro/net/scratch1/christensen/MolecH/ShockTube/st_glass/glass_highd_codetest/"

dirs = ['glass_highd_noH2_metal/', 'glass_highd_noH2cool/','glass_highd_noH2cool_noMcool/', 'glass_highd_noMcool/', 'glass_highd_1.28.09/']
dirs = ['glass_highd_noH2_metal_highT/','glass_highd_noH2cool_highT/','glass_highd_2.4.10/','glass_highd_noH2_z0/','glass_highd_noH2cool_z0/','glass_highd_z0/']
dirs = ['glass_highd_1.3/','glass_highd_1.2/']
dirs = ['H2/again/','noH2/again/']
dirs = ['glass_highd_1.3/','glass_highd_codetest/']
dirs = ['zsol0.01_rho1_H2/','zsol0.01_rho8_H2/','zsol0.1_rho64_H2/','zsol_rho1_H2/','zsol_rho8_H2/','zsol0.01_rho64_H2/','zsol0.1_rho1_H2/','zsol0.1_rho8_H2/','zsol_rho64_H2/']
dirs = ['zsol0.01_rho1_H2/','zsol0.01_rho8_H2/','zsol0.01_rho64_H2/', $
        'zsol0.1_rho1_H2/','zsol0.1_rho8_H2/','zsol0.1_rho64_H2/',    $
        'zsol_rho1_H2/','zsol_rho8_H2/','zsol0.1_rho64_H2/']
dirs = ['zsol0.01_rho1/','zsol0.01_rho8/','zsol0.01_rho64/', $
        'zsol0.1_rho1/','zsol0.1_rho8/','zsol0.1_rho64/',    $
        'zsol_rho1/','zsol_rho8/','zsol0.1_rho64/']
dirs = ['zsol0.01_rho1/','zsol0.01_rho1_H2/','zsol0.1_rho1/','zsol0.1_rho1_H2/']
dirs = ['zsol0.01_rho1_H2/','zsol0.01_rho8_H2/','zsol0.01_rho64_H2/', $
        'zsol0.1_rho1_H2/','zsol0.1_rho8_H2/','zsol0.1_rho64_H2/',    $
        'zsol_rho1_H2/','zsol_rho8_H2/','zsol_rho64_H2/', $
        'zsol0.01_rho1/','zsol0.01_rho8/','zsol0.01_rho64/', $
        'zsol0.1_rho1/','zsol0.1_rho8/','zsol0.1_rho64/',    $
        'zsol_rho1/','zsol_rho8/','zsol_rho64/']
;dirs = ['zsol0.01_rho8_H2/','zsol0.1_rho64_H2/','zsol_rho1_H2/','zsol_rho8_H2/','zsol0.01_rho64_H2/','zsol0.1_rho8_H2/','zsol_rho64_H2/']
dirs = root + dirs

key = ['no H2, z_sol', 'no H2 cool, z_sol', 'H2, z_sol', 'no H2, z0', 'no H2 cool, z0','H2, z0']
key = ['After H2','Before H2']
key = ['1.3','New Code']
key = ['zsol0.01_rho1_H2','zsol0.01_rho8_H2','zsol0.01_rho64_H2', $
       'zsol0.1_rho1_H2','zsol0.1_rho8_H2','zsol0.1_rho64_H2',    $
       'zsol_rho1_H2','zsol_rho8_H2','zsol0.1_rho64_H2']
key = ['zsol0.01_rho1','zsol0.01_rho8','zsol0.01_rho64', $
       'zsol0.1_rho1','zsol0.1_rho8','zsol0.1_rho64',    $
       'zsol_rho1','zsol_rho8','zsol0.1_rho64']
key = ['zsol0.01_rho1/','zsol0.01_rho1_H2/','zsol0.1_rho1/','zsol0.1_rho1_H2/']
key = ['zsol0.01_rho1_H2','zsol0.01_rho8_H2','zsol0.01_rho64_H2', $
       'zsol0.1_rho1_H2','zsol0.1_rho8_H2','zsol0.1_rho64_H2',    $
       'zsol_rho1_H2','zsol_rho8_H2','zsol_rho64_H2', $
       'zsol0.01_rho1','zsol0.01_rho8','zsol0.01_rho64', $
       'zsol0.1_rho1','zsol0.1_rho8','zsol0.1_rho64',    $
       'zsol_rho1','zsol_rho8','zsol_rho64']
;key = ['zsol0.01_rho8_H2','zsol0.1_rho64_H2','zsol_rho1_H2','zsol_rho8_H2','zsol0.01_rho64_H2','zsol0.1_rho8_H2','zsol_rho64_H2']

outplot = '/astro/net/scratch1/christensen/MolecH/ShockTube/st_glass/glass_highd_1.3/tarfile/cooling_curve.eps'
outplot = '/astro/net/scratch1/christensen/MolecH/ShockTube/st_glass/glass_highd_codetest/cooling_curve.eps'
outplot = '/astro/net/scratch1/christensen/MolecH/ShockTube/st_glass/glass_highd_codetest/cooling_curve.eps'

sym = [-2,-1]
sym = [3,-2,-1,3,-2,-1,-4,-2,-1]
;sym = [-2,-1,-3,-2,-1,-2,-1]
sym = [1,2,4,5]
sym = fltarr(N_ELEMENTS(key))
color = [fltarr(3)+240,fltarr(3)+100,fltarr(3)+50,fltarr(3)+240,fltarr(3)+100,fltarr(3)+50]
linestyle = [REFORM(TRANSPOSE([[fltarr(3)+1],[fltarr(3)+2],[fltarr(3)]]),9),REFORM(TRANSPOSE([[fltarr(3)+1],[fltarr(3)+2],[fltarr(3)]]),9)]
thick = [fltarr(9)+2,fltarr(9)+1]
compCooling,dirs,key = key,sym = sym,color = color,linestyle = linestyle,thick = thick;,outplot = outplot

END
