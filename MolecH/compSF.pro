;This program compares the 
;1) phase diagram where the stars are formed
;2) the stellar and gas profile
;3) the star formation history

; .r /astro/users/christensen/code/idl4tipsy/rstarlog.pro

PRO compSF_master
loadct,39
dir = '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/'
dir = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/'


;tfiles = dir + ['h516.cosmo25cmb.1536g3HBWK/steps/h516.cosmo25cmb.1536g3HBWK.00240.dir/h516.cosmo25cmb.1536g3HBWK.00240.halo.1','h516.cosmo25cmb.1536g6MbwK/steps/h516.cosmo25cmb.1536g6MbwK.00240.dir/h516.cosmo25cmb.1536g6MbwK.00240.halo.1']
;tfiles = dir + ['h516.cosmo25cmb.1536g3HBWK/steps/h516.cosmo25cmb.1536g3HBWK.00512.dir/h516.cosmo25cmb.1536g3HBWK.00512.halo.1',$
;        'h516.cosmo25cmb.1536g3HBWK/steps_noH2SF/h516.cosmo25cmb.1536g3HBWK_noH2SF.00512.dir/h516.cosmo25cmb.1536g3HBWK_noH2SF.00512.halo.1',$
;        'h516.cosmo25cmb.1536g6MbwK/steps/h516.cosmo25cmb.1536g6MbwK.00512.dir/h516.cosmo25cmb.1536g6MbwK.00512.halo.1']
;tfiles = dir + ['h516.cosmo25cmb.1536g3HBWK/steps/h516.cosmo25cmb.1536g3HBWK.00324.dir/h516.cosmo25cmb.1536g3HBWK.00324.halo.1',$
;        'h516.cosmo25cmb.1536g3HBWK/steps_noH2SF/h516.cosmo25cmb.1536g3HBWK_noH2SF.00324.dir/h516.cosmo25cmb.1536g3HBWK_noH2SF.00324.halo.1',$
;        'h516.cosmo25cmb.1536g6MbwK/steps/h516.cosmo25cmb.1536g6MbwK.00328.dir/h516.cosmo25cmb.1536g6MbwK.00328.halo.1']
;tfiles = dir + ['h516.cosmo25cmb.1536g14HBWK/steps/h516.cosmo25cmb.1536g14HBWK.00512.dir/h516.cosmo25cmb.1536g14HBWK.00512.halo.1', $
;               'h516.cosmo25cmb.1536g1MBWK/steps/h516.cosmo25cmb.1536g1MBWK.00512.dir/h516.cosmo25cmb.1536g1MBWK.00512.halo.1']
tfiles = dir + ['h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.00492.dir/h516.cosmo25cmb.3072g1MBWK.00492.halo.1',$
                'h516.cosmo25cmb.2304g/h516.cosmo25cmb.2304g14HBWK/steps/h516.cosmo25cmb.2304g14HBWK.00512.dir/h516.cosmo25cmb.2304g14HBWK.00512.halo.1']

;tfiles = dir + ['h516.cosmo25cmb.1536g3HBWK/steps/h516.cosmo25cmb.1536g3HBWK.00324.dir/h516.cosmo25cmb.1536g3HBWK.00324.halo.1','h516.cosmo25cmb.1536g6MbwK/steps/h516.cosmo25cmb.1536g6MbwK.00328.dir/h516.cosmo25cmb.1536g6MbwK.00328.halo.1']
;slfiles = dir + ['h516.cosmo25cmb.1536g3HBWK/h516.cosmo25cmb.1536g3HBWK.starlog','h516.cosmo25cmb.1536g6MbwK/h516.cosmo25cmb.1536g6MbwK.starlog']
;slfiles = dir + ['h516.cosmo25cmb.1536g3HBWK/h516.cosmo25cmb.1536g3HBWK','h516.cosmo25cmb.1536g3HBWK/h516.cosmo25cmb.1536g3HBWK_noH2SF','h516.cosmo25cmb.1536g6MbwK/h516.cosmo25cmb.1536g6MbwK.starlog']
;slfiles = dir + ['h516.cosmo25cmb.1536g14HBWK/h516.cosmo25cmb.1536g14HBWK.starlog','h516.cosmo25cmb.1536g1MBWK/h516.cosmo25cmb.1536g1MBWK.starlog']
slfiles =  dir + ['h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/h516.cosmo25cmb.3072g1MBWK.starlog',$
                  'h516.cosmo25cmb.2304g/h516.cosmo25cmb.2304g14HBWK/h516.cosmo25cmb.2304g14HBWK.starlog']

;dMsolUnit = [2.310e15,2.310e15 ]
;dKpcUnit = [25000.0,25000.0]
;key = ['H2','No H2']

;dMsolUnit       =  [2.310e15,2.310e15,2.310e15]
;dKpcUnit        = [25000.,25000.,25000.]
;key = ['Standard','No H2 SF','Fabio']

;dMsolUnit       =  [2.310e15, 2.310e15,2.310e15,2.310e15]
;dKpcUnit        = [25000., 25000.,25000.,25000.]
;key = ['new H2, lr','new H2, high T','old H2','no H2']
;color = [240,60]

dMsolUnit       =  [2.310e15, 2.310e15]
dKpcUnit        =  [25000., 25000.]
key             = ['no '+textoidl('H_2'), textoidl('H_2')]
color = [0,0]

compSF,tfiles,slfiles,dMsolUnit,dKpcUnit,color = color,key = key,outplot = '~/h516.cosmo25cmb.paper',thick = [2,5]
stop

prefix = '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/'
base0 = "h516.cosmo25cmb.1536g3HBWK"
base1 = "h516.cosmo25cmb.1536g6MbwK"
base2 = "h516.cosmo25cmb.1536g6HBWK.jeans.prev"
base3 = "h516.cosmo25cmb.1536g3HBWK" 
steps = [12,24,36,48,60,72,84,96,108,120,132,180,192,204,216,228,240,252,264,276,288,300,312,324,336,348,360,372,384,396,408,420,432,444,456,464]
step0 = ['00195','00240','00324','00384','00408','00456','00480','00512']
step1 = ['00192','00240','00328','00384','00406','00455','00480','00512'];192
step2 = ['00144','00192','00276','00336','00360','00408','00432','00464'];144
step3 = ['00192','00240','00324','00384','00408','00456','00480','00512']

nsteps = N_ELEMENTS(step0) - 1

cd,prefix
start = 0
dt = 1
slfile0 = "h516.cosmo25cmb.1536g3HBWK/steps/h516.cosmo25cmb.1536g3HBWK"
slfile1 =  "h516.cosmo25cmb.1536g6MbwK/steps/h516.cosmo25cmb.1536g6MbwK"
slfile2 = "h516.cosmo25cmb.1536g6HBWK/Jeans_oldLW/h516.cosmo25cmb.1536g6HBWK.jeans.prev"
slfile3 = "h516.cosmo25cmb.1536g3HBWK/steps_noH2SF/h516.cosmo25cmb.1536g3HBWK_noH2SF"
slfile = ["h516.cosmo25cmb.1536g6HBWK/Jeans_oldLW/h516.cosmo25cmb.1536g6HBWK.jeans.prev"]
slfile = [slfile1,slfile3,slfile0,slfile2]
key = 'fracLim'
key = ['no H2','H2','H2 SF','MC SF']
dMsolUnit       =  (fltarr(N_ELEMENTS(slfile))+ 1)*2.310e15
dKpcUnit        =  (fltarr(N_ELEMENTS(slfile))+ 1)*25000.


FOR i = start/dt, nsteps DO BEGIN 
;    step = i*dt
;    step = steps[i]
;    if (step lt 10) THEN step = '0000'+STRTRIM(step,2) ELSE BEGIN
;        if (step lt 100) THEN step = '000'+STRTRIM(step,2) ELSE step = '00'+STRTRIM(step,2)
;    ENDELSE
;    filename0 = base+"/Jeans_oldLW/steps/"+base+"."+step+".dir/"+base + "." + step+'.halo.1';+STRTRIM(halo,2)

    filename0 = base0 + "/steps/"+base0+"."+step0[i]+".dir/"+base0 + "." + step0[i]+'.halo.1'
    filename1 = base1 + "/steps/"+base1+"."+step1[i]+".dir/"+base1 + "." + step1[i]+'.halo.1'
    filename2 = "h516.cosmo25cmb.1536g6HBWK/Jeans_oldLW/steps/"+base2+"."+step2[i]+".dir/"+base2 + "." + step2[i]+'.halo.1';+STRTRIM(halo,2)
    filename3 = base3 + "/steps_noH2SF/"+base3+"_noH2SF."+step3[i]+".dir/"+base0 + "_noH2SF." + step3[i]+'.halo.1'
    filename = [filename1,filename3,filename0,filename2]
    compSF,filename,slfile,dMsolUnit,dKpcUnit,key =key,color = [30,90,140,240],outplot = 'stellarprof.'+step3[i];,redshift = redshift[i];,outfile = outfile,/color
;    stop
ENDFOR
END

PRO compSF,tfiles,slfiles,dMsolUnit,dKpcUnit,key = key,outplot = outplot,color = color,thick = thick
!Y.STYLE = 1
!X.STYLE = 1
!P.THICK = 3.5
if (KEYWORD_SET(outplot)) then begin
    set_plot,'ps' 
    !P.CHARTHICK=4
    !X.THICK=4
    !Y.THICK=4
    !P.charsize=1.0
    !x.charsize=1.5
    !y.charsize=1.5
    l_charsize = 0.75
    !X.MARGIN = [12,3]
    !Y.MARGIN = [6,2]
    device,filename = outplot+'_sprof.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 12,xoffset =  2,yoffset =  2 
endif else begin
    set_plot,'x'
    !P.CHARTHICK=1.5
    !X.THICK=1.5
    !Y.THICK=1.5
    !p.charsize=1.0
    !x.charsize=1.5
    !y.charsize=1.5  
    l_charsize = 1.0
    !X.MARGIN = [10,3]
    !Y.MARGIN = [4,2]
    window,0,xsize = 600,ysize = 600
endelse

hubble = 73.0
cm_per_kpc = 3.08568021d21
gm_per_msol = 1.98892d33
amu_per_gm = 6.022142d23
gm_per_H = 1.673534d-24
H_per_gm = 5.9753790e+23
dDelta = 0.003
sec_per_year = 31556926
molec_weight = (0.76*1 + 0.24*4.0)
IF NOT KEYWORD_SET(color) THEN color = (findgen(N_ELEMENTS(slfiles)) + 1)/N_ELEMENTS(slfiles)*240
IF NOT KEYWORD_SET(thick) THEN thick = fltarr(N_ELEMENTS(slfiles)) + 1
;set_plot,'x'
;!p.multi = [0,1,2]
!p.multi = 0
FOR i = 0, N_ELEMENTS(slfiles) - 1 DO BEGIN
    timeunit = SQRT((dKpcUnit[i]*3.086d21)^3/(6.67d-8*dMsolUnit[i]*1.99d33))/(3600.*24.*365.24)
    rtipsy,tfiles[i],h,g,d,s
    a = h.time
    dKpcUnit_a = dKpcUnit[i]*a
    dens_convert =  dMsolUnit[i] * gm_per_msol * 5.9753790e+23/dKpcUnit[i]^3/cm_per_kpc^3
    g.x = g.x * dKpcUnit[i]
    g.y = g.y * dKpcUnit[i]
    g.z = g.z * dKpcUnit[i]
    g.mass = g.mass*dMsolUnit[i]
    s.x = s.x * dKpcUnit[i]
    s.y = s.y * dKpcUnit[i]
    s.z = s.z * dKpcUnit[i]
    s.mass = s.mass*dMsolUnit[i]
    s.tform = s.tform*timeunit
    dens = g.dens*dens_convert/a^3
     IF i eq 0 then begin
         data = rstarlog(slfiles[i]) 
;         stop
     endif else data = rstarlog(slfiles[i],/molecularH)
    data = data[where(data.timeform*timeunit le MAX(s.tform))]
    rhoform = data.rhoform*dens_convert
;    IF (i EQ 0) THEN plot,dens,g.tempg,psym = 3,/xlog,/ylog,xtitle = 'Density',ytitle = 'Temperature',xrange = [1e-7,1e4],yrange = [1e2,1e7],xstyle = 1,ystyle = 1
;    oplot,dens,g.tempg,psym = 3,color = color[i]
;    oplot,rhoform,data.tempform,psym = 1,color = color[i]
ENDFOR
IF NOT KEYWORD_SET(key) then key=slfiles
;legend,key,psym = (fltarr(N_ELEMENTS(slfiles)) + 1.0),color =color,/right,charsize = l_charsize

FOR i = 0, N_ELEMENTS(slfiles) - 1 DO BEGIN
    rtipsy,tfiles[i],h,g,d,s
    g.x = g.x * dKpcUnit[i]
    g.y = g.y * dKpcUnit[i]
    g.z = g.z * dKpcUnit[i]
    g.mass = g.mass*dMsolUnit[i]
    s.x = s.x * dKpcUnit[i]
    s.y = s.y * dKpcUnit[i]
    s.z = s.z * dKpcUnit[i]
    s.mass = s.mass*dMsolUnit[i]
    s.tform = s.tform*SQRT((dKpcUnit[i]*3.086d21)^3/(6.67d-8*dMsolUnit[i]*1.99d33))/(3600.*24.*365.24)
    syoung = s[where(s.tform gt MAX(s.tform) - 1e9)]
;    smed = s[where(s.tform gt MAX(s.tform) - 4e9)]
;    sold = s[where(s.tform gt MAX(s.tform) - 6e9)]    
    prof_s = prof(s,'star',MAX(s.tform),nbins = 50,rmax = 25)
    prof_g = prof(g,'gas',MAX(s.tform),nbins = 50,rmax = 25)
    prof_syoung = prof(syoung,'star',MAX(s.tform),nbins = 50,rmax = 25)
;    prof_smed = prof(smed,'star',MAX(s.tform),nbins = 50,rmax = 25)
;    prof_sold = prof(sold,'star',MAX(s.tform),nbins = 50,rmax = 25)
    IF (i EQ 0) THEN plot,prof_s.rbins,prof_s.rho,/ylog,xtitle = 'Radius [kpc]',ytitle = textoidl('Stellar Density [M')+sunsymbol()+textoidl(' kpc^{-2}]'),yrange = [1e2,1e8],xrange = [0,8],thick = thick[i]
    oplot,prof_s.rbins,prof_s.rho,color = color[i],thick = thick[i]
;    oplot,prof_g.rbins,prof_g.rho,color = color[i],linestyle = 2
    oplot,prof_syoung.rbins,prof_syoung.rho,linestyle = 3,color = color[i],thick = thick[i];color = color[i]-20,thick = 0.5
;    oplot,prof_smed.rbins,prof_smed.rho,color = color[i]-20,thick = 0.5,linestyle = 4
;    oplot,prof_sold.rbins,prof_sold.rho,color = color[i]-20,thick = 0.5,linestyle = 5
    print,tfiles[i]
    print,'Stellar Mass: ',TOTAL(s.mass)
    print,'Gas Mass:     ',TOTAL(g.mass)
    print,''
ENDFOR
legend,['All Stars',textoidl('Age < 10^9 years')],linestyle = [0,3],/right,charsize = l_charsize
legend,[key],linestyle = (fltarr(N_ELEMENTS(tfiles))),color = color,thick = thick,charsize = l_charsize
;legend,['Stars: ','Gas'],linestyle = [0,2],/bottom

IF keyword_set(outplot) then begin
    device,/close
    device,filename = outplot+'_shist.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 12,xoffset =  2,yoffset =  2 
endif

FOR i = 0, N_ELEMENTS(slfiles) - 1 DO BEGIN
    rtipsy,tfiles[i],h,g,d,s
 ;   g.x = g.x * dKpcUnit[i]
 ;   g.y = g.y * dKpcUnit[i]
 ;   g.z = g.z * dKpcUnit[i]
 ;   g.mass = g.mass*dMsolUnit[i]
 ;   s.x = s.x * dKpcUnit[i]
 ;   s.y = s.y * dKpcUnit[i]
 ;   s.z = s.z * dKpcUnit[i]
 ;   s.mass = s.mass*dMsolUnit[i]
 ;   s.tform = s.tform*SQRT((dKpcUnit[i]*3.086d21)^3/(6.67d-8*dMsolUnit[i]*1.99d33))/(3600.*24.*365.24)
    IF (i eq 0) THEN sfr,s,massunit = dMsolUnit[i],timeunit = timeunit,binsize = 1e8,yrange = [0,MAX(s.tform)],xstyle = 1,xrange = [0,13],thick = thick[0]
    sfr,s,/overplot,color = color[i],massunit = dMsolUnit[i],timeunit = timeunit,binsize = 1e8,xstyle = 1,thick = thick[i]
;    print,TOTAL(s[where(s.tform*timeunit ge 1.3634376e+10 AND s.tform*timeunit le 1.3734376e+10)].mass)/TOTAL(s.mass)/1e8
;    print,TOTAL(s.mass)*dMsolUnit[i]
ENDFOR
legend,key,linestyle = (fltarr(N_ELEMENTS(slfiles))),color = color,/left,thick = thick
IF keyword_set(outplot) then  device,/close else stop
;wait,1
close,/all
END
