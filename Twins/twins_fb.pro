pro twins_fb,dir,files,pfiles,outplot = outplot,linestyle = linestyle, color = color, psym = psym, symsize = symsize, thick = thick,scale = scale,keys = keys,vline = vline
formatplot,outplot = outplot
;window,xoffset = 2,yoffset = 2
;color = [50,245]
color = [30,245]
highz_time = 1.93996e+09/1e9 ;dtime = 5.002289e-02 ;72, z = 3.43597, a = 0.22542965
;highz_time = 2.6905329E+00; 100

ZSOLAR  =  0.0130215
filebase = 'masstrack.dat'
n = N_ELEMENTS(dir)
IF KEYWORD_SET(outplot) THEN BEGIN 
    fgcolor = 0
    bgcolor = 255
ENDIF ELSE BEGIN
    fgcolor = 255
    bgcolor = 0
ENDELSE
IF KEYWORD_SET(color) THEN BEGIN
    loadct,39
    if color[0] eq 1 then  colors = (findgen(n) + 1)*240/n else colors = color
    IF NOT KEYWORD_SET(thick) THEN thick = fltarr(n) + 2;1
    IF NOT KEYWORD_SET(linestyle) THEN linestyle = fltarr(n) 
    IF NOT KEYWORD_SET(psym) THEN psym =-1*( fltarr(n) + 4)
    IF NOT KEYWORD_SET(symsize) THEN symsize = fltarr(n) + 2
ENDIF ELSE BEGIN
    loadct,0    
    colors = fltarr(n) + fgcolor
    IF NOT KEYWORD_SET(thick) THEN thick = fltarr(n) + 2;1 ;findgen(n)*2
    IF NOT KEYWORD_SET(linestyle) THEN linestyle = findgen(n)*2   
    IF NOT KEYWORD_SET(psym) THEN psym = -1*(fltarr(n) + 4)
    IF NOT KEYWORD_SET(symsize) THEN symsize = fltarr(n) + 2
ENDELSE
IF NOT KEYWORD_SET(scale) THEN scale = fltarr(n) + 1.0

ageUniverse = 13.7346*1e9 ;wmap3_lookback(100)
z = REVERSE(findgen(100)*10.0/100.0)
t = ageUniverse - wmap3_lookback(z)
ticktime_in_z = findgen(7)*2*1e9 ;redshift major plot
tickred_in_z = spline(t,z,ticktime_in_z )

tickred_in_t = REVERSE(findgen(5))
ticktime_in_t = (ageUniverse - wmap3_lookback(tickred_in_t))/1e9

;------------------ FB ---------------------
if (KEYWORD_SET(outplot)) then begin
    device,filename = outplot + '_fb.eps',/color,bits_per_pixel= 8,/times,ysize = 16,xsize = 18,xoffset =  2,yoffset =  2 
endif else begin
    window,9,xsize = 712,ysize = 600
endelse
binsize = 5e8
nbins = ROUND(14.0/(binsize/1e9))

FOR i = 0, n - 1 DO BEGIN
    cd,dir[i]
    readcol,filebase,halo,time,z,mtot,mgas,mstar,mdark,metal,mcoldg,mH2,H2frac,Pressure,r25,SFR
    fbase = strmid(files[i],0,strlen(files[i]) - 7)
    rtipsy,fbase,h,g,d,s,/justhead
    readarr,fbase + '.iord',h,iord,/ascii,part = 'gas'
    readarr,fbase + '.amiga.grp',h,grp,/ascii,part = 'gas'
    readarr,fbase + '.coolontime',h,coolon,/ascii,part = 'gas'

;------------------ Outflow -------------------------
    outflowm_all = mrdfits(dir[i] + '../grp1.rvir.mass_at_outflow.fits')
    outflowiord_all = mrdfits(dir[i] + '../grp1.rvir.outflow_iord.fits')
    outflowz_all = mrdfits(dir[i] + '../grp1.rvir.outflow_z.fits')
    
;------------------ Select for the gas that leaves the galaxy ----
    outflowm = outflowm_all[where(outflowz_all ne 99.0)]
    outflowiord = outflowiord_all[where(outflowz_all ne 99.0)]
    outflowz = outflowz_all[where(outflowz_all ne 99.0)]

;------------------ Select for gas that is heated by supernovae -----
    match,iord,outflowiord,ind_iord,ind_iordout
    coolonout = where(coolon[ind_iord] ne 0)
    outflow_m = outflowm[ind_iordout[coolonout]]
    outflow_iord = outflowiord[ind_iordout[coolonout]]
    outflow_z = outflowz[ind_iordout[coolonout]]
    outflow_t = (ageUniverse - wmap3_lookback(outflow_z))/1e9

;    match,iord,outflowiord,iordind,ofiordind
;    iord0 = iord[iordind];Particles once accreted that still exist
;    grp0 = grp[iordind]
;    outflowiord0 = outflowiord[ofiordind]
;    outflowm0 =    outflowm[ofiordind]
;    outflowt0 =    outflowt[ofiordind]
;    outflowm00 =   outflowm0[where(grp0 ne 1)]
;    outflowt00 =   outflowt0[where(grp0 ne 1)]
    
;    outflowmbin = weighted_histogram(outflowt00,weight = outflowm00,min = 0.000001, max = 14.0,nbins = nbins,locations = outflowtbin)
    outflowmbin = weighted_histogram(outflow_t,weight = outflow_m,min = 0.000001, max = 14.0,nbins = nbins,locations = outflowtbin)
    cumarr=fltarr(nbins,nbins)
    ii=lindgen(nbins*nbins)
    cumarr(where(ii mod nbins le ii/nbins)) = 1
    cumoutflow=reform(outflowmbin # cumarr)
    ii=0
    cumarr=0

; ------------------- Inflow -----------------------------
    inflow_z = mrdfits(dir[i] + '../grp1.accrz.fits')
    inflow_m = mrdfits(dir[i] + '../grp1.mass_at_accr.fits')
    inflow_t = (ageUniverse - wmap3_lookback(inflowz))/1e9     
    inflowmbin  = weighted_histogram(inflow_t, weight = inflow_m, min = 0.000001, max = 14.0,nbins = nbins,locations = inflowtbin)
    cumarr=fltarr(nbins,nbins)
    ii=lindgen(nbins*nbins)
    cumarr(where(ii mod nbins le ii/nbins)) = 1
    cuminflow =reform(inflowmbin  # cumarr)
    ii=0
    cumarr=0
  
;-------------------- Bin Outflow and Inflows according to the same
;                     steps in the data file
    a_outflowm = fltarr(N_ELEMENTS(mtot))
    a_outflowm_cum = fltarr(N_ELEMENTS(mtot))
    a_inflowm = fltarr(N_ELEMENTS(mtot))
    a_inflowm_cum = fltarr(N_ELEMENTS(mtot))
    a_time = [0,time]

    FOR j = 0, N_ELEMENTS(mtot) - 1 DO BEGIN 
        t0 = a_time[j] 
        t1 = a_time[j + 1]  
        IF (where(outflow_t GT t0 AND outflow_t LE t1))[0] NE -1 THEN a_outflowm[j]     = total(outflow_m[where(outflow_t GT t0 AND outflow_t LE t1)]) 
        IF (where(outflow_t LE t1))[0]                     NE -1 THEN a_outflowm_cum[j] = total(outflow_m[where(outflow_t LE t1)]) 
        IF (where(inflow_t  GT t0 AND inflow_t  LE t1))[0] NE -1 THEN a_inflowm[j]      = total( inflow_m[where(inflow_t  GT t0 AND  inflow_t LE t1)]) 
        IF (where(inflow_t  LE t1))[0]                     NE -1 THEN a_inflowm_cum[j]  = total( inflow_m[where(inflow_t  LE t1)]) 
    ENDFOR
    IF i eq 0 THEN a_array0 = [[time],[a_inflowm],[a_inflowm_cum],[a_outflowm],[a_outflowm_cum]] ELSE a_array1 = [[time],[a_inflowm],[a_inflowm_cum],[a_outflowm],[a_outflowm_cum]]

    IF i EQ 0 THEN BEGIN
        plot,inflowtbin,cuminflow*scale[i],ytitle = 'Gas  Mass [M' + sunsymbol() + ']',xstyle =9,xrange=[0,ageUniverse/1e9],xtitle = 'Time [Gyr]',yrange = [-2.5e10,5e10],/nodata,ystyle = 9,xmargin = [17,17],ymargin = [6,6]
        axis,ystyle = 1,yaxis = 1,ytitle = textoidl('8 X M_{star}') + '[M' + sunsymbol() + ']'
        axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'z_from_time_axes',xtitle = 'z',xtickv = ticktime_in_t,xticks = N_ELEMENTS(ticktime_in_t) - 1        
        oplot,[0,14],[0,0],thick = 2
    ENDIF
    oplot, inflowtbin,             (cuminflow)*scale[i],linestyle = 2,thick = 2,color = colors[i] ;gas mass in
    oplot,outflowtbin,       (-1.0*cumoutflow)*scale[i],linestyle = 3,thick = 2,color = colors[i];gas mass out
    oplot,outflowtbin,(cuminflow - cumoutflow)*scale[i],linestyle = 0,thick = 2,color = colors[i]
;    oplot,time,mgas*scale[i],linestyle = 2,thick = 2,color = colors[i]
;    oplot,time,mstar*scale[i],linestyle = 2,thick = 2,color = colors[i]
;    oplot,time,mstar*scale[i] + mgas*scale[i],linestyle = 0,thick = 2,color = colors[i]
;    oplot,time,cuminflow*scale[i] - (mstar*scale[i] + mgas*scale[i]),linestyle = 1,thick = 2,color = colors[i]
    IF i EQ 0 THEN outflowratio0 = [[outflowtbin],[cumoutflow/cuminflow],[outflowmbin]] ELSE outflowratio1 = [[outflowtbin],[cumoutflow/cuminflow],[outflowmbin]]
ENDFOR
legend,keys,color = colors,linestyle = [0,0],/bottom
legend,['Gas Accreted','Accreted Gas Still in Halo','Gas Lost'],linestyle = [2,0,3]

;------------------ Outflow Ratio ---------------------
if (KEYWORD_SET(outplot)) then begin
    device,/close
    device,filename = outplot + '_fbratio.eps',/color,bits_per_pixel= 8,/times,ysize = 12,xsize = 18,xoffset =  2,yoffset =  2 
endif else begin
    stop
    window,10,xsize = 712,ysize = 392
endelse

FOR i = 0, n - 1 DO BEGIN
    cd,dir[i]
    readcol,filebase,halo,time,z,mtot,mgas,mstar,mdark,metal,mcoldg,mH2,H2frac,Pressure,r25,SFR
    IF i EQ 0 THEN BEGIN
        plot,time,a_array0[*,3]/mgas,ytitle = 'Gas Mass Loss/Gas Mass',xstyle =9,xrange=[0,ageUniverse/1e9],xtitle = 'Time [Gyr]',yrange = [0,0.8],/nodata,xmargin = [17,17],ymargin = [6,6]
        axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'z_from_time_axes',xtitle = 'z',xtickv = ticktime_in_t,xticks = N_ELEMENTS(ticktime_in_t) - 1        
        oplot,time,a_array0[*,3]/mgas,thick = 2,color = colors[i]
;        oplot,time,a_array0[*,3]/mcoldg,thick = 2,color = colors[i]
    ENDIF ELSE oplot,time,a_array1[*,3]/mgas,thick = 2,color = colors[i]
;oplot,time,a_array1[*,3]/mcoldg,thick = 2,color = colors[i]
ENDFOR
legend,keys,color = colors,linestyle = [0,0],/top,/right

; Mass of gas lost versus mass of gas accreted so far
;FOR i = 0, n - 1 DO BEGIN
;    IF i EQ 0 THEN BEGIN
;        plot,outflowratio0[*,0],outflowratio0[*,1],ytitle = 'Gas Mass Loss/Gas Mass Accreted',xstyle =9,xrange=[0,ageUniverse/1e9],xtitle = 'Time [Gyr]',yrange = [0,0.8],/nodata,xmargin = [17,17],ymargin = [6,6]
;        axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'z_from_time_axes',xtitle = 'z',xtickv = ticktime_in_t,xticks = N_ELEMENTS(ticktime_in_t) - 1        
;        oplot,outflowratio0[*,0],outflowratio0[*,1],thick = 2,color = colors[i]
;    ENDIF ELSE oplot,outflowratio1[*,0],outflowratio1[*,1],thick = 2,color = colors[i]
;ENDFOR
;legend,keys,color = colors,linestyle = [0,0],/bottom,/right

;------------------ FB efficiency ---------------------
if (KEYWORD_SET(outplot)) then begin
    device,/close
    device,filename = outplot + '_fbeff.eps',/color,bits_per_pixel= 8,/times,ysize = 12,xsize = 18,xoffset =  2,yoffset =  2 
endif else begin
    stop
    window,12,xsize = 712,ysize = 392
endelse

FOR i = 0, n - 1 DO BEGIN
    cd,dir[i]
    readcol,filebase,halo,time,z,mtot,mgas,mstar,mdark,metal,mcoldg,mH2,H2frac,Pressure,r25,SFR
    dtime = (time - [0,time[0:N_ELEMENTS(time) - 1]])*1e9
    IF i EQ 0 THEN BEGIN
        plot,time,a_array0[*,3]/dtime/SFR,ytitle = 'Gas Mass Lost Per Year/SFR',xstyle =9,xrange=[0,ageUniverse/1e9],xtitle = 'Time [Gyr]',yrange = [0,100],/nodata,xmargin = [17,17],ymargin = [6,6]
        axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'z_from_time_axes',xtitle = 'z',xtickv = ticktime_in_t,xticks = N_ELEMENTS(ticktime_in_t) - 1        
        oplot,time,a_array0[*,3]/dtime/SFR,thick = 2,color = colors[i],psym = 10
    ENDIF ELSE oplot,time,a_array1[*,3]/dtime/SFR,thick = 2,color = colors[i],psym = 10
ENDFOR
legend,keys,color = colors,linestyle = [0,0],/top,/right
if (KEYWORD_SET(outplot)) then device,/close ELSE stop
end

pro twins_fb_master

dir   = ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/steps/',$
         '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/',$
         '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072gs1MbwK/steps/',$
         '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/',$
         '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/',$
         '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/']

dir   = ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/',$
         '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/',$
         '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/']

dir   = ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/',$
;         '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/',$
         '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/']
files = ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.00512.dir/h603.cosmo50cmb.3072g14HBWK.00512.halo.1',$
;         '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/',$
         '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.00512.dir/h516.cosmo25cmb.3072g14HBWK.00512.halo.1']
scale = [1,8]
pfiles = ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/h603.cosmo50cmb.3072g14HBWK.param',$
;         '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/',$
         '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/h516.cosmo25cmb.3072g14HBWK.param']
maxdistance_photo = 6

;dir   = ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/']
;files = ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.00512.dir/h516.cosmo25cmb.3072g14HBWK.00512.halo.1']
;scale = [1]
;pfiles = ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/h516.cosmo25cmb.3072g14HBWK.param']


;densRadius,[],[50000.,25000],[1.84793e16,2.310e15],maxdistance = maxdistance_photo ;,outplot = outplot
;mass_evol,dir,files,pfiles,scale = scale,linestyle = [0,0],psym = [14,14],/color,outplot = '~/plots/twins/h603_h516_evol.talk'
twins_fb,dir,files,pfiles,scale = scale,linestyle = [0,0],psym = [14,14],/color,keys = ['Spiral'+textoidl(' Galaxy'),'Dwarf Galaxy'],outplot = '~/plots/twins/h603_h516_evol_flash'
end
