function time_from_z_axes,axis,index,value
ageUniverse = 13.7346*1e9 ;wmap3_lookback(10000)
time = (ageUniverse - wmap3_lookback(value))/1e9
return,strtrim(ROUND(time),2)
end

function z_from_time_axes,axis,index,value
ageUniverse = 13.7346*1e9 ;wmap3_lookback(10000)
z = REVERSE(findgen(100)*10.0/100.0)
t = ageUniverse - wmap3_lookback(z)
redshift = spline(t,z,value*1e9)
return,strtrim(ROUND(redshift),2)
end

;function read_masstrack
;filebase = 'masstrack.dat'
;readcol,filebase,halo,time,z,mtot,mgas,mstar,mdark,metal,mcoldg,mH2,H2frac,Pressure,r25,SFR
;data = 
;end

;--------------------------------------------------------------------
PRO fb_plots,dir,files,pfiles,outplot = outplot,linestyle = linestyle, color = color, thick = thick,keys = keys, ctables = ctables
highz_time = 1.93996e+09/1e9 ;dtime = 5.002289e-02 ;72, z = 3.43597, a = 0.22542965

filebase = 'masstrack.dat'
n = N_ELEMENTS(dir)

ageUniverse = 13.7346*1e9 ;wmap3_lookback(100)
z = REVERSE(findgen(100)*10.0/100.0)
t = ageUniverse - wmap3_lookback(z)
ticktime_in_z = findgen(7)*2*1e9 ;redshift major plot
tickred_in_z = spline(t,z,ticktime_in_z )
tickred_in_t = REVERSE(findgen(5))
ticktime_in_t = (ageUniverse - wmap3_lookback(tickred_in_t))/1e9

binsize = 5e8
nbins = ROUND(14.0/(binsize/1e9))

pinout = PTRARR(n, /ALLOCATE_HEAP)
pdata = PTRARR(n, /ALLOCATE_HEAP)

;----------------------- READ DATA -------------------------
FOR i = 0, n - 1 DO BEGIN
    cd,dir[i]

;------------------ Outflow -------------------------
    outflowz_all = mrdfits(dir[i] + '../grp1.outflow_z.fits')
    outflowm_all = mrdfits(dir[i] + '../grp1.mass_at_outflow.fits')
    outflowiord_all = mrdfits(dir[i] + '../grp1.outflow_iord.fits')
    
    outflowz = outflowz_all[where(outflowz_all ne 99.0)]
    outflowm = outflowm_all[where(outflowz_all ne 99.0)]
    outflowiord = outflowiord_all[where(outflowz_all ne 99.0)]
    outflowt = (ageUniverse - wmap3_lookback(outflowz))/1e9

    fbase = strmid(files[i],0,strlen(files[i]) - 7)
    rtipsy,fbase,h,g,d,s,/justhead
    readarr,fbase + '.iord',h,iord,/ascii,part = 'gas'
    readarr,fbase + '.amiga.grp',h,grp,/ascii,part = 'gas'

    match,iord,outflowiord,iordind,ofiordind
    iord0 = iord[iordind];Particles once accreted that still exist
    grp0 = grp[iordind]
    outflowiord0 = outflowiord[ofiordind]
    outflowm0 =    outflowm[ofiordind]
    outflowt0 =    outflowt[ofiordind]
    outflowm00 =   outflowm0[where(grp0 ne 1)]
    outflowt00 =   outflowt0[where(grp0 ne 1)]
    
    outflowmbin = weighted_histogram(outflowt00,weight = outflowm00,min = 0.000001, max = 14.0,nbins = nbins,locations = outflowtbin)
    cumarr=fltarr(nbins,nbins)
    ii=lindgen(nbins*nbins)
    cumarr(where(ii mod nbins le ii/nbins)) = 1
    cumoutflow=reform(outflowmbin # cumarr)
    ii=0
    cumarr=0

; ------------------- Inflow -----------------------------
    inflowz = mrdfits(dir[i] + '../grp1.accrz.fits')
    inflowm = mrdfits(dir[i] + '../grp1.mass_at_accr.fits')
    inflowt = (ageUniverse - wmap3_lookback(inflowz))/1e9     
    inflowmbin = weighted_histogram(inflowt,weight = inflowm,min = 0.000001, max = 14.0,nbins = nbins,locations = outflowtbin)
    cumarr=fltarr(nbins,nbins)
    ii=lindgen(nbins*nbins)
    cumarr(where(ii mod nbins le ii/nbins)) = 1
    cuminflow=reform(inflowmbin # cumarr)
    ii=0
    cumarr=0
 
;-------------------- Write inflow and outflow structures
    inout = {time:0d0,inflow:0d0,inflowcum:0d0,outflow:0d0,outflowcum:0d0}
    inout = replicate(inout,N_ELEMENTS(outflowtbin))
    inout.time = outflowtbin
    inout.inflow = inflowmbin
    inout.inflowcum = cuminflow
    inout.outflow = outflowmbin
    inout.outflowcum = cumoutflow
    *pinout[i] = inout

;-------------------- Bin Outflow and Inflows according to the same
;                     steps in the data file
    readcol,filebase,halo,time,z,mtot,mgas,mstar,mdark,metal,mcoldg,mH2,H2frac,Pressure,r25,SFR
    a_outflowm = fltarr(N_ELEMENTS(mtot))
    a_outflowm_cum = fltarr(N_ELEMENTS(mtot))
    a_inflowm = fltarr(N_ELEMENTS(mtot))
    a_inflowm_cum = fltarr(N_ELEMENTS(mtot))
    a_time = [0,time]

    FOR j = 0, N_ELEMENTS(mtot) - 1 DO BEGIN $
&        t0 = a_time[j] $
&        t1 = a_time[j + 1] $ 
&        IF (where(outflowt00 GT t0 and outflowt00 LE t1))[0] NE -1 THEN a_outflowm[j] = total(outflowm00[where(outflowt00 GT t0 and outflowt00 LE t1)]) $
&        IF (where(outflowt00 LE t1))[0] NE -1 THEN a_outflowm_cum[j] = total(outflowm00[where(outflowt00 LE t1)]) $
&        IF (where(inflowt GT t0 and inflowt LE t1))[0] NE -1 THEN a_inflowm[j] = total(inflowm[where(inflowt GT t0 and inflowt LE t1)]) $
&        IF (where(inflowt LE t1))[0] NE -1 THEN a_inflowm_cum[j] = total(inflowm[where(inflowt LE t1)]) $
&    ENDFOR

;---------------------- Write Data to Structure --------------------
    data = {halo:0,time:0d0,z:0d0,mtot:0d0,mgas:0d0,mstar:0d0,mdark:0d0,metal:0d0,mcoldg:0d0,mH2:0d0,H2frac:0d0,Pressure:0d0,r25:0d0,SFR:0d0,inflow:0d0,inflowcum:0d0,outflow:0d0,outflowcum:0d0}
    data = replicate(data,N_ELEMENTS(halo))
    data.halo = halo
    data.time = time
    data.z = z
    data.mtot = mtot
    data.mgas = mgas
    data.mstar = mstar
    data.mdark = mdark
    data.metal = metal
    data.mcoldg = mcoldg
    data.mH2 = mH2
    data.H2frac = data.H2frac
    data.Pressure = Pressure
    data.r25 = r25
    data.SFR = SFR
    data.inflow = a_inflowm
    data.inflowcum = a_inflowm_cum
    data.outflow = a_outflowm
    data.outflowcum = a_outflowm_cum
    *pdata[i] = data
ENDFOR

;----------------------------------- Plots --------------------
;IF KEYWORD_SET(outplot) THEN set_plot,'ps' ELSE set_plot,'x'
formatplot,outplot = outplot,/thick
IF KEYWORD_SET(outplot) THEN BEGIN 
    fgcolor = 0
    bgcolor = 255
ENDIF ELSE BEGIN
    fgcolor = 255
    bgcolor = 0
ENDELSE
IF KEYWORD_SET(color) THEN BEGIN
    loadct,39
    IF color[0] EQ 1 THEN  colors = (findgen(n) + 1)*240/n ELSE colors = color
    IF NOT KEYWORD_SET(thick) THEN thick = fltarr(n) + 2;1
    IF NOT KEYWORD_SET(linestyle) THEN linestyle = fltarr(n) 
ENDIF ELSE BEGIN
    loadct,0    
    colors = fltarr(n) + fgcolor
    IF NOT KEYWORD_SET(thick) THEN thick = fltarr(n) + 2;1 ;findgen(n)*2
    IF NOT KEYWORD_SET(linestyle) THEN linestyle = findgen(n)*2   
ENDELSE
IF (KEYWORD_SET(outplot)) THEN BEGIN
    device,filename = outplot + '_fb.eps',/color,bits_per_pixel= 8,/times,ysize = 16,xsize = 18,xoffset =  2,yoffset =  2 
ENDIF ELSE BEGIN
    window,0,xsize = 712,ysize = 600
ENDELSE

FOR i = 0, n - 1 DO BEGIN
    yrange = [-2.5e9,4e9]
    yrange = [-3e10,6e10]
    loadct,ctables[i]

    IF i EQ 0 THEN BEGIN
        plot,(*pinout[0]).time,(*pinout[0]).inflowcum,ytitle = 'Gas  Mass [M' + sunsymbol() + ']',xstyle =9,xrange=[0,ageUniverse/1e9],xtitle = 'Time [Gyr]',yrange = yrange,/nodata,ystyle = 1 ;,xmargin = [17,17],ymargin = [6,6],ystyle = 1
        axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'z_from_time_axes',xtitle = 'z',xtickv = ticktime_in_t,xticks = N_ELEMENTS(ticktime_in_t) - 1        
        oplot,[0,14],[0,0],thick = thick[i]
    ENDIF
    oplot,(*pinout[i]).time,(*pinout[i]).inflowcum,linestyle = 2,thick = thick[i],color = colors[i] ;gas mass in
    oplot,(*pinout[i]).time,-1.0*(*pinout[i]).outflowcum,linestyle = 0,thick = thick[i],color = colors[i];gas mass out
    oplot,(*pinout[i]).time,((*pinout[i]).inflowcum - (*pinout[i]).outflowcum),linestyle = 3,thick = thick[i],color = colors[i]
;    oplot,(*pdata[i]).time,(*pdata[i]).mgas,linestyle = 2,thick = 2,color = colors[i]
;    oplot,(*pdata[i]).time,(*pdata[i]).mstar,linestyle = 2,thick = 2,color = colors[i]
;    oplot,(*pdata[i]).time,(*pdata[i]).mstar + (*data[i]).mgas,linestyle = 0,thick = 2,color = colors[i]
;    oplot,(*pdata[i]).time,(*pdata[i]).inflowcum - ((*pdata[i]).mstar + (*pdata[i]).mgas),linestyle = 1,thick = 2,color = colors[i]
ENDFOR
legend,keys,color = colors,linestyle = [0,0],/bottom
legend,['Gas Accreted','Accreted Gas Still in Halo','Gas Lost'],linestyle = [2,3,0]

;------------------ Outflow Ratio ---------------------
if (KEYWORD_SET(outplot)) then begin
    device,/close
    stop
    device,filename = outplot + '_fbratio.eps',/color,bits_per_pixel= 8,/times,ysize = 12,xsize = 18,xoffset =  2,yoffset =  2 
endif else begin
    stop
    window,1,xsize = 712,ysize = 392
endelse

yrange = [0,0.8]
yrange = [0,0.5]
FOR i = 0, n - 1 DO BEGIN
    loadct,ctables[i]
    IF i EQ 0 THEN BEGIN
        plot,(*pdata[i]).time,(*pdata[i]).outflow/(*pdata[i]).mgas,ytitle = 'Gas Mass Loss/Gas Mass',xstyle =9,xrange=[0,ageUniverse/1e9],xtitle = 'Time [Gyr]',yrange = yrange,/nodata,xmargin = [17,17],ymargin = [6,6],ystyle = 1
        axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'z_from_time_axes',xtitle = 'z',xtickv = ticktime_in_t,xticks = N_ELEMENTS(ticktime_in_t) - 1        
    ENDIF 
    oplot,(*pdata[i]).time,(*pdata[i]).outflow/(*pdata[i]).mgas,thick = thick[i],color = colors[i]
;   oplot,(*pdata[i]).time,(*pdata[i]).outflow/(*pdata[i]).mcoldg,thick = 2,color = colors[i]
ENDFOR
legend,keys,color = colors,linestyle = [0,0],/top,/right

;------------------ FB efficiency ---------------------
if (KEYWORD_SET(outplot)) then begin
    device,/close
    device,filename = outplot + '_fbeff.eps',/color,bits_per_pixel= 8,/times,ysize = 12,xsize = 18,xoffset =  2,yoffset =  2 
endif else begin
    stop
    window,2,xsize = 712,ysize = 392
endelse

yrange = [0,100]
yrange = [0,15]
FOR i = 0, n - 1 DO BEGIN
    loadct,ctables[i]
    dtime = ((*pdata[i]).time - [0,(*pdata[i])[0:N_ELEMENTS((*pdata[i]).time) - 1].time])*1e9
    IF i EQ 0 THEN BEGIN
        plot,(*pdata[i]).time,(*pdata[i]).outflow/dtime/(*pdata[i]).SFR,ytitle = 'Gas Mass Lost Per Year/SFR',xstyle =9,xrange=[0,ageUniverse/1e9],xtitle = 'Time [Gyr]',yrange = yrange,/nodata,xmargin = [17,17],ymargin = [6,6],ystyle = 1
        axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'z_from_time_axes',xtitle = 'z',xtickv = ticktime_in_t,xticks = N_ELEMENTS(ticktime_in_t) - 1        
    ENDIF
    oplot,(*pdata[i]).time,(*pdata[i]).outflow/dtime/(*pdata[i]).SFR,thick = thick[i],color = colors[i],psym = 10
ENDFOR
legend,keys,color = colors,linestyle = [0,0],/top,/right
if (KEYWORD_SET(outplot)) then device,/close ELSE stop

END


PRO fb_plots_master,outplot = outplot
IF KEYWORD_SET(outplot) THEN fgcolor = 0 ELSE fgcolor = 255
IF KEYWORD_SET(outplot) THEN thicks = [5,8] ELSE thicks = [1,3]

dir   = ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/',$
         '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/']
files = ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.00512.dir/h516.cosmo25cmb.3072g1MBWK.00512.halo.1',$
         '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.00512.dir/h516.cosmo25cmb.3072g14HBWK.00512.halo.1']
pfiles = ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/h516.cosmo25cmb.3072g1MBWK.param',$
         '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/h516.cosmo25cmb.3072g14HBWK.param']
keys = ['DnoH2','DH2']


dir   = ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/steps/',$
         '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/']
files = ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/steps/h603.cosmo50cmb.3072gs1MbwK.00512.dir/h603.cosmo50cmb.3072gs1MbwK.00512.halo.1',$
         '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.00512.dir/h603.cosmo50cmb.3072g14HBWK.00512.halo.1']
pfiles = ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/h603.cosmo50cmb.3072gs1MbwK.param',$
          '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/h603.cosmo50cmb.3072g14HBWK.param']
keys = ['warm','H2']

colors = [fgcolor,254]
ctables = [0,39]

colors = [80,254]
ctables = [39,39]
thick = [6,6]

fb_plots,dir,files,pfiles,color = colors,keys = keys,outplot = outplot,ctables = ctables,thick = thick
END
