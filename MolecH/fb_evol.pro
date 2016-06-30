

pro fb_evol,dir,files,pfiles,outplot = outplot,linestyle = linestyle, color = color, psym = psym, symsize = symsize, thick = thick,scale = scale,keys = keys,vline = vline
formatplot,outplot = outplot
;------------------ FB ---------------------
if (KEYWORD_SET(outplot)) then begin
    device,/close
    device,filename = outplot + '_fb.eps',/color,bits_per_pixel= 8,/times,ysize = 16,xsize = 18,xoffset =  2,yoffset =  2 
endif else begin
;    stop
    window,9,xsize = 712,ysize = 600
endelse
binsize = 5e8
nbins = ROUND(14.0/(binsize/1e9))

FOR i = 0, n - 1 DO BEGIN
    cd,dir[i]
    readcol,filebase,halo,time,z,mtot,mgas,mstar,mdark,metal,mcoldg,mH2,H2frac,Pressure,r25,SFR
    rtipsy,file,h,g,d,s,/justhead
    readarr,file_coolon,h,coolon,part = 'gas',/ascii
    readarr,file_iord,h,iord,part = 'gas',/ascii
;------------------ Outflow -------------------------
    outflowz_all = mrdfits(dir[i] + '../grp1.rvir.outflow_z.fits')
    outflowm_all = mrdfits(dir[i] + '../grp1.rvir.mass_at_outflow.fits')
    outflowiord_all = mrdfits(dir[i] + '../grp1.rvir.outflow_iord.fits')
    
    outflowz = outflowz_all[where(outflowz_all ne 99.0)]
    outflowm = outflowm_all[where(outflowz_all ne 99.0)]
    outflowiord = outflowiord_all[where(outflowz_all ne 99.0)]
    match,iord,outflowiord,ind_iord,ind_iordout
    coolonout = where(coolon[ind_iord] ne 0)

    outflow_z = outflowz[ind_iordout[coolonout]]
    outflow_iord = outflowiord[ind_iordout[coolonout]]
    outflow_mass = outflowm[ind_iordout[coolonout]]

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

  
;-------------------- Bin Outflow and Inflows according to the same
;                     steps in the data file
    a_outflowm = fltarr(N_ELEMENTS(mtot))
    a_outflowm_cum = fltarr(N_ELEMENTS(mtot))
    a_inflowm = fltarr(N_ELEMENTS(mtot))
    a_inflowm_cum = fltarr(N_ELEMENTS(mtot))
    a_time = [0,time]

;IDL> for i=1,10 do begin $
;IDL> & print,i $
;IDL> & endfor

    FOR j = 0, N_ELEMENTS(mtot) - 1 DO BEGIN 
        t0 = a_time[j] 
        t1 = a_time[j + 1]  
        IF (where(outflowt00 GT t0 and outflowt00 LE t1))[0] NE -1 THEN a_outflowm[j] = total(outflowm00[where(outflowt00 GT t0 and outflowt00 LE t1)]) 
        IF (where(outflowt00 LE t1))[0] NE -1 THEN a_outflowm_cum[j] = total(outflowm00[where(outflowt00 LE t1)]) 
        IF (where(inflowt GT t0 and inflowt LE t1))[0] NE -1 THEN a_inflowm[j] = total(inflowm[where(inflowt GT t0 and inflowt LE t1)]) 
        IF (where(inflowt LE t1))[0] NE -1 THEN a_inflowm_cum[j] = total(inflowm[where(inflowt LE t1)]) 
    ENDFOR
    IF i eq 0 THEN a_array0 = [[time],[a_inflowm],[a_inflowm_cum],[a_outflowm],[a_outflowm_cum]] ELSE a_array1 = [[time],[a_inflowm],[a_inflowm_cum],[a_outflowm],[a_outflowm_cum]]

    IF i EQ 0 THEN BEGIN
        plot,outflowtbin,cuminflow*scale[i],ytitle = 'Gas  Mass [M' + sunsymbol() + ']',xstyle =9,xrange=[0,ageUniverse/1e9],xtitle = 'Time [Gyr]',yrange = [-2.5e10,5e10],/nodata,ystyle = 9,xmargin = [17,17],ymargin = [6,6]
        axis,ystyle = 1,yaxis = 1,ytitle = textoidl('8 X M_{star}') + '[M' + sunsymbol() + ']'
        axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'z_from_time_axes',xtitle = 'z',xtickv = ticktime_in_t,xticks = N_ELEMENTS(ticktime_in_t) - 1        
        oplot,[0,14],[0,0],thick = 2
    ENDIF
    oplot,outflowtbin,cuminflow*scale[i],linestyle = 2,thick = 2,color = colors[i] ;gas mass in
    oplot,outflowtbin,-1.0*cumoutflow*scale[i],linestyle = 3,thick = 2,color = colors[i];gas mass out
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

pro mass_evol_master

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
mass_evol,dir,files,pfiles,scale = scale,linestyle = [0,0],psym = [14,14],/color,keys = ['Spiral'+textoidl(' Galaxy'),'Dwarf Galaxy'],outplot = '~/plots/twins/h603_h516_evol_flash'
end
