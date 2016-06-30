PRO times_ejected,dirs,files,halo = halo,finalstep = finalstep,colors = colors,outplot = outplot,expelled = expelled,accrdisk = accrdisk,outflow = outflow
zmax = 3.2
n = n_elements(dirs)
IF NOT keyword_set(halo) THEN halo = strarr(n) + '1'
IF NOT keyword_set(finalstep) THEN finalstep = '00512'

formatplot,outplot = outplot,thick = formatthick

IF keyword_set(outplot) THEN BEGIN
    fgcolor = 0 
    bgcolor = 255
    xsize = 18
    ysize = 12
    thick = 14
    symsize = 2
ENDIF ELSE BEGIN
    fgcolor = 255
    bgcolor = 0
    xsize = 800
    ysize = 500
    thick = 4
    symsize = 3
ENDELSE
IF keyword_set(colors) THEN BEGIN
    loadct,39
    IF colors[0] eq 1 THEN  colors = (findgen(n) + 1)*240/n else colors = colors
    IF NOT keyword_set(ctables) THEN ctables = 39 + fltarr(n)
    IF NOT keyword_set(linestyles) THEN linestyles = fltarr(n) ;REVERSE(findgen(n)*2)
ENDIF ELSE BEGIN
    loadct,0    
    colors = (findgen(n) + 1)*fgcolor;  fltarr(n) + 5
    IF NOT keyword_set(ctables) THEN ctables = findgen(n)
    IF NOT keyword_set(linestyles) THEN linestyles = REVERSE(findgen(n)*2)  
 ENDELSE

colors = ['indian red','orange','Coral','Firebrick','Deep Pink',$
         'thistle','violet','medium orchid','purple',$
          'navy','blue','royal blue','sky blue','cyan',$
          'sea green','olive drab','lime green','gold','red']
;colors = ['cornflower blue','sea green','orange red','tomato','dark orchid',$
;          'slate grey','blue','gold','dark green','hot pink',$
;          'plum','brown','violet','dodger blue','goldenrod',$
;          'forest green','red','blue violet','gold','navy']
maxtimes = n_elements(colors) + 1
maxtimes = 15
;colors = (findgen(maxtimes) + 1)*254/maxtimes
nejected = fltarr(n,maxtimes)

smass = fltarr(n)
vmass = fltarr(n)

FOR i = 0, n - 1 DO BEGIN
   stat = read_stat_struc_amiga(dirs[i] + files[i] + '.' + finalstep + '/' + files[i] + '.' +  finalstep + '.amiga.stat')
   ind = (where(stat.group eq halo[i]))[0]
   vmass[i] = stat[ind].m_tot
   smass[i] = stat[ind].m_star
   IF keyword_set(accrdisk) THEN BEGIN
      IF keyword_set(outflow) THEN statei = mrdfits(dirs[i] + '/grp' + halo[i] + '.reaccrdiskheat_iord.fits',0) ELSE statei = mrdfits(dirs[i] + '/grp' + halo[i] + '.reaccrdisk_iord.fits',0) 
   ENDIF ELSE BEGIN
      If keyword_set(expelled) THEN statei = mrdfits(dirs[i] + '/grp' + halo[i] + '.reexpell_iord.fits',0) $ 
      ELSE BEGIN
         IF keyword_set(outflow) THEN statei = mrdfits(dirs[i] + '/grp' + halo[i] + '.reheat_iord.fits',0) ELSE  statei = mrdfits(dirs[i] + '/grp' + halo[i] + '.reeject_iord.fits',0)
      ENDELSE
   ENDELSE
   IF keyword_set(accrdisk) THEN BEGIN
       If keyword_set(expelled) THEN BEGIN
           statez = mrdfits(dirs[i] + '/grp' + halo[i] + '.reaccrdisk_z.fits',0)
           expelli = mrdfits(dirs[i] + '/grp' + halo[i] + '.reexpell_iord.fits',0)
           expellz = mrdfits(dirs[i] + '/grp' + halo[i] + '.reexpell_z.fits',0)               
           match2,expelli,statei,ind1,ind2
           statei = statei[where(ind2 NE -1)]
           statez = statez[where(ind2 NE -1)]
           statei1 = [0]
           FOR j = 0, n_elements(expelli) - 1 DO BEGIN
               ind = where(statei EQ expelli[j])
               IF (where(statez[ind] LT expellz[j]))[0] NE -1 THEN statei1 = [statei1,expelli[j]]
           ENDFOR
           statei = [expelli[uniq(expelli,sort(expelli))],statei1[1:n_elements(statei1) - 1]]
       ENDIF ELSE BEGIN
          IF keyword_set(outflow) THEN BEGIN
             ejected = mrdfits(dirs[i] + '/grp' + halo[i] + '.reheat_iord.fits',0)
             ejectedz = mrdfits(dirs[i] + '/grp' + halo[i] + '.reheat_z.fits',0)
             ejected = ejected[where(ejectedz LE zmax)]
             match2,ejected,statei,ind0,ind1
             statei = statei[where(ind1 NE -1)] ;Select for those particles that are ever ejected from the disk
          ENDIF ELSE BEGIN
             ejected = mrdfits(dirs[i] + '/grp' + halo[i] + '.reeject_iord.fits',0)
             ejectedz = mrdfits(dirs[i] + '/grp' + halo[i] + '.reeject_z.fits',0)
             ejected = ejected[where(ejectedz LE zmax)]
             match2,ejected,statei,ind0,ind1
             statei = statei[where(ind1 NE -1)] ;Select for those particles that are ever ejected from the disk
           ENDELSE
       ENDELSE
   ENDIF
   hist = histogram(statei)
   hist2 = histogram(hist,min = 1,max = 100)
;   stop
   nejected[i,0:maxtimes - 1] = float(hist2[0:maxtimes - 1])/float(total(hist2))
;   print,hist2
;   print,nejected[i,*]/total(nejected[i,*])
ENDFOR

nejected_scale = nejected
;FOR i = 0, n - 1 DO nejected_scale[i,*] = nejected[i,*]/total(nejected[i,*])
sorted  = sort(smass)
nejected_scale = nejected_scale[sorted,*]
smass = smass[sorted]
vmass = vmass[sorted]
baselines = fltarr(n); + 0.8; + 0.9

IF keyword_set(accrdisk) THEN BEGIN
    IF keyword_set(expelled) THEN BEGIN
        ytitle = textoidl('% Expelled Particles Accreted N Times') 
        outext = '_timesexpelled_accr.eps'
;        baselines = baselines + 0.95
    ENDIF ELSE BEGIN
       IF keyword_set(outflow) THEN BEGIN
        ytitle = textoidl('% Heated Particles Accreted N Times')
        outext = '_timesheated_accr.eps'
       ENDIF ELSE BEGIN
        ytitle = textoidl('% Ejected Particles Accreted N Times')
        outext = '_timesejected_accr.eps' 
     ENDELSE
    ENDELSE
ENDIF ELSE BEGIN
    IF keyword_set(expelled) THEN BEGIN
        ytitle = textoidl('% Particles Expelled N Times') 
        outext = '_timesexpelled.eps'
    ENDIF ELSE BEGIN
       IF keyword_set(outflow) THEN BEGIN
          ytitle = textoidl('% Particles Heated N Times')
          outext = '_timesheated.eps'
       ENDIF ELSE BEGIN
          ytitle = textpod('% Particles Ejected N Times') 
          outext = '_timesejected.eps' 
       ENDELSE
    ENDELSE    
ENDELSE

xrange = [3e9,1e12]
IF keyword_set(outplot) THEN  device,filename = outplot + outext,/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize

;Bar Plot
;cgBarPlot,(fltarr(n) + 1) - baselines,ytitle = ytitle,xtitle = 'Log Virial Mass [M' + sunsymbol() + ']',colors = colors[0],yrange = [baselines[0],1.01],baselines = baselines,barcoords = alog10(vmass),barwidth = 0.1,xrange = xrange;,barnames = cgnumber_formatter(alog10(vmass),decimal = 1);,yrange = ;'cornflower blue'
;axis,xaxis = 0,color = 'black',xrange = xrange
;plotv = fltarr(n)
;FOR j = 1, maxtimes  DO BEGIN
;   FOR i = 0,n - 1 DO plotv[i] = total(nejected_scale[i,0:maxtimes - j])
;   cgBarPlot,plotv - baselines,colors = colors[j],/overplot,baselines = baselines,barcoords = alog10(vmass),barwidth = 0.2
;ENDFOR

nejected_scale_total = nejected_scale
FOR j = 1, maxtimes - 1  DO nejected_scale_total[*,maxtimes - j - 1] = nejected_scale[*,maxtimes - j - 1] + nejected_scale_total[*,maxtimes - j] 

plot,vmass,1 - nejected_scale_total[*,maxtimes - 1],ytitle = ytitle,xtitle = 'Log Virial Mass [M' + sunsymbol() + ']',yrange = [0,105],xrange = xrange,psym = symcat(18),/nodata,/xlog,ycharsize = 1
loadct,0
FOR i = 0,n - 1 DO BEGIN
   oplot,[vmass[i],vmass[i]],[0,100],color = 100,thick = thick
   FOR j = 1, maxtimes - 1 DO oplot,[vmass[i],vmass[i]],[0,100*(1 - nejected_scale_total[i,maxtimes - j])],color = cgcolor(colors[j - 1 - (maxtimes - n_elements(colors) - 1)]),thick = thick
ENDFOR
FOR j = 1, maxtimes - 1  DO BEGIN
   oplot,vmass,100*(1 - nejected_scale_total[*,maxtimes - j]),psym = symcat(18),symsize = symsize,color = cgcolor(colors[j - 1 - (maxtimes - n_elements(colors) - 1)])
   oplot,vmass,100*(1 - nejected_scale_total[*,maxtimes - j]),psym = symcat(sym_outline(18)),symsize = symsize
ENDFOR

IF keyword_set(outplot) THEN device, /close ELSE stop
END
