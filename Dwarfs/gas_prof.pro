pro gas_prof,infile,outfile,densmin,densmax,SECONDFILE = secondfile,TITLE = TITLE
;This plots a series of stellar profiles showing the evolution in time
;gas_prof,'../1E5R/10M/o10M_1.00300','10M_gprof',1e2,1e8,secondfile=
;'../1E5R/10M_original/o10M_1.00300',title = '1e10 Solar Mass' 

msol = 2.362e5
SYS_AMUCC = 0.0096302124  ;  conversion between system units and amu/cc
timeunit = 1e9
nframes = 20
maxtime = 3e9
nbins = 25
rmax = 10
rmin = 0
loadct,39

rtipsy,infile,h,g,d,s
IF (KEYWORD_SET(secondfile)) THEN rtipsy,secondfile,h2,g2,d2,s2
IF NOT (KEYWORD_SET(TITLE)) THEN title = infile

time = (findgen(nframes)+1.)*maxtime/nframes
files = sindgen(nframes)

ind = WHERE(ALOG10(g.dens) gt 3.0)
;g_dens = g[ind]
g_prof = prof(g, 'gas', h.time, nbins = nbins, rmax = rmax)
;g_prof_dens = prof(g_dens, 'gas', h.time, nbins = nbins, rmax = rmax)
;    fitpar =  [0,0,0,0,0]    
;g_prof.rho = s_prof.rho * msol  ; Msol/kpc^2

IF (KEYWORD_SET(secondfile)) THEN BEGIN
    ind = WHERE(ALOG10(g2.dens) gt 3.0)   
    g_dens2 = g2[ind]
    g_prof2 = prof(g2, 'gas', h2.time, nbins = nbins, rmax = rmax)
    g_prof_dens2 = prof(g_dens2, 'gas', h2.time, nbins = nbins, rmax = rmax)
;    fitpar =  [0,0,0,0,0]    
;    s2_prof.rho = s2_prof.rho * msol ; Msol/kpc^2
ENDIF

set_plot,'x'
plot,g_prof.rbins,g_prof.rho,/ylog,title = title,xtitle = 'Radius',ytitle = 'GAs Density',color = 240,linestyle = 0
;oplot,g_prof_dens.rbins,g_pro_dens.rho,color = 240, linestyle = 1
IF (KEYWORD_SET(secondfile)) THEN BEGIN
    oplot,g_prof2.rbins,g_prof2.rho,color = 50,linestyle = 0
    oplot,g_prof_dens2.rbins,g_prof_dens2.rho,color = 50,linestyle = 1
    legend,['Higher SF Threshold - All','Original -- All','Original -- Dense'],linestyle=[0,0,1],color = [240,50,50]
ENDIF   

set_plot,'ps'
device,filename='gasprof_10M.eps',/color,bits_per_pixel=8 
plot,g_prof.rbins,g_prof.rho,/ylog,title = title,xtitle = 'Radius',ytitle = 'GAs Density',color = 240,linestyle = 0
;oplot,g_prof_dens.rbins,g_pro_dens.rho,color = 240, linestyle = 1
IF (KEYWORD_SET(secondfile)) THEN BEGIN
    oplot,g_prof2.rbins,g_prof2.rho,color = 50,linestyle = 0
    oplot,g_prof_dens2.rbins,g_prof_dens2.rho,color = 50,linestyle = 1
    legend,['Higher SF Threshold - All','Original -- All','Original -- Dense'],linestyle=[0,0,1],color = [240,50,50]
ENDIF
device,/close
    set_plot,'x'

END
