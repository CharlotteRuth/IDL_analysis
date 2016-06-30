PRO VIEWGAS_DIST
!P.CHARSIZE = 1.25
!P.thick = 1.5
!X.Charsize = 1
!Y.Charsize = 1
!X.style = 1

msol_unit = 2.362d5
msol = 1.989d33 ;solar mass in g
kpc_unit = 1
kpc = 3.0857d21 ;kpc in cm
mH = 1.673534d-24
molec_mass = 1.72000 ;average mass of an atom in the universe (0.76 + 0.24*4)
grav = 6.67259d-8 ;in cgs
SYS_AMUCC = 0.0096302124  ;  conversion between system units and amu/cc
loadct,0
dens_conv = msol*msol_unit/(kpc_unit*kpc)^3/mH*molec_mass ;convers system density into atoms per cm^3
set_plot,'x'
set_plot,'ps'
device,filename = '/astro/net/scratch1/christensen/DwarfResearch/results/gas_distro.eps',/color,bits_per_pixel=8,/times
!P.FONT = 0

filenames = ['/1E5R/10M_k/o10M_1.00300','/1E6R/10M_k/o10M_1.00300','/1E6R/12M_k/o12M_1.00300']
;filenames = ['/1E5R/10M_k/o10M_1.00300','/1E5R/10M_k/o10M_1.00300','/1E5R/12M_k/o12M_1.00300']
!P.MULTI = [0, 3, 1, 0, 0]
ytitles = [textoidl('R_{vir}'),' ',' ']
yticknames = [['','','','',''],[' ',' ',' ',' ',' '],[' ',' ',' ',' ',' ']]
xticknames = [['-4','-2','0','2','4'],[' ','-2','0','2','4'],[' ','-2','0','2','4']]
rvirs = [44L,44L,206L]
scale = rvirs*rvirs/rvirs[0]/rvirs[0]
label = [textoidl('10^{10} M')+sunsymbol()+', '+textoidl('10^{5}')+' DM particles',textoidl('10^{10} M')+sunsymbol()+', '+textoidl('10^{6}')+' DM particles',textoidl('10^{12} M')+sunsymbol()+', '+textoidl('10^{6}')+' DM particles']
!X.Style = 1
!Y.Style = 1
FOR ct = 0,2 DO BEGIN
    base = '/astro/net/scratch1/christensen/DwarfResearch/ResMassTests'
    filename = filenames[ct]
;    window,ct

    rtipsy, base+filename, h, g, d, s
    nlevels = 21
    maxdist = 5.4
    xbinsize = 0.2
    ybinsize = 0.2
    xmax = 4.0;200
    xmin = -4.0;-200
    ymax = 8.0;400
    ymin = -8.0;-400
    threshold = alog10(mean(g.mass)*5/xbinsize/ybinsize) ;maxdist/nlevels
    threshold = 0.25
    maxdist = 5.25
    levels = (findgen(nlevels)*(maxdist-threshold)/(nlevels - 1.0) + threshold);*scale[ct]
;    threshold = MIN(levels)
    print,levels
    print,TOTAL(g.mass)
;    stop
    multiplot
    contour_plus,g.x/rvirs[ct],g.z/rvirs[ct],weight = g.mass/scale[ct],xmin = xmin, xmax = xmax, ymin = ymin, ymax  = ymax, xbinsize = xbinsize, ybinsize = ybinsize, threshold = threshold,levels = levels,/loglevel, ytitle = ytitles[ct],ytickname = yticknames[*,ct], xtickname = xticknames[*,ct]
;    stop
ENDFOR
xyouts,[2900,7700,12400],[11000,11000,11000],label,charsize = 1.0 ,/device
multiplot,/reset
;stop
device,/close
END
