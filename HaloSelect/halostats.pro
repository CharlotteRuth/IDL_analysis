pro halostats
formatplot,outplot = 1

;haloindfiles = 'BOTTOM_QUARTILE_STELLAR_MASS_SELECTION'
haloindfiles = 'TOP_QUARTILE_STELLAR_MASS_SELECTION'
filebase = 'cosmo50cmb.256g2MbwK.00512'
kpcunit = 50000.0
msununit = 18.48e15
timeunit=SQRT((kpcunit*3.086d21)^3/(6.67d-8*msununit*1.99d33))/(3600.*24.*365.24)

readcol,haloindfiles,haloind
amigastat = read_stat_struc_amiga(filebase+'.amiga.stat')
rtipsy,filebase,h,g,d,s
z = 1/h.time - 1
readarr,filebase+'.amiga.grp',h,grp,part = 'star',/ascii
close,/all
;openw,1,filebase+'_halodat_blue.dat'
openw,1,filebase+'_halodat_red.dat'
printf,1,'Group','Mass','Stellar Mass','Star Mass/Mass','Gas Mass/Bar. Mass','V_peak','xc','yc','zc',format='(9A19)'
openw,2,'halosplot.macro'
printf,2,'halosplot'
printf,2,''
printf,2,'openb '+filebase
printf,2,'loads 1'
printf,2,'loadall'
printf,2,'redshift 0 c 73 0.24 0.76 p 1'
printf,2,'readarray '+filebase+'.amiga.grp'
FOR i =0, N_ELEMENTS(haloind) - 1 DO BEGIN
    j = where(haloind[i] eq amigastat.group)
    center = [amigastat[j].xc, amigastat[j].yc, amigastat[j].zc]/kpcunit*1000 - 0.5
    radius = amigastat[j].rvir/kpcunit*1000
    printf,1,amigastat[j].group,amigastat[j].m_tot,amigastat[j].m_star,amigastat[j].m_star/amigastat[j].m_tot,amigastat[j].m_gas/(amigastat[j].m_star +amigastat[j].m_gas),amigastat[j].vc,amigastat[j].xc, amigastat[j].yc, amigastat[j].zc,format='(I19,2E19,6F19)'
    device,filename = 'SFH_' + strtrim(haloind[i],2) + '.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 12,xoffset =  2,yoffset =  2
    sfr,s[where(grp eq haloind[i])],massunit = msununit,timeunit = timeunit,xrange = [0,14],sarray = sarray, tarray  = tarray, mintime = 0, maxtime = 14.0, binsize = 1e8,title = 'Halo ID: ' + strtrim(haloind[i],2)
    device,/close
    printf,2,'markarray 0 ' + strtrim(haloind[i] - 0.5,2) + ' ' + strtrim(haloind[i] + 0.5,2)
    printf,2,'setsphere 1 ' + strtrim(center[0],2) + ' ' + strtrim(center[1],2) + ' ' + strtrim(center[2],2) + ' 0.01'
    printf,2,'abox 1'
    printf,2,'viewarray all default 140 800 clip'
    printf,2,'hard dump ' + filebase + '.' + strtrim(haloind[i],2) + '.xwd'
ENDFOR
printf,2,'closeb'
printf,2,'end'
printf,2,''
close,1
close,2
END
