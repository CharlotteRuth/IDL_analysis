;filename = '/astro/net/scratch2/christensen/MolecH/11M/Disk_Iso_1e5/largeStar/MW_disk.00003'

PRO selectBar, filename
limit = 0.00004
height = 5e-6/2

rtipsy,filename,h,g,d,s
g = g[where(ABS(g.z) lt height)]
ind_in = where(SQRT(g.x*g.x + g.y*g.y) lt limit*SQRT(2.0))
gin = g[ind_in]
gin.dens = alog10(gin.dens)
minc = MIN(gin.dens)
maxc = MAX(gin.dens)
color = FIX((gin.dens - minc)/(maxc - minc)*200 + 40)

regionBAR1 = selectregion_tipsy(gin.x,gin.y,color = color,XRANGEVEC=[-1.0*limit,limit],YRANGEVEC =[-1.0*limit,limit])
regionBAR = regionBAR1
WHILE (1) DO BEGIN
    print,'Would you like to select another bar? [y or n]'
    i = GET_KBRD()
    IF i eq 'n' THEN BREAK
    regionBAR1 = selectregion_tipsy(gin.x,gin.y,color = color,XRANGEVEC=[-1.0*limit,limit],YRANGEVEC =[-1.0*limit,limit],/noplot)
    regionBAR = [regionBAR,regionBAR1]
END

regionINNER1 = selectregion_tipsy(gin.x,gin.y,color = color,XRANGEVEC=[-1.0*limit,limit],YRANGEVEC =[-1.0*limit,limit])
regionINNER = regionINNER1
WHILE (1) DO BEGIN
    print,'Would you like to select another inner region? [y or n]'
    i = GET_KBRD()
    IF i eq 'n' THEN BREAK
    regionINNER1 = selectregion_tipsy(gin.x,gin.y,color = color,XRANGEVEC=[-1.0*limit,limit],YRANGEVEC =[-1.0*limit,limit],/noplot)
   regionINNER = [regionINNER,regionINNER1]
END 

plot,gin.x,gin.y,psym = 3,xrange = [-1.0*limit,limit],yrange = [-1.0*limit,limit]
FOR i = 0LL, N_ELEMENTS(gin.x) - 1 DO oplot,[gin[i].x,gin[i].x],[gin[i].y,gin[i].y],psym = 3,color = color[i]
oplot,gin[regionINNER].x,gin[regionINNER].y,color = 240,psym = 3



stop
END
