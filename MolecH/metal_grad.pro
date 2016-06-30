pro metal_grad
loadct = 37

set_plot,'ps'
device,filename = '/astro/users/christensen/code/MolecH/metal_grad.ps'


filename = "/astro/net/scratch1/abrooks/FABIO/MW1.1024g1bwK/MW1.1024g1bwK.00512/MW1.1024g1bwK.00512.1.std"
dKpcUnit        = 2.85714e4
dMsolUnit       = 3.171e15
height = 4.8e-5/2.0
radius = 1e-3/2.0
radius = 50.0/dKpcUnit/5.0 ;The 5 is scaling the viral radii -- kinda
rtipsy,filename,h,g,d,s

g_r = SQRT(g.x*g.x + g.y*g.y)
disk = g[where(g.z le height AND g.z ge -1.0*height AND g_r le radius)]
rad  = SQRT(disk.x*disk.x + disk.y*disk.y)*dKpcUnit
plot,rad,disk.zmetal,psym = 3,xtitle = 'Radius',ytitle = 'Z_metal'
oplot,rad,disk.zmetal,psym = 3,color = 50
nbins = 20.
bin_size = MAX(rad)/nbins
rad_array = (findgen(nbins + 1))*bin_size
plot_rad_array = (findgen(nbins) + 0.5)*bin_size
ave_zmetal = fltarr(nbins)
stderr_zmetal = fltarr(nbins)
m_stars = fltarr(nbins)

for i = 0,nbins-1 do begin
    ave_zmetal[i] = MEAN(disk[where((rad ge rad_array[i]) AND rad lt rad_array[i + 1])].zmetal) 
    stderr_zmetal[i] = stddev(disk[where(rad ge rad_array[i] AND rad lt rad_array[i + 1])].zmetal)
    m_stars[i] = ALOG10(TOTAL(disk[where(rad ge rad_array[i] AND rad lt rad_array[i + 1])].mass)/(rad_array[i + 1]*rad_array[i + 1] - rad_array[i]*rad_array[i])*dMsolUnit)
endfor

normalize = 0.025/MAX(m_stars)
oplot,plot_rad_array,ave_zmetal,psym = 1
oploterr,plot_rad_array,ave_zmetal,stderr_zmetal
oplot,plot_rad_array,m_stars*normalize,linestyle = 2,color = 50

;-----------------------------------------------------------------
filename = "/astro/net/scratch2/christensen/MolecH/11M/Disk_Iso_1e5/MW_disk"
dKpcUnit        = 1e5
dMsolUnit       = 1.36e17

rtipsy,filename + '.std',h,g,d,s
g.zmetal = 0.025*(fltarr(N_ELEMENTS(g.zmetal)) + 1)
wtipsy,filename + '_solar.std',h,g,d,s,/standard

rad  = SQRT(g.x*g.x + g.y*g.y)*dKpcUnit/5.0 ;The 5 is scaling the viral radii -- kinda
for i = 0L,nbins-1 do begin
    ind = where(rad ge rad_array[i] AND rad lt rad_array[i + 1])
    for i2 = 0L, N_ELEMENTS(ind) - 1 do begin
      ;  g[ind[i2]].zmetal = ave_zmetal[i] + stderr_zmetal[i]*RANDOMN(i*nbins+i2)
      ;  if g[ind[i2]].zmetal lt 0 then g[ind[i2]].zmetal = 0
    endfor
    m_stars[i] = ALOG10(TOTAL(g[where(rad ge rad_array[i] AND rad lt rad_array[i + 1])].mass)/(rad_array[i + 1]*rad_array[i + 1] - rad_array[i]*rad_array[i])*dMsolUnit)
endfor

;normalize = 0.025/MAX(m_stars)
;oplot,plot_rad_array,m_stars*normalize,linestyle = 3,color = 240
;oplot,rad,g.zmetal,psym = 3,color = 240

;wtipsy,filename + '_zgrad.std',h,g,d,s,/standard
;print,minmax(g.zmetal)
;print,mean(g.zmetal)
;device,/close
end
