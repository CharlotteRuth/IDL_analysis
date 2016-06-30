pro length_distro,out
loadct,39
stat_file = '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g1MBWK/steps/h516.cosmo25cmb.1536g1MBWK.00512.dir/h516.cosmo25cmb.1536g1MBWK.00512.amiga.stat'

input_halos = READ_STAT_STRUC_AMIGA(stat_file)
lengths = alog10(input_halos.rvir)
mass = input_halos.m_tot/0.73

wstars = where(input_halos.m_star ne 0)
lengths_wstars = alog10(input_halos[wstars].rvir)
mass_wstars = input_halos[wstars].m_tot/0.73


db = (MAX(lengths) - MIN(lengths))/50.
lengths_y = histogram(lengths,binsize = db, min = MIN(lengths))
lengths_x = findgen(N_elements(lengths_y))*db + MIN(lengths)

lengths_y_wstars = histogram(lengths[wstars],binsize = db, min = MIN(lengths))

min_mass = MIN(alog(mass))
db = (MAX(alog(mass)) - min_mass)/50.
mass_y = histogram(alog(mass),binsize = db, min = min_mass)
mass_x = exp(findgen(N_elements(mass_y))*db + min_mass)

mass_y_wstars = histogram(alog(mass[wstars]),binsize = db, min = min_mass)

;LocalGroup volume = 677912.0*0.000273762
;h516 volume
volume = 8.090540e-04*15625


set_plot,'x'
set_plot,'ps'
device,filename = 'PowerSpectrum_length'+out+'.eps',/color,bits_per_pixel=8,/times
plot,lengths_x,lengths_y/volume,xtitle = 'R_vir',ytitle = 'dn/log10(M) [Mpc^3]',psym = 10,/xlog,/ylog,yrange = [1e-2,10],xrange = [0.3,2],charsize = 2,xstyle = 1,ystyle = 1
oplot,lengths_x,lengths_y_wstars/volume,psym = 10,color = 240
device,/close

set_plot,'ps'
device,filename = 'PowerSpectrum_mass'+out+'.eps',/color,bits_per_pixel=8,/times
plot,mass_x,mass_y/volume,xtitle = 'Log(M_vir) [Solar Mass]',ytitle = 'dn(m)/ln(M) [Mpc^-3]',psym = 10,/ylog,/xlog,yrange = [1e-2,10],charsize = 2,ystyle = 1,xstyle = 1,xrange = [1e7,1e11]
oplot,mass_x,mass_y_wstars/volume,psym = 10,color = 240
device,/close

set_plot,'ps'
device,filename = 'PowerSpectrum_mass2'+out+'.eps',/color,bits_per_pixel=8,/times
plot,mass_x,mass_y/volume,xtitle = 'Log(M_vir) [Solar Mass]',ytitle = 'dn(m)/ln(M) [Mpc^-3]',psym = 10,/ylog,/xlog,yrange = [1e-10,1],charsize = 2,ystyle = 1,xstyle = 1,xrange = [1e10,1e16]
oplot,mass_x,mass_y_wstars/volume,psym = 10,color = 240
device,/close
end
