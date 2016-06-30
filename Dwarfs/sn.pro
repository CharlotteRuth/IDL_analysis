;Make SNII.out file this way:
;grep 'feedback 0:'
;/astro/net/scratch1/christensen/DwarfResearch/ResSpecTests/1E5R/10M/200pc/out.116421
;| awk '{print$3,$4,$5}' >
;/astro/net/scratch1/christensen/DwarfResearch/ResSpecTests/1E5R/10M/200pc/SNII.out
;grep 'feedback 0:' /astro/net/scratch1/christensen/DwarfResearch/ResSpecTests/1E5R/10M/200pc/out.116437 | awk '{print$3,$4,$5}' >> /astro/net/scratch1/christensen/DwarfResearch/ResSpecTests/1E5R/10M/200pc/SNII.out
;snova,MASSUNIT=2.362d38,ERGUNIT=7.15279d53,DESN=4d50
pro snova,MASSUNIT=massunit,ERGUNIT=ergunit,DESN=desn

;  if ( keyword_set(massunit) EQ 0) then massunit=4.6980189d28
  if ( keyword_set(massunit) EQ 0) then massunit=2.362d38
  if ( keyword_set(ergunit) EQ 0) then ergunit=7.15279d53
  if ( keyword_set(desn) EQ 0) then desn=4.0d50

res = fltarr(3)
;print, 'What happens in each 1e4 M_sol of stars'
;print, "Massunit:",massunit,"  Ergunit:",ergunit,"  dESN:",desn
; stmass = read_ascii_array('onestar.00320.massform')  ;mass of star when it forms
; smass=total(stmass)		;total of newly formed mass in stars
   readcol, '/astro/net/scratch1/christensen/DwarfResearch/ResSpecTests/1E5R/10M/200pc/SNII.out',format = 'F3,F3,F3',mass, energy, metals
;  print, 'Total energy output by all SNII over all timesteps
;  considered: ', total(energy)
  nsII = ergunit/desn*energy*1d4
  print, 'Number of SNII:',(ergunit/desn)*total(energy)*1d4   ;*1d4/(smass*massunit)
  print, 'Total mass ejected by SNII:  ', total(mass)  ;,smass*massunit,'  
  print, 'Mass of SNII Ejecta [M_sol]: ',total(mass)*massunit 
  plot,3.0*lindgen(N_ELEMENTS(nsII))/N_ELEMENTS(nsII),nsII,xtitle = 'Time', ytitle = "Number of Supernovas",xrange=[0,3]
  set_plot,'ps'
  device,filename="/astro/net/scratch1/christensen/DwarfResearch/ResSpecTests/1E5R/10M/200pc/supernovas.eps"
  plot,3.0*lindgen(N_ELEMENTS(nsII))/N_ELEMENTS(nsII),nsII,xtitle = 'Time', ytitle = "Number of Supernovas",xrange=[0,3]
  print, 'Minimum Number of Supernovas per time step: ',MIN(nsII[Where(nsII ne 0)])
  device,/close
  set_plot,'x'
  stop
;   mass = mass*massunit
;   res[0] = total(mass)

;   rdfloat, 'SNIa.out',mass, energy, metals
;  print, 'Number of SNIa: ', (ergunit/desn)*total(energy)  ;*1d4/(smass*massunit)
;  print, 'Mass of SNIa Ejecta [M_sol]:', total(mass)*massunit 
;res[1] = total(mass)*massunit

;   rdfloat,'winds.out',mass, energy, metals
;  print, 'Mass of Ejecta in Winds [M_sol]:', total(mass)*massunit 
;res[2] = total(mass)*massunit

END
