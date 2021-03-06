;Make SNII.out file this way:
;grep 'feedback 0:'
;/astro/net/scratch1/christensen/DwarfResearch/ResSpecTests/1E5R/10M/200pc/out.116421
;| awk '{print$3,$4,$5}' >
;/astro/net/scratch1/christensen/DwarfResearch/ResSpecTests/1E5R/10M/200pc/SNII.out
;grep 'feedback 0:' /astro/net/scratch1/christensen/DwarfResearch/ResSpecTests/1E5R/10M/200pc/out.116437 | awk '{print$3,$4,$5}' >> /astro/net/scratch1/christensen/DwarfResearch/ResSpecTests/1E5R/10M/200pc/SNII.out

function sn,MASSUNIT=massunit,ERGUNIT=ergunit,DESN=desn

  if ( keyword_set(massunit) EQ 0) then massunit=2.362d5
  if ( keyword_set(ergunit) EQ 0) then ergunit=1.52251d15
  if ( keyword_set(desn) EQ 0) then desn=4.0d50

res = fltarr(3)
;print, 'What happens in each 1e4 M_sol of stars'
;print, "Massunit:",massunit,"  Ergunit:",ergunit,"  dESN:",desn
; stmass = read_ascii_array('onestar.00320.massform')  ;mass of star when it forms
; smass=total(stmass)		;total of newly formed mass in stars
   rdfloat, '/astro/net/scratch1/christensen/DwarfResearch/ResSpecTests/1E5R/10M/200pc/SNII.out',mass, energy, metals
;  print, 'Total energy output by all SNII over all timesteps
;  considered: ', total(energy)
  nsII = ergunit/desn*energy
  print, 'Number of SNII:',(ergunit/desn)*total(energy)   ;*1d4/(smass*massunit)
  print, 'Total mass ejected by SNII:  ', total(mass)  ;,smass*massunit,'  
  print, 'Mass of SNII Ejecta [M_sol]: ',total(mass)*massunit 
  plot,indgen(N_ELEMENTS(nsII),nsII
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

return, res 

END
