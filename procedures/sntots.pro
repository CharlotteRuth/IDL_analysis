pro sntots,MASSUNIT=massunit,ERGUNIT=ergunit,DESN=desn

  if ( keyword_set(massunit) EQ 0) then massunit=13.5985d16
  if ( keyword_set(ergunit) EQ 0) then ergunit=1.592d67
  if ( keyword_set(desn) EQ 0) then desn=1d48

print, 'What happens in each 1e4 M_sol of stars'
print, "Massunit:",massunit,"  Ergunit:",ergunit,"  dESN:",desn
 stmass = read_ascii_array('mw_105k.00200.massform')
  smass=total(stmass)
  rdfloat,'SNII.out',mass, energy, metals
  print, total(energy),'  Number of SNII:',ergunit/desn*total(energy)*1d4/(smass*massunit)
  print, total(mass),smass*massunit,'  Mass of SNII Ejecta [M_sol]:',total(mass)/smass*1d4

  rdfloat,'SNIa.out',mass, energy, metals
  print, 'Number of SNIa:',ergunit/desn*total(energy)*1d4/(smass*massunit)
  print, 'Mass of SNIa Ejecta [M_sol]:',total(mass)/smass*1d4
END
