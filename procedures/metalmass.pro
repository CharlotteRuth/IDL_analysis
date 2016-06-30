pro metalmass,base,orig=orig

  m=read_ascii_array(base+'.massform')
  rtipsy,base,h,g,d,s
  o=read_ascii_array(base+'.OxMassFrac')
  fe=read_ascii_array(base+'.FeMassFrac')
  if(keyword_set(d)) then pm=[g.mass,d.mass,s.mass] else pm=[g.mass,s.mass]
  
  fem = total(fe*pm)
  om = total(o*pm)
  print,"Oxygen present:",om
  print,"Iron present:",fem
  if(keyword_set(orig)) then begin
    rtipsy,orig,oh,og,od,os,/justhead
    gd=oh.ngas-h.ngas
    mform=total(m[oh.n-gd:h.n-1])
  endif else mform = total(m)
  print,"Stellar mass formed:",mform
  print,"Fe produced by 1e4 M_sol stars:",1e4*fem/mform
  print,"Ox produced by 1e4 M_sol stars:",1e4*om/mform

END
