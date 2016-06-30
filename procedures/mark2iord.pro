function mark2iord, markfile, iordfile, h
  rdfloat, markfile, marks, skipline=1,/long
;  iords = rbarray(iordfile)
  readarr,iordfile,h,iords,type = 'long',/ascii
  iordarr = iords[marks]
  return, iordarr
end
