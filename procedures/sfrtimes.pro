FUNCTION sfrtimes,tipsbase

  restore,'/net/mega-2/stinson/MW/sfrtimetemp.sav'
  b=read_ascii(tipsbase+'.sfrtimes',template=temp)
  return,b
END
