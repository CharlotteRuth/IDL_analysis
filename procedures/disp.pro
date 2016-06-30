pro disp,DIR=dir,LEGEND=legend,XLABEL=xlabel,YLABEL=ylabel,_EXTRA=_extra

restore,'/net/grads-1/stinson/trace/tipsyprofile.sav'
 b=read_ascii(dir+'/mw_105k.prof',template=temp)
 d=read_ascii(dir+'/mw_21k.prof',template=temp)
 f=read_ascii(dir+'/mw_525k.prof',template=temp)
 if(keyword_set(xlabel)) THEN plot,d.r*1e5,d.sigr*2418,xtit='                               r [kpc]',linestyle=1,_EXTRA=_extra ELSE $
 if(keyword_set(ylabel)) THEN plot,d.r*1e5,d.sigr*2418,ytit=textoidl(' \sigma')+'!lr!n [km s!u-1!n]                                        ',linestyle=1,_EXTRA=_extra ELSE $
 plot,d.r*1e5,d.sigr*2418,linestyle=1,_EXTRA=_extra
 oplot,b.r*1e5,b.sigr*2418,line=2
 oplot,f.r*1e5,f.sigr*2418
 if(keyword_set(legend)) then legend,['11k g+s particles','55k g+s particles','275k g+s particles'],linestyle=[1,2,0],/bottom,/right,charsize=0.8

END
