
;Make TWO graphs as .ps files (see tutorial 3)
;First graph: luminosity versus age
;Second graph: luminosity versus metallicity

;Label, try log axes, try different symbols, try different ranges (see
;IDL help on plot command)

;Write a sentence description of what is happening in each

pro starformation

;laroy chase 2/5/10
;read in sfrhist.fits and plot stellar data


filename="/astro/users/lacowboy/sunrise/MW1.256g1bwdK.00512.1/sfrhist.fits"
data=mrdfits(filename,6,header)
help,data,/struct
print,header
print,data[0].age
print,data[0:100].age

plot,data.age,data.mass,psym = 3,yrange = [1.7e6,2.4e6],xtitle = 'Age [yr]',ytitle='Mass [Msun]'

plot,data.age,data.l_bol,yrange = [3.4e32,6e35],xrange = [1e6,1e11],ytitle = 'l_bol [W]',xtitle = 'Age [yr]',psym = 3,/ylog,/xlog



set_plot,'ps'
device,filename = 'starformation.ps'
plot,data.age,data.l_bol,yrange = [3.4e32,6e35],xrange = [1e6,1e11],ytitle = 'l_bol [W]',xtitle = 'Age [yr]',psym = 3,/ylog,/xlog
device,/close
set_plot,'x'


set_plot,'ps'
device,filename = 'starformation.eps'
plot,data.metallicity,data.l_bol,xtitle = 'metallicity',ytitle = 'l_bol[W]',psym = 3,/ylog
device,/close
set_plot,'x'

stop
loadct,39  ;allow you to plot in color
ind = where(data.age gt 1e6 AND data.age lt 1e8)
oplot,data[ind].metallicity,data[ind].l_bol,psym = 3,color = 240
ind = where(data.age gt 1e10 AND data.age lt 1e12)
oplot,data[ind].metallicity,data[ind].l_bol,psym = 3,color = 60

stop
data=mrdfits(filename,6,header)

end

