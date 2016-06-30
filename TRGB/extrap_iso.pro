function extrap_iso,filename,elements
;9/10/07 Charlotte Christensen
;This program take one of Leo's isochrones and return a file with the
;magnitudes and colors of the RGB, complete with extapolation beyond the tip.

;called from iso.pro, finish_tips.pro, iso_segment.pro and plottip.pro
READCOL,filename,z,logage,m_ini,m_act,f475mag,f606mag,f814mag,FORMAT='F,F,F,F,X,X,X,X,X,F,X,X,F,X,X,X,X,F',/SILENT
padova=REPLICATE({z:0.0, logage:0.0, m_ini:0.0, m_act:0.0,umag:0.0,bmag:0.0,vmag:0.0,rmag:0.0,imag:0.0,jmag:0.0,hmag:0.0,kmag:0.0,f475mag:0.0,f606mag:0.0,f814mag:0.0,j2mass:0.0,h2mass:0.0,k2mass:0.0,w300mag:0.0,w439mag:0.0,w555mag:0.0,w606mag:0.0,w814mag:0.0},N_ELEMENTS(z))
length = N_ELEMENTS(z)
padova.z=z
padova.logage=logage 
padova.m_ini=m_ini   
padova.m_act=m_act   
padova.f475mag=f475mag  
padova.f606mag=f606mag  
padova.f814mag=f814mag  

additional_elements = elements - length
padova_exten = REPLICATE({z:0.0, logage:0.0, m_ini:0.0, m_act:0.0,umag:0.0,bmag:0.0,vmag:0.0,rmag:0.0,imag:0.0,jmag:0.0,hmag:0.0,kmag:0.0,f475mag:0.0,f606mag:0.0,f814mag:0.0,j2mass:0.0,h2mass:0.0,k2mass:0.0,w300mag:0.0,w439mag:0.0,w555mag:0.0,w606mag:0.0,w814mag:0.0},elements)
padova_exten.z = z[0]
padova_exten.logage=logage[0]

f475last = [padova[length - 3:length - 1].f475mag] ;b
f606last = [padova[length - 3:length - 1].f606mag] ;v
f814last = [padova[length - 3:length - 1].f814mag] ;i

color1 = padova.f606mag - padova.f814mag
color1last = [color1[length - 3:length - 1]] ;v-i
fit = poly_fit(color1last,f814last,1)
color1_added =  4.0*findgen(additional_elements)/FLOAT(additional_elements) + color1[length-1]
f814_added = fit[1]*color1_added  + fit[0]
padova_exten.f814mag = [padova.f814mag,f814_added];extending f814 magnitude

color2 = padova.f475mag - padova.f606mag
color2last = [color2[length - 3:length - 1]] ;b-v
fit = poly_fit(color2last,f606last,1)
color2_added =  4.0*findgen(additional_elements)/FLOAT(additional_elements) + color2[length-1]
f606_added = fit[1]*color2_added + fit[0]
padova_exten.f606mag = [padova.f606mag,f606_added];extending f606 magnitude
padova_exten.f475mag = [padova.f475mag,color2_added + f606_added]
padova_exten.f606mag = [padova.f606mag,color1_added + f814_added]
return, padova_exten
END
