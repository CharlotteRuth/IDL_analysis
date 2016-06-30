pro west,file,hdu,alpha,Q,m,scale=scale,bands=bands, imscale = imscale, image_path=image_path,one_each=one_each,noerase=noerase

sz=size(hdu)  
if(sz[0] gt 0) then begin
    start=hdu[0]
    stop=hdu[1]
end else begin
    start=hdu
    stop=hdu
end

if not keyword_set (image_path) then begin
    if not keyword_set(noerase) then $
      erase
    xyouts,4,8,file
end

for v=start,stop do begin
image = mrdfits (file, v)

Sz = size (image)
xs = sz [1]
ys = sz [2]
if not keyword_set(scale) then scale=[.0065,.0055,.0065]*0.01
nonlinearity = 2.0
origin=[0.,0.,0.]
if not keyword_set (bands) then bands = [2,3,4,5,6]
u_band = image (*,*, bands [0]) ; default is SDSS i band
g_band = image (*,*, bands [1]) ; default is SDSS r band
r_band = image (*,*, bands [2]) ; default is SDSS g band
i_band = image (*,*, bands [3]) ; default is SDSS z band
z_band = image (*,*, bands [4]) ; default is SDSS u band


;r = ((i_band-1000)+(z_band-1000))
;g = ((r_band-1000))
;b = ((g_band-1000)+(u_band-1000))

r = ((z_band))
g = ((r_band))
b = ((u_band))

nx=Sz[1]
ny=Sz[2]
img= dblarr(nx,ny,3)
img[*,*,0]= r
img[*,*,1]= g
img[*,*,2]= b


satind=where(img[*,*,0] ge 8000. or img[*,*,1] ge 8000. or img[*,*,2] ge 8000.,count)

; scale and set colors
img = nw_scale_rgb(img,scales=scale)
img = nw_arcsinh_fit(img,nonlinearity=nonlinearity)
;img = nw_fit_to_box(img,origin=origin)
;img = nw_float_to_byte(img)

;return,im

;****************************

print, max(img)

b=bytscl(img,min=0,max=1)

if not keyword_set (imscale) then imscale = 1.
imsize=3
margin=.5
lmargin=.5
bmargin=.5
imnum=v-start
;tv,b,v-start,true=3,xsize=1,/inches,/order
if not keyword_set (image_path) then begin
    if not keyword_set(one_each) then begin
        xpos=(lmargin+(imnum mod 2)*(imsize+margin))*imscale
        ypos=(bmargin+(imnum/2)*(imsize+margin))*imscale
    end else begin
        xpos=lmargin*imscale
        ypos=bmargin*imscale
    end
    tv,b,xpos,ypos,true=3,xsize=imsize,/inches,/order
end else begin
    file_name = image_path + '/'+ file + "-" + $
                string (imnum, format = '(i3.3)')+".jpg"
    print, "Saving image " + file_name
    write_jpeg, file_name, b, true = 3,qual=99
end

end
end

pro west_all,base,view,alpha,Q,m,scale=scale,bands=bands, imscale = imscale, save_images= save_images,image_path=image_path,start_hdu=start_hdu,noerase=noerase

;spawn,"ls "+base+"_???.fits",files
files=["broadband_1.fits"]
start_file=files [0]
if keyword_set(start_hdu) then $
  start_hdu = find_HDU(start_file, start_hdu) $
else $
  start_hdu = find_HDU(start_file, "CAMERA0-BROADBAND")

print, 'Generating color images starting with HDU ', start_HDU
if keyword_set (save_images) then begin
    if not keyword_set(image_path) then $
      image_path = "./color/"
    spawn, "mkdir "+image_path
end

for i =0,n_elements(files) - 1 do begin
    west,files[i],view+start_hdu,alpha,Q,m, scale = scale, $
               bands = bands, imscale =imscale, image_path = image_path, $
               noerase=noerase
end
end
