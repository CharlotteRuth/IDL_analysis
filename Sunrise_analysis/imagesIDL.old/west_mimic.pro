pro west,file,hdu,alpha,Q,m,scale=scale,imscale = imscale, image_path=image_path,one_each=one_each,noerase=noerase

print,''
print,'Chose 3 bands, in form x,y,z for red green and blue from:' 
print,'0 = GALEX FUV'
print,'1 = GALEX NUV'  
print,'2 = SDSS u band'
print,'3 = SDSS g band'
print,'4 = SDSS r band'
print,'5 = SDSS i band'
print,'6 = SDSS z band'
print,'7 = HST ACS F435'
print,'8 = HST ACS F606'
print,'9  = HST ACS F775'
print,'10 = HST ACS F850'
print,'11 = IRAC1 SIRTF'   
print,'12 = IRAC2 SIRTF'   
print,'13 = MIPS24 SIRTF'   
print,'14 = MIPS70 SIRTF'  
print,'15 = MIPS160 SIRTF'   

;bands = intarr(3)
;print,'bands ='
;read,bands

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
if not keyword_set(scale) then scale=[1.5,0.9,.9]
if not keyword_set (bands) then bands = [13,12,11]
red = image (*,*, bands [0]) ; GALEX FUV
green = image (*,*, bands [1]) ; GALEX NUV  
blue =  image (*,*, bands [2]) ; SDSS u band
 
;r = (green)*scale[0]
;g = (4*green+blue)*scale[1]
;b = (blue)*scale[2]

r = (red)*scale[0]
g = (green)*scale[1]
b = (blue)*scale[2]


ii = (r+g+b)/3+1e-20
print,max(ii)

input= alpha*Q*(ii-m)

rr = r*alog(input+(input*input+1)^0.5)/(Q*ii)
gg = g*alog(input+(input*input+1)^0.5)/(Q*ii)
bb = b*alog(input+(input*input+1)^0.5)/(Q*ii)

img = [[[rr] ], [[gg] ], [[bb] ]]

print, max(img)

b=bytscl(img,min=0,max=1)

if not keyword_set (imscale) then imscale = 1.0
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

pro west_all,base,view,alpha,Q,m,scale=scale, imscale = imscale, save_images= save_images,image_path=image_path,start_hdu=start_hdu,noerase=noerase

spawn,"ls "+base+"_???.fits",files
files[0]="broadband.fits"
start_file=files[0]
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
               imscale =imscale, image_path = image_path, $
               noerase=noerase
end
end
