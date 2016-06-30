function tp,file

restore,'/net/grads-1/stinson/trace/tipsyprofile.sav'
p=read_ascii(file,template=temp)

return,p
END
