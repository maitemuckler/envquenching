library(SciServer)

Authentication_loginName = 'maitemuckler';
Authentication_loginPassword = '**Luamor27**'
token1 = Authentication.login(Authentication_loginName, Authentication_loginPassword)

#get an image cutout
SkyServer_DataRelease = "DR16"
img = SkyServer.getJpegImgCutout(ra=177.33785, dec=-0.38458351, 
                                 width=512, height=512, scale=0.1, opt = "GLSQM",
                                 dataRelease=SkyServer_DataRelease,
                                 query="SELECT TOP 100 p.objID, p.ra, p.dec, p.r FROM fGetObjFromRectEq(197.6,18.4,197.7,18.5) n, PhotoPrimary p WHERE n.objID=p.objID")

#png('file.png') 
plot(0:1, 0:1, type = "n")
rasterImage(img, 0, 0, 1, 1)
dev.off()


# wget --output-document=Images/1_image.jpg 'http://skyserver.sdss.org/dr12/SkyserverWS/ImgCutout/getjpeg?ra=226.783793&dec=-2.662714&scale=0.100000&width=512&height=512'
# wget --output-document=Images/1_spec.png 'http://skyserver.sdss.org/dr14/en/get/specById.ashx?ID=1038091849438357504'


system("wget --output-document=image.jpg 'http://skyserver.sdss.org/dr12/SkyserverWS/ImgCutout/getjpeg?ra=226.783793&dec=-2.662714&scale=0.100000&width=512&height=512'")
system("wget --output-document=spec.png 'http://skyserver.sdss.org/dr14/en/get/specById.ashx?ID=1038091849438357504'")
