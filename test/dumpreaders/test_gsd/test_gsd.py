import dumpreader
import gsd
import gsd.hoomd

gsdf = gsd.hoomd.open('init.gsd','rb')
print(gsdf[0].particles.position)


