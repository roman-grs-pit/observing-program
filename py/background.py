import numpy as np
from matplotlib import pyplot as plt

#d = np.loadtxt('back_testout.txt').transpose()
d = np.loadtxt('background_2027all_daystep10_1.4.txt').transpose()

ra1 = 60.0
dec1 = -15.14
ra2 = 30.511950 
dec2 = -51.66
ra3 = 89.488050
dec3 = -51.66
ra4 = 39.316320
dec4 = -84.86

s1 = d[0] == ra1
s1 &= d[1] == dec1

s2 = d[0] == ra2
s2 &= d[1] == dec2

s3 = d[0] == ra3
s3 &= d[1] == dec3

s4 = d[0] == ra4
s4 &= d[1] == dec4

plt.plot(d[-5][s1],d[-8][s1],label='ra,dec='+str(ra1)+','+str(dec1))
plt.plot(d[-5][s2],d[-8][s2],label='ra,dec='+str(ra2)+','+str(dec2))
plt.plot(d[-5][s3],d[-8][s3],label='ra,dec='+str(ra3)+','+str(dec3))
plt.plot(d[-5][s4],d[-8][s4],label='ra,dec='+str(ra4)+','+str(dec4))
plt.xlabel('day in 2027')
plt.ylabel('total background (MJy/sr)')
plt.legend()
plt.savefig('Testback_4corn.png')
plt.show()
plt.plot(d[-5][s4],d[-8][s4],label='ra,dec='+str(ra4)+','+str(dec4))
plt.xlabel('day in 2027')
plt.ylabel('total background (MJy/sr)')
plt.legend()
plt.savefig('Testback_60n50.png')
plt.show()

plt.plot(d[-5][s1],d[-12][s1],label='ra,dec='+str(ra1)+','+str(dec1))
plt.plot(d[-5][s2],d[-12][s2],label='ra,dec='+str(ra2)+','+str(dec2))
plt.plot(d[-5][s3],d[-12][s3],label='ra,dec='+str(ra3)+','+str(dec3))
plt.plot(d[-5][s4],d[-12][s4],label='ra,dec='+str(ra4)+','+str(dec4))
plt.xlabel('day in 2027')
plt.ylabel('zodi background (MJy/sr)')
plt.legend()
plt.savefig('Testzodiback_4corn.png')
plt.show()
plt.plot(d[-5][s4],d[-12][s4],label='ra,dec='+str(ra4)+','+str(dec4))
plt.xlabel('day in 2027')
plt.ylabel('total background (MJy/sr)')
plt.legend()
plt.savefig('Testzodiback_60n50.png')
plt.show()
