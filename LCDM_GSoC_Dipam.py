# GSoC 2014 Solution to Warm-Up Question
# Submitted By -- Dipam Changede

from matplotlib.colors import LogNorm as lm
from astropy.io import fits as pyfits
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from numpy import mgrid, sum
from pylab import figure, cm
import numpy
import math

#calculating 2nd order moments Ixx, Ixy, Iyy
def moments2e(tbdata,centx,centy,thresh):
  moments = {}
  moments['Ixx']=0
  moments['Iyy']=0
  moments['Ixy']=0
  for i  in range (0,313):
     for j in range (0,313):
         if tbdata[i][j] >= (thresh):
             moments['Ixx'] =  moments['Ixx'] + ((j-centx)**2*tbdata[i][j])
             moments['Iyy'] =  moments['Iyy'] + ((i-centy)**2*tbdata[i][j])
             moments['Ixy'] =  moments['Ixy'] + ((j-centx)*(i-centy)*tbdata[i][j])

  return moments

#Reading given FITS file
hdulist = pyfits.open('F:\M86.fits')
tbdata = hdulist[0].data
meann= numpy.mean(tbdata) # mean of all values of the image
stdn= numpy.std(tbdata) #standard deviation
threshold = meann+stdn # threshold set as mean + one std dev
weights=0
weig_sum_x = 0
weig_sum_y = 0

#loop for calculating weighted sum of all eligible pixels, i.e pixels whose intensity above the threshold
for i  in range (0,313):
    for j in range (0,313):
        if tbdata[i][j] >= (threshold):
            weig_sum_x = weig_sum_x + (tbdata[i][j]*j)
            weig_sum_y = weig_sum_y + (tbdata[i][j]*i)
            weights = weights + tbdata[i][j]
            
# location of x and y coordinate of centroid rounded to nearest int
centerx = round(float(weig_sum_x)/float(weights)) 
centery = round(float(weig_sum_y)/float(weights)) 
plt.plot(centerx,centery,'ro') #plotting Centroid on image
moments = moments2e(tbdata,centerx,centery,threshold) #calling moments2e function
angle = (math.atan(2 * moments['Ixy'] / (moments['Iyy'] - moments['Ixx']))) / 2 #calculating angle 
angledeg = math.degrees(angle)
angle2=(angle+1.57079633)

#plotting moment axis on the image
linelength = 15
linelength1 = 20
x = [centerx-linelength*math.cos(angle), centerx+linelength*math.cos(angle)]
y = [centery-linelength*math.sin(angle), centery+linelength*math.sin(angle)]
x1 = [centerx-linelength1*math.cos(angle2), centerx+linelength1*math.cos(angle2)]
y1 = [centery-linelength1*math.sin(angle2), centery+linelength1*math.sin(angle2)]
plt.plot(x,y,c='yellow')
plt.plot(x1,y1,c='green')
imgplot = plt.imshow(tbdata, norm=lm(meann + 0.5 * stdn, numpy.max(tbdata), clip='True'),cmap=cm.gray, origin="lower")

#customizing the output axis
plt.xlim(60,180)   
plt.xticks([60,80,100,120,140,160,180])
plt.ylim(100,210)
plt.yticks([100,120,140,160,180,200])

# axis label and title of th output image
plt.title("GSoC 2014 Solution")
plt.xlabel("X")
plt.ylabel("Y")

#image saved in the current directory with name "GSoC_2014_Sol_Dipam.png"
plt.savefig("GSoC_2014_Sol_Dipam.png",bbox_inches='tight')
plt.show()
