#!/usr/bin/python
# -*- coding: utf-8 -*-

import math

N = 2
deginpix = 0.000366

class Star():
  def __init__(self):
    self.x1 = 0.0
    self.y1 = 0.0
    self.x2 = 0.0
    self.y2 = 0.0
    self.mag1 = 0.0
    self.mag2 = 0.0
    self.magerr1 = 0.0
    self.magerr2 = 0.0
    self.flux1 = 0.0
    self.flux2 = 0.0
    self.dist = 0.0
    self.ang = 0.0
    self.fwhm1 = 0.0
    self.fwhm2 = 0.0
  def calcdist(self):
    distdegr = math.sqrt((self.x1 - self.x2)**2 + (self.y1 - self.y2)**2)
    self.dist = distdegr/deginpix
  def calcang(self):
    ang = math.atan((self.y1 - self.y2) / (self.x1 - self.x2))
    if self.x1 - self.x2 < 0:
      ang = ang + math.pi
    elif self.y1 - self.y2 < 0 and self.x1 - self.x2 > 0:
      ang = 2*math.pi + ang
    self.ang = math.degrees(ang)

def parseline(line):
  x = float(line.split()[1])
  y = float(line.split()[2])
  mag = float(line.split()[3])
  magerr = float(line.split()[4])
  fwhm = float(line.split()[5])
  fwhm = fwhm/deginpix
  flux = 10**(-0.4*mag)
  return (x,y,mag,magerr,flux,fwhm)

def getdata(n):
  st = Star()
  fop = open("mag000"+str(n)+".dat",'r')
  line = fop.readline()
  st.x1,st.y1,st.mag1,st.magerr1,st.flux1,st.fwhm1 = parseline(line)
  line = fop.readline()
  st.x2,st.y2,st.mag2,st.magerr2,st.flux2,st.fwhm2 = parseline(line)
  fop.close()
  return st

def process(n):
  st = getdata(n)
  mag = math.log10(0.5*(st.flux1 + st.flux2))
  pol = 100*(st.flux1 - st.flux2)/(st.flux1 + st.flux2)
  flerr1 = (10**(-0.4*st.mag1) - 10**(-0.4*(st.mag1+st.magerr1)))/st.flux1
  flerr2 = (10**(-0.4*st.mag2) - 10**(-0.4*(st.mag2+st.magerr2)))/st.flux2
  flerr = 100*0.5*(flerr1 + flerr2)                                             #Check this
  st.calcdist()
  st.calcang()
  return st.x1,st.y1,st.x2,st.y2,mag,pol,flerr,st.dist,st.ang,st.fwhm1,st.fwhm2

def out(stars):
  fop = open("stand.dat",'w')
  for output in stars:
    fop.write(str(output[0])+" ")
    fop.write(str(output[1])+" ")
    fop.write(str(output[2])+" ")
    fop.write(str(output[3])+" ")
    fop.write(str(output[4])+" ")
    fop.write(str(output[5])+" ")
    fop.write(str(output[6])+" ")
    fop.write(str(output[7])+" ")
    fop.write(str(output[8])+" ")
    fop.write(str(output[9])+" ")
    fop.write(str(output[10])+"\n")
  fop.close()

def calcmean(stars):
  pmean = 0
  fwhm1mean = 0
  fwhm2mean = 0
  sigma = 0
  for star in stars:
    pmean = pmean + star[5]
    fwhm1mean = fwhm1mean + star[9]
    fwhm2mean = fwhm2mean + star[10]
  pmean = pmean/len(stars)
  fwhm1mean = fwhm1mean/len(stars)
  fwhm2mean = fwhm2mean/len(stars)
  for star in stars:
    sigma = sigma + (pmean - star[5])**2 #Maybe wrong in polpar4.bas str.35
  sigma = math.sqrt(sigma/len(stars))
  sigma = sigma/math.sqrt(len(stars))
  fop = open("1.dat",'w')
  fop.write(str(pmean)+" ")
  fop.write(str(sigma)+" ")
  fop.write(str(fwhm1mean)+" ")
  fop.write(str(fwhm2mean)+" ")
  fop.close()
  return pmean,sigma,fwhm1mean,fwhm2mean

def main():
  stars = []
  for i in range(1,N+1):
    stars.append(process(i))
  out(stars)
  starsparmean = calcmean(stars)
  

if __name__ == "__main__":
  main()
