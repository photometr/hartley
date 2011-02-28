#!/usr/bin/python
# -*- coding: utf-8 -*-

import math

N = 5
date = "19"
tetacam = 0
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

def getdata(fname):
  st = Star()
  fop = open(fname,'r')
  line = fop.readline()
  st.x1,st.y1,st.mag1,st.magerr1,st.flux1,st.fwhm1 = parseline(line)
  line = fop.readline()
  st.x2,st.y2,st.mag2,st.magerr2,st.flux2,st.fwhm2 = parseline(line)
  fop.close()
  return st

def process(fname):
  st = getdata(fname)
  mag = math.log10(0.5*(st.flux1 + st.flux2))
  pol = 100*(st.flux1 - st.flux2)/(st.flux1 + st.flux2)
  flerr1 = (10**(-0.4*st.mag1) - 10**(-0.4*(st.mag1+st.magerr1)))/st.flux1
  flerr2 = (10**(-0.4*st.mag2) - 10**(-0.4*(st.mag2+st.magerr2)))/st.flux2
  flerr = 100*0.5*(flerr1 + flerr2)                                             #Check this
  st.calcdist()
  st.calcang()
  return st.x1,st.y1,st.x2,st.y2,mag,pol,flerr,st.dist,st.ang,st.fwhm1,st.fwhm2

def out(stars,filt):
  fop = open(filt+".dat",'w')
  fop.write("x1  y1  x2    y2   mag   pol   flerr    dist    ang       fwhm1      fwhm2\n")
  for output in stars:
    fop.write(str(output[0])+" ")#x1
    fop.write(str(output[1])+" ")#y1
    fop.write(str(output[2])+" ")#x2
    fop.write(str(output[3])+" ")#y2
    fop.write(str(output[4])+" ")#mag
    fop.write(str(output[5])+" ")#pol
    fop.write(str(output[6])+" ")#flerr
    fop.write(str(output[7])+" ")#dist
    fop.write(str(output[8])+" ")#ang
    fop.write(str(output[9])+" ")#fwhm1
    fop.write(str(output[10])+"\n")#fwhm2
  fop.close()

def calcmean(stars,comet,filt):
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
  flerr = comet[0][6]
  pol = comet[0][5]
  sigma1 = math.sqrt(sigma**2 + flerr**2)
  pol = pol - pmean

  fop = open(filt+"1.dat",'w')
  fop.write(str(comet[0][4])+" ")
  fop.write(str(pol)+" ")
  fop.write(str(sigma)+" ")
  fop.write(str(sigma1)+" ")
  fop.write(str(fwhm1mean)+" ")
  fop.write(str(fwhm2mean)+" ")
  fop.close()
  return comet[0][4],pol,sigma,sigma1,fwhm1mean,fwhm2mean

def calcPol(comparmean):
  mQ = comparmean[0][0]
  mU = comparmean[1][0]
  Q = comparmean[0][1]
  U = comparmean[1][1]
  sigmaQ = comparmean[0][2]
  sigmaU = comparmean[1][2]
  meanFWHM = comparmean[0][4]+comparmean[0][5]
  meanFWHM = meanFWHM + comparmean[1][4]+comparmean[1][5]
  meanFWHM = meanFWHM/4
  P = math.sqrt(Q**2 + U**2)
  sigmaP = math.sqrt(sigmaQ**2 + sigmaU**2)
  sigmat = 28.7 * sigmaP / P
  teta = math.atan(U/Q)
  if Q < 0 and U > 0: teta = teta + math.pi
  if Q < 0 and U < 0: teta = teta + math.pi
  if Q > 0 and U < 0: teta = teta + 2 * math.pi
  teta = math.degrees(teta)/2 - tetacam
  if teta < 0: teta = teta + 180
  fop = open("res"+str(date)+".dat",'w')
  fop.write(str(0.5*(mQ+mU))+" ")
  fop.write(str(P)+" ")
  fop.write(str(sigmaP)+" ")
  fop.write(str(teta)+" ")
  fop.write(str(sigmat)+" ")
  fop.write(str(meanFWHM)+"\n")
  fop.close()
  return 0

def main():
  stars = {}
  comet = {}
  comparmean = []
  for filt in ["Q","U"]:
    stars[filt] = []
    comet[filt] = []
    for i in range(1,N+1):
      stars[filt].append(process("stars"+date+filt+"rs.fts000"+str(i)+".dat"))
    out(stars[filt],filt)
    comet[filt].append(process("hartley"+date+filt+"rs.fts.dat"))
    out(comet[filt],filt+"comet")
    comparmean.append(calcmean(stars[filt],comet[filt],filt))
  calcPol(comparmean)


if __name__ == "__main__":
  main()
