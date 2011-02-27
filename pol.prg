DEFINE/PARAMETER P1 + F "Enter day number:"
DEFINE/MAXPAR 1
DEFINE/LOCAL name/C*15/1/4
DEFINE/LOCAL mean/R/1/1 0.0
DEFINE/LOCAL stdev/R/1/1 0.0
DEFINE/LOCAL mfwhmx/R/1/1 0.0
DEFINE/LOCAL mfwhmy/R/1/1 0.0
DEFINE/LOCAL mfwhm/R/1/1 0.0
DEFINE/LOCAL arrfwhm/R/1/2 0.0
DEFINE/LOCAL step/R/1/2 0.0,0.0
DEFINE/LOCAL apert/R/1/1 0.0
DEFINE/LOCAL Fsiz/R/1/1 0.0
DEFINE/LOCAL n/i/1/1 0
DEFINE/LOCAL i/i/1/1 0
DEFINE/LOCAL k/i/1/1 0

COMPUTE/KEYWORD name(1) = "stars" // P1 // "Qrs.fts"
COMPUTE/KEYWORD name(2) = "stars" // P1 // "Urs.fts"
COMPUTE/KEYWORD name(3) = "hartley" // P1 // "Qrs.fts"
COMPUTE/KEYWORD name(4) = "hartley" // P1 // "Urs.fts"

DO i = 1 4
  !loading image
  CLEAR/CHAN 2
  INDISK/FITS {name({i})} tmpim
  WRITE/OUT {name({i})}
  STATIS/IMAG tmpim ? ? ? ? stattab.tbl
  mean = m$value(stattab,:MEAN,@1)
  stdev = m$value(stattab,:STD,@1)
  WRITE/OUT mean = {mean}
  WRITE/KEYW lcut/R/1/1 1.0
  lcut = {mean} - 1 * {stdev}}
  WRITE/KEYW hcut/R/1/1 1.0
  hcut = {mean} + 1 * {stdev}
  LOAD/IMA tmpim scale=2,2 cuts={lcut},{hcut}
  
  !calculate fwhm
  IF i .LT. 3 THEN !field stars
    !getting good stars
    CENTER/GAUS ? fwhm.tbl
    READ/TABL   fwhm.tbl

    !calculate mean fwhm
    STATISTICS/TABLE fwhm.tbl :XFWHM
    mfwhmx = {OUTPUTR(3)}
    WRITE/OUT mean fwhm x = {mfwhmx}

    STATISTICS/TABLE fwhm.tbl :YFWHM
    mfwhmy = {OUTPUTR(3)}
    WRITE/OUT mean fwhm y = {mfwhmy}

    mfwhm = ({mfwhmy} + {mfwhmx})*0.5
    WRITE/OUT mean FWHM = {mfwhm}
    arrfwhm({i}) = {mfwhm}
  ELSE !comet
    k = {i} - 2
    mfwhm = {arrfwhm({k})}
    WRITE/OUT mean FWHM = {mfwhm}
  ENDIF

  !calculate optimal aperture size
  COPY/DKEY tmpim STEP step
  apert = M$ABS(1.55*{mfwhm}/{step} + 1)
  WRITE/OUT aperture = {apert}
  Fsiz = {apert} * 2 + 1
  WRITE/OUT Fsiz = {Fsiz}

  WRITE/OUT "Now mark pairs for the measurment"
  IF i .LT. 3 THEN !field stars
    DO n = 1 2
      CENTER/GAUS ? polapar.tbl
      MAGNIT/CIRC tmpim,polapar.tbl mag.tbl @{Fsiz},@4,@5

      !FWHM calculation
      CREATE/COLU mag.tbl FWHM "pix" E12.5 R*4
      CREATE/COLU mag.tbl XFWHM "pix" E12.5 R*4
      CREATE/COLU mag.tbl YFWHM "pix" E12.5 R*4
      COPY/TT     polapar.tbl :XFWHM mag.tbl :XFWHM
      COPY/TT     polapar.tbl :YFWHM mag.tbl :YFWHM
      COMPUT/TABL mag.tbl :FWHM = 0.5 * (:XFWHM + :YFWHM)

      !output into ascii
      ASSIGN/PRIN FILE {name({i})}{n}.dat
      PRINT/TABL  mag.tbl :XCEN :YCEN :MAGNITUDE :MAG_SIGMA :FWHM N 200
    ENDDO
  ELSE !comet
    CENTER/GAUS ? polapar.tbl
    MAGNIT/CIRC tmpim,polapar.tbl mag.tbl @{Fsiz},@4,@5

    !FWHM calculation
    CREATE/COLU mag.tbl FWHM "pix" E12.5 R*4
    CREATE/COLU mag.tbl XFWHM "pix" E12.5 R*4
    CREATE/COLU mag.tbl YFWHM "pix" E12.5 R*4
    COPY/TT     polapar.tbl :XFWHM mag.tbl :XFWHM
    COPY/TT     polapar.tbl :YFWHM mag.tbl :YFWHM
    COMPUT/TABL mag.tbl :FWHM = 0.5 * (:XFWHM + :YFWHM)

    !output into ascii
    ASSIGN/PRIN FILE {name({i})}.dat
    PRINT/TABL  mag.tbl :XCEN :YCEN :MAGNITUDE :MAG_SIGMA :FWHM N 200
  ENDIF
ENDDO




!DEFINE/LOCAL n/i/1/1 0
!DEFINE/LOCAL k/i/1/1 0
!DEFINE/LOCAL l/i/1/1 0
!CENTER/GAUS ? tmptab.tbl
!DO n = 1 12
!  k = {n} * 2
!  l = (24-({n}*2))
!  MAGNIT/CIRC qqq.bdf,tmptab.tbl mag{n}.tbl @{k},@{l},@10
!  COMPUTE/TABLE mag{n}.tbl :SB = :MAGNITUDE/(3.14159265*({k}*0.5*0.395034)**2)
!ENDDO
!MAGNIT/CIRC qqq.bdf,tmptab.tbl mag1.tbl @2,@22,@10
!MAGNIT/CIRC qqq.bdf,tmptab.tbl mag12.tbl @24,@0,@10

!COPY/TABLE mag0001.tbl magall.tbl
!DO n = 2 12
!  merge/table magall.tbl mag{n}.tbl magalltmp.tbl
!  COPY/TABLE magalltmp.tbl magall.tbl
!ENDDO


!ASSIGN/PRIN FILE mag.dat
!PRINT/TABL  magall.tbl :MAGNITUDE :MAG_SIGMA :SB
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!crea/disp
!load/imag LINE.bdf scale=2,2
!DEFINE/LOCAL nofiter/i/1/1 10
!DEFINE/LOCAL n/i/1/1 0
!DEFINE/LOCAL nol/I/1/10 127,135,145,153,160,170,175,180,190,198
!WRITE/OUT nofiter = {nofiter}

!DO n = 1 nofiter

!WRITE/OUT nol = {nol({n})} 
!WRITE/OUT n = {n}
!EXTR/RTRAC ? {nol({n})}.bdf PLOT
!CENTER/GAUSS {nol({n})}.bdf {nol({n})}.tbl
!crea/column {nol({n})}.tbl :Y R*4 "px" E12.5
!write/table {nol({n})}.tbl :Y @1 {nol({n})}
!ENDDO

!merge/table 0127.tbl 0135.tbl 0145.tbl 0153.tbl 0160.tbl 0170.tbl 0175.tbl ALLLINES1.tbl
!merge/table 0175.tbl 0180.tbl 0190.tbl 0198.tbl ALLLINES2.tbl
!merge/table ALLLINES1.tbl ALLLINES2.tbl ALLLINES.tbl
!COMPUTE/TABLE ALLLINES.tbl :V = 299792458*(:XCEN - 6562.82)/6562.82
!PLOT/TABLE ALLLINES.tbl :Y :V
!COMPUTE/BARYCORR u07s.bdf 16,58,30 58,56,44 41,44166667 43,6533333
!COMPUTE/TABLE ALLLINES.tbl :Vhel = :V - 20
!COMPUTE/TABLE ALLLINES.tbl :Vgal = :Vhel + 238.1
!COMPUTE/TABLE ALLLINES.tbl :Arcsec = :Y + 0.409
!COMPUTE/TABLE ALLLINES.tbl :Arcsec = :Arcsec * 0.4
!PLOT/TABLE ALLLINES.tbl :Arcsec :Vhel
