DEFINE/PARAMETER P1 + F "Enter day number:"
DEFINE/MAXPAR 1
DEFINE/LOCAL name/C*20/1/4
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
    DO n = 1 5
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

