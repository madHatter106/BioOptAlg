; docformat:'rst'

 
;+
; A group of procedures and functions for transferring absorption and scattering coefficient of pure seawater from given files to an IDL structure, and 
; for querying this structure to obtain the pure seawater IOPs at specific wavelengths.
;
; :version: 3
; 
; :Author:
;
;   Frédéric Mélin (European Commission/JRC/IES/WRES)
;   
;   Gert Sclep (European Commission/JRC/IES/WRES)
;
; :Copyright:
; 
;   All rights reserved
;-   

;+
; :Author:
;   
;   Frédéric Mélin 
;
; :version: 1
;   
; :History:
; 
;   Created by Frédéric Mélin, 10/2010
;   
; :Categories:
; 
;   IOP, aw, bw, absorption, scattering, pure water
;   
; :Description:
; 
; Read absorption and scattering of pure seawater from tabulated values.
; 
; :Examples:
;
;   Example syntax::
;     filename = 'water_spectra.dat'
;     swater = ReadWater(filename)  
;   
; :Pre:
; 
; Indicated file needs to contain data in three columns sepeated by a space. The first column should contain the wavelength (at chosen precision),
; the second column the absorption in m-1 and the third column the scattering in m-1.
; 
; Standard, the absorption is taken from Pope & Fry (1997) and Kou et al. (1993), and the pure water scattering from 
; Smith & Baker (1981), according to::
;   - Kou, L., D. Labrie, P. Chylek, 1993, "Refractive indices of water and ice in the 0.65-2.5 m spectral range," Appl. Opt.,32, 3531-3540 (1993).
;   - Pope, R. M. and E. S. Fry, 1997, "Absorption spectrum (380-700 nm) of pure water. II. Integrating cavity measurements," Appl. Opt.,36, 8710-8723.
;   - Smith, R.C. and K.S. Baker, 1981, "Optical properties of the clearest natural waters (200-800 nm)," Appl. Opt.,20, 177-184.
; 
; See file 'water_spectra.dat' (taken from SeaDAS 6.4 distribution, distributed along with this code)
;   
; :Requires:
; 
;   IDL 7.1
;   
; :Returns:
;  
;   Structure containing following fields::
;   
;     lambda: wavelength in nm (type=double)
;     aw: pure water absorption in m-1 (type=double)
;     bw: pure water scattering in m-1 (type=double)
;     
; :Params:
; 
;   filename: in, required, type=string
;     input file (ex: 'water_spectra.dat', as taken from SeaDAS 6.4)
;         
;-
FUNCTION readWater,filename

aw = dblarr(3000)
bw = dblarr(3000)
lambda = dblarr(3000)

CLOSE,1
OPENR,1,filename

stringline = ' '

FOR i=0,29-1 DO readf,1,stringline
il=0
WHILE NOT eof(1) DO BEGIN

    stringline=' '
    
    readf,1,stringline
    
    IF ( stringline NE '' ) THEN BEGIN
    
       strarray = str_sep( strcompress(stringline),' ')
     
       stringarray = strarray
       ee=0
       FOR e=0,n_elements(strarray)-1 DO BEGIN
          IF (strlen(strcompress(strarray(e),/remove_all)) NE 0) THEN BEGIN
              stringarray(ee) = strarray(e)
              ee = ee +1
          ENDIF
       ENDFOR
       stringarray = stringarray(0:ee-1)

    ENDIF ELSE BEGIN
         print,"ERROR"
         STOP
    ENDELSE
    
    lambda[il] = double(stringarray[0])
    aw[il] = double(stringarray[1])
    bw[il] = double(stringarray[2])

    il = il + 1
         
ENDWHILE

CLOSE,1

NBANDS = il

aw = aw(0:NBANDS-1)
bw = bw(0:NBANDS-1)
lambda = lambda(0:NBANDS-1)

swater =  { lambda:lambda, aw:aw, bw:bw }

RETURN,swater

END

;+
; :Author:
;   
;   Frédéric Mélin 
;
; :version: 1
;   
; :History:
; 
;   Created by Frédéric Mélin, 10/2010
;   
; :Categories:
; 
;   IOP, bw, scattering, pure water
;   
; :Description:
; 
;   Read scattering coefficients of pure seawater from tabulated values.
;   
; :Examples:
;   
;   Example syntax::
;   filename = 'bw_spectrum_s35t22.dat'
;   swater = ReadbWater(filename)
;
; :Pre:
;   
;   Indicated file needs to contain data in two columns sepeated by a space. The first column should contain the wavelength (at chosen precision) and the second column the scattering in m-1.
;   
;   Standard, the scattering values are calculated for a given temperature and salinity according to Zhang, 2009::
;     X., Hu, L., He, M.-X.: Scattering by pure seawater: Effect of salinity. Optics Express, 17, 5698-5710, 2009.
;     
;   See file bw_spectrum_s35t22.dat (distributed along with this code).
; 
; :Requires:
;   
;   IDL 7.1
; 
; :Returns:
; 
;   Structure containing following fields::
;   
;     lambda: wavelength in nm (type=double)
;     bw: pure water scattering in m-1(type=double)
;
; :Params:
; 
;   filename: in, required, type=string
;     input file (ex: 'bw_spectrum_s35t22.dat', pure water scattering for temperature 22 deg C and salinity of 35, calculated according to Zhang, 2009 (see below))
;         
;-
FUNCTION readbWater,filename
  
NBANDS=2250
NMAX = 2500
bw = DBLARR(NMAX)
lambda = DBLARR(NMAX)

close,1
openr,1,filename
il = 0

stringline = ' '
ip=-1
WHILE ( ip LT 0 ) DO BEGIN

   readf,1,stringline
   ip = STRPOS(stringline,'end_header')

ENDWHILE
WHILE NOT EOF(1) DO BEGIN

    stringline=' '
    
    readf,1,stringline
    
    IF ( stringline NE '' ) THEN BEGIN
    
       strarray = STR_SEP( STRCOMPRESS(stringline),' ')
     
       stringarray = strarray
       ee=0
       FOR e=0,N_ELEMENTS(strarray)-1 DO BEGIN
           IF (STRLEN(STRCOMPRESS(strarray(e),/REMOVE_ALL)) NE 0) THEN BEGIN
              stringarray(ee) = strarray(e)
              ee = ee +1
          ENDIF
       ENDFOR
       stringarray = stringarray(0:ee-1)

    ENDIF ELSE BEGIN
         print,"ERROR"
         stop
    ENDELSE
    
    lambda[il] = DOUBLE(stringarray[0])
    bw[il] = DOUBLE(stringarray[1])

    il = il + 1
         
ENDWHILE

close,1
IF ( il NE NBANDS ) THEN STOP

bw = bw(0:NBANDS-1)
lambda_w = lambda(0:NBANDS-1)

sdata = { lambda:lambda_w, bw:bw }

RETURN, sdata

END


;+
; :Author:
;   
;   Frédéric Mélin 
;
; :version: 1
;   
; :History:
; 
;   Created by Frédéric Mélin
;   
; :Categories:
;   
;   IOP, pure water
;   
; :Description:
; 
;   Read pure water spectral properties, absorption and scattering coefficient, into an IDL structure.
;   
;   When choosing to mimick SeaDAS, see `mimick_seadas`, the pure water absorption and scattering are read from the 'water_spectra.dat' of the SeaDAS 6.4 distribution.
;   This means the absorption is taken from Pope & Fry (1997) and Kou et al. (1993), and the pure water scattering from from Smith & Baker (1981), according to::
;   - Kou, L., D. Labrie, P. Chylek, 1993, "Refractive indices of water and ice in the 0.65-2.5 m spectral range," Appl. Opt.,32, 3531-3540 (1993).
;   - Pope, R. M. and E. S. Fry, 1997, "Absorption spectrum (380-700 nm) of pure water. II. Integrating cavity measurements," Appl. Opt.,36, 8710-8723.
;   - Smith, R.C. and K.S. Baker, 1981, "Optical properties of the clearest natural waters (200-800 nm)," Appl. Opt.,20, 177-184.
;   
;   When opting not to mimick SeaDAS, the scattering is taken from file 'bw_spectrum_s35t22.dat'. This means that pure watter scattering is calculated for a temperature
;   of 22 degrees Celcius and a salinity of 35, according to Zhang, 2009::
;     X., Hu, L., He, M.-X.: Scattering by pure seawater: Effect of salinity. Optics Express, 17, 5698-5710, 2009.
; 
; :Pre:
;   Files 'bw_spectrum_s35t22.dat' and 'water_spectra.dat' need to be present in the working directory.
;   
; :Requires:
; 
;   IDL 7.1
;   
; :Params:
;   
;   idir: in, type=string
;			
;     directory containing pure water IOPs
;
;   sdata: out, type=structure
;   
;     Structure containing following fields::
;   
;       lambda: wavelength in nm from 200.0 to 2449.0 (type=double)
;     
;       aw: pure water absorption in m-1 for wavelengths in lambda (type=double)
;     
;       bw: pure water scattering in m-1 for wavelengths in lambda (type=double) 
;
; :Keywords:
; 
;   seadas_mimick:  in, optional, type=int
;     if set to 1, the pure water scattering coefficients from Smith & Baker are not replaced by those of Zhang (for T 22 deg C and salinity 35)
;-
PRO readSpectralDataPureWater,idir,sdata, seadas_mimick=seadas_mimick

  ; Builds pure water IOPs
  filename = STRCOMPRESS(idir+'/'+'bw_spectrum_s35t22.dat',/REMOVE_ALL)
  b=ReadbWater(filename)  
  filename = STRCOMPRESS(idir+'/'+'water_spectra.dat',/REMOVE_ALL)
  sdata=ReadWater(filename)  
  ll = ABS(b.lambda-sdata.lambda)
  IF ( MAX(ll) GT 0. ) THEN BEGIN
     print,'problem with lambda' & STOP
  ENDIF  
  IF (NOT(KEYWORD_SET(seadas_mimick))) THEN sdata.bw(*) = b.bw(*) ; substitute with Zhang et al. values

END

;+
; :Author:
;   
;   Frédéric Mélin 
;   
;   Gert Sclep
;
; :version: 3
;    
; :History:
; 
;   Created by Frédéric Mélin
;   
;   Modified by Gert Sclep, 07/2010:: 
;       
;       Slight modification of routine, making the function identical to the aw_spectra and bbw_spectra functions of SeaDAS 6.1 update 1
;       
;       Now the step size between the wavelengths in the input file (read by ReadSpectralData) is also taken into account.
;       
;       Instead of half bandwidth now the entire bandwidth has to be given as input.
;       
;   Modified by Gert Sclep, 06/07/2012::
;   
;      Bug fix in calculation of wave_index   
;   
; :Categories:
;   
;   IOP, pure water
;   
; :Description:
; 
;   Get pure water spectral properties (aw and bw), for a given wavelength, using a structure of previously read in absorption and scattering coefficients.
;   
;   If the given wavelength is within the structure's range of the wavelengths, the structure's wavelength closest to the given wavelength and lower or equal to the given wavelength is taken as a reference.
;   
;   Either the IOPs at that reference wavelength are returned, if no bandwidth averaging is requested. If bandwidth averaging is requested however, the average IOPs are returned of the wavelength range defined as::
;   
;     from: reference wavelength - the downwards rounded half of the bandwidth 
;     to: the reference wavelength + the downwards rounded half of the bandwidth
;   
; :Examples:
;   
;   Example syntax::
;     readSpectralData, sdata
;     st = GetSpectralData(sdata, 400., aw, bw)
;     if (st EQ 1) then print, aw, bw
;     
; :Pre:
; 
;   Structure containing wavelengths, pure water absorption and scattering coefficients has to be previously created, according to readSpectralData.
;   
; :Requires:
; 
;   IDL 7.1
;   
; :Returns:
; 
;   -1 when requested wavelength was out of range, 0 when the corresponding aw and bw could be retrieved
;   
; :Params:
; 
;   sdata: in, required, type=structure
;     Structure according to output of readSpectralData containing following fields::
;   
;      lambda: wavelength in nm (type=double)
;      
;      aw: pure water absorption in m-1 for wavelengths in lambda (type=double)
;      
;      bw: pure water scattering in m-1 for wavelengths in lambda (type=double)
;
;   wave: in, required, type=float
;     wavelength for which pure water a and b need to be retrieved
;
;   aw: out, type=float
;     pure water absorption coefficient for the given wavelength
;     
;   bw: out, type=float
;     water scattering coefficient for the given wavelength
;
; :Keywords:
;   bandwidth: in, optional, type=int, default=10
;     size of window around chosen wavelength (or more exactly: around wavelength in sdata, closest to the chosen wavelength and lower or
;     equal to chosen wavelength) over which the a and b values will be averaged (from minus downwards rounded half of bandwidth to
;     plus downwards rounded half of bandwidth)
;     
;     if no bandwidth averaging is desired, bandwidth should be set to 0
;- 
FUNCTION getSpectralDataPureWater, sdata, wave, aw, bw, bandwidth=bandwidth

    routine  = 'GetSpectralData'

    IF N_ELEMENTS(bandwidth) EQ 1 THEN width=bandwidth ELSE width=10
  

; ----------------------------------------------------------------
;   Tests on wave
; ----------------------------------------------------------------

    step_size = sdata.lambda[1]-sdata.lambda[0]
    num_waves = N_ELEMENTS(sdata.lambda)

    wave_index = FIX((wave - sdata.lambda[0])/step_size, TYPE=3)
    
    IF (wave_index LT 0 OR wave_index GT num_waves ) THEN BEGIN

        print, routine+' - ERROR : wave out of range ', wave
        RETURN, -1

    ENDIF


; ----------------------------------------------------------------
;   Get interval and calculate averages
; ----------------------------------------------------------------

    imin = MAX([wave_index - width/2./step_size,    0])
    imax = MIN([wave_index + width/2./step_size, num_waves])    
    
    x  = sdata.aw(imin:imax)
    aw = MEAN(x)

    y  = sdata.bw(imin:imax)
    bw = MEAN(y)
    
    RETURN, 0

END


