; NAME:
;       SpecPro
; 
; PURPOSE: 
;       Interactive widget display of astronomical spectra and
;       associated stamp images / spectral energy distribution. 
;
;       Note that the routine can be used to view a subset of 
;       this data.
;
; CALLING SEQUENCE:
;       specpro [,slitnumber], [,/histogram], [,/space], [/OH], [,/small], [,/basic]
;     
; REQUIRED INPUTS:
;       None. If specpro is run with no inputs, the window widgets
;       are realized and a spectrum can subsequently be passed to the 
;       routine by selecting a slit number from within the interface.
;
;       The routine searches for five files for the currently selected
;       slit number: 
; 
;       1D spectrum file (fits/ascii)
;       2D spectrum file (fits)
;       Stamp image file (fits)
;       Photometry file  (ascii)
;       Information file (ascii)
; 
;       None of these files is necessary for the program to run. Refer
;       to the website http://specpro.caltech.edu for a description of
;       the expected data formats.
;
; OPTIONAL INPUTS:
;       slitnumber: A valid slit from the current working directory. 
;       The files are assumed to be organized by slit number on a
;       multislit mask, but any sequential numbering scheme can be
;       used to differentiate them.
;
;       /histogram: If given, displays the 1D specrum as a histogram.
;
;       /space: Use this keyword if you are viewing a spectrum taken
;       in space, e.g. WFC3 grism spectra. This will prevent the sky
;       bands from being displayed.
;
;       /OH: If given, overlays the positions of prominent near-IR OH emission 
;       features on the 1D.
;
;       /small: If the keyword /small is given, a smaller version of specpro
;       is started, more suitable for laptops or smaller 
;       monitors.
;    
;       /basic: If the keyword /basic is given, a version of the interface
;       that displays only the 1D/2D spectrum is launched. 
;         
;
; USAGE NOTES:
;       Please refer to http://specpro.caltech.edu for a brief
;       tutorial. A user manual is provided there as well.
;       Most of the interactive features are self-explanatory. 
;       
; OUTPUTS:
;       The output fields in the lower part of the GUI can be used to 
;       save the redshift, confidence, your initials, and notes in a 
;       formatted ascii file. By default this file is called maskname_zinfo.dat,
;       where maskname is the name of the mask being examined, but the user 
;       can choose a different output location. 
;
;       The output file contains the following fields:
;  
;       mask_name  slit_number  RA  DEC  source_name  redshift  confidence  initials  notes
;  
;       Double clicking any plot window allows the user to save a bitmap (.tif)
;       file of the image.  The button 'Save Spec1D' allows the user to save the 
;       lambda and flux values of the currently displayed 1D plot to
;       an ascii file, for plotting in other programs.
; 
; PROGRAMS USED:
;
;       This program uses a modified version of zfind.pro, originally
;       written by D. Schlegel and modified by Michael C. Cooper.       
;       We have also adapted the DEEP2 pipeline routine extract1d.pro 
;       for use with SpecPro, originally written by
;       Michael Cooper.
    
;       We use screenread.pro, an example program distributed
;       by Liam Gumley with his excellent book "Practical IDL
;       Programming".
;
;       A number of routines written by others are made use of in the
;       automated cross-correlation against spectral templates. These routines
;       are in a separate folder called "external", and a README file there
;       contains references for each routine. The "external" directory
;       should be added to the top of the IDL_PATH environment variable if you wish to
;       use these versions of the routines (which are tested and work
;       with SpecPro). If you have your own versions installed feel
;       free to use them instead.
;
;
; Written by Dan Masters (danmasters1@gmail.com), with contributions
; from Peter Capak. 
;
; Version 3.0
;
; -Modifications to zfindspec.pro to make cross-correlation more
; robust for spectra from different instruments, including WFC3.
;
; -Adding "/space" keyword for spectra taken on space-based
; telescopes. Prevents sky lines from being displayed. 
;
; -Adding "/hist" keyword to allow 1D spectrum to be displayed as a 
;  histogram.
;
; -Adding "/OH" keyword to overlay positions of near-IR OH sky features on 1D
;
; -Changed "Auto-update OFF" button to "Auto-z OFF". The meaning and
; functionality should be more clear.
;
; -Display microns in the infrared.
;
; -Adding buttons "Alt 1D file" and "Alt 2D file". The purpose of
; these is to search for other 1D/2D files with the same slit number
; as currently selected, letting the user quickly jump between them
; (useful, for example, if different orients for a WFC3
; grism spectrum exist, or if there is a calibrated and uncalibrated
; version of the same spectrum).
;
; -Double clicking on 1D spectrum now returns the wavelength (both in
; observed and rest frame)
;
; -Now allows 1-D data to be stored in ascii format
;
; -1D spectra saved to ascii in engineering format

;------------------------------------------------------
;Function to get template spectrum array from fits file
function loadtemplate, filename
   fits_read, filename, array_temp, head 
   crval1=sxpar(head,'CRVAL1')
   cd=sxpar(head,'CDELT1')
   crpix1=sxpar(head,'CRPIX1')
   dim_temp=size(array_temp)
     
   xtemp = findgen(dim_temp(1))
   template=findgen(2,dim_temp(1))
   xtemp=crval1+cd*(xtemp+1-crpix1)
   template[0,*] = xtemp
   template[1,*] = array_temp
   return,template

end

;------------------------------------------------------
function allen_k,wave
  ;used for SED fitting

  Rv=3.1

  ;wavelength vector
  wa=[1000.0,1110.0,1250.0,1430.0,1670.0,2000.0,2220.0,2500.0,2850.0, $
      3330.0,3650.0,4000.0,4400.0,5000.0,5530.0,6700.0,9000.0,10000.0,$ 
      20000.0,100000.0]

  ;k vector
  ka=[4.20,3.70,3.30,3.00,2.70,2.80,2.90,2.30,1.97,1.69, $
      1.58,1.45,1.32,1.13,1.00,0.74,0.46,0.38,0.11,0.00]

  k = interpol(ka,wa,wave,/spline)

  return,k*Rv
end
  
;------------------------------------------------------
function prevot_k,wave
  ;used for SED fitting

  Rv=2.72

;wavelength vector
 wa=[1275,1330,1385,1435,1490,1545,1595,1647,1700,1755,1810,1860,1910,2000,2115, $
     2220,2335,2445,2550,2665,2778,2890,2995,3105,3704,4255,5291,12500,16500, $
     22000,24000,26000,28000,30000,32000,34000,36000,38000,40000]

 ;ee vector
 ee=[13.54,12.52,11.51,10.8,9.84,9.28,9.06,8.49,8.01,7.71,7.17,6.9,6.76,6.38,5.85, $
     5.3,4.53,4.24,3.91,3.49,3.15,3,2.65,2.29,1.81,1,0,-2.02,-2.36,-2.47,-2.51,-2.55, $
     -2.59,-2.63,-2.67,-2.71,-2.75,-2.79,-2.83]


  k=make_array(n_elements(wave),VALUE=0,/DOUBLE)

  short=where(wave LE 34500,scnt)

  if scnt GT 0 then k(short) = interpol(ee,wa,wave(short))+Rv
  
  long=where(wave GT 34500,lcnt)

  if lcnt GT 0 then k(long)=0

  return,k
end
  
;-------------------------------------------------------
function calz_k,wave
  ;used for SED fitting

  ;take input input wave vector in angstroms

  Rv=4.05

  w = wave*1.0e-4; convert to microns
  iwave = 1.0/w ; invert

  k=make_array(n_elements(wave),VALUE=0,/DOUBLE)

  short=where(w LT 0.63,scnt)
  medium = where(w GE 0.63,mcnt)


  if scnt gt 0 then k(short)=(2.659*(- 2.156 + (1.509*iwave(short)) - (0.198*iwave(short)^2) + (0.011*iwave(short)^3)))+Rv
		
  if mcnt gt 0 then k(medium)=(2.659*(-1.857+(1.040*iwave(medium))))+Rv

  bad = where(k lt 0,bcnt)
  if bcnt GT 0 then k(bad)=0

  return,k
end

;------------------------------------------------------
function seaton_k,wave
  ;used for SED fitting
  
   Rv=3.1

   w = wave*1.0e-4; convert to microns
   iwave = 1.0/w ; invert
   
   ee=make_array(n_elements(wave),VALUE=0,/DOUBLE)

   lo=4.595
   gamma=1.051
   c1=-0.38
   c2=0.74
   c3=3.96
   c4=0.26


   short = where(iwave GE 5.9 ,scnt)

   if scnt GT 0 then ee(short) = c1 + c2*iwave(short) + c3/( (iwave(short) - (lo^2)*w(short))^2  + $
                                 gamma^2)+ c4*(0.539*((iwave(short) - 5.9)^2)+0.0564*((iwave(short)-5.9)^3))

 
   medium = where(iwave LT 5.9 and iwave GE 2.74,mcnt)

   if mcnt GT 0 then ee(medium) = c1 + c2*iwave(medium) + c3/( (iwave(medium) - (lo^2)*w(medium))^2  + gamma^2)

   long = where(iwave LT 2.74 ,lcnt)

   if lcnt GT 0 then ee(long) = allen_k(wave(long))-Rv
   

   k=ee+Rv

   return,k
end	

;------------------------------------------------------
function fitzpatrick_k,wave
  ;used for SED fitting
  
   Rv=3.1

   w = wave*1.0e-4; convert to microns
   iwave = 1.0/w ; invert
   
   ee=make_array(n_elements(wave),VALUE=0,/DOUBLE)

   lo=4.608
   gamma=0.994
   c1=-0.69
   c2=0.89
   c3=2.55
   c4=0.50


   short = where(iwave GE 5.9 ,scnt)

   if scnt GT 0 then ee(short) = c1 + c2*iwave(short) + c3/( (iwave(short) - (lo^2)*w(short))^2  + gamma^2)+ $
                                 c4*(0.539*((iwave(short) - 5.9)^2)+0.0564*((iwave(short)-5.9)^3))

 
   medium = where(iwave LT 5.9 and iwave GE 3.0,mcnt)

   if mcnt GT 0 then ee(medium) = c1 + c2*iwave(medium) + c3/( (iwave(medium) - (lo^2)*w(medium))^2  + gamma^2)

   long = where(iwave LT 3.0 ,lcnt)

   if lcnt GT 0 then ee(long) = allen_k(wave(long))-Rv
   

   k=ee+Rv

   return,k
end

;------------------------------------------------------	
function apply_dust,wave,flux,ebv,law
  ;used for SED fitting

 ;apply a dust curve to a flux vector
 ;fist vector should be wavelength,second flux, third ebv value
 ;dust laws are as follows
 ;0 = MW Allen 1976  
 ;1 = MW Seaton 1979 
 ;2 = LMC Fitzpatrick 1986
 ;3 = SMC Prevot et al. 1984
 ;4 = SB Calzetti et al. 2000

 if law EQ 0 then a = allen_k(wave)*ebv
 if law EQ 1 then a = seaton_k(wave)*ebv
 if law EQ 2 then a = FITZPATRICK_K(wave)*ebv
 if law EQ 3 then a = prevot_k(wave)*ebv
 if law EQ 4 then a = calz_k(wave)*ebv

 if ((law LT 0) OR (law GT 4)) then begin
   print,'Invalid dust law, no correction applied!'
   a = make_array(n_elements(wave),VALUE=0,/DOUBLE)
 endif

 return,(10^(-0.4*a)*flux)
end

;------------------------------------------------------
function dust,wave,law
 ;used for SED fitting

 ;apply a dust curve to a flux vector
 ;fist vector should be wavelength,second flux, third ebv value
 ;dust laws are as follows
 ;0 = MW Allen 1976  
 ;1 = MW Seaton 1979 
 ;2 = LMC Fitzpatrick 1986
 ;3 = SMC Prevot et al. 1984
 ;4 = SB Calzetti et al. 2000

 if law EQ 0 then a = allen_k(wave)
 if law EQ 1 then a = seaton_k(wave)
 if law EQ 2 then a = FITZPATRICK_K(wave)
 if law EQ 3 then a = prevot_k(wave)
 if law EQ 4 then a = calz_k(wave)

 if ((law LT 0) OR (law GT 4)) then begin
   print,'Invalid dust law, no correction applied!'
   a = make_array(n_elements(wave),VALUE=0,/DOUBLE)
 endif

 return,a
end

;------------------------------------------------------
FUNCTION madau_teff1,x,z   ;PART OF EXTINCTION DUE TO H-CLOUDS (Madau)
  ;used for SED fitting

  teff1=0.0

  lambda=1216.0
  zlambda=lambda*(1.+z)
  IF (x LT zlambda) THEN teff1=0.0037*(x/lambda)^3.46
     

  lambda=1026.0
  zlambda=lambda*(1.+z)
  IF (x LT zlambda) THEN teff1=teff1+0.00177*(x/lambda)^3.46
     
  lambda=973.0
  zlambda=lambda*(1.+z)
  IF (x LT zlambda) THEN teff1=teff1+0.00106*(x/lambda)^3.46
      
  lambda=950.0
  zlambda=lambda*(1.+z)
  IF (x LT zlambda) THEN teff1=teff1+0.000584*(x/lambda)^3.46
    
  lambda=938.1
  zlambda=lambda*(1.+z)
  IF (x LT zlambda) THEN teff1=teff1+0.00044*(x/lambda)^3.46
     
  lambda=931.0
  zlambda=lambda*(1.+z)
  IF(x LT zlambda) THEN teff1=teff1+0.00040*(x/lambda)^3.46

  lambda=926.5
  zlambda=lambda*(1.+z)
  IF (x LT zlambda) THEN teff1=teff1+0.00037*(x/lambda)^3.46

  lambda=923.4
  zlambda=lambda*(1.+z)
  IF (x LT zlambda) THEN teff1=teff1+0.00035*(x/lambda)^3.46

  lambda=921.2
  zlambda=lambda*(1.+z)
  IF (x LT zlambda) THEN teff1=teff1+0.00033*(x/lambda)^3.46

  lambda=919.6
  zlambda=lambda*(1.+z)
  IF (x LT zlambda) THEN teff1=teff1+0.00032*(x/lambda)^3.46

  lambda=918.4
  zlambda=lambda*(1.+z)
  IF (x LT zlambda) THEN teff1=teff1+0.00031*(x/lambda)^3.46


  return,teff1
end

;------------------------------------------------------
FUNCTION madau_teff2,x,z  ;PART OF EXTINCTION DUE TO H-CLOUDS (Madau)
  ;used for SED fitting

  zlambda=912.0*(1.0 + z)
  teff2=0.0
  IF (x LT zlambda) THEN begin
    xc=x/912.0
    xem=1.0 + z
    teff2=0.25*xc^3*(xem^0.46 - xc^0.46)
    teff2=teff2 + 9.4*xc^1.5*(xem^0.18 - xc^0.18) 
    teff2=teff2 - 0.7*xc^3*(xc^(-1.32) - xem^(-1.32))
    teff2=teff2 - 0.023*(xem^1.68 - xc^1.68)
  ENDIF

  IF (teff2 LT 0.0) THEN teff2=0.0

  return, teff2
END

;------------------------------------------------------
FUNCTION madau_teff,x,y
  ;used for SED fitting

  tf= madau_teff1(x,y) + madau_teff2(x,y)
  teff=exp(-tf)

  return,teff
END
     
;------------------------------------------------------
FUNCTION madau,wave,z ;CALCULATES EXTINCTION DUE TO H-CLOUDS
  ;used for SED fitting

  att = make_array(n_elements(wave),/double,value=0)  
  
  for I=0,n_elements(wave)-1 do att(I)=madau_teff(wave(I)*(1.0+z),z) 

  return,att
END

;------------------------------------------------------
FUNCTION apply_madau,wave,flux,z ;CALCULATES EXTINCTION DUE TO H-CLOUDS
  ;used for SED fitting
  att = make_array(n_elements(wave),/double,value=0)  
  
  for I=0,n_elements(wave)-1 do att(I)=madau_teff(wave(I)*(1.0+z),z) 

  return,(att*flux)
END

;------------------------------------------------------
function read_seds,dir,list,lmin,lmax,Nlowres,Nhighres
  ;used for SED fitting

  readcol, dir+list, F='A,A,A',names,types,files,/SILENT

  ;declare the structure of the sed array
  sedstr = {wave:dblarr(Nlowres),flux:dblarr(Nlowres),hrwave:dblarr(Nhighres),hrflux:dblarr(Nhighres),name:'',type:'',file:''}
  outdata= replicate(sedstr,N_elements(names))

  ;make the sed array
  lrwave = make_array(Nlowres,/double,/index)
  dll=(alog10(lmax)-alog10(lmin))/Nlowres
  lrwave=10^(lrwave*dll + alog10(lmin))

  hrwave = make_array(Nhighres,/double,/index)
  dll=(alog10(lmax)-alog10(lmin))/Nhighres
  hrwave=10^(hrwave*dll + alog10(lmin))

  for I=0,n_elements(files)-1 do begin
    readcol, dir+strtrim(files(I),2), F='F,F',swave,sflux,/SILENT
    sflux=sflux*swave*swave
    nrm=median(sflux)
    sflux=sflux/nrm
        
    outdata[I].wave = lrwave
    outdata[I].flux = interpol(sflux,swave,lrwave)
    bad=where(lrwave LT min(swave) or lrwave GT max(swave),bcnt)
    if (bcnt gt 0) then outdata[I].flux(bad)=-1

    outdata[I].hrwave = hrwave
    outdata[I].hrflux = interpol(sflux,swave,hrwave)
    bad=where(hrwave LT min(swave) or hrwave GT max(swave),bcnt)
    if (bcnt gt 0) then outdata[I].hrflux(bad)=-1

    outdata[I].name = names[I]
    outdata[I].type = types[I]
    outdata[I].file = files[I]

  endfor

  return,outdata
end

;------------------------------------------------------
function fit_sed,sed,ebvmax,debv,rlaws,z,lambda,dlambda,phot,err,TYPE=type

  ;setup the photometry
  goodphot = where(phot GT 0,ngood)
  measphot = where(phot GT 0 and phot LT 90,nmeas)
  badphot = where(phot LT 0,nbad)
  limits = where(phot GE 90,nlim)

  flux = make_array(n_elements(phot),/double,value=0)
  fluxerr = make_array(n_elements(phot),/double,value=0)

  ;convert mags to fluxes
  if ngood GT 0 then begin
     flux(goodphot) = 10^(-0.4*(phot(goodphot)-8.9))
     fluxerr(goodphot) = err(goodphot)*flux(goodphot)/1.085736205
  endif

  ;put in limits
  if nlim GT 0 then begin
    flux(limits) = (10^(-0.4*(err(limits)-8.9)))
    fluxerr(limits) = flux(limits)
  endif

  ;set bad mags to negative fluxes
  if nbad GT 0 then begin
    flux(badphot) = -99
    fluxerr(badphot) = -99
  endif
    
 
  ;get wavelength ranges
  wmin=lambda-dlambda/2.0
  wmax=lambda+dlambda/2.0

  ;figure out if there is data longward of the 1216 break
  uvcon=where(wmin(measphot)/(1.0+z) GT 912 ,uvcnt)

  stars = where(sed.type eq 'STAR')
  Nstars = N_elements(stars)

  if keyword_set(TYPE) then begin

    if TYPE eq 'STAR' then begin 
       z=0.0
       other = where(sed.type ne 'STAR')
    endif else begin
       if TYPE eq 'NONE' then begin
          other = where(sed.type ne 'STAR')
       endif else begin
          other = where(sed.type eq TYPE)
       endelse
    endelse 

  endif else begin
    other = where(sed.type ne 'STAR')
  endelse

  Nsed=N_elements(other)
  SEDsize=N_elements(sed[0].flux)

  zsed = make_array(Nsed,SEDsize,/double,value=0)
  wave = sed[0].wave

  ;calculate the madau correction
  madau_corr = madau(wave,z)

  ;apply madau correctoion
  for I=0,Nsed-1 do begin
    index=other[I]
    zsed(I,*)=sed[index].flux*madau_corr
  endfor

  Nebv=uint(ebvmax/debv)+1
  Nred=n_elements(rlaws)
  ;make array of obscuration curves
  obs=make_array(Nred,SEDsize,/double,value=0)

  FOR R=0,Nred-1 do obs(R,*)= dust(wave,rlaws[R])
  
  ;find the best fit SED

  chimin=9e99
  bestsed=0
  bestebv=0
  bestobs=0
  bestNF=0

  ;set up the current working SED
  current_sed = make_array(SEDsize,/double,value=0)

  if z GT 0.0 then begin
  
  FOR S=0,Nsed-1 do begin  ; loop over SEDs
    FOR R=0,Nebv-1 do begin ;loop over obscurations
      ;apply redening
      ebv = R*debv

      FOR RT=0,Nred-1 do begin ;loop over obscuration laws
        
        current_sed=zsed(S,*)*10^(-0.4*obs(RT,*)*ebv)

	;get template photometry
        tphot=make_array(n_elements(phot),/double,value=0)

        FOR P=0,n_elements(phot)-1 do begin
          tindex = where(wave*(1.0+z) GE wmin(P) and wave*(1.0+z) LE wmax(P),tcnt)
         
          if tcnt GT 0 then begin
             if (min(current_sed(tindex)) LT 0) then begin
               tphot(P)=-99
             endif else begin
               tphot(P)=mean(current_sed(tindex))
             endelse
          endif else begin
	       tphot(P) = interpol(current_sed,wave*(1.0+z),lambda(P))
               if tphot(P) LT 0 then tphot(P) = -99
          endelse

        ENDFOR


        ;calculate normalization factor   
        nf=0
        nfwt=0

        if uvcnt GT 0 then begin 

          FOR P=0,nmeas-1 do begin 
            if ((tphot(measphot(P)) GT 0) and (wmin(measphot(P))/(1.0+z) GT 1216) ) then begin
	       nf=nf + (flux(measphot(P))/tphot(measphot(P)))/(fluxerr(measphot(P))^2)
	       nfwt = nfwt + 1.0/(fluxerr(measphot(P))^2)
            endif
          ENDFOR           
        endif else begin
          FOR P=0,nmeas-1 do begin 
            if ((tphot(measphot(P)) GT 0) and (wmin(measphot(P))/(1.0+z) GT 912) ) then begin
	       nf=nf + (flux(measphot(P))/tphot(measphot(P)))/(fluxerr(measphot(P))^2)
	       nfwt = nfwt + 1.0/(fluxerr(measphot(P))^2)
            endif
          ENDFOR        
        endelse

        nf=nf/nfwt 

        ;do the fit

        chisq=0;
  
        FOR P=0,ngood-1 do begin
          ;print,lambda(goodphot(P)),flux(goodphot(P)),nf*tphot(goodphot(P))

          if (tphot(goodphot(P)) GT 0) then begin
	     chisq=chisq + (flux(goodphot(P))- nf*tphot(goodphot(P)))^2/(fluxerr(goodphot(P))^2)
          endif
          ;print,' '

        ENDFOR

        if (chisq LT chimin) then begin
          bestsed=other(S)
          bestebv=ebv
          bestobs=rlaws(RT)
          bestNF=nf
          chimin=chisq
        endif

      ENDFOR ;end loop over Obscuration laws
    ENDFOR ;end loop over obscurations
  ENDFOR ;end loop over sed types
   
  endif else begin 

    chistarmin=9e99  

    FOR S=0,Nstars-1 do begin  ; loop over the stars
      index = stars(S)

      current_sed=sed[index].flux
      tphot=make_array(n_elements(phot),/double,value=0)

      FOR P=0,n_elements(phot)-1 do begin
         tindex = where(wave GE wmin(P) and wave LE wmax(P),tcnt)
         
         if tcnt GT 0 then begin
           if (min(current_sed(tindex)) LT 0) then begin
             tphot(P)=-99
           endif else begin
             tphot(P)=mean(current_sed(tindex))
           endelse
         endif else begin
           tphot(P) = interpol(current_sed,wave,lambda(P))
           if tphot(P) LT 0 then tphot(P) = -99
         endelse

       ENDFOR

       ;calculate normalization factor   
       nf=0
       nfwt=0

       FOR P=0,nmeas-1 do begin 

          if (tphot(measphot(P)) GT 0) then begin
	    nf=nf + (flux(measphot(P))/tphot(measphot(P)))/(fluxerr(measphot(P))^2)
	    nfwt = nfwt + 1.0/(fluxerr(measphot(P))^2)
          endif

       ENDFOR           

       nf=nf/nfwt 
       ;do the fit

       chisq=0;
  
       FOR P=0,ngood-1 do begin

          if (tphot(goodphot(P)) GT 0) then begin
	    chisq=chisq + (flux(goodphot(P))- nf*tphot(goodphot(P)))^2/(fluxerr(goodphot(P))^2)
          endif

       ENDFOR

       if (chisq LT chistarmin) then begin
         Sbestsed=index
         Sbestebv=0.0
         Sbestobs=-1
         SbestNF=nf
         chistarmin=chisq
       endif
    ENDFOR ;end loop over star types

      bestsed=Sbestsed
      bestebv=Sbestebv
      bestobs=Sbestobs
      bestNF=SbestNF
      chimin=chistarmin

  endelse
  

  bfsed=apply_madau(sed[bestsed].hrwave,apply_dust(sed[bestsed].hrwave,sed[bestsed].hrflux,bestebv,bestobs),z)*bestNF
  
  outdata = {wave:sed[bestsed].hrwave*(1.0+z),flux:bfsed,z:z,chisq:chimin,ebv:bestebv,OBSlaw:bestobs,normfac:bestNF,name:sed[bestsed].name, type:sed[bestsed].type,file:sed[bestsed].file}

  return,outdata  
end

;------------------------------------------------------
;Function to return cleaned-up, scaled stamp image
function stampimage, flux, ivar 
  replaceval = max(flux)

  ;Check for and replace nan
  idxfluxnan = where(finite(flux) eq 0, count1)
  idxivarnan = where(finite(ivar) eq 0, count2)
  idxivarne0 = where(ivar ne 0, count3)
  idxivargood = where(ivar ne 0 and finite(ivar) ne 0, count4)

  if count1 gt 0 then flux[idxfluxnan] = replaceval

  if count3 eq 0 then begin
     image = flux
  endif else if count4 eq 0 then begin
    flux[*,*] = 1000
    image = flux
  endif else begin
    avg_error = mean(sqrt(1.0/(ivar(idxivargood))))
    image = flux / avg_error
  endelse

  return, image

end 

;------------------------------------------------------
;Function used by auto_find_z to determine appropriate
;range of z to check for a given template
function getzrange, event, info, temp
  specpro_get_state, event, info
  template = loadtemplate(temp)

  minspec = info.xrange[0]
  maxspec = info.xrange[1]
 
  templambda = template[0,*]
  mintemp = min(templambda)
  maxtemp = max(templambda)
 
  if maxtemp ge ((maxspec-minspec)/4. + minspec) then begin
    ;covers at least a quarter of spectrum at z = 0.0
     zmin = 0.0
  endif else begin
     ;force zmin such that template cover at least a half
     zmin = ((maxspec-minspec)/2. + minspec) / maxtemp - 1.0
     ;just set to zero if close
     if zmin le 0.05 then zmin = 0.0
  endelse

  zmax = (maxspec - (maxspec-minspec)/2.) / mintemp - 1

  if zmax lt 0.0 then zmax = 0.0
  return, [zmin, zmax]
end

;------------------------------------------------------
;Function used by zfindspec
function findpix, lambda, lambda0
  diff = abs(lambda - lambda0)
  minval = min(diff, minsub)
  minsub = (minsub - 1) > 0
  return, minsub
end

;------------------------------------------------------
;function to create bitmap image of plots (from Gumley)
FUNCTION SCREENREAD, X0, Y0, NX, NY, DEPTH=DEPTH

;- Check arguments
if (n_elements(x0) eq 0) then x0 = 0
if (n_elements(y0) eq 0) then y0 = 0
if (n_elements(nx) eq 0) then nx = !d.x_vsize - x0
if (n_elements(ny) eq 0) then ny = !d.y_vsize - y0

;- Check for TVRD capable device
tvrd_true = !d.flags and 128
if (tvrd_true eq 0) then message, $
  'TVRD is not supported on this device: ' + !d.name

;- On devices which support windows, check for open window
win_true = !d.flags and 256
if (win_true gt 0) and (!d.window lt 0) then message, $
  'No graphics window are open'

;- Get IDL version number
version = float(!version.release)

;- Get display depth
depth = 8
if (win_true gt 0) then begin
  if (version ge 5.1) then begin
    device, get_visual_depth=depth
  endif else begin
    if (!d.n_colors gt 256) then depth = 24
  endelse
endif

;- Set decomposed color mode on 24-bit displays
if (depth gt 8) then begin
  entry_decomposed = 0
  if (version gt 5.1) then $
    device, get_decomposed=entry_decomposed
  device, decomposed=1
endif

;- Get the contents of the window
if (depth gt 8) then true = 1 else true = 0
image = tvrd(x0, y0, nx, ny, order=0, true=true)

;- Restore decomposed color mode on 24-bit displays
if (depth gt 8) then device, decomposed=entry_decomposed

;- Return result to caller
return, image

END

;------------------------------------------------------
;Function to do histogram clipping (from Gumley)
function imclip, image, PERCENT=PERCENT

  ;Skipping argument checks

  if (n_elements(percent) eq 0) then percent = 2.0

  ;Get image max and min
  min_value = min(image, max=max_value)

  ;Compute histogram 
  nbins = 100
  binsize = float(max_value-min_value) / float(nbins)
  hist = histogram(float(image),binsize=binsize)
  bins = lindgen(nbins+1)*binsize + min_value

  ;Compute normalized cumulative sum
  sum = fltarr(n_elements(hist))
  sum[0]=hist[0]
  for i = 1L,n_elements(hist)-1L do $
     sum[i] = sum[i-1]+hist[i]
  sum=100.0*(sum/float(n_elements(image)))

  ;Find and return range
  range = [min_value, max_value]
  index = where((sum ge percent) and $
  (sum le (100.0-percent)),count)
  if (count ge 2) then $
     range = [bins[index[0]],bins[index[count-1]]]
  return, range

end

;-----------------------------------------------------------------------
pro get_serendip_radec, event, info, serendipra=serendipra, serendipdec=serendipdec
  ;Retrieve info structure
  specpro_get_state, event, info

  newextractpos = info.extractpos ;where serendip was extracted
  if info.missinginfofile ne 1 then begin
     ;Get the info for this source
     fmt = 'A,D'
     readcol,info.infofile,F=fmt,fields,values, /silent
     DECidx = where(fields eq 'DEC' or fields eq 'dec')
     RAidx = where(fields eq 'RA' or fields eq 'ra')
     slitPAidx = where(fields eq 'slitPA' or fields eq 'slitpa')  
     slitpixidx = where(fields eq 'pixscale2D' or fields eq 'pixscale2d')
     oldextractposidx = where(fields eq 'extractpos')
     if DECidx ne -1 and slitPAidx ne -1 then begin 
        slitPA = values[slitpaidx]
        if slitpixidx ne -1 then begin
           pixscale = values[slitpixidx]
        endif else begin ;default to DEIMOS pixel scale
           pixscale = 0.1185 ;arcsec / pixel
        endelse
        oldextractpos = values[oldextractposidx]
        pixel_distance = newextractpos - oldextractpos
        old_declination = values[DECidx] ;declination of main source
        old_ra = values[RAidx]
        ;Calculate the RA / DEC of the serendip
        serendipdec = old_declination + (pixel_distance * pixscale * cos(slitpa*!pi/180)) / 3600.
        serendipra =  old_ra + (pixel_distance * pixscale * sin(slitpa*!pi/180) / cos(old_declination*!pi/180)) / 3600.
     endif else begin
        serendipra=-99
        serendipdec=-99
     endelse
  endif else begin
        serendipra=-99
        serendipdec=-99
  endelse

end ;get_serendip_radec

;------------------------------------------------------
pro overplot_slit, stamps, sz, info
  ;utility for stamp update, if slit overplot is on.
  pixscale = stamps.pixscale

  ;First make a grid of the pixel -> RA, DEC mapping for stamp
  pixRA = (indgen(sz[1])-sz[1]/2.+0.5)*pixscale/3600. + stamps.RA
  pixRA = reverse(pixRA)
  pixDEC = (indgen(sz[2])-sz[2]/2.+0.5)*pixscale/3600. + stamps.DEC
  maxRA = max(pixRA)
  minRA = min(pixRA)
  maxDEC = max(pixDEC)
  minDEC = min(pixDEC)
  slitpa = info.slitPA
  slitwid = info.slitwid
 
  ;get the pixel positions of the corners (these might be outside the frame)
  corner1xpix = 50.5 - (info.corner1_ra - stamps.RA)*3600/pixscale
  corner1ypix = (info.corner1_dec - stamps.DEC)*3600/pixscale + 49.5
  corner2xpix = 50.5 - (info.corner2_ra - stamps.RA)*3600/pixscale
  corner2ypix = (info.corner2_dec - stamps.DEC)*3600/pixscale + 49.5
  corner3xpix = 50.5 - (info.corner3_ra - stamps.RA)*3600/pixscale 
  corner3ypix = (info.corner3_dec - stamps.DEC)*3600/pixscale + 49.5
  corner4xpix = 50.5 - (info.corner4_ra - stamps.RA)*3600/pixscale
  corner4ypix = (info.corner4_dec - stamps.DEC)*3600/pixscale + 49.5
  
  ;get x,y pixel positions for corner 1
  ;first check if it falls outside of the frame
  if info.corner1_ra lt minRA $
     or info.corner1_ra gt maxRA $
     or info.corner1_dec lt minDEC $
     or info.corner1_dec gt maxDEC then begin

     ;Here we treat the case of a slit extending outside of 
     ;frame.
     
     ;get the slope of the line connecting it to its partner, which 
     ;is corner 4 in this case
     slope14 = (corner4ypix - corner1ypix)/(corner4xpix - corner1xpix)
     ;solve for where this line intersects the frame edges
     lefty = slope14*(0-corner1xpix)+corner1ypix
     righty = slope14*(99-corner1xpix)+corner1ypix
     bottomx = (0-corner1ypix)/slope14 + corner1xpix
     topx = (99-corner1ypix)/slope14 + corner1xpix
     intersectlist = [[0, lefty],$
                      [99, righty],$
                      [bottomx, 0],$
                      [topx, 99]]

     ;remove bad answers
     good = [1, 1, 1, 1]
     for j = 0,n_elements(intersectlist[0,*])-1 do begin
        if max(intersectlist[*,j]) gt 99 or min(intersectlist[*,j]) lt 0 then begin
           good[j] = 0
        endif
     endfor

     intlist = intersectlist[*,where(good eq 1)]
     
     ;find minimum separation
     minsep = (intlist[0,0]-corner1xpix)^2+(intlist[1,0]-corner1ypix)^2
     minidx = 0
     for jj = 0,n_elements(intlist[0,*])-1 do begin
        sep = (intlist[0,jj]-corner1xpix)^2+(intlist[1,jj]-corner1ypix)^2
        if sep lt minsep then minidx = jj
     endfor
        
     corner1_xpix = round(intlist[0,minidx])
     corner1_ypix = round(intlist[1,minidx])
     if corner1_ypix lt 0 then begin 
        corner1_ypix = 0
     endif else if corner1_ypix gt 99 then begin
        corner1_ypix = 99
     endif
     if corner1_xpix lt 0 then begin 
        corner1_xpix = 0
     endif else if corner1_xpix gt 99 then begin
        corner1_xpix = 99
     endif

     corner1out = 1

  endif else begin
     diff1_x = abs(pixRA - info.corner1_ra) 
     diff1_y = abs(pixDEC - info.corner1_dec)
     corner1_xpix = where(diff1_x eq min(diff1_x))
     corner1_ypix = where(diff1_y eq min(diff1_y))
  
     corner1_xpix = corner1_xpix[0]
     corner1_ypix = corner1_ypix[0]
     
     corner1out = 0
  
  endelse

  ;get x,y pixel positions for corner 2
  ;first check if it falls outside of the frame
  if info.corner2_ra lt minRA $
     or info.corner2_ra gt maxRA $
     or info.corner2_dec lt minDEC $
     or info.corner2_dec gt maxDEC then begin

     ;Here we treat the case of a slit extending outside of 
     ;frame.
  
     ;get the slope of the line connecting it to its partner, which 
     ;is corner 3 in this case
     slope23 = (corner3ypix - corner2ypix)/(corner3xpix - corner2xpix)
     ;solve for where this line intersects the frame edges
     lefty = slope23*(0-corner2xpix)+corner2ypix
     righty = slope23*(99-corner2xpix)+corner2ypix
     bottomx = (0-corner2ypix)/slope23 + corner2xpix
     topx = (99-corner2ypix)/slope23 + corner2xpix
     intersectlist = [[0, lefty],$
                      [99, righty],$
                      [bottomx, 0],$
                      [topx, 99]]

     ;remove bad answers
     good = [1, 1, 1, 1]
     for j = 0,n_elements(intersectlist[0,*])-1 do begin
        if max(intersectlist[*,j]) gt 99 or min(intersectlist[*,j]) lt 0 then begin
           good[j] = 0
        endif
     endfor

     intlist = intersectlist[*,where(good eq 1)]
     
     ;find minimum separation
     minsep = (intlist[0,0]-corner2xpix)^2+(intlist[1,0]-corner2ypix)^2
     minidx = 0
     for jj = 0,n_elements(intlist[0,*])-1 do begin
        sep = (intlist[0,jj]-corner2xpix)^2+(intlist[1,jj]-corner2ypix)^2
        if sep lt minsep then minidx = jj
     endfor
        
     corner2_xpix = round(intlist[0,minidx])
     corner2_ypix = round(intlist[1,minidx])
     if corner2_ypix lt 0 then begin 
        corner2_ypix = 0
     endif else if corner2_ypix gt 99 then begin
        corner2_ypix = 99
     endif
     if corner2_xpix lt 0 then begin 
        corner2_xpix = 0
     endif else if corner2_xpix gt 99 then begin
        corner2_xpix = 99
     endif

     corner2out = 1

  endif else begin
     diff2_x = abs(pixRA - info.corner2_ra) 
     diff2_y = abs(pixDEC - info.corner2_dec)
     corner2_xpix = where(diff2_x eq min(diff2_x))
     corner2_ypix = where(diff2_y eq min(diff2_y))

     corner2_xpix = corner2_xpix[0]
     corner2_ypix = corner2_ypix[0] 

     corner2out = 0

  endelse 

  ;get x,y pixel positions for corner 3
  ;first check if it falls outside of the frame
  if info.corner3_ra lt minRA $
     or info.corner3_ra gt maxRA $
     or info.corner3_dec lt minDEC $
     or info.corner3_dec gt maxDEC then begin

     ;Here we treat the case of a slit extending outside of 
     ;frame.
     
     ;get the slope of the line connecting it to its partner, which 
     ;is corner 2 in this case
     slope32 = (corner2ypix - corner3ypix)/(corner2xpix - corner3xpix)
     ;solve for where this line intersects the frame edges
     lefty = slope32*(0-corner3xpix)+corner3ypix
     righty = slope32*(99-corner3xpix)+corner3ypix
     bottomx = (0-corner3ypix)/slope32 + corner3xpix
     topx = (99-corner3ypix)/slope32 + corner3xpix
     intersectlist = [[0, lefty],$
                      [99, righty],$
                      [bottomx, 0],$
                      [topx, 99]]

     ;remove bad answers
     good = [1, 1, 1, 1]
     for j = 0,n_elements(intersectlist[0,*])-1 do begin
        if max(intersectlist[*,j]) gt 99 or min(intersectlist[*,j]) lt 0 then begin
           good[j] = 0
        endif
     endfor

     intlist = intersectlist[*,where(good eq 1)]
     
     ;find minimum separation
     minsep = (intlist[0,0]-corner3xpix)^2+(intlist[1,0]-corner3ypix)^2
     minidx = 0
     for jj = 0,n_elements(intlist[0,*])-1 do begin
        sep = (intlist[0,jj]-corner3xpix)^2+(intlist[1,jj]-corner3ypix)^2
        if sep lt minsep then minidx = jj
     endfor
       
     corner3_xpix = round(intlist[0,minidx])
     corner3_ypix = round(intlist[1,minidx])
     if corner3_ypix lt 0 then begin 
        corner3_ypix = 0
     endif else if corner3_ypix gt 99 then begin
        corner3_ypix = 99
     endif
     if corner3_xpix lt 0 then begin 
        corner3_xpix = 0
     endif else if corner3_xpix gt 99 then begin
        corner3_xpix = 99
     endif
        
     corner3out = 1

  endif else begin
     diff3_x = abs(pixRA - info.corner3_ra) 
     diff3_y = abs(pixDEC - info.corner3_dec)
     corner3_xpix = where(diff3_x eq min(diff3_x))
     corner3_ypix = where(diff3_y eq min(diff3_y))

     corner3_xpix = corner3_xpix[0]
     corner3_ypix = corner3_ypix[0]

     corner3out = 0  
  endelse

  ;get x,y pixel positions for corner 4
  ;first check if it falls outside of the frame
  if info.corner4_ra lt minRA $
     or info.corner4_ra gt maxRA $
     or info.corner4_dec lt minDEC $
     or info.corner4_dec gt maxDEC then begin

     ;Here we treat the case of a slit extending outside of 
     ;frame.
      
     ;get the slope of the line connecting it to its partner, which 
     ;is corner 1 in this case
     slope14 = (corner4ypix - corner1ypix)/(corner4xpix - corner1xpix)
     ;solve for where this line intersects the frame edges
     lefty = slope14*(0-corner4xpix)+corner4ypix
     righty = slope14*(99-corner4xpix)+corner4ypix
     bottomx = (0-corner4ypix)/slope14 + corner4xpix
     topx = (99-corner4ypix)/slope14 + corner4xpix
     intersectlist = [[0, lefty],$
                      [99, righty],$
                      [bottomx, 0],$
                      [topx, 99]]

     ;remove bad answers
     good = [1, 1, 1, 1]
     for j = 0,n_elements(intersectlist[0,*])-1 do begin
        if max(intersectlist[*,j]) gt 99 or min(intersectlist[*,j]) lt 0 then begin
           good[j] = 0
        endif
     endfor

     intlist = intersectlist[*,where(good eq 1)]
     
     ;find minimum separation
     minsep = (intlist[0,0]-corner4xpix)^2+(intlist[1,0]-corner4ypix)^2
     minidx = 0
     for jj = 0,n_elements(intlist[0,*])-1 do begin
        sep = (intlist[0,jj]-corner4xpix)^2+(intlist[1,jj]-corner4ypix)^2
        if sep lt minsep then minidx = jj
     endfor
        
     corner4_xpix = round(intlist[0,minidx])
     corner4_ypix = round(intlist[1,minidx])
     if corner4_ypix lt 0 then begin 
        corner4_ypix = 0
     endif else if corner4_ypix gt 99 then begin
        corner4_ypix = 99
     endif
     if corner4_xpix lt 0 then begin 
        corner4_xpix = 0
     endif else if corner4_xpix gt 99 then begin
        corner4_xpix = 99
     endif
       
     corner4out = 1

  endif else begin
     diff1_x = abs(pixRA - info.corner4_ra) 
     diff1_y = abs(pixDEC - info.corner4_dec)
     corner4_xpix = where(diff1_x eq min(diff1_x))
     corner4_ypix = where(diff1_y eq min(diff1_y))
  
     corner4_xpix = corner4_xpix[0]
     corner4_ypix = corner4_ypix[0]
     
     corner4out = 0
  
  endelse

  ;Connect the corners to plot the slit
  if corner1out eq 0 or corner2out eq 0 then $
  plots, [corner1_xpix,corner2_xpix], [corner1_ypix,corner2_ypix], thick=1.0, linestyle=0, color='00ff00'XL, /device
  plots, [corner2_xpix,corner3_xpix], [corner2_ypix,corner3_ypix], thick=1.0, linestyle=0, color='00ff00'XL, /device
  if corner3out eq 0 or corner4out eq 0 then $
  plots, [corner3_xpix,corner4_xpix], [corner3_ypix,corner4_ypix], thick=1.0, linestyle=0, color='00ff00'XL, /device
  plots, [corner4_xpix,corner1_xpix], [corner4_ypix,corner1_ypix], thick=1.0, linestyle=0, color='00ff00'XL, /device

end

;------------------------------------------------------
pro select_slit, event, info, prevnext = prevnext, init=init
  ;This routine is triggered when the user selects a slit number.
  ;Finds all other needed files and updates plots.

  ;Called on startup if a filenumber is specified by user.
  
  ;Get state information
  if n_elements(init) eq 0 then $
     specpro_get_state, event, info
 
  ;Look at the current output fields.  If they are non-empty, 
  ;save the output to the default output file.
  if n_elements(prevnext) eq 0 and info.outputdatasaved eq 0 then begin
    widget_control, info.zoutput_id, get_value=zval
    widget_control, info.zconfidence_id, get_value=conf
    widget_control, info.initials_id, get_value=initials
    if zval ne '' and conf ne '' and initials ne '' then  $
        save_output_to_file, event, info, /default
  endif

  ;reset outputdatasaved flag
  info.outputdatasaved = 0

  ;reset reextracted flag 
  info.reextracted = 0

  ;reset RAclicked, DECclicked
  info.RAclicked = 0.0
  info.DECclicked = 0.0
  
  ;reset lambdamin, lambdamax
  info.lambdamin = 0
  info.lambdamax = 0 ;these are set in makeplot

  ;Get the current text in the field.
  valid=0
  on_ioerror, bad_entry
  if n_elements(prevnext) eq 0 and n_elements(init) eq 0 then begin
     widget_control, event.id, GET_VALUE=slitno
     ;Get rid of possible blank spaces
     remchar,slitno,' '
     ;Convert this value to an integer. If not possible, print an error.
     slitno_int = fix(slitno)
     ;If a negative value is entered, set to 0
     if slitno_int lt 0 then slitno = 0
     valid=1
     bad_entry: if ~valid then print, 'Slit entered not valid.' 
    ;Update slitnumber field of info
    info.slitnumber = slitno_int
  endif else begin
     slitno_int = info.slitnumber
     slitno = string(slitno_int)
     remchar,slitno,' '
  endelse

  ;Create a string with the slit number, find relevant files using it
  slitlen = strlen(slitno)
  
  case slitlen of
     1: slit = '00'+slitno
     2: slit = '0'+slitno
     3: slit = slitno
     4: slit = slitno
     else: message, 'Invalid slit number entered'
  endcase
  
  ;Find spec1d file. 
  info.missing1dfile=0
  info.manualselect1d=0
  info.manualselect2d=0
  
  spec1dfilematch = file_search('spec1d.*.'+slit+'.*.fits')

  if n_elements(spec1dfilematch) gt 1 then begin
     print, 'Warning: more than one spec1d file found.'
     filters = 'spec1d.*.'+slit+'.*.fits'
     spec1dfile = dialog_pickfile(title='Mulitple 1D files found. Please select one.',FILTER=filters,file=spec1dfilematch[0])
  endif else begin
     spec1dfile = spec1dfilematch[0]
  endelse     
     
  if spec1dfile eq '' then begin 
     ;search for ascii 1D file
     spec1dfilematch = file_search('spec1d.*.'+slit+'.*.*')
     spec1dfile = spec1dfilematch[0]
     if spec1dfile eq '' then begin
        info.missing1dfile=1
        print, 'Missing spec1d file'
     endif 
  endif
 
  info.spec1dfile = spec1dfile
  if spec1dfile ne '' then begin
     if strpos(spec1dfile,'.fits') ne -1 then begin
        sp = mrdfits(info.spec1dfile,1)
        ptr_free, info.spec1Dptr
        info.spec1Dptr = ptr_new(sp)
     endif else begin
        readcol, info.spec1dfile, lambda, flux, error, f='d,d,d'
        ivar = 1 / error^2
        sp = {flux:flux, lambda:lambda, ivar:ivar}
        ptr_free, info.spec1Dptr
        info.spec1Dptr = ptr_new(sp)
     endelse
  endif

  ;Find spec2d file
  info.missing2dfile=0
  spec2dfilematch = file_search('spec2d.*.'+slit+'.*.fits')
  if n_elements(spec2dfilematch) gt 1 then begin
     print, 'Warning: more than one spec1d file found.'
     filters = 'spec2d.*.'+slit+'.*.fits'
     spec2dfile = dialog_pickfile(title='Mulitple 2D files found. Please select one.',FILTER=filters,file=spec2dfilematch[0])
  endif else begin
     spec2dfile = spec2dfilematch[0]
  endelse   

  if spec2dfile eq '' then begin 
     print, 'Missing spec2d file'
     info.missing2dfile = 1
  endif

  info.spec2Dfile = spec2dfile
  if spec2dfile ne '' then begin
       ;Do the mrdfits now, use pointers to reference them later
       spec2d = mrdfits(info.spec2dfile,1)
       ptr_free, info.spec2Dptr
       info.spec2Dptr = ptr_new(spec2d)
  endif

  ;Find stamp file
  info.missingstampfile = 0
  stampfile = file_search('stamps.*.'+slit+'.*.fits')
  
  ;check whether two or more files found
  if n_elements(stampfile) gt 1 then begin
     print, 'Warning: more than one stamp file found. Using first match.'
     stampfile = stampfile[0]
  endif 
  if stampfile eq '' then begin 
     info.missingstampfile=1
     print, 'Missing stamps file'
  endif
  info.stampfile = stampfile
  if stampfile ne '' then begin
    stamps = mrdfits(info.stampfile,1)
    ptr_free, info.stampsptr
    info.stampsptr = ptr_new(stamps) 
  endif  ;have stamps

  ;Find photfile
  info.missingphotfile = 0
  photfile = file_search('phot.*.'+slit+'.*.dat')
  ;check whether two or more files found
  if n_elements(photfile) gt 1 then begin
     print, 'Warning: more than one phot file found. Using first match.'
     photfile = photfile[0]
  endif 
  if photfile eq '' then begin 
     info.missingphotfile=1
     print, 'Missing photometry file'
  endif
  info.photfile = photfile

  ;Find infofile
  missinginfofile = 0
  infofile = file_search('info.*.'+slit+'.*.dat')
  ;check whether two or more files found
  if n_elements(infofile) gt 1 then begin
     print, 'Warning: more than one info file found.'
     filters = 'info.*.'+slit+'.*.dat'
     infofile = dialog_pickfile(title='Mulitple info files found. Please select one.',FILTER=filters,file=infofile[0])
  endif 
  if infofile eq '' then begin 
     info.missinginfofile=1
     print = 'Missing info file'
  endif
  info.infofile = infofile

  ;Get the photometric redshift for this source
  if info.missinginfofile eq 0 then begin  
     fmt = 'A,D'
     readcol,infofile[0],F=fmt,fields,values, /silent
     zpos = where(fields eq 'zphot')
     extractposidx = where(fields eq 'extractpos')
     extractwidthidx = where(fields eq 'extractwidth')
     RApos = where(fields eq 'RA')
     DECpos = where(fields eq 'DEC')
     info.targetRA = float(values(RApos))
     info.targetDEC = float(values(DECpos)) 
     zphot = values(zpos)
    
     ;Use photz as initial redshift guess
     info.redshift = zphot

     ;set the objpos (vertical pixel position of 1d extraction) and
     ;the pixel width of extraction
     if extractposidx ne -1 then begin
        info.extractpos = values(extractposidx)
        info.extractwidth = values(extractwidthidx)
     endif else begin
        info.extractpos = 0
        info.extractwidth = 10
     endelse

     if info.usingbasicversion ne 1 then begin
        ;Get the corners of the slit in RA, DEC, assuming we have this info
        ;If these values are incorrect (or all zeros) it won't hurt anything.
        slitRApos = where(fields eq 'slitRA' or fields eq 'slitra')
        slitDECpos = where(fields eq 'slitDEC' or fields eq 'slitdec')
        slitlenpos = where(fields eq 'slitlen')
        slitwidpos = where(fields eq 'slitwid')
        slitPApos = where(fields eq 'slitPA' or fields eq 'slitpa')        
        slitRA = values[slitRApos]
        slitDEC = values[slitDECpos] 
        slitlen = values[slitlenpos]
        slitwid = values[slitwidpos]
        slitPA = values[slitpapos]
        if slitRApos ne -1 and slitDECpos ne -1 and slitlenpos ne -1 and slitwidpos ne -1 and slitPApos ne -1 then begin
           info.slitPA = slitPA
           info.slitwid = slitwid
        
           ;get slit corners in RA, DEC, assuming the slit
           ;is aligned with RA. Use slitPA to rotate as needed.
           info.corner1_ra =     slitRA+0.5*(slitlen/3600.)
           info.corner2_ra =     slitRA+0.5*(slitlen/3600.)
           info.corner3_ra =     slitRA-0.5*(slitlen/3600.)
           info.corner4_ra =     slitRA-0.5*(slitlen/3600.)  
           info.corner1_dec =    slitDEC-0.5*(slitwid/3600.)
           info.corner2_dec =    slitDEC+0.5*(slitwid/3600.)
           info.corner3_dec =    slitDEC+0.5*(slitwid/3600.)
           info.corner4_dec =    slitDEC-0.5*(slitwid/3600.)
   
           if (slitPA ne 90.0) and (slitPA ne 180.0) then begin
              ;rotate corners
              ;define rotation matrix
              angle = 90.0 - slitPA
              rotmatrix = [[cos(angle*!pi/180.0), -sin(angle*!pi/180.0)], [sin(angle*!pi/180.0), cos(angle*!pi/180.0)]]
              
              corner1vec = [info.corner1_ra-slitRA, info.corner1_dec-slitDEC]       
              corner2vec = [info.corner2_ra-slitRA, info.corner2_dec-slitDEC] 
              corner3vec = [info.corner3_ra-slitRA, info.corner3_dec-slitDEC] 
              corner4vec = [info.corner4_ra-slitRA, info.corner4_dec-slitDEC] 
              
                                ;get rotated corners
              corner1rot = rotmatrix ## corner1vec
              corner2rot = rotmatrix ## corner2vec
              corner3rot = rotmatrix ## corner3vec
              corner4rot = rotmatrix ## corner4vec

              ;Update the corners
              info.corner1_ra =  corner1rot[0] + slitRA
              info.corner2_ra =  corner2rot[0] + slitRA
              info.corner3_ra =  corner3rot[0] + slitRA
              info.corner4_ra =  corner4rot[0] + slitRA
              info.corner1_dec = corner1rot[1] + slitDEC
              info.corner2_dec = corner2rot[1] + slitDEC
              info.corner3_dec = corner3rot[1] + slitDEC
              info.corner4_dec = corner4rot[1] + slitDEC
           endif
        endif else begin
           ;info file does not have all the
           ;necessary info to plot slit
           print, 'Warning: incomplete slit information in info file.'
           info.corner1_ra =  0
           info.corner2_ra =  0
           info.corner3_ra =  0
           info.corner4_ra =  0
           info.corner1_dec = 0
           info.corner2_dec = 0
           info.corner3_dec = 0
           info.corner4_dec = 0
        endelse
     endif

  endif else begin                         ;if missinginfofile eq 0   
     if spec2dfile ne '' then begin
           ;Do the mrdfits now, use pointers to reference them later
           sz = size(spec2d.flux)
           info.extractpos = fix(sz[2]/2.0)
           info.extractwidth = 0
     endif
  endelse

  info.autozflag = 0 ;indicates that auto-z has not yet been performed
  info.corr_results = fltarr(3,6) ;reset the auto-z results to zero
  info.reextractflag = 0  ;this flag goes to one if the user extracts at new slit position
  
  ;Update the info structure. If this is initializing, no event is triggering
  ;this routine, so update info by hand.
  if n_elements(init) eq 0 then begin
     specpro_set_state,event,info
  endif else begin
    ;Update the displayed slit number
     widget_control, info.enter_slit_id, set_value=slitno
     ;Create pointer to info structure
     infoptr = ptr_new(info)
     ;Store point to info structure in top level base
     widget_control, info.tlb_id, set_uvalue=infoptr
     *infoptr = info 
  endelse

  ;Update plots
  spec1d_plot_update,event,info, missing=info.missing1dfile
  spec2d_plot_update,event,info, missing=info.missing2dfile
  text_update,event,info, missing=info.missinginfofile
  if info.usingbasicversion ne 1 then begin
     phot_plot_update,event,info, missing=info.missingphotfile
     stamp_update,event,info, missing=info.missingstampfile
  endif

  ;force the overlaid template to none, to begin
  fakeevent = create_struct(event,'index',0)
  spec_template_update, fakeevent, info
  widget_control, info.spectemp_list_id, set_droplist_select = 0

  ;force the selected auto-z value to none, to begin
  autoz_solution_update,fakeevent,info
  widget_control, info.autoz_list_id, set_droplist_select = 0 

  ;reset the 1d zoom, if necessary
  if info.zoom ne 0 then begin
     reset_zoom, fakeevent, info
  endif
  
  if missinginfofile eq 0 then redshift_update,event,info,init=1

  ;clear output fields
  widget_control, info.zoutput_id, set_value=''
  widget_control, info.zconfidence_id, set_value=''
  widget_control, info.notes_id, set_value=''
  widget_control, info.initials_id, set_value=''

end ;select_slit

;-----------------------------------------------------------------------
pro pick_alt1d_file, event, info
  specpro_get_state, event, info

  slitno_int = info.slitnumber
  slitno = string(slitno_int)
  remchar,slitno,' '
  slitlen = strlen(slitno) 
  case slitlen of
     1: slit = '00'+slitno
     2: slit = '0'+slitno
     3: slit = slitno
     4: slit = slitno
     else: message, 'Invalid slit number entered'
  endcase

  filters = 'spec1d.*.'+slit+'.*.fits'
  spec1dfile = dialog_pickfile(title='Select 1D file:',FILTER=filters)    
     
  if spec1dfile eq '' then begin 
     print, 'No file selected!'
     return
  endif
 
  info.spec1dfile = spec1dfile
  sp = mrdfits(info.spec1dfile,1)
  ptr_free, info.spec1Dptr
  info.spec1Dptr = ptr_new(sp)
  specpro_set_state, event, info
  spec1d_plot_update,event,info

end

;-----------------------------------------------------------------------
pro pick_alt2d_file, event, info
  specpro_get_state, event, info

  slitno_int = info.slitnumber
  slitno = string(slitno_int)
  remchar,slitno,' '
  slitlen = strlen(slitno) 
  case slitlen of
     1: slit = '00'+slitno
     2: slit = '0'+slitno
     3: slit = slitno
     4: slit = slitno
     else: message, 'Invalid slit number entered'
  endcase

  filters = 'spec2d.*.'+slit+'.*.fits'
  spec2dfile = dialog_pickfile(title='Select 2D file:',FILTER=filters)    
     
  if spec2dfile eq '' then begin 
     print, 'No file selected!'
     return
  endif
 
  info.spec2dfile = spec2dfile
  sp = mrdfits(info.spec2dfile,1)
  ptr_free, info.spec2Dptr
  info.spec2Dptr = ptr_new(sp)
  specpro_set_state, event, info
  spec2d_plot_update,event,info

end
;-----------------------------------------------------------------------
pro spec1d_plot_update, event, info, zoomevent=zoomevent, missing=missing, $
                        autoz_update=autoz_update

  ;Retrieve info structure
  specpro_get_state, event, info
  widget_control, info.spec1d_id, get_value=winid
  wset,winid

  ;Error code if spec1d not found
  if n_elements(missing) ne 0 then begin
      if missing eq 1 then begin
         erase
         return
      endif
  endif 

  if n_elements(zoomevent) ne 0 then begin
      makeplot,info.spec1dfile,info.rebin, info.smooth, info.redshift, info.linetemplates, $
      info.zoom, info.x1Dpress, info.x1Drelease, info.y1Dpress, $
      info.y1Drelease, event, info, zoomevent=zoomevent
  endif else if keyword_set(autoz_update) then begin
      makeplot,info.spec1dfile,info.rebin, info.smooth, info.redshift, info.linetemplates, $
      info.zoom, info.x1Dpress, info.x1Drelease, info.y1Dpress, $
      info.y1Drelease, event, info, /autoz_update
  endif else begin
      makeplot,info.spec1dfile,info.rebin, info.smooth, info.redshift, info.linetemplates, $
      info.zoom, info.x1Dpress, info.x1Drelease, info.y1Dpress, $
      info.y1Drelease, event, info
  endelse
 
end

;-----------------------------------------------------------------------
pro spec2d_plot_update, event, info, missing=missing, reextract=reextract, $
                             zoombox=zoombox
  ;Retrieve info structure
  specpro_get_state, event, info
  widget_control, info.spec2d_id, get_value=winid
  wset,winid
  erase

  ;Error code if spec2d file not found
  if n_elements(missing) ne 0 then begin
      if missing eq 1 then return
  endif 

  if n_elements(reextract) ne 0 then begin
    make2Dplot, info.redshift, info.linetemplates, $
                  info.showspecpos, info.extractpos, event, info, reextract=1
  endif else if n_elements(zoombox) ne 0 then begin 
    make2Dplot, info.redshift, info.linetemplates, $
                  info.showspecpos, info.extractpos, event, info, zoombox=1
  endif else begin
    make2Dplot, info.redshift, info.linetemplates, $
                  info.showspecpos, info.extractpos, event, info
  endelse 

end

;-----------------------------------------------------------------------
pro phot_plot_update, event, info, missing=missing
  ;Retrieve info structure
  specpro_get_state, event, info
  ;Get current redshift
  z = info.redshift
  widget_control, info.sed_id, get_value=winid
  wset,winid
  
  ;Error code if phot file not found
  if n_elements(missing) ne 0 then begin
      if missing eq 1 then begin
        erase
        return
      endif
  endif

  if info.sed_fit_flag eq 0 then begin
     photplot,info.photfile,info.redshift, event, info
  endif else begin
     photplot,info.photfile,info.redshift, event, info, sed=info.sed_templates
  endelse
  
end

;-----------------------------------------------------------------------
pro stamp_update, event, info, missing=missing
  ;Retrieve info structure
  specpro_get_state, event, info
  
  if info.usingsmallversion ne 1 and info.usingbasicversion ne 1 then begin
   widget_control, info.u_stamp_id, get_value=winidu
   widget_control, info.b_stamp_id, get_value=winidb
   widget_control, info.v_stamp_id, get_value=winidv
   widget_control, info.r_stamp_id, get_value=winidr
   widget_control, info.i_stamp_id, get_value=winidi
   widget_control, info.z_stamp_id, get_value=winidz
   widget_control, info.j_stamp_id, get_value=winidj
   widget_control, info.g_stamp_id, get_value=winidg
   widget_control, info.h_stamp_id, get_value=winidh
   widget_control, info.k_stamp_id, get_value=winidk
   widget_control, info.optional_stamp_id, get_value=winidoption

   ;Error code if stamp file not found
   if n_elements(missing) ne 0 then begin
      if missing eq 1 then begin
        wset,winidu
        erase
        wset,winidb
        erase
        wset,winidv
        erase
        wset,winidr
        erase
        wset,winidi
        erase
        wset,winidz
        erase
        wset,winidj
        erase    
        wset,winidg
        erase
        wset,winidh
        erase
        wset,winidk
        erase
        wset,winidoption
        erase
       return
     endif
   endif

   stamps = *(info.stampsptr)

   addcontrast = info.stampcontrast * 0.2

  ;the names of the window ids here come from previous hard-coded settings. 
  stampwinids = [winidu, winidb, winidg, winidv, winidr, winidi, winidz, winidj, winidh, winidk, winidoption]
  labelids = ['u_label_id','b_label_id','g_label_id','v_label_id','r_label_id','i_label_id','z_label_id','j_label_id','h_label_id','k_label_id']

  ;initialize a list that will be the option stamps
  optionlist = ''

  for i = 0, n_elements(stamps)-1 do begin
     if i le 9 then begin
        wset, stampwinids[i]
        erase

        ;set the name of the stamp image
        ;this complicated stuff is just formatting
        stampname = strtrim(stamps[i].name)
        namelen = strlen(stampname)
        tempstring = '            '
        templen = strlen(tempstring)
        if namelen lt 12 then begin
           strput,tempstring,stampname,templen/2-namelen/2
           stampname = tempstring
       endif

        stringcommand = 'labid = info.'+labelids[i]
        void = execute(stringcommand)
        widget_control, labid, set_value=stampname

        ;Check that the image is 100x100 good pixels.
        ;if it's smaller, try to make it fit the window
        ;also remove NaN values; replace with zero.
        thisflux = stamps[i].flux
        thisivar = stamps[i].ivar
        idxnan = where(finite(thisflux) eq 0)
        if n_elements(idxnan) gt 1 then begin
           thisflux[idxnan] = 0
           thisflux[idxnan] = max(thisflux)
        endif else if idxnan ne -1 then begin
           thisflux[idxnan] = 0
           thisflux[idxnan] = max(thisflux)
        endif

        idxgood = where(thisflux gt -1e8)
        if n_elements(idxgood) eq 1e4 then begin
           ;have a full 100x100 pixel stamp
           tv,bytscl(stampimage(thisflux, thisivar), min=-10+addcontrast,max=10)
           rescale = 0
        endif else begin
           ;have something smaller; make sure it's square, if not return
           side = sqrt(n_elements(idxgood))
           if side eq floor(side) then begin
              fl = thisflux[0:fix(side)-1,0:fix(side)-1]
              iv = thisivar[0:fix(side)-1,0:fix(side)-1]
              image = stampimage(fl, iv)
              fitimage = congrid(image, 100, 100)
              tv,bytscl(fitimage, min=-10+addcontrast,max=10)
              rescale = 1
           endif else begin
              continue
           endelse
        endelse 

        if info.drawslit ne 0 then begin
           sz = size(stamps[i].flux)
           overplot_slit, stamps[i], sz, info
           ;show pixel scale as well
           if rescale eq 0 then begin
              arcsec = 100 * stamps[i].pixscale
              outval = string(fix(arcsec))
              outval = strcompress(outval, /remove_all)
              plots, [0,100], [90,90], thick=1.0, linestyle=0, color='ffff00'XL, /device
              xyouts,45,92, outval+'"', charsize = 0.8, color='ffff00'XL, /device
           endif else begin
              arcsec = side * stamps[i].pixscale
              outval = string(fix(arcsec))
              outval = strcompress(outval, /remove_all)
              plots, [0,100], [90,90], thick=1.0, linestyle=0, color='ffff00'XL, /device
              xyouts,45,92, outval+'"', charsize = 0.8, color='ffff00'XL, /device
           endelse
        endif    

     endif else begin
        ;we are now grabbing the overflow stamp images and collecting them
        info.optionstamplist[i-10] = strtrim(stamps[i].name) 
     endelse
  
  endfor

  if n_elements(stamps) ge 11 then begin
     widget_control, info.stampoptionlist_id, set_value = info.optionstamplist
     
     idx = where(info.optionstamplist eq info.optionstamp)
     if idx eq -1 then begin
        info.optionstamp = info.optionstamplist[0]
        widget_control, info.stampoptionlist_id, set_list_select = 0     
     endif else begin
        widget_control, info.stampoptionlist_id, set_list_select = idx
     endelse
  endif

  specpro_set_state, event, info  

  ;update the option stamp window with the appropriate image, if we
  ;have more than 10 stamps
  if n_elements(stamps) ge 11 then begin
     widget_control, info.optional_label, set_value='Other: '+info.optionstamp
     option_stamp_update, event, info
  endif


 ;--------------  
 ;if using small version

 endif else if info.usingbasicversion ne 1 then begin 
   widget_control, info.u_stamp_id, get_value=winidu
   widget_control, info.b_stamp_id, get_value=winidb
   widget_control, info.v_stamp_id, get_value=winidv
   widget_control, info.r_stamp_id, get_value=winidr
   widget_control, info.i_stamp_id, get_value=winidi
   widget_control, info.optional_stamp_id, get_value=winidoption

   ;Error code if stamp file not found
   if n_elements(missing) ne 0 then begin
      if missing eq 1 then begin
        wset,winidu
        erase
        wset,winidb
        erase
        wset,winidv
        erase
        wset,winidr
        erase
        wset,winidi
        erase
        wset,winidoption
        erase
       return
     endif
   endif

   stamps = *(info.stampsptr)

   addcontrast = info.stampcontrast * 0.2

  ;the names of the window ids here come from previous hard-coded settings. 
  stampwinids = [winidu, winidb, winidv, winidr, winidi, winidoption]
  labelids = ['u_label_id','b_label_id','v_label_id','r_label_id','i_label_id']

  for i = 0, n_elements(stamps)-1 do begin
     if i le 4 then begin
        wset, stampwinids[i]
        erase

        ;set the name of the stamp image
        ;this complicated stuff is just formatting
        stampname = strtrim(stamps[i].name)
        namelen = strlen(stampname)
        tempstring = '            '
        templen = strlen(tempstring)
        if namelen lt 12 then begin
           strput,tempstring,stampname,templen/2-namelen/2
           stampname = tempstring
        endif

        stringcommand = 'labid = info.'+labelids[i]
        void = execute(stringcommand)
        widget_control, labid, set_value=stampname

        ;Check that the image is 100x100 good pixels.
        ;if it's smaller, try to make it fit the window
        thisflux = stamps[i].flux
        thisivar = stamps[i].ivar
        idxnan = where(finite(thisflux) eq 0)
        if n_elements(idxnan) gt 1 then begin
           thisflux[idxnan] = 0
           thisflux[idxnan] = max(thisflux)
        endif else if idxnan ne -1 then begin
           thisflux[idxnan] = 0
           thisflux[idxnan] = max(thisflux)
        endif
        
        idxgood = where(thisflux gt -1e6)
        if n_elements(idxgood) eq 1e4 then begin
           ;have a full 100x100 pixel stamp
           tv,bytscl(stampimage(thisflux, thisivar), min=-10+addcontrast,max=10)
           rescale = 0
        endif else begin
           ;have something smaller; make sure it's square, if not return
           side = sqrt(n_elements(idxgood))
           if side eq floor(side) then begin
              fl = thisflux[0:fix(side)-1,0:fix(side)-1]
              iv = thisivar[0:fix(side)-1,0:fix(side)-1]
              image = stampimage(fl, iv)
              fitimage = congrid(image, 100, 100)
              tv,bytscl(fitimage, min=-10+addcontrast,max=10)
              rescale = 1
           endif else begin
              return
           endelse
        endelse

        if info.drawslit ne 0 then begin
           sz = size(stamps[i].flux)
           overplot_slit, stamps[i], sz, info
           ;show pixel scale as well
           if rescale eq 0 then begin
              arcsec = 100 * stamps[i].pixscale
              outval = string(fix(arcsec))
              outval = strcompress(outval, /remove_all)
              plots, [0,100], [90,90], thick=1.0, linestyle=0, color='ffff00'XL, /device
              xyouts,45,92, outval+'"', charsize = 0.8, color='ffff00'XL, /device
           endif else begin
              arcsec = side * stamps[i].pixscale
              outval = string(fix(arcsec))
              outval = strcompress(outval, /remove_all)
              plots, [0,100], [90,90], thick=1.0, linestyle=0, color='ffff00'XL, /device
              xyouts,45,92, outval+'"', charsize = 0.8, color='ffff00'XL, /device
           endelse
        endif     

     endif else begin
        ;we are now grabbing the overflow stamp images and collecting them
        info.optionstamplist[i-5] = strtrim(stamps[i].name)         
     endelse

  endfor

  if n_elements(stamps) ge 6 then begin
     widget_control, info.stampoptionlist_id, set_value = info.optionstamplist
  
     idx = where(info.optionstamplist eq info.optionstamp)
     if idx eq -1 then begin
        info.optionstamp = info.optionstamplist[0]
        widget_control, info.stampoptionlist_id, set_list_select = 0     
     endif else begin
        widget_control, info.stampoptionlist_id, set_list_select = idx
     endelse
  endif

  specpro_set_state, event, info  

  ;update the option stamp window with the appropriate image
  if n_elements(stamps) ge 6 then begin
     widget_control, info.optional_label, set_value='Other: '+strtrim(info.optionstamp)
     option_stamp_update, event, info
  endif

endif ;using small version


end

;-----------------------------------------------------------------------
pro text_update, event, info, missing=missing
  ;Retrieve info structure
  specpro_get_state, event, info
  
  ;Error code if info file is missing
  if n_elements(missing) ne 0 then begin
     if missing eq 1 then begin
       widget_control, info.target_name_id,    set_value='Source Name:  '
       widget_control, info.target_id,         set_value='Source ID:    '
       widget_control, info.target_ra_id,      set_value='RA:           '
       widget_control, info.target_dec_id,     set_value='DEC:          '
       widget_control, info.target_zphot_id,   set_value ='zphot:       '
       widget_control, info.target_zpdf_id,    set_value='zpdf:         '
       widget_control, info.target_zpdf_up_id, set_value='zpdf-up:      '
       widget_control, info.target_zpdf_low_id,set_value='zpdf-low:     '
       return
     endif
  endif

  ;Get the info for this source
  fmt = 'A,D'
  readcol,info.infofile,F=fmt,fields,values, /silent
  zpos = where(fields eq 'zphot')
  zqsopos = where(fields eq 'zqso')
  zpdfpos = where(fields eq 'zpdf')
  zpdfuppos = where(fields eq 'zpdf-up')
  zpdflowpos = where(fields eq 'zpdf-low')
  idpos = where(fields eq 'ID')
  RApos = where(fields eq 'RA')
  DECpos = where(fields eq 'DEC')
  
  if zpos ne -1 then begin 
     zphot = string(values(zpos))
     remchar,zphot,' '
  endif else begin
     zphot = ''
  endelse
  if zpdfpos ne -1 then begin
     zpdf = string(values(zpdfpos))
     remchar,zpdf,' '
  endif else begin
     zpdf = ''
  endelse
  if zpdfuppos ne -1 then begin
     zpdfup = string(values(zpdfuppos))
     remchar,zpdfup,' '
  endif else begin
     zpdfup = ''
  endelse
  if zpdflowpos ne -1 then begin
     zpdflow = string(values(zpdflowpos))
     remchar,zpdflow,' '
  endif else begin
     zpdflow = ''
  endelse
  if RApos ne -1 then begin
     ra = string(double(values(RApos)))
     remchar, ra, ' '
  endif else begin
     ra = ''
  endelse
  if DECpos ne -1 then begin
     dec = string(double(values(DECpos)))
     remchar, dec, ' '  
  endif else begin
     dec = ''
  endelse
  ;getting the id from the infofile can be somewhat tricky
  if idpos ne -1 then begin
     id = double(values(idpos))
     id = string(id)
     temp = strsplit(id, '.', /extract)
     id = temp[0]
     id = strcompress(id,  /remove_all)
  endif else begin
     id = ''
  endelse
  ;get the name from the spec1d file name
  name = info.spec1dfile
  if strlen(name) gt 1 then begin
     temp = strsplit(name, '.', /extract)
     name = temp[n_elements(temp)-2]
     name = strcompress(name, /remove_all)
  endif else begin
     name = '???'
  endelse
  
  widget_control, info.target_name_id,     set_value='Source Name:     '+name[0]
  widget_control, info.target_id,          set_value='Source ID:       '+id[0]
  widget_control, info.target_ra_id,       set_value='RA:              '+ra[0]
  widget_control, info.target_dec_id,      set_value='DEC:             '+dec[0]
  widget_control, info.target_zphot_id,    set_value ='zphot:           '+zphot[0]
  widget_control, info.target_zpdf_id,     set_value='zpdf:            '+zpdf[0]
  widget_control, info.target_zpdf_up_id,  set_value='zpdf-up:         '+zpdfup[0]
  widget_control, info.target_zpdf_low_id, set_value='zpdf-low:        '+zpdflow[0]


end

;-----------------------------------------------------------------------
pro redshift_update, event, info, init=init
  ;Get the current text in the redshift field
  if n_elements(init) eq 0 then begin
    valid=0
    on_ioerror, bad_entry
    widget_control, event.id, GET_VALUE=ztest
    ;Convert this value to a float. If not possible, print an error.
    redshift = float(ztest)
    valid=1
    bad_entry: if ~valid then print, 'Redshift entered not valid.'
  
    ;Update the plots with the appropriate redshift.
    if valid eq 1 then begin 
      ;Update the info structure 
      specpro_get_state, event, info
      info.redshift = redshift
      specpro_set_state, event, info

      ;Update sed plot
      if info.usingbasicversion ne 1 then begin
         widget_control, info.sed_id, get_value=winid
         wset,winid
         phot_plot_update, event, info, missing = info.missingphotfile
      endif
      ;Update spec1d plot
      spec1d_plot_update,event,info,missing=info.missing1dfile
      ;Update spec2d plot
       spec2d_plot_update, event, info, missing = info.missing2dfile
    endif

  endif else begin
     ;Set the redshift field to the photz estimate, which is used as 
     ;as the initial guess. 
     initredshift = string(info.redshift)
     remchar,initredshift,' '
     ;Update the displayed redshift
     widget_control, info.zin, set_value=initredshift
  endelse

end

;-----------------------------------------------------------------------
pro correlation_plot_update, event, info
  specpro_get_state,event,info
  
  ;if "None" is selected, then just stay put
  if event.index eq 0 then begin
     spec1d_plot_update, event, info,missing=info.missing1dfile
     return
  endif

  ;update the redshift
  idx = info.corr_number ;of top six results
  info.redshift = info.corr_results[0,idx]
  specpro_set_state, event, info
  ;Update the redshift field
  redshiftstring = string(info.redshift)
  remchar,redshiftstring,' '
  ;Update the displayed redshift
  widget_control, info.zin, set_value=redshiftstring

  ;Update plots
  if info.usingbasicversion ne 1 then begin
     widget_control, info.sed_id, get_value=winid
     wset,winid
     phot_plot_update, event, info, missing = info.missingphotfile
  endif
  spec1d_plot_update,event,info, /autoz_update, missing = info.missing1dfile
  spec2d_plot_update, event, info, missing = info.missing2dfile

end

;-----------------------------------------------------------------------
pro line_template_update, event, info
  ;Retrieve info structure
  specpro_get_state, event, info
  templates = info.linetemplates

  ;Get the button name
  widget_control, event.id, get_value=template_name  
  dir_em = getenv('SPVIEW') 

  ;Put in error code in case it comes with a final slash
  len = strlen(dir_em)
  lastslash = strpos(dir_em, '/', /reverse_search)
  if lastslash eq len-1 then begin
     dir_em = strmid(dir_em, 0, len-1)
  endif

  ;probably not the most efficient way of doing this, but it works
  case template_name of
   'Emission': begin
       if event.select eq 1 then begin
          info.linetemplates[0] = dir_em + '/' + 'emlines.dat'
       endif else begin
          info.linetemplates[0] = '0'
       endelse
    end     
   'High-z': begin
       if event.select eq 1 then begin
          info.linetemplates[1] = dir_em + '/' + 'highzlines.dat'
       endif else begin
          info.linetemplates[1] = '0'
       endelse
    end     
   'Seyfert emission':  begin
       if event.select eq 1 then begin
          info.linetemplates[2] = dir_em + '/' + 'seyfertemissionlines.dat'
       endif else begin
          info.linetemplates[2] = '0'
       endelse
    end  
   'QSO absorption':  begin
       if event.select eq 1 then begin
          info.linetemplates[3] = dir_em + '/' + 'qsoabsorptionlines.dat'
       endif else begin
          info.linetemplates[3] = '0'
       endelse
    end     
   'QSO emission':  begin
       if event.select eq 1 then begin
          info.linetemplates[4] = dir_em + '/' + 'qsoemissionlines.dat'        
       endif else begin
          info.linetemplates[4] = '0'
       endelse
     end     
    'QSO emission (subset)':  begin
       if event.select eq 1 then begin
          info.linetemplates[5] = dir_em + '/' + 'prominentqsoemlines.dat'        
       endif else begin
          info.linetemplates[5] = '0'
       endelse
     end    
   'Elliptical':  begin
       if event.select eq 1 then begin
          info.linetemplates[6] = dir_em + '/' + 'ellipticalabsorptionlines.dat'
       endif else begin
          info.linetemplates[6] = '0'
       endelse
    end     
   'Elliptical (subset)':  begin
       if event.select eq 1 then begin
          info.linetemplates[7] = dir_em + '/' + 'prominentellipticalabsorptionlines.dat'
       endif else begin
          info.linetemplates[7] = '0'
       endelse
    end   
     else: message, 'template name not recognized'
  endcase

  ;Update the info structure
  specpro_set_state, event, info
   
  ;redo 1D plot
  spec1d_plot_update,event,info,missing=info.missing1dfile
  ;remake 2D plot
  spec2d_plot_update,event,info, missing = info.missing2dfile

end

;-----------------------------------------------------------------------
pro spec_template_update, event, info
  ;Retrieve info structure
  specpro_get_state, event, info
  
  dir = getenv('SPVIEW')

  info.spectempidx = event.index
  
  case event.index of 
    0: begin
      ;None selected
      if ptr_valid(info.spectemplateptr) eq 0 then return
      ptr_free, info.spectemplateptr
      info.spectemplateptr = ptr_new()
      ;clear any previous auto-z results
      info.corr_results = fltarr(3,6)
      info.autozflag = 0 
      widget_control, info.autoz_list_id, set_value = 'None'
      widget_control, info.autoz_list_id, set_droplist_select = 0 
      specpro_set_state, event, info
      spec1d_plot_update, event, info, missing=info.missing1dfile
      specpro_set_state, event, info
      ;Update SED fit
      catch, err
      if err ne 0 then begin
          catch, /cancel
          print, 'Caught error when accessing sed_fit_flag, skipping phot_plot_update.'
          return
      endif
      if tag_exist(info, 'sed_fit_flag') and info.sed_fit_flag then begin
         phot_plot_update,event,info,missing=info.missingphotfile
      endif
      return
     end
    1: begin
      template = loadtemplate(dir+'/template5_olf.fits')
      ;adjust the flux scale
      template[1,*] = template[1,*]*template[0,*]^2
     end
    2: begin
      template = loadtemplate(dir+'/template0_olf.fits')
      ;adjust the flux scale
      template[1,*] = template[1,*]*template[0,*]^2
     end
    3: begin
      template = loadtemplate(dir+'/template1_olf.fits')
      ;adjust the flux scale
      template[1,*] = template[1,*]*template[0,*]^2
      end
    4: begin
      template = loadtemplate(dir+'/template2_olf.fits')
      ;adjust the flux scale
      template[1,*] = template[1,*]*template[0,*]^2
      end
    5: begin
      template = loadtemplate(dir+'/template3_olf.fits')
      ;adjust the flux scale
      template[1,*] = template[1,*]*template[0,*]^2
      end
    6: begin
      template = loadtemplate(dir+'/template4_olf.fits')
      ;adjust the flux scale
      template[1,*] = template[1,*]*template[0,*]^2
     end
    7: begin
      template = loadtemplate(dir+'/sdss_qso.fits')
      ;adjust the flux scale
      template[1,*] = template[1,*]*template[0,*]^2
     end
    8: begin
      template = loadtemplate(dir+'/gal001vel0.fits')
    end
    9: begin
      template = loadtemplate(dir+'/gal015vel0.fits')
    end
    10: begin
      template = loadtemplate(dir+'/gal025vel0.fits')
    end
    11: begin
      ;LBG shapley
      template = loadtemplate(dir+'/lbg.fits')
    end
    ;new bals
    12: begin
      fmt = 'D, D'
      readcol, dir+'/sdss_lobal.spec', F=fmt, lambda, flux
      template = findgen(2,n_elements(lambda))
      template[0,*] = lambda
      template[1,*] = flux
    end   
    13: begin
      fmt = 'D, D'
      readcol, dir+'/sdss_hibal.spec', F=fmt, lambda, flux
      template = findgen(2,n_elements(lambda))
      template[0,*] = lambda
      template[1,*] = flux
    end 
    14: begin
      fmt = 'D, D'
      readcol, dir+'/uka0v.dat', F=fmt, lambda, flux
      template = findgen(2,n_elements(lambda))
      template[0,*] = lambda
      template[1,*] = flux
    end
    15: begin
      fmt = 'D, D'
      readcol, dir+'/ukf0v.dat', F=fmt, lambda, flux
      template = findgen(2,n_elements(lambda))
      template[0,*] = lambda
      template[1,*] = flux
     end
    16: begin
      fmt = 'D, D'
      readcol, dir+'/ukg0v.dat', F=fmt, lambda, flux
      template = findgen(2,n_elements(lambda))
      template[0,*] = lambda
      template[1,*] = flux
    end
    17: begin
      fmt = 'D, D'
      readcol, dir+'/ukk0v.dat', F=fmt, lambda, flux
      template = findgen(2,n_elements(lambda))
      template[0,*] = lambda
      template[1,*] = flux
    end
    18: begin
      fmt = 'D, D'
      readcol, dir+'/ukm0v.dat', F=fmt, lambda, flux 
      template = findgen(2,n_elements(lambda))
      template[0,*] = lambda
      template[1,*] = flux  
    end
    19: begin
      fmt = 'D, D'
      readcol, dir+'/ukm6v.dat', F=fmt, lambda, flux 
      template = findgen(2,n_elements(lambda))
      template[0,*] = lambda
      template[1,*] = flux  
    end
 endcase

  ;Free memory from previous template pointer
  if ptr_valid(info.spectemplateptr) ne 0 then ptr_free, info.spectemplateptr
  
   info.spectemplateptr = ptr_new(template)
  ;Update info structure
   specpro_set_state, event, info

  if info.auto_update_flag eq 1 then begin
     ;calculate the auto-z results, as long as this is a galaxy template
     if event.index le 13 then begin
        if info.missing1dfile ne 1 then begin
           auto_find_z, event, info
        endif
        ;set droplist to show first solution
        widget_control, info.autoz_list_id, set_droplist_select = 1 
        info.corr_number = 0
        specpro_set_state,event,info
        if info.missing1dfile ne 1 then begin
           correlation_plot_update,event,info
        endif else if info.missing1dfile ne 1 then begin
           spec1d_plot_update,event,info,missing=info.missing1dfile
        endif
        if info.usingbasicversion ne 1 then begin
           phot_plot_update,event,info,missing=info.missingphotfile
        endif
        return
     end
  endif

  if info.usingbasicversion ne 1 then begin
     ;Update SED fit
     if info.sed_fit_flag then begin
        phot_plot_update,event,info,missing=info.missingphotfile
     endif
  endif

  if info.missing1dfile ne 1 then begin
     spec1d_plot_update,event,info,missing=info.missing1dfile
  endif

end

;-----------------------------------------------------------------------
pro show_sky_update, event, info
  ;Retrieve info structure
  specpro_get_state, event, info
  if event.select eq 1 then info.showsky = 1
  if event.select eq 0 then info.showsky = 0

  ;Update the info structure
  specpro_set_state, event, info
   
  ;redo 1D plot
  spec1d_plot_update,event,info,missing=info.missing1dfile

end

;-----------------------------------------------------------------------
pro option_stamp_update, event, info
  ;Retrieve info structure
  specpro_get_state, event, info

  ;Get extension of optional stamp
  option = info.optionstamp

  ;stamps = mrdfits(info.stampfile,1)
  stamps = *(info.stampsptr)
  idx = where(strtrim(stamps.name) eq option)

  widget_control, info.optional_stamp_id, get_value=winidoption
  wset,winidoption
  erase
  
  addcontrast = info.stampcontrast * 0.2
  stopdraw = 0
  thisflux = stamps[idx].flux
  thisivar = stamps[idx].ivar
  idxnan = where(finite(thisflux) eq 0)
  if n_elements(idxnan) gt 1 then begin
     thisflux[idxnan] = 0
     thisflux[idxnan] = max(thisflux)
  endif else if idxnan ne -1 then begin
     thisflux[idxnan] = 0
     thisflux[idxnan] = max(thisflux)
  endif
 
  idxgood = where(thisflux gt -1e8)
  if n_elements(idxgood) eq 1e4 then begin
     ;have a full 100x100 pixel stamp
     tv,bytscl(stampimage(thisflux, thisivar), min=-10+addcontrast,max=10)
     rescale = 0
  endif else begin
     ;have something smaller; make sure it's square, if not return
     side = sqrt(n_elements(idxgood))
     if side eq floor(side) then begin
        fl = thisflux[0:fix(side)-1,0:fix(side)-1]
        iv = thisivar[0:fix(side)-1,0:fix(side)-1]
        image = stampimage(fl, iv)
        fitimage = congrid(image, 100, 100)
        tv,bytscl(fitimage, min=-10+addcontrast,max=10)
        rescale = 1
     endif else begin
        return
     endelse
  endelse     
 

  if info.drawslit ne 0 then begin
     sz = size(stamps[idx].flux)
     overplot_slit, stamps[idx], sz, info
     ;show pixel scale as well
     if rescale eq 0 then begin
        arcsec = 100 * stamps[idx].pixscale
        outval = string(fix(arcsec))
        outval = strcompress(outval, /remove_all)
        plots, [0,100], [90,90], thick=1.0, linestyle=0, color='ffff00'XL, /device
        xyouts,45,92, outval+'"', charsize = 0.8, color='ffff00'XL, /device
     endif else begin
        arcsec = side * stamps[idx].pixscale
        outval = string(fix(arcsec))
        outval = strcompress(outval, /remove_all)
        plots, [0,100], [90,90], thick=1.0, linestyle=0, color='ffff00'XL, /device
        xyouts,45,92, outval+'"', charsize = 0.8, color='ffff00'XL, /device
     endelse
  endif   


end ;option_stamp_update

;-----------------------------------------------------------------------
pro rebin_update, event, info
  ;Retrieve info structure
  specpro_get_state, event, info
  ;Get selected rebin factor
  widget_control,event.id,get_value=vals
  rebin = fix(vals(event.index))
  info.rebin = rebin
  specpro_set_state, event, info
  ;update the 1D plot
  spec1d_plot_update,event,info, missing=info.missing1dfile

end ;rebin_update

;-----------------------------------------------------------------------
pro smooth_update, event, info
  ;Retrieve info structure
  specpro_get_state, event, info
  ;Get selected smooth factor
  widget_control,event.id,get_value=vals
  smooth = fix(vals(event.index))
  info.smooth = smooth
  specpro_set_state, event, info
  ;update the 1D plot
  spec1d_plot_update,event,info, missing=info.missing1dfile

end ;smooth_update

;-----------------------------------------------------------------------
pro next_slit, event, info
  specpro_get_state, event, info
 
  if info.outputdatasaved eq 0 then begin
    ;if output fields are filled in, save to default file
    widget_control, info.zoutput_id, get_value=zval
    widget_control, info.zconfidence_id, get_value=conf
    widget_control, info.initials_id, get_value=initials
    if zval ne '' and conf ne '' and initials ne '' then  $
       save_output_to_file, event, info, /default
  endif
  
  ;Set updated slit number
  info.slitnumber = info.slitnumber + 1
  specpro_set_state, event, info

  slitno = string(info.slitnumber)
  remchar,slitno,' '
  ;Update the displayed slit number
   widget_control, info.enter_slit_id, set_value=slitno
  
  ;Call routine to update
  select_slit, event, info, prevnext=1

end ;next_slit

;-----------------------------------------------------------------------
pro previous_slit, event, info
  specpro_get_state, event, info

  if info.outputdatasaved eq 0 then begin
    ;if output fields are filled in, save to default file
    widget_control, info.zoutput_id, get_value=zval
    widget_control, info.zconfidence_id, get_value=conf
    widget_control, info.initials_id, get_value=initials
    if zval ne '' and conf ne '' and initials ne '' then  $
       save_output_to_file, event, info, /default
  endif

  ;Set updated slit number
  info.slitnumber = info.slitnumber - 1
  specpro_set_state, event, info

  ;Update the displayed slit number
   slitno = string(info.slitnumber)
   remchar,slitno,' '
   widget_control, info.enter_slit_id, set_value=slitno

  ;Call routine to update
  select_slit, event, info, prevnext=1

end ;previous_slit

;-----------------------------------------------------------------------
pro specpro_get_state, event, info
  ;Get pointer to info structure
  widget_control, event.top, get_uvalue=infoptr
  info = *infoptr 

end

;-----------------------------------------------------------------------
pro specpro_set_state, event, info
  ;Set pointer to info structure
  widget_control, event.top, get_uvalue=infoptr
  *infoptr = info 
end

;-----------------------------------------------------------------------
pro quit_application, event, info
  specpro_get_state, event, info
  widget_control, event.top, get_uvalue=infoptr
  ;Deallocate pointers
  ptr_free, info.spec2Dptr 
  ptr_free, info.spec1Dptr
  ptr_free, info.stampsptr
  ptr_free, info.spectemplateptr
  ptr_free, infoptr
  ;Close GUI
  widget_control, event.top, /destroy
end

;-----------------------------------------------------------------------
pro set_option_stamp, event, info
  ;Get current info structure
  specpro_get_state, event, info
  idx = event.index
  val = info.optionstamplist[idx]
  ;Reset which optional stamp is viewed
  info.optionstamp = val
  specpro_set_state, event, info
  ;Call the routine to redraw the option stamp
  option_stamp_update, event, info
  ;Change its label
  widget_control, info.optional_label, set_value='Other: '+val
  
end

;-----------------------------------------------------------------------
pro increasez, event, info
  ;Increase z by .0005
  specpro_get_state, event, info
  info.redshift = info.redshift+.0005
  specpro_set_state, event, info
  ;Update the redshift field
  redshiftstring = string(info.redshift)
  remchar,redshiftstring,' '
  ;Update the displayed redshift
  widget_control, info.zin, set_value=redshiftstring
  
  ;Remake plots
  ;Update sed plot
  if info.usingbasicversion ne 1 then begin
     widget_control, info.sed_id, get_value=winid
     wset,winid
     phot_plot_update,event, info, missing = info.missingphotfile
  endif
  ;Update spec1d plot
  spec1d_plot_update,event,info, missing = info.missing1dfile
  ;Update spec2d plot
  spec2D_plot_update,event,info, missing = info.missing2dfile

end
;-----------------------------------------------------------------------
pro increasez_coarse, event, info
  ;Increase z by .05
  specpro_get_state, event, info
  info.redshift = info.redshift+.02
  specpro_set_state, event, info
  ;Update the redshift field
  redshiftstring = string(info.redshift)
  remchar,redshiftstring,' '
  ;Update the displayed redshift
  widget_control, info.zin, set_value=redshiftstring
  
  ;Remake plots
  ;Update sed plot
  if info.usingbasicversion ne 1 then begin
     widget_control, info.sed_id, get_value=winid
     wset,winid
     phot_plot_update,event, info, missing = info.missingphotfile
  endif
  ;Update spec1d plot
  spec1d_plot_update,event,info, missing = info.missing1dfile
  ;Update spec2d plot
  spec2D_plot_update,event,info, missing = info.missing2dfile

end
;-----------------------------------------------------------------------
pro decreasez, event, info
  ;Decrease z by .0005
  specpro_get_state, event, info
  info.redshift = info.redshift-.0005
  specpro_set_state, event, info
  ;Update the redshift field
  redshiftstring = string(info.redshift)
  remchar,redshiftstring,' '
  ;Update the displayed redshift
  widget_control, info.zin, set_value=redshiftstring
  
  ;Remake plots
  ;Update sed plot
  if info.usingbasicversion ne 1 then begin
     widget_control, info.sed_id, get_value=winid
     wset,winid
     phot_plot_update,event, info, missing = info.missingphotfile
  endif
  ;Update spec1d plot
  spec1d_plot_update,event,info, missing = info.missing1dfile
  ;Update spec2d plot
  spec2D_plot_update,event,info, missing = info.missing2dfile

end

;-----------------------------------------------------------------------
pro decreasez_coarse, event, info
  ;Decrease z by .05
  specpro_get_state, event, info
  info.redshift = info.redshift-.02
  specpro_set_state, event, info
  ;Update the redshift field
  redshiftstring = string(info.redshift)
  remchar,redshiftstring,' '
  ;Update the displayed redshift
  widget_control, info.zin, set_value=redshiftstring
  
  ;Remake plots
  ;Update sed plot
  if info.usingbasicversion ne 1 then begin
     widget_control, info.sed_id, get_value=winid
     wset,winid
     phot_plot_update,event, info, missing = info.missingphotfile
  endif
  ;Update spec1d plot
  spec1d_plot_update,event,info, missing = info.missing1dfile
  ;Update spec2d plot
  spec2D_plot_update,event,info, missing = info.missing2dfile

end

;-----------------------------------------------------------------------
pro display_extract_pos, event, info
  ;Controls whether extraction position is displayed on 2D plot and
  ;slit is drawn on stamps
  specpro_get_state, event, info
  info.showspecpos = event.select
  info.drawslit = event.select
  specpro_set_state, event, info

  ;Update 2D plot
  spec2D_plot_update,event,info, missing = info.missing2dfile
 
  ;Update stamp images
  stamp_update, event, info, missing=info.missingstampfile

end

;-----------------------------------------------------------------------
pro zoom_spec1D_event, event, info
  specpro_get_state,event,info
  if info.missing1dfile ne 1 then begin
     ;Detects zoom request via mouse click/drag

     ;a double click returns the wavelength of the current position
     if event.clicks eq 2 then begin
        xrange = info.xrange
        if event.x gt 60 and event.x lt 1080 then begin
           lambdaclicked = ((event.x-60.)/(1080.-60.))*(xrange[1]-xrange[0])+xrange[0]
        endif
        ;print result
        print, 'You clicked on '+strcompress(string(lambdaclicked,F='(D8.2)'),/remove_all) + $
               ' angstroms (rest-frame '+strcompress(string(lambdaclicked/(1+info.redshift),F='(D8.2)'),/remove_all)+' angstroms)'
        return
     endif

     if event.press gt 0 then begin
        ;Update the x, y pixel coordinates of button press
        info.x1Dpress = event.x
        info.y1Dpress = event.y
        specpro_set_state,event,info
     endif else if event.release gt 0 then begin
       ;Update the zoom for 1D plot
        info.x1Drelease = event.x
        info.y1Drelease = event.y
        info.zoom = 1
        specpro_set_state,event,info
        spec1D_plot_update,event,info,zoomevent=1,missing=info.missing1dfile
     endif
  endif

end

;-----------------------------------------------------------------------
pro reset_zoom, event, info
  specpro_get_state,event,info
  info.zoom = 0
  specpro_set_state,event,info
  spec1D_plot_update,event,info,missing=info.missing1dfile
end

;-----------------------------------------------------------------------
pro contrast_2d_update, event, info
  specpro_get_state,event,info
  info.contrast2D = event.value
  specpro_set_state,event,info
  spec2D_plot_update,event,inf, missing = info.missing2dfile
end

;-----------------------------------------------------------------------
pro sigma_update, event, info
  specpro_get_state,event,info
  info.upsigma = -event.index + 10
  specpro_set_state,event,info
  spec2D_plot_update,event,info, missing = info.missing2dfile
end

;-----------------------------------------------------------------------
pro contrast_stamp_update, event, info
  specpro_get_state,event,info
  info.stampcontrast = event.value
  specpro_set_state,event,info
  stamp_update,event,info, missing=info.missingstampfile
end

;-----------------------------------------------------------------------
pro update_template_scale, event, info
  specpro_get_state,event,info
  info.templatescale = event.index + 1.0
  specpro_set_state,event,info
  spec1D_plot_update,event,info,missing=info.missing1dfile
end

;-----------------------------------------------------------------------
pro like_redshift, event, info
  specpro_get_state,event,info
  zstring = string(info.redshift, format = '(F5.3)') 
  widget_control, info.zoutput_id, set_value=zstring
end

;-----------------------------------------------------------------------
pro no_redshift, event, info
  specpro_get_state,event,info 
  widget_control, info.zoutput_id, set_value='-99'
  widget_control, info.zconfidence_id, set_value='0'
end

;-----------------------------------------------------------------------
pro save_image_to_file, event, info
  if event.clicks eq 2 then begin
    specpro_get_state,event,info
    widget_control, event.id, get_value=wid
    wset, wid
    image = screenread(depth=depth)

    imagefile = dialog_pickfile(/WRITE, title='Save as', file='image.tif')
    if imagefile ne '' then begin
      write_tiff, imagefile, reverse(image,3),1
    endif
    
  endif  
end

;-----------------------------------------------------------------------
pro show_serendip_posradec, event, info, serendipra, serendipdec
   ;This routine plots the position of an extracted serendip on the 
   ;stamp images

     if info.usingsmallversion ne 1 then begin
       widget_control, info.u_stamp_id, get_value=winidu
       widget_control, info.b_stamp_id, get_value=winidb
       widget_control, info.v_stamp_id, get_value=winidv
       widget_control, info.r_stamp_id, get_value=winidr
       widget_control, info.i_stamp_id, get_value=winidi
       widget_control, info.z_stamp_id, get_value=winidz
       widget_control, info.j_stamp_id, get_value=winidj
       widget_control, info.g_stamp_id, get_value=winidg
       widget_control, info.h_stamp_id, get_value=winidh
       widget_control, info.k_stamp_id, get_value=winidk
       windowvec = [winidu,winidb,winidv,winidr,winidi,winidz,winidj,winidg,winidh,winidk]
     endif else begin
        widget_control, info.u_stamp_id, get_value=winidu
        widget_control, info.b_stamp_id, get_value=winidb
        widget_control, info.v_stamp_id, get_value=winidv
        widget_control, info.r_stamp_id, get_value=winidr
        widget_control, info.i_stamp_id, get_value=winidi
        windowvec = [winidu,winidb,winidv,winidr,winidi]
     endelse      

    stamps = *(info.stampsptr)
    centerRA = stamps[0].RA
    centerDEC = stamps[0].DEC
    ;Calculate offset
    offset_ra = serendipra - centerRA
    offset_dec = serendipdec - centerDEC
    if info.usingsmallversion ne 1 then begin
       for i = 0, 9 do begin
          stampwid = windowvec[i]
          case stampwid of
             winidu:scale=stamps[0].pixscale
             winidb:scale=stamps[1].pixscale
             winidv:scale=stamps[2].pixscale
             winidr:scale=stamps[3].pixscale
             winidi:scale=stamps[4].pixscale
             winidz:scale=stamps[5].pixscale
             winidj:scale=stamps[6].pixscale
             winidg:scale=stamps[7].pixscale
             winidh:scale=stamps[8].pixscale
             winidk:scale=stamps[9].pixscale
             else:return
          endcase
         ;Plot the position of serendip
         serendipx = 50 - fix(offset_ra*3600./scale)
         serendipy = fix(offset_dec*3600./scale) + 50
         wset,stampwid
         xyouts, serendipx, serendipy, 'o', charsize=0.8, color='ff0000'XL, /device
      endfor

    endif else begin
       stampwid = windowvec[i]
       case stampwid of
          winidu:scale=stamps[0].pixscale
          winidb:scale=stamps[1].pixscale
          winidv:scale=stamps[2].pixscale
          winidr:scale=stamps[3].pixscale
          winidi:scale=stamps[4].pixscale
          else:return
       endcase
         ;Plot the position of serendip
         serendipx = 50 - fix(offset_ra*3600./scale)
         serendipy = fix(offset_dec*3600./scale) + 50
         wset,stampwid
         xyouts, serendipx, serendipy, 'o', charsize=0.8, color='ff0000'XL, /device
    endelse

 end ;show_serendip_posradec

;-----------------------------------------------------------------------
pro stamp_button_event, event, info
   specpro_get_state,event,info

   ;a double click prompts user to save image to file
   if event.clicks eq 2 then begin
      save_image_to_file, event, info
      return
   endif else if event.clicks eq 1 then begin
      
      if info.usingsmallversion ne 1 then begin
         widget_control, info.u_stamp_id, get_value=winidu
         widget_control, info.b_stamp_id, get_value=winidb
         widget_control, info.v_stamp_id, get_value=winidv
         widget_control, info.r_stamp_id, get_value=winidr
         widget_control, info.i_stamp_id, get_value=winidi
         widget_control, info.z_stamp_id, get_value=winidz
         widget_control, info.j_stamp_id, get_value=winidj
         widget_control, info.g_stamp_id, get_value=winidg
         widget_control, info.h_stamp_id, get_value=winidh
         widget_control, info.k_stamp_id, get_value=winidk
       endif else begin
          widget_control, info.u_stamp_id, get_value=winidu
          widget_control, info.b_stamp_id, get_value=winidb
          widget_control, info.v_stamp_id, get_value=winidv
          widget_control, info.r_stamp_id, get_value=winidr
          widget_control, info.i_stamp_id, get_value=winidi
       endelse      

      ;Calculate and display the coordinates
      widget_control, event.id, get_value=stampwid
      wset, stampwid
      if info.usingsmallversion ne 1 then begin
         stamps = *(info.stampsptr)
         case stampwid of
            winidu:scale=stamps[0].pixscale
            winidb:scale=stamps[1].pixscale
            winidv:scale=stamps[2].pixscale
            winidr:scale=stamps[3].pixscale
            winidi:scale=stamps[4].pixscale
            winidz:scale=stamps[5].pixscale
            winidj:scale=stamps[6].pixscale
            winidg:scale=stamps[7].pixscale
            winidh:scale=stamps[8].pixscale
            winidk:scale=stamps[9].pixscale
            else:return
         endcase
      endif else begin
         stamps = *(info.stampsptr)
         case stampwid of
            winidu:scale=stamps[0].pixscale
            winidb:scale=stamps[1].pixscale
            winidv:scale=stamps[2].pixscale
            winidr:scale=stamps[3].pixscale
            winidi:scale=stamps[4].pixscale
            else:return
         endcase
      endelse

      ;Calculate the RA/DEC of clicked position
      centerRA = stamps[0].RA
      centerDEC = stamps[0].DEC
      xoff = event.x - 50
      yoff = event.y - 50
      RAclicked = centerRA - xoff*scale/3600.
      DECclicked = centerDEC + yoff*scale/3600.
      print, 'You clicked: ', RAclicked, DECclicked
      info.RAclicked = RAclicked
      info.DECclicked = DECclicked
      xyouts, event.x-1, event.y-1, '*', charsize=0.8, color='ff0000'XL, /device
      specpro_set_state,event,info
   endif
end

;-----------------------------------------------------------------------
pro spec2d_button_event, event, info
     
   specpro_get_state,event,info

  if info.missing2dfile ne 1 then begin
    ;Detects zoom request via mouse click/drag

    ;a double click prompts user to save image to file
     if event.clicks eq 2 then begin
        save_image_to_file, event, info
        return
     endif
  
     ;A press / release along a horizontal line triggers a re-extraction of 1d spec.
     if event.press gt 0 then begin
     ;Update the x, y pixel coordinates of button press
        info.x2Dpress = event.x
        info.y2Dpress = event.y
        specpro_set_state,event,info
     endif else if event.release gt 0 then begin
       ;Update the zoom for 1D plot
        info.x2Drelease = event.x
        info.y2Drelease = event.y
        if abs(info.y2Drelease - info.y2Dpress) le 4. and info.x2Dpress ne info.x2Drelease then begin
           info.reextractpix = (info.y2Dpress + info.y2Drelease) / 2.0
           info.reextractflag = 1
           specpro_set_state,event,info
           spec2d_plot_update,event,info,reextract=1     
        endif else if info.x2Dpress ne info.x2Drelease then begin 
           ;pop up a plot showing unbinned 2D of box
           specpro_set_state,event,info
           spec2d_plot_update,event,info,zoombox=1
        endif
        
     endif
  endif
    
end
;-----------------------------------------------------------------------
pro create_reextracted_infofile, event, info, serendipra, serendipdec
  ;if 1-D has been reextracted and saved, this routine creates a 
  ;corresponding infofile with the correct extractpos
 
  ;get current info file
  if info.infofile ne '' then begin
     readcol, info.infofile, f='a,d', name, value
     spec1dsp = strsplit(info.spec1dfile,'/',/extract)
     spec1dsplit = strsplit(spec1dsp[n_elements(spec1dsp)-1],'.',/extract)
     newinfofilename = 'info.'+spec1dsplit[1]+'.'+spec1dsplit[2]+'.'+spec1dsplit[3]+'.dat'
     openw, lun,newinfofilename, /get_lun, /append
     idxextractpos = where(name eq 'extractpos')
     value[idxextractpos] = info.extractpos
     idxname = where(name eq 'name' or name eq 'Name')
     if idxname ne -1 then $
        value[idxname] = value[idxname]+'-serendip'
     if serendipra ne -99 then begin
        idxra = where(name eq 'RA' or name eq 'ra')
        value[idxra] = serendipra
        idxdec = where(name eq 'DEC' or name eq 'dec')
        value[idxdec] = serendipdec    
     endif
     for i = 0, n_elements(name)-1 do begin
        printf, lun, f='(a15, d)', name[i], value[i]
     endfor
     close,lun
  endif
    
end

;-----------------------------------------------------------------------
pro phot_button_event, event, info
  ;lets user save .ps of image
  specpro_get_state,event,info

  if info.missingphotfile ne 1 then begin

     ;a double click prompts user to save image to file
     if event.clicks eq 2 then begin
        save_image_to_file, event, info
        return
     endif
  endif
end

;-----------------------------------------------------------------------
pro save_1d_to_ascii, event, info
  ;save current 1d spec (maybe binned, smoothed, and zoomed) to ascii
  specpro_get_state,event,info
  xrange = info.xrange
  yrange = info.yrange
  rebin = info.rebin
  sm = info.smooth

  sp = *(info.spec1Dptr)

  numpts = n_elements(sp.flux)
  
  while (numpts mod rebin) ne 0 do begin
     numpts = numpts + 1
  end

  spec = dblarr(numpts)
  ivar = dblarr(numpts)
  lambda = dblarr(numpts)
  
  spec[0:n_elements(sp.flux)-1] = sp.flux
  ivar[0:n_elements(sp.flux)-1] = sp.ivar
  lambda[0:n_elements(sp.flux)-1] = sp.lambda
 
  ;Get number of pixels in binned image
  binsize = numpts / rebin

  ;Perform the desired binning
  specbin = rebin(spec*ivar,binsize)/rebin(ivar,binsize)
  wavebin = rebin(lambda*ivar,binsize)/rebin(ivar,binsize)
  ivarbin = rebin(ivar,binsize)
 
  ;check the scale factors to make sure smoothing doesn't crash
  if sm GE 3 then begin
    flux = ivarsmooth(specbin,ivarbin,sm)
  endif else begin
    flux = specbin 
  endelse
  
  ;this renaming is relic from older code
  allspec = flux
  allivar = ivarbin
  alllambda = wavebin

  idxfinite = where(finite(alllambda) eq 1)
  alllambda = alllambda[idxfinite]
  allspec = allspec[idxfinite]
  allivar = allivar[idxfinite]

  ;Check in case we switched to plotting microns
  if xrange[1] lt min(alllambda) then begin
     xrange[0] *= 1e4
     xrange[1] *= 1e4
  endif

  idxkeep = where(alllambda gt xrange[0] and alllambda lt xrange[1])
  
  outfile = dialog_pickfile(/WRITE, title='Save spec1D to ascii', file='spec1d.dat')

  if outfile ne '' then begin
    openw, lun, outfile, /get_lun
    for i=0,n_elements(idxkeep)-1 do begin
       printf, lun, alllambda(idxkeep[i]), allspec(idxkeep[i]), format = '(d9.3, 2x, e)'
    endfor
    free_lun, lun
    print,'Spec1d saved to file.'
  endif else begin
    print, 'Action canceled. Data not saved.'
  endelse
   
end

;-----------------------------------------------------------------------
pro save_output_to_file, event, info, default=default
  specpro_get_state,event,info
  widget_control, info.zoutput_id, get_value=z
  widget_control, info.zconfidence_id, get_value=conf
  widget_control, info.notes_id, get_value=notes
  widget_control, info.initials_id, get_value=initials

  z = float(z)
  confidence = fix(conf)
  remchar, initials, ' '

  maskname = info.stampmaskname
  slitnum = info.slitnumber
  slitno = string(slitnum)
  remchar, slitno, ' '
  ;get source name from whichever files exist
  infiles = [info.spec1dfile, info.spec2dfile, info.stampfile, $
               info.photfile]
  idxfiles = where(infiles ne '')
  name = infiles[idxfiles[0]]
  temp = strsplit(name, '.', /extract)
  name = temp[n_elements(temp)-2]
  sourcename = strcompress(name, /remove_all)

  slitlen = strlen(slitno)
  
  case slitlen of
     1: slit = '00'+slitno
     2: slit = '0'+slitno
     3: slit = slitno
     4: slit = slitno
  endcase
  
  if info.outputfilesaved eq 0 then begin ;nothing yet saved by user
     outfile = dialog_pickfile(/WRITE, title='Select file for writing', file=info.outfilename)
     info.outfilename = outfile
     info.outputfilesaved = 1
  endif else begin
     outfile = dialog_pickfile(/WRITE, title='Select file for writing', file=info.outfilename)
  endelse
  
  ;If the user has reextracted the spectrum, and clicked on one of the stamp
  ;images to get the serendip coords, the clicked coords are saved to file.
  if info.reextracted eq 1 and info.RAclicked ne 0.0 then begin
     ra = info.RAclicked
     dec = info.DECclicked
  endif else begin
     ra = info.targetRA
     dec = info.targetDEC  
  endelse

  if outfile ne '' then begin
    openw, lun, outfile, /get_lun, /append
    printf, lun, maskname, slit, ra, dec, sourcename, z, conf, initials, notes, $
            format = '(a, 2x, a4, 1x, d13.8, d13.8, 2x, a15, 1x, d10.4, 1x, f4.1, 1x, a5, 2x, a)'
    free_lun, lun
    print,'File updated'
  endif else begin
    print, 'Action canceled. Data not saved.'
  endelse

  info.outputdatasaved = 1
  specpro_set_state, event, info

end

;-----------------------------------------------------------------------
pro autoz_solution_update, event, info
  ;determines which solution to show
  specpro_get_state,event,info

  if info.autozflag eq 0 and event.index ne 0 then begin ;and info.reextractflag eq 0 then begin
    print, 'First select a galaxy template for autocorrelation!'
    widget_control, info.autoz_list_id, set_droplist_select = 0 
    return
  endif ;else if event.index ne 0 and info.reextractflag eq 1 then begin
    ;print, 'Auto-z not currently working for reextracted 1d spectra.'
    ;widget_control, info.autoz_list_id, set_droplist_select = 0 
    ;return
  ;end
    
  info.corr_number = event.index - 1
  specpro_set_state,event,info

  correlation_plot_update,event,info

end

;-----------------------------------------------------------------------
;pro select_1d_file, event, info
;  ;Select 1D file manually (no naming convention required)
;  specpro_get_state,event,info
;  spec1dfile = dialog_pickfile(title='Select 1D Spectrum to View')
;  
;  ;check that it is in the proper format
;  if strpos(spec1dfile, '.fits') eq -1 then begin
;     print, 'Spectrum must be in fits file!'
;     return
;  endif
;
;  spec1d = mrdfits(spec1dfile, 1)
;  exts = tag_names(spec1d)
;  idxflux = where(exts eq 'flux' or exts eq 'FLUX')
;  idxlambda = where(exts eq 'lambda' or exts eq 'LAMBDA')
;  idxivar = where(exts eq 'ivar' or exts eq 'IVAR')
;  
;  if idxflux eq -1 or idxlambda eq -1 or idxivar eq -1 then begin
;     print, 'Spec1D fits file must contain extensions FLUX, IVAR, and LAMBDA.'
;     return
;  endif

;  ptr_free, info.spec1Dptr
;  info.spec1Dptr = ptr_new(spec1d)
;  info.spec1dfile = spec1dfile
;  info.manualselect1d = 1
;  info.missingphotfile = 1
 
;  specpro_set_state, event, info

;  spec1d_plot_update,event,info

;end

;-----------------------------------------------------------------------
pro auto_update_set, event, info
  ;Set the flag for whether to automatically calculate z when
  ;template is selected
  specpro_get_state,event,info
  info.auto_update_flag = abs(event.select - 1)
  specpro_set_state, event, info
end
;-----------------------------------------------------------------------
pro fit_sed_update, event, info
  ;Set the flag for fitting / not fitting SED
  specpro_get_state,event,info
  info.sed_fit_flag = event.select
  specpro_set_state, event, info
  phot_plot_update,event,info,missing=info.missingphotfile

end
;-----------------------------------------------------------------------
;pro select_2d_file, event, info
;  ;Select 2D file manually (no naming convention required)
;  specpro_get_state,event,info
;  spec2dfile = dialog_pickfile(title='Select 2D Spectrum to View')
;  
;  ;check that it is in the proper format
;  if strpos(spec2dfile, '.fits') eq -1 then begin
;     print, 'Spectrum must be in fits file!'
;     return
;  endif
;
;  spec2d = mrdfits(spec2dfile, 1)
;  exts = tag_names(spec2d)
;  idxflux = where(exts eq 'flux' or exts eq 'FLUX')
;  idxlambda = where(exts eq 'lambda' or exts eq 'LAMBDA')
;  idxivar = where(exts eq 'ivar' or exts eq 'IVAR')
  
;  if idxflux eq -1 or idxlambda eq -1 or idxivar eq -1 then begin
;     print, 'Spec2D fits file must contain extensions FLUX, IVAR, and LAMBDA'
;     return
;  endif

;  ptr_free, info.spec2Dptr
;  info.spec2Dptr = ptr_new(spec2d)
;  info.spec2dfile = spec2dfile
;  info.manualselect2d = 1
;  info.missingphotfile = 1
; 
;  specpro_set_state, event, info
;
;  spec2d_plot_update,event,info

;end

;-----------------------------------------------------------------------
pro recenter_update, event, info
  ;computes the best z in a small window around current redshift (z-0.05, z+0.05)
  specpro_get_state,event,info

  if info.autozflag eq 0 then begin ;and info.reextractflag eq 0 then begin
    print, 'Galaxy template must be selected for autocorrelation!'
    return
  ;endif else if info.reextractflag eq 1 then begin
  ;  print, 'Auto-z not currently working for reextracted 1d spectra.'
  ;  return
  endif

  current_z = info.redshift 
  minz = current_z - 0.05
  maxz = current_z + 0.05
  auto_find_z, event, info, minz=minz, maxz=maxz
  ;set droplist to show first solution
  widget_control, info.autoz_list_id, set_droplist_select = 1 
  info.corr_number = 0
  specpro_set_state,event,info
  fakeevent = create_struct(event,'index',1)
  correlation_plot_update,fakeevent,info    

end

;-----------------------------------------------------------------------
pro auto_find_z, event, info, minz=minz, maxz=maxz
  specpro_get_state,event,info
  
  dir = getenv('SPVIEW')
  
    ;in the following set min / max z based on coverage of template
    case info.spectempidx of 
      0: begin
         ;no template selected
         return
      end
      1: begin
        ;VVDS LBG
        template = dir+'/template5_olf.fits'

      end
      2: begin
        ;VVDS elliptical
        template = dir+'/template0_olf.fits'

      end
      3: begin
        ;VVDS S0
        template = dir+'/template1_olf.fits'

      end
      4: begin
        ;VVDS Early Spiral
        template = dir+'/template2_olf.fits'

      end
      5: begin
        ;VVDS Spiral
        template = dir+'/template3_olf.fits'
    
      end
      6: begin
        ;VVDS Starburst
        template = dir+'/template4_olf.fits'
 
      end
      7: begin
        ;SSDS QSO
        template = dir+'/sdss_qso.fits'
 
       end
      8: begin
        ;Red galaxy
        template = dir+'/gal001vel0.fits'
 
      end
      9: begin
        ;Green galaxy
        template = dir+'/gal015vel0.fits'
 
      end
      10: begin
        ;Blue galaxy
        template = dir+'/gal025vel0.fits'
  
       end
      11: begin
        ;LBG shapley
        template = dir+'/lbg.fits'
  
      end
      12: begin
        ;LoBAL
        template = dir+'/first_lobals.fits'
  
      end
      13: begin
        ;HiBAL
        template = dir+'/first_hibals.fits'
   
      end

   endcase
 

   if n_elements(minz) ne 0 then begin
      ;recenter update, force the z bounds
      zmin = minz
      zmax = maxz
   endif else if n_elements(minz) eq 0 and info.usephotzforautoz eq 1 then begin 
      ;set z range from photz
      fmt = 'A,D'
      readcol,info.infofile,F=fmt,fields,values, /silent
      zpdfuppos = where(fields eq 'zpdf-up')
      zpdflowpos = where(fields eq 'zpdf-low') 
      zpdfup = values(zpdfuppos)
      zpdflow = values(zpdflowpos)
   
      if zpdfup[0] gt zpdflow[0] and abs(zpdfup[0]) lt 7 and abs(zpdflow[0]) lt 7 then begin
        zmin = zpdflow[0]
        zmax = zpdfup[0]
      endif else begin
        print, 'Not computing auto-z due to bad zpdf_up / zpdf_low'
        return
      endelse
   endif else begin
       zrange = getzrange(event, info, template)
       zmin = zrange[0]
       zmax = zrange[1]
   endelse

    pspace = 1;controls how many pixels to move for each correlation
    width = 3*pspace
    nfind = 6 ;number of answers to be found by zcompute 
    result = zfindspec(template, event, info, /linear_lambda, zmin=zmin, zmax=zmax, $
       nfind=nfind, width=width, pspace = pspace)
    
    if size(result,/tn) ne 'STRUCT' then return

    zstring = strarr(nfind+1)
    zstring[0] = 'None'
    
    for i = 0, nfind-1 do begin
      ;update the stored auto-z vals
      outvec = [result[i].z, result[i].z_err, result[i].chi2]
      info.corr_results[0:2,i] = outvec 
      zstring[i+1] = string(result[i].z, format = '(F5.3)')   
    endfor
  
    ;Force the droplist to show the z solutions
    widget_control, info.autoz_list_id, set_value = zstring
    info.autozflag = 1 ;tells routine that auto-z has been calculated
 
    specpro_set_state,event,info

end

;-----------------------------------------------------------------------
pro makeplot, file, rebin, sm, z, linetemplates, zoom, xin1, $
              xin2, yin1, yin2, event, info, zoomevent=zoomevent, $
              autoz_update=autoz_update

;Make a 1-D binned plot of the spectrum.  

if info.manualselect2d eq 1 and info.manualselect1d eq 0 then return

specpro_get_state,event,info
xrange = info.xrange
yrange = info.yrange
sp = *(info.spec1Dptr)

numpts = n_elements(sp.flux)

while (numpts mod rebin) ne 0 do begin
   numpts = numpts + 1
end

spec = dblarr(numpts)
ivar = dblarr(numpts)
lambda = dblarr(numpts)

spec[0:n_elements(sp.flux)-1] = sp.flux
ivar[0:n_elements(sp.flux)-1] = sp.ivar
lambda[0:n_elements(sp.flux)-1] = sp.lambda

;Get number of pixels in binned image
binsize = numpts / rebin

;Perform the desired binning

;first make points where ivar = 0 to be finite but
;small
ivarzeros = where(ivar eq 0)
ivarnonzeros = where(ivar ne 0)
if n_elements(ivarnonzeros) gt 1 then begin
   minnonzero = min(ivar[ivarnonzeros])
endif else begin
   minnonzero = 0
endelse
if n_elements(ivarzeros) gt 1 and minnonzero ne 0 then begin
   ivar[ivarzeros] = minnonzero / 100.
endif

specbin = rebin(spec*ivar,binsize)/rebin(ivar,binsize)
wavebin = rebin(lambda*ivar,binsize)/rebin(ivar,binsize)
ivarbin = rebin(ivar,binsize)

;check the scale factors to make sure smoothing doesn't crash
if sm GE 3 then begin
  flux = ivarsmooth(specbin,ivarbin,sm)
endif else begin
  flux = specbin 
endelse

;this renaming is relic from older code
allspec = flux
allivar = ivarbin
alllambda = wavebin

;set the x scale
minlambda=min(sp.lambda)
maxlambda=max(sp.lambda)

;store min/max values for the phot-z plotting
info.lambdamin = minlambda
info.lambdamax = maxlambda

xsize=maxlambda-minlambda
xmin = minlambda-xsize*.05
xmax = maxlambda+xsize*.05

;shortest line viewable, for this redshift
lmin = xmin/(1.0+z)
;longest line viewable, for this redshift
lmax = xmax/(1.0+z)

;set the scale
;rms=max(allivar)
;minflux=-5*rms;
good=where(allspec gt -1e30 and allspec lt 1e30)
if n_elements(good) eq 1 then begin
  print,'Bad spec1d'
  erase
  return
endif
maxflux=max(allspec(good))
minflux=min(allspec(good))
allspec = allspec(good)
alllambda = alllambda(good) 

;Scale for micron display
if mean(alllambda) ge 1e4 then begin
   sc=1e-4
endif else begin
   sc=1
endelse

;define some plot parameters
ymargin=(maxflux-minflux)*0.05
ymax=maxflux+ymargin*2
ymin=minflux-ymargin*2
ylmax= ymax - (ymax-ymin)*.1
ylmin= ymin + (ymax-ymin)*.1
ylabelpos = ylmax+(ymax-ymin)*.05

xtick = 500

if info.usingsmallversion ne 1 then begin
  ;if this is a zoom update, force min and max values
  if zoom ne 0 and n_elements(zoomevent) ne 0 then begin
     ;the lower corner of plot is (60,50), upper is (1080,250)
     ;this may change; probably want to find way to automatically find this.
     if min([xin1,xin2]) gt 60 then begin
       xmin = ((min([xin1,xin2])-60.)/(1080.-60.))*(xrange[1]-xrange[0])+xrange[0]
     endif else begin 
       xmin = xrange[0]
    endelse
     if max([xin1,xin2]) lt 1080 then begin
       xmax = ((max([xin1,xin2])-60.)/(1080.-60.))*(xrange[1]-xrange[0])+xrange[0]
     endif else begin
       xmax = xrange[1]
    endelse
     if min([yin1,yin2]) gt 50 then begin
       ymin = ((min([yin1,yin2])-50.)/(250.-50.))*(yrange[1]-yrange[0])+yrange[0]
     endif else begin 
       ymin = yrange[0]
     endelse
     if max([yin1,yin2]) lt 250 then begin
       ymax = ((max([yin1,yin2])-50.)/(250.-50.))*(yrange[1]-yrange[0])+yrange[0]
     endif else begin
       ymax = yrange[1]
     endelse
     ;Check to see if we are showing far too much in the y-dimension. If so, scale
     ;appropriately.
     xzoomregion = where(alllambda gt xmin and alllambda lt xmax)
     ymaxzoom = max(allspec[xzoomregion])
     yminzoom = min(allspec[xzoomregion])
     yzoomrange = ymaxzoom - yminzoom
   
     ;if zoom region too big, reset
     if yminzoom - ymin gt yzoomrange/2. and ymax - ymaxzoom gt yzoomrange/2. then begin 
         yzoommargin=yzoomrange*0.05
         ymax=ymaxzoom+yzoommargin*2
         ymin=yminzoom-yzoommargin*2
     endif
          
     ylabelpos = ymax - (ymax-ymin)*0.1
     xtick = round((xmax-xmin)/100.)*10
   endif else if zoom ne 0 and n_elements(zoomevent) eq 0 then begin
     xmin = xrange[0]
     xmax = xrange[1]
     ymin = yrange[0]
     ymax = yrange[1]
     ylabelpos = ymax - (ymax-ymin)*0.1
     xtick = round((xmax-xmin)/100.)*10
   end
endif else begin  ;using small version
     if zoom ne 0 and n_elements(zoomevent) ne 0 then begin
     ;the lower corner of plot is (60,50), upper is (1000,190)
     ;this may change; probably want to find way to automatically find this.
     if min([xin1,xin2]) gt 60 then begin
       xmin = ((min([xin1,xin2])-60.)/(1000.-60.))*(xrange[1]-xrange[0])+xrange[0]
     endif else begin 
       xmin = xrange[0]
    endelse
     if max([xin1,xin2]) lt 1080 then begin
       xmax = ((max([xin1,xin2])-60.)/(1000.-60.))*(xrange[1]-xrange[0])+xrange[0]
     endif else begin
       xmax = xrange[1]
    endelse
     if min([yin1,yin2]) gt 50 then begin
       ymin = ((min([yin1,yin2])-50.)/(190.-50.))*(yrange[1]-yrange[0])+yrange[0]
     endif else begin 
       ymin = yrange[0]
     endelse
     if max([yin1,yin2]) lt 250 then begin
       ymax = ((max([yin1,yin2])-50.)/(190.-50.))*(yrange[1]-yrange[0])+yrange[0]
     endif else begin
       ymax = yrange[1]
     endelse

     ;Check to see if we are showing far too much in the y-dimension. If so, scale
     ;appropriately.
     xzoomregion = where(alllambda*sc gt xmin and alllambda*sc lt xmax)
     ymaxzoom = max(allspec[xzoomregion])
     yminzoom = min(allspec[xzoomregion])
     yzoomrange = ymaxzoom - yminzoom
    
     ;if zoom region too big, reset
     if yminzoom - ymin gt yzoomrange/2. and ymax - ymaxzoom gt yzoomrange/2. then begin 
         yzoommargin=yzoomrange*0.05
         ymax=ymaxzoom+yzoommargin*2
         ymin=yminzoom-yzoommargin*2
     endif

     ylabelpos = ymax - (ymax-ymin)*0.1
     xtick = round((xmax-xmin)/100.)*10
   endif else if zoom ne 0 and n_elements(zoomevent) eq 0 then begin
     xmin = xrange[0]
     xmax = xrange[1]
     ymin = yrange[0]
     ymax = yrange[1]
     ylabelpos = ymax - (ymax-ymin)*0.1
     xtick = round((xmax-xmin)/100.)*10
  end
 endelse 

;Check in case we switched to plotting microns
if xmax lt min(alllambda) then begin
   xmin *= 1e4
   xmax *= 1e4
endif

;ignore clicks with no motion
if xmin eq xmax then return
if ymin eq ymax then return

;plot spectrum
if info.hist eq 1 then begin
   if mean(allspec) ge 1e5 or mean(allspec) le 1e-5  then begin
      if mean(alllambda) lt 1e4 then begin  
         plot,alllambda, allspec, xrange=[xmin,xmax], yrange=[ymin,ymax], thick=1, xstyle=9, $
              XTITLE='Observed Wavelength (Angstroms)', YTITLE='Flux', xtickinterval=xtick, $
              xticklen=0.06, ymargin=[5,5], psym=10, ytickformat='(E9.1)'
      endif else begin
         plot,alllambda/1e4, allspec, xrange=[xmin/1e4,xmax/1e4], yrange=[ymin,ymax], thick=1, xstyle=9, $
              XTITLE='Observed Wavelength (Microns)', YTITLE='Flux', xtickinterval=xtick/1e4, $
              xticklen=0.06, ymargin=[5,5], psym=10, ytickformat='(E9.1)'
      endelse   
   endif else begin
      if mean(alllambda) lt 1e4 then begin  
         plot,alllambda, allspec, xrange=[xmin,xmax], yrange=[ymin,ymax], thick=1, xstyle=9, $
              XTITLE='Observed Wavelength (Angstroms)', YTITLE='Flux', xtickinterval=xtick, $
              xticklen=0.06, ymargin=[5,5], psym=10
      endif else begin
         plot,alllambda/1e4, allspec, xrange=[xmin/1e4,xmax/1e4], yrange=[ymin,ymax], thick=1, xstyle=9, $
              XTITLE='Observed Wavelength (Microns)', YTITLE='Flux', xtickinterval=xtick/1e4, $
              xticklen=0.06, ymargin=[5,5], psym=10
      endelse
   endelse
  
endif else begin ;not plotting as histogram
   if mean(allspec) ge 1e5 or mean(allspec) le 1e-5  then begin 
      if mean(alllambda) lt 1e4 then begin  
         plot,alllambda, allspec, xrange=[xmin,xmax], yrange=[ymin,ymax], thick=1, xstyle=9, $
              XTITLE='Observed Wavelength (Angstroms)', YTITLE='Flux', xtickinterval=xtick, $
              xticklen=0.06, ymargin=[5,5], ytickformat='(E9.1)'
      endif else begin
         plot,alllambda/1e4, allspec, xrange=[xmin/1e4,xmax/1e4], yrange=[ymin,ymax], thick=1, xstyle=9, $
              XTITLE='Observed Wavelength (Microns)', YTITLE='Flux', xtickinterval=xtick/1e4, $
              xticklen=0.06, ymargin=[5,5], ytickformat='(E9.1)'
      endelse 
   endif else begin
      if mean(alllambda) lt 1e4 then begin  
         plot,alllambda, allspec, xrange=[xmin,xmax], yrange=[ymin,ymax], thick=1, xstyle=9, $
              XTITLE='Observed Wavelength (Angstroms)', YTITLE='Flux', xtickinterval=xtick, $
              xticklen=0.06, ymargin=[5,5]
      endif else begin
         plot,alllambda/1e4, allspec, xrange=[xmin/1e4,xmax/1e4], yrange=[ymin,ymax], thick=1, xstyle=9, $
              XTITLE='Observed Wavelength (Microns)', YTITLE='Flux', xtickinterval=xtick/1e4, $
              xticklen=0.06, ymargin=[5,5]
      endelse   
   endelse
endelse

info.xrange = !x.crange
info.yrange = !y.crange

;if xrange in microns, then change to angstroms
if max(info.xrange le 100) then begin
   info.xrange *= 1e4
endif

specpro_set_state,event,info

if info.showsky eq 1 then begin
   print, '[L#3853] oplot,wavebin*sc,1.0/sqrt(ivarbin*rebin*sm)'
   oplot,wavebin*sc,1.0/sqrt(ivarbin*rebin*sm),color='00FC7C'XL  
end

if info.space ne 1 then begin
     border = (ymax-ymin)*.05
     ;atmospheric absorbtion bands
     if xmin lt 7708. and xmax gt 7586 then begin
        abandmin=7586.0*sc
        abandmax=7708.0*sc
        polyfill,[abandmin,abandmax,abandmax,abandmin],[ymin,ymin,ymax-border,ymax-border],/LINE_FILL,spacing=1,color='ff0000'XL
        xyouts, abandmin,ylabelpos, '!6A Band!X',charsize=1
     endif

     if xmin lt 6945. and xmax gt 6864. then begin
        bbandmin=6864.0*sc
        bbandmax=6945.0*sc
        polyfill,[bbandmin,bbandmax,bbandmax,bbandmin],[ymin,ymin,ymax-border,ymax-border],/LINE_FILL,spacing=1,color='ff0000'XL
        xyouts, bbandmin, ylabelpos, '!6B Band!X',charsize=1
     endif

     if xmin lt 9700. and xmax gt 9200. then begin     
        wbandmin=9200.0*sc
        wbandmax=9700.0*sc
        polyfill,[wbandmin,wbandmax,wbandmax,wbandmin],[ymin,ymin,ymax-border,ymax-border],/LINE_FILL,spacing=1,color='ff0000'XL
        xyouts, wbandmin, ylabelpos, '!6Water!X',charsize=1
     endif

     if xmin lt 11600. and xmax gt 11000. then begin     
        wbandmin=11000.0*sc
        wbandmax=11600.0*sc
        polyfill,[wbandmin,wbandmax,wbandmax,wbandmin],[ymin,ymin,ymax-border,ymax-border],/LINE_FILL,spacing=1,color='ff0000'XL
        xyouts, wbandmin, ylabelpos, '!6Water!X',charsize=1
     endif

     if xmin lt 12800. and xmax gt 12500. then begin     
        O2min=12500.0*sc
        O2max=12800.0*sc
        polyfill,[O2min,O2max,O2max,O2min],[ymin,ymin,ymax-border,ymax-border],/LINE_FILL,spacing=1,color='ff0000'XL
        xyouts, O2min, ylabelpos, '!6O2!X',charsize=1
     endif

     if xmin lt 14500. and xmax gt 13400. then begin     
        wbandmin=13400.0*sc
        wbandmax=14500.0*sc
        polyfill,[wbandmin,wbandmax,wbandmax,wbandmin],[ymin,ymin,ymax-border,ymax-border],/LINE_FILL,spacing=1,color='ff0000'XL
        xyouts, wbandmin, ylabelpos, '!6Water!X',charsize=1
     endif

     if xmin lt 15800. and xmax gt 15700. then begin     
        CO2min=15700.0*sc
        CO2max=15800.0*sc
        polyfill,[CO2min,CO2max,CO2max,CO2min],[ymin,ymin,ymax-border,ymax-border],/LINE_FILL,spacing=1,color='ff0000'XL
        xyouts, CO2min, ylabelpos, '!6CO2!X',charsize=1
     endif

     if xmin lt 16200. and xmax gt 16000. then begin     
        CO2min=16000.0*sc
        CO2max=16200.0*sc
        polyfill,[CO2min,CO2max,CO2max,CO2min],[ymin,ymin,ymax-border,ymax-border],/LINE_FILL,spacing=1,color='ff0000'XL
        xyouts, CO2min, ylabelpos, '!6CO2!X',charsize=1
     endif

     if xmin lt 19700. and xmax gt 17900. then begin     
        wbandmin=17900.0*sc
        wbandmax=19700.0*sc
        polyfill,[wbandmin,wbandmax,wbandmax,wbandmin],[ymin,ymin,ymax-border,ymax-border],/LINE_FILL,spacing=1,color='ff0000'XL
        xyouts, wbandmin, ylabelpos, '!6Water!X',charsize=1
     endif

     if xmin lt 20300. and xmax gt 20000. then begin     
        CO2min=20000.0*sc
        CO2max=20300.0*sc
        polyfill,[CO2min,CO2max,CO2max,CO2min],[ymin,ymin,ymax-border,ymax-border],/LINE_FILL,spacing=1,color='ff0000'XL
        xyouts, CO2min, ylabelpos, '!6CO2!X',charsize=1
     endif

     if xmin lt 20700. and xmax gt 20400. then begin     
        CO2min=20400.0*sc
        CO2max=20700.0*sc
        polyfill,[CO2min,CO2max,CO2max,CO2min],[ymin,ymin,ymax-border,ymax-border],/LINE_FILL,spacing=1,color='ff0000'XL
        xyouts, CO2min, ylabelpos, '!6CO2!X',charsize=1
     endif
endif
   
;Display near-IR OH lines if requested
if info.OH eq 1 then begin
     border = (ymax-ymin)*.05
     ;atmospheric absorbtion bands
     if xmax gt 11770. then begin
        polyfill,[11740.*sc,11770.*sc,11770.*sc,11740.*sc],[ymin,ymin,ymax-border,ymax-border],color='0000ff'XL, spacing=.1
        ;xyouts, 11750*sc, ylabelpos, '!6OH!X',charsize=1
     endif

     if xmax gt 10820. then begin
        polyfill,[10820.*sc,10850.*sc,10850.*sc,10820.*sc],[ymin,ymin,ymax-border,ymax-border],color='0000ff'XL, spacing=.1
        ;xyouts, 10830*sc, ylabelpos, '!6OH!X',charsize=1
     endif

     if xmax gt 11300. then begin
        polyfill,[11300.*sc,11330.*sc,11330.*sc,11300.*sc],[ymin,ymin,ymax-border,ymax-border],color='0000ff'XL, spacing=.1
        ;xyouts, 11310*sc, ylabelpos, '!6OH!X',charsize=1
     endif

     if xmax gt 12000. then begin
        polyfill,[12000.*sc,12030.*sc,12030.*sc,12000.*sc],[ymin,ymin,ymax-border,ymax-border],color='0000ff'XL, spacing=.1
        ;xyouts, 12010*sc, ylabelpos, '!6OH!X',charsize=1
     endif

     if xmax gt 12110. then begin
        polyfill,[12110.*sc,12140.*sc,12140.*sc,12110.*sc],[ymin,ymin,ymax-border,ymax-border],color='0000ff'XL, spacing=.1
        ;xyouts, 12120*sc, ylabelpos, '!6OH!X',charsize=1
     endif

     if xmax gt 12700. then begin
        polyfill,[12700.*sc,12730.*sc,12730.*sc,12700.*sc],[ymin,ymin,ymax-border,ymax-border],color='0000ff'XL, spacing=.1
        ;xyouts, 12710*sc, ylabelpos, '!6OH!X',charsize=1
     endif

     if xmax gt 12760. then begin
        polyfill,[12760.*sc,12810.*sc,12810.*sc,12760.*sc],[ymin,ymin,ymax-border,ymax-border],color='0000ff'XL, spacing=.1
        ;xyouts, 12770*sc, ylabelpos, '!6OH!X',charsize=1
     endif

     if xmax gt 14370. then begin
        polyfill,[14370.*sc,14400.*sc,14400.*sc,14370.*sc],[ymin,ymin,ymax-border,ymax-border],color='0000ff'XL, spacing=.1
        ;xyouts, 14380*sc, ylabelpos, '!6OH!X',charsize=1
     endif

     if xmax gt 14790. then begin
        polyfill,[14790.*sc,14810.*sc,14810.*sc,14790.*sc],[ymin,ymin,ymax-border,ymax-border],color='0000ff'XL, spacing=.1
        ;xyouts, 14800*sc, ylabelpos, '!6OH!X',charsize=1
     endif

     if xmax gt 15070. then begin
        polyfill,[15070.*sc,15100.*sc,15100.*sc,15070.*sc],[ymin,ymin,ymax-border,ymax-border],color='0000ff'XL, spacing=.1
        ;xyouts, 15080*sc, ylabelpos, '!6OH!X',charsize=1
     endif

     if xmax gt 15500. then begin
        polyfill,[15500.*sc,15530.*sc,15530.*sc,15500.*sc],[ymin,ymin,ymax-border,ymax-border],color='0000ff'XL, spacing=.1
        ;xyouts, 15510*sc, ylabelpos, '!6OH!X',charsize=1
     endif

     if xmax gt 15820. then begin
        polyfill,[15820.*sc,15850.*sc,15850.*sc,15820.*sc],[ymin,ymin,ymax-border,ymax-border],color='0000ff'XL, spacing=.1
        ;xyouts, 15830*sc, ylabelpos, '!6OH!X',charsize=1
     endif

     if xmax gt 16370. then begin
        polyfill,[16370.*sc,16400.*sc,16400.*sc,16370.*sc],[ymin,ymin,ymax-border,ymax-border],color='0000ff'XL, spacing=.1
        ;xyouts, 16380*sc, ylabelpos, '!6OH!X',charsize=1
     endif

     if xmax gt 16700. then begin
        polyfill,[16700.*sc,16750.*sc,16750.*sc,16700.*sc],[ymin,ymin,ymax-border,ymax-border],color='0000ff'XL, spacing=.1
        ;xyouts, 16720*sc, ylabelpos, '!6OH!X',charsize=1
     endif
    
  endif ;info.OH eq 1

;Plot the lines from template(s), if provided
idxlinetemps = where(linetemplates ne '0', count)

if count ne 0 then begin
   for temps = 0,count-1 do begin
      fmt = 'D,A'
      tempnum = idxlinetemps[temps]
      readcol,linetemplates[tempnum],F=fmt,linelambda,linename, /silent
      idxr = where(linename eq 'r', count)
      if count gt 0 then linename(idxr) = ''
      case tempnum of
         0: color = '00ff00'XL ;green
         1: color = '00FF5A'XL ;orange
         2: color = '0000ff'XL ;red
         3: color = 'ff7084'XL ;light blue
         4: color = '4782FF'XL ;sienna
         5: color = '4782FF'XL
         6: color = '00FFFF'XL ;yellow
         7: color = '00FFFF'XL
      endcase
      for i = 0,n_elements(linelambda)-1 do begin
         ;overplot the lines
         if linelambda[i]*(1.0+z) gt xmin and linelambda[i]*(1.0+z) lt xmax then begin
           oplot,[linelambda[i]*(1.0+z)*sc,linelambda[i]*(1.0+z)*sc],[ylmin,ylmax],color=color, linestyle=1
           xyouts, linelambda[i]*(1.0+z)*sc, ylabelpos, '!6'+linename[i]+'!X', orientatio=90, charsize=1.5
         endif
      endfor
   endfor
endif

;Plot template galaxy spectrum
if ptr_valid(info.spectemplateptr) ne 0 then begin
   template = *(info.spectemplateptr)
   templambda = template[0,*]
   tempspec = template[1,*]
   ;see if template and spectra overlap

   if info.spectempidx le 13 then begin  ;using galaxy template
      templmin=min(templambda)*(1.0+z)
      templmax=max(templambda)*(1.0+z)
  
      ovl=where(alllambda gt templmin and alllambda lt templmax and abs(allspec) lt 1e30 ,ovlcnt)

      if ovlcnt gt 0 then begin
        ;calculate the variance weighted renormalization for the template
        temp_spl = spline(templambda*(1.0+z),tempspec,alllambda(ovl),/DOUBLE)
        ;normfac=total( (allspec(ovl)/temp_spl)*allivar(ovl))/total(allivar(ovl))
        
        temp = (allspec(ovl)/temp_spl)*allivar(ovl)
        good = where(finite(temp), count)
        if count gt 0 then begin
           good_ovl = ovl[good]
           normfac = total(temp[good]) / total(allivar(good_ovl))
        endif else normfac = 1.0
        
        tempnorm = tempspec*normfac ;normalize template
        tempzp = mean(allspec(ovl)) ;get spec mean
        scaled_template = (tempnorm-tempzp)*info.templatescale + tempzp
        oplot,templambda*(1.0+z)*sc,scaled_template,color='ff00ff'XL 
      endif
   endif else begin 
        templmin=min(templambda)
        templmax=max(templambda)
  
      ovl=where(alllambda gt templmin and alllambda lt templmax and allspec lt 1e10 ,ovlcnt)

      if ovlcnt gt 0 then begin
        ;calculate the variance weighted renormalization for the template
        temp_spl = spline(templambda,tempspec,alllambda(ovl),/DOUBLE)
        normfac=total( (allspec(ovl)/temp_spl)*allivar(ovl))/total(allivar(ovl))
        tempnorm = tempspec*normfac ;normalize template
        tempzp = mean(allspec(ovl)) ;get spec mean
        scaled_template = (tempnorm-tempzp)*info.templatescale + tempzp
        print, '[L#4080] oplot,templambda*sc,scaled_template'
        oplot,templambda*sc,scaled_template,color='ff00ff'XL 
     endif
   endelse
   ; DEBUG
   print, min(templambda*(1.0+z)), max(templambda*(1.0+z))
   print, min(scaled_template), max(scaled_template)
   print, sc
endif

;if this call is due to auto-z, put the values used on the plot
if keyword_set(autoz_update) then begin 
   zinfoposy = ymin + (ymax-ymin)*0.02
   zinfoposx = xmax - (xmax-xmin)*0.45
   indx = info.corr_number
   z = string(info.corr_results[0,indx])
   z_err = string(info.corr_results[1,indx])
   chi2 = string(info.corr_results[2,indx])
   remchar, z, ' '
   remchar, z_err, ' '
   remchar, chi2, ' '
   zinfo = 'z = ' + z + ', z_err = '+ z_err + ', chi2 = ' + chi2
   xyouts, zinfoposx, zinfoposy, zinfo, charsize=1.5
endif

;if the user has manually selected 1D file, put the file name on the plot
if info.manualselect1d eq 1 then begin 
  nameposy = ymin + (ymax-ymin)*0.01
  nameposx = xmax - (xmax-xmin)*0.9
  namesplit = strsplit(info.spec1dfile,'/',/extract)
  spec1dfilename = namesplit[n_elements(namesplit)-1]
  xyouts, nameposx*sc, nameposy, 'Viewing: '+spec1dfilename, charsize=1.5
endif

;plot the pixel coords on the top of the window
;make array of pixel co-ords
pixstepsize = (numpts / 1000.) * 100
numpixsteps = ceil(float(numpts) / pixstepsize)
pix=intarr(numpixsteps)
for i = 0, n_elements(pix)-1 do begin
   pix[i] = i*pixstepsize
endfor
;find corresponding lambdas
xtick=dblarr(numpixsteps)
lam = sp.lambda
xtick = lam[pix]
;xtick=sp.lambda(pix)

pixlabel=string(pix)

;remove some lables to make it more readable
blanklabel=pixlabel
blanklabel(*)=' '
;pixlable(clearlable)=' '
xtickoff=xtick-xsize/150.0

AXIS,XAXIS=1, XTICKS=33,XTICKV=xtick,XTITLE='Pixel Position',XTICKLEN=0.06,XTICKN=blanklabel
AXIS,XAXIS=1, XTICKS=34,XTICKV=xtickoff,XTITLE='Pixel Position',XTICKLAYOUT=1,XTICKN=pixlabel

end ;makeplot

;------------------------------------------------------
pro make2Dplot, z, linetemplates, showspecpos, specpos, $
                event, info, reextract=reextract, zoombox=zoombox

  if info.manualselect1d eq 1 and info.manualselect2d eq 0 then return

  specpro_get_state,event,info   

  spec2d = *(info.spec2Dptr)

  ;Bin the spectrum
  specsize = size(spec2d.flux)
  xpix = specsize[1]
  ypix = specsize[2]
  
  ;First determine appropriate bin factors for x, y
  xbincands = [1., 2., 3., 4., 5., 6., 7., 8., 9., 10.]
  ybincands = [1., 2., 3., 4., 5.]
  xresult = xpix / xbincands
  yresult = ypix / ybincands

  ;the following should work with existing deimos data, 
  ;and allow greater flexibility for other spectra  
  xgoodidx = where(xresult le 1024.)
  ygoodidx = where(yresult le 80.)
  xrebin = xbincands[xgoodidx[0]]
  yrebin = ybincands[ygoodidx[0]]    

  while (xpix mod xrebin) ne 0 do begin
     xpix = xpix+1
  end
  while (ypix mod yrebin) ne 0 do begin
     ypix = ypix+1
  end

  flux = dblarr(xpix,ypix)
  flux(0:specsize[1]-1,0:specsize[2]-1)= spec2d.flux
  ivar = dblarr(xpix,ypix)
  ivar(0:specsize[1]-1,0:specsize[2]-1)= spec2d.ivar
  lambda_map = dblarr(xpix,ypix)
  lambda_map(0:specsize[1]-1,0:specsize[2]-1)= spec2d.lambda

  ivar = double(ivar)
  flux = double(flux)
  lambda_map = double(lambda_map)

  ;replace bad pix
  goodpix = where(ivar gt 0)
  if n_elements(goodpix) ne 1 then begin
     minpix = min(ivar(goodpix))
  endif else begin
     minpix = 0.0
  endelse
  badpix = where(ivar le minpix)
  ivar(badpix) = minpix

  newxpix = xpix / xrebin
  newypix = ypix / yrebin
  binnedspec = rebin(flux*ivar,newxpix,newypix)/ $
        sqrt(rebin(ivar,newxpix,newypix))
  ;interpolate wavelength solution
  lambda_map = congrid(lambda_map, newxpix, newypix)

  ;increase the x pixels for small 2D spectra 
  if xpix lt 800 then begin
     binnedspec = congrid(binnedspec,800,newypix)
     lambda_map = congrid(lambda_map,800,newypix)
     newxpix = 800
  endif

  ;increase the y pixels for small 2D images
  if ypix lt 40 then begin
     binnedspec = congrid(binnedspec,newxpix,40)
     lambda_map = congrid(lambda_map,newxpix,40)
     newypix = 40
     yrebin = ypix / 40.
  endif

  ;Find the max flux, for setting the brightness
  ;First check for and remove -NaN values
  missing = -1000
  idx = where(finite(binnedspec) eq 0, count)
  if (count gt 0) then binnedspec[idx] = missing
  maxflux = max(binnedspec)
  idx2 = where(binnedspec gt -1000)

  if n_elements(idx2) eq 1 then begin
     if idx2 eq -1 then return
  endif

  range = median(binnedspec(idx2))+[-3,info.upsigma]*stddev(binnedspec(idx2))
  
  fullrange = range[1]-range[0]
  contrastmod = fullrange * (info.contrast2D / 200.)  
  image = bytscl(binnedspec,min=range[0]+contrastmod,max=range[1],top=(!d.table_size-1))

  ;position = [0.1,0.1,0.9,0.9]
  ;xsize = (position[2]-position[0])*!d.x_vsize
  ;ysize = (position[3]-position[1])*!d.y_vsize
  xstart = (!d.x_vsize - newxpix)/2.0
  ystart = (!d.y_vsize - newypix)/2.0
  ypixtot = !d.y_vsize

  ;this code only called on a user reextraction event, triggered by a sliding click/release on 2D plot.
  if n_elements(reextract) ne 0 then begin
     specpos = (info.reextractpix - ystart)*yrebin
     new_spec_box = extract1d_gen(info.spec2dfile, specpos, info.extractwidth, /boxsprof)
     new_spec_opt = extract1d_gen(info.spec2dfile, specpos, info.extractwidth, /horne)
     dim_spec = size(new_spec_box.spec)
 
     new_spec1d_box = {flux:dblarr(dim_spec[1]), ivar:dblarr(dim_spec[1]), lambda:dblarr(dim_spec[1])}
     new_spec1d_opt = {flux:dblarr(dim_spec[1]), ivar:dblarr(dim_spec[1]), lambda:dblarr(dim_spec[1])}

     new_spec1d_box.flux = new_spec_box.spec
     new_spec1d_box.lambda = new_spec_box.lambda
     new_spec1d_box.ivar = new_spec_box.ivar

     new_spec1d_opt.flux = new_spec_opt.spec
     new_spec1d_opt.lambda = new_spec_opt.lambda
     new_spec1d_opt.ivar = new_spec_opt.ivar

     ;here change format so it matches new one
     ptr_free, info.spec1Dptr
     info.spec1Dptr = ptr_new(new_spec1d_box)
     ;update position of extraction
     info.extractpos = specpos
     specpro_set_state,event,info
     spec1d_plot_update,event,info
     widget_control, info.spec2d_id, get_value=winid
     wset,winid
  endif

  ;this code only called on a user zoom event, triggered by a click/release box on 2D plot.
  if n_elements(zoombox) ne 0 then begin
     minxpix = min([info.x2Dpress, info.x2Drelease])
     maxxpix = max([info.x2Dpress, info.x2Drelease])
     minypix = min([info.y2Dpress, info.y2Drelease])
     maxypix = max([info.y2Dpress, info.y2Drelease])

     minx = (minxpix - xstart)*xrebin
     maxx = (maxxpix - xstart)*xrebin
     miny = (minypix - ystart)*yrebin
     maxy = (maxypix - ystart)*yrebin    
  
     specsize = size(spec2d.flux)
     xnumpix = specsize[1]
     ynumpix  = specsize[2]
     ;check if box went off of 2D
     if minx le 0 then minx = 0
     if maxx ge xnumpix then maxx = xnumpix-1
     if miny le 0 then miny = 0
     if maxy ge ynumpix then maxy = ynumpix-1

     allspec = spec2d.flux
     
     idxgood = where(allspec gt -1000)     
     range = median(allspec(idxgood))+[-3,info.upsigma]*stddev(allspec(idxgood))
     fullrange = range[1]-range[0]
     contrastmod = fullrange * (info.contrast2D / 200.)  
     zoomimage = bytscl(allspec,min=range[0]+contrastmod,max=range[1],top=(!d.table_size-1))
     window, xsize = round(1.5*(maxx-minx)), ysize = round(1.5*(maxy-miny)), $
                 title = 'Unbinned zoom region'
     tv, zoomimage[minx:maxx, miny:maxy], round(.25*(maxx-minx)), round(.25*(maxy-miny))
     widget_control, info.spec2d_id, get_value = wid
     wset, wid
  endif

  ;display the spectrum
  tv, image, xstart, ystart
 
  ;get lambda from 2D wavelength map
  lambda = lambda_map[*,fix(specpos*float(newypix/ypix))]

  if info.space ne 1 then begin
     ;display A, B and water regions 
     aidxmin = fix(float(where(abs(lambda-7586.) eq min(abs(lambda-7586.)))))
     aidxmin = aidxmin[0]
     aidxmax = fix(float(where(abs(lambda-7703.) eq min(abs(lambda-7703.)))))
     aidxmax = aidxmax[0]
     bidxmin = fix(float(where(abs(lambda-6864.) eq min(abs(lambda-6864.)))))
     bidxmin = bidxmin[0]
     bidxmax = fix(float(where(abs(lambda-6945.) eq min(abs(lambda-6945.)))))
     bidxmax = bidxmax[0]
     wateridxmin = fix(float(where(abs(lambda-9200.) eq min(abs(lambda-9200.)))))
     wateridxmin = wateridxmin[0]
     wateridxmax = fix(float(where(abs(lambda-9700.) eq min(abs(lambda-9700.)))))
     wateridxmax = wateridxmax[0]
     water2idxmin = fix(float(where(abs(lambda-11000.) eq min(abs(lambda-11000.)))))
     water2idxmin = water2idxmin[0]
     water2idxmax = fix(float(where(abs(lambda-11600.) eq min(abs(lambda-11600.)))))
     water2idxmax = water2idxmax[0]
     water3idxmin = fix(float(where(abs(lambda-13400.) eq min(abs(lambda-13400.)))))
     water3idxmin = water3idxmin[0]
     water3idxmax = fix(float(where(abs(lambda-14500.) eq min(abs(lambda-14500.)))))
     water3idxmax = water3idxmax[0]
     water4idxmin = fix(float(where(abs(lambda-17900.) eq min(abs(lambda-17900.)))))
     water4idxmin = water4idxmin[0]
     water4idxmax = fix(float(where(abs(lambda-19700.) eq min(abs(lambda-19700.)))))
     water4idxmax = water4idxmax[0]
     O2idxmin = fix(float(where(abs(lambda-12500.) eq min(abs(lambda-12500.)))))
     O2idxmin = O2idxmin[0]
     O2idxmax = fix(float(where(abs(lambda-12800.) eq min(abs(lambda-12800.)))))
     O2idxmax = O2idxmax[0]
     CO2idxmin = fix(float(where(abs(lambda-15700.) eq min(abs(lambda-15700.)))))
     CO2idxmin = CO2idxmin[0]
     CO2idxmax = fix(float(where(abs(lambda-15800.) eq min(abs(lambda-15800.)))))
     CO2idxmax = CO2idxmax[0]
     CO22idxmin = fix(float(where(abs(lambda-16000.) eq min(abs(lambda-16000.)))))
     CO22idxmin = CO22idxmin[0]
     CO22idxmax = fix(float(where(abs(lambda-16200.) eq min(abs(lambda-16200.)))))
     CO22idxmax = CO22idxmax[0]
     CO23idxmin = fix(float(where(abs(lambda-20000.) eq min(abs(lambda-20000.)))))
     CO23idxmin = CO23idxmin[0]
     CO23idxmax = fix(float(where(abs(lambda-20300.) eq min(abs(lambda-20300.)))))
     CO23idxmax = CO23idxmax[0]
     CO24idxmin = fix(float(where(abs(lambda-20400.) eq min(abs(lambda-20400.)))))
     CO24idxmin = CO24idxmin[0]
     CO24idxmax = fix(float(where(abs(lambda-20700.) eq min(abs(lambda-20700.)))))
     CO24idxmax = CO24idxmax[0]

     aavg = fix(float(aidxmin+aidxmax)/2.0)
     bavg = fix(float(bidxmin+bidxmax)/2.0)
     wateravg = fix(float(wateridxmin+wateridxmax)/2.0)
     water2avg = fix(float(water2idxmin+water2idxmax)/2.0)
     water3avg = fix(float(water3idxmin+water3idxmax)/2.0)
     water4avg = fix(float(water4idxmin+water4idxmax)/2.0)
     O2avg = fix(float(O2idxmin+O2idxmax)/2.0)
     CO2avg = fix(float(CO2idxmin+CO2idxmax)/2.0)
     CO22avg = fix(float(CO22idxmin+CO22idxmax)/2.0)
     CO23avg = fix(float(CO23idxmin+CO23idxmax)/2.0)
     CO24avg = fix(float(CO24idxmin+CO24idxmax)/2.0)

     if aidxmin ne 0 and aidxmin ne n_elements(lambda)-1 then begin
        plots,[xstart+aidxmin,xstart+aidxmax],[ystart+newypix,ystart+newypix],thick=.5,linestyle=0,color='0000ff'XL,/device 
        plots,[xstart+aidxmin,xstart+aidxmax],[ystart,ystart],thick=.5,linestyle=0,color='0000ff'XL,/device
        xyouts,[xstart+aavg,xstart+aavg],ystart+fix(newypix*1.1), 'A Band', alignment=0.5,charsize=1.2,/device
     endif

     if bidxmin ne 0 and bidxmin ne n_elements(lambda)-1 then begin
        plots,[xstart+bidxmin,xstart+bidxmax],[ystart+newypix,ystart+newypix],thick=.5,linestyle=0,color='0000ff'XL,/device
        plots,[xstart+bidxmin,xstart+bidxmax],[ystart,ystart],thick=.5,linestyle=0,color='0000ff'XL,/device
        xyouts,[xstart+bavg,xstart+bavg],ystart+fix(newypix*1.1), 'B Band', alignment=0.5,charsize=1.2,/device
     endif

     if wateridxmin ne 0 and wateridxmin ne n_elements(lambda)-1 then begin
        plots,[xstart+wateridxmin,xstart+wateridxmax],[ystart+newypix,ystart+newypix],thick=.5,linestyle=0,color='0000ff'XL,/device
        plots,[xstart+wateridxmin,xstart+wateridxmax],[ystart,ystart],thick=.5,linestyle=0,color='0000ff'XL,/device
        xyouts,[xstart+wateravg,xstart+wateravg],ystart+fix(newypix*1.1), 'Water', alignment=0.5,charsize=1.2,/device
     endif

     if water2idxmin ne 0 and water2idxmin ne n_elements(lambda)-1 then begin
        plots,[xstart+water2idxmin,xstart+water2idxmax],[ystart+newypix,ystart+newypix],thick=.5,linestyle=0,color='0000ff'XL,/device
        plots,[xstart+water2idxmin,xstart+water2idxmax],[ystart,ystart],thick=.5,linestyle=0,color='0000ff'XL,/device
        xyouts,[xstart+water2avg,xstart+water2avg],ystart+fix(newypix*1.1), 'Water', alignment=0.5,charsize=1.2,/device
     endif
 
     if water3idxmin ne 0 and water3idxmin ne n_elements(lambda)-1 then begin
        plots,[xstart+water3idxmin,xstart+water3idxmax],[ystart+newypix,ystart+newypix],thick=.5,linestyle=0,color='0000ff'XL,/device
        plots,[xstart+water3idxmin,xstart+water3idxmax],[ystart,ystart],thick=.5,linestyle=0,color='0000ff'XL,/device
        xyouts,[xstart+water3avg,xstart+water3avg],ystart+fix(newypix*1.1), 'Water', alignment=0.5,charsize=1.2,/device
     endif

     if water4idxmin ne 0 and water4idxmin ne n_elements(lambda)-1 then begin
        plots,[xstart+water4idxmin,xstart+water4idxmax],[ystart+newypix,ystart+newypix],thick=.5,linestyle=0,color='0000ff'XL,/device
        plots,[xstart+water4idxmin,xstart+water4idxmax],[ystart,ystart],thick=.5,linestyle=0,color='0000ff'XL,/device
        xyouts,[xstart+water4avg,xstart+water4avg],ystart+fix(newypix*1.1), 'Water', alignment=0.5,charsize=1.2,/device
     endif

     if O2idxmin ne 0 and O2idxmin ne n_elements(lambda)-1 then begin
        plots,[xstart+O2idxmin,xstart+O2idxmax],[ystart+newypix,ystart+newypix],thick=.5,linestyle=0,color='0000ff'XL,/device
        plots,[xstart+O2idxmin,xstart+O2idxmax],[ystart,ystart],thick=.5,linestyle=0,color='0000ff'XL,/device
        xyouts,[xstart+O2avg,xstart+O2avg],ystart+fix(newypix*1.1), 'O2', alignment=0.5,charsize=1.2,/device
     endif

     if CO2idxmin ne 0 and CO2idxmin ne n_elements(lambda)-1 then begin
        plots,[xstart+CO2idxmin,xstart+CO2idxmax],[ystart+newypix,ystart+newypix],thick=.5,linestyle=0,color='0000ff'XL,/device
        plots,[xstart+CO2idxmin,xstart+CO2idxmax],[ystart,ystart],thick=.5,linestyle=0,color='0000ff'XL,/device
        xyouts,[xstart+CO2avg,xstart+CO2avg],ystart+fix(newypix*1.1), 'CO2', alignment=0.5,charsize=1.2,/device
     endif

     if CO22idxmin ne 0 and CO22idxmin ne n_elements(lambda)-1 then begin
        plots,[xstart+CO22idxmin,xstart+CO22idxmax],[ystart+newypix,ystart+newypix],thick=.5,linestyle=0,color='0000ff'XL,/device
        plots,[xstart+CO22idxmin,xstart+CO22idxmax],[ystart,ystart],thick=.5,linestyle=0,color='0000ff'XL,/device
        xyouts,[xstart+CO22avg,xstart+CO22avg],ystart+fix(newypix*1.1), 'CO2', alignment=0.5,charsize=1.2,/device
     endif

     if CO23idxmin ne 0 and CO23idxmin ne n_elements(lambda)-1 then begin
        plots,[xstart+CO23idxmin,xstart+CO23idxmax],[ystart+newypix,ystart+newypix],thick=.5,linestyle=0,color='0000ff'XL,/device
        plots,[xstart+CO23idxmin,xstart+CO23idxmax],[ystart,ystart],thick=.5,linestyle=0,color='0000ff'XL,/device
        xyouts,[xstart+CO23avg,xstart+CO23avg],ystart+fix(newypix*1.1), 'CO2', alignment=0.5,charsize=1.2,/device
     endif

     if CO24idxmin ne 0 and CO24idxmin ne n_elements(lambda)-1 then begin
        plots,[xstart+CO24idxmin,xstart+CO24idxmax],[ystart+newypix,ystart+newypix],thick=.5,linestyle=0,color='0000ff'XL,/device
        plots,[xstart+CO24idxmin,xstart+CO24idxmax],[ystart,ystart],thick=.5,linestyle=0,color='0000ff'XL,/device
        xyouts,[xstart+CO24avg,xstart+CO24avg],ystart+fix(newypix*1.1), 'CO2', alignment=0.5,charsize=1.2,/device
     endif

  endif

  ;Plot the lines from template(s), if provided
  idx = where(linetemplates ne '0', count)
  minlambda=min(lambda)
  maxlambda=max(lambda)

  ;following parameter indicate width of extraction in pixels 
  extractwidth = info.extractwidth

  if count ne 0 then begin
     for temps = 0,count-1 do begin
        fmt = 'D,A'
        tempnum = idx[temps]
        readcol,linetemplates[tempnum],F=fmt,linelambda,linename,/silent
        idxr = where(linename eq 'r', count)
        if count gt 0 then linename(idxr) = ''
        case tempnum of
           0: color = '00ff00'XL ;green
           1: color = '00FF5A'XL ;orange
           2: color = '0000ff'XL ;red
           3: color = 'ff7084'XL ;light blue
           4: color = '4287FF'XL ;sienna
           5: color = '4287FF'XL
           6: color = '00FFFF'XL ;yellow
           7: color = '00FFFF'XL
        endcase
        for i = 0,n_elements(linelambda)-1 do begin
           ;overplot the lines
           if linelambda[i]*(1.0+z) gt minlambda and linelambda[i]*(1.0+z) lt maxlambda then begin
             ;get pixel destination for this wavelength
             if extractwidth ne 0 then begin
                pix = fix(float(where(abs(lambda-linelambda[i]*(1.0+z)) eq min(abs(lambda-linelambda[i]*(1.0+z))))))
                plots,[xstart+pix,xstart+pix],[ystart,ystart+fix(float(specpos-extractwidth/2.)/yrebin)],thick=2.0, $
                      linestyle=0,color=color,/device
                plots,[xstart+pix,xstart+pix],[ystart+fix(float(specpos+extractwidth/2.)/yrebin),ystart+newypix], $
                      thick=2.0,linestyle=0,color=color,/device
             endif else begin
                pix = fix(float(where(abs(lambda-linelambda[i]*(1.0+z)) eq min(abs(lambda-linelambda[i]*(1.0+z))))))
                plots,[xstart+pix,xstart+pix],[ystart,ystart+fix(float(specpos-extractwidth/2.)/yrebin)],thick=2.0, $
                      linestyle=1,color=color,/device
                plots,[xstart+pix,xstart+pix],[ystart+fix(float(specpos+extractwidth/2.)/yrebin),ystart+newypix], $
                      thick=2.0,linestyle=1,color=color,/device
             endelse
   
             free_ypix = (ypixtot - newypix)/2.
             if free_ypix ge 30 then begin
                xyouts,[xstart+pix,xstart+pix],ystart+newypix+15, '!6'+linename[i]+'!X', $
                       alignment=0.5,charsize=1.1, orientation=90, /device
             endif else begin
                xyouts,[xstart+pix,xstart+pix],ystart+newypix+(free_ypix-15), '!6'+linename[i]+'!X', $
                       alignment=0.5,charsize=1.1, orientation=90, /device
             endelse

           endif
        endfor
     endfor
  endif

  ;Indicate the position of the extraction
  if showspecpos ne 0 then begin
     top = ystart+fix(float(fix(specpos+extractwidth/2.))/yrebin)
     bottom = ystart+fix(float(fix(specpos-extractwidth/2.))/yrebin)
     plots,[xstart,xstart+fix(newxpix*.05)],[top,top],thick=0.6,linestyle=0,color='ffffff'XL,/device
     plots,[xstart,xstart+fix(newxpix*.05)],[bottom,bottom],thick=0.6,linestyle=0,color='ffffff'XL,/device
     plots,[xstart+newxpix/2,xstart+newxpix/2+fix(newxpix*.05)],[top,top],thick=0.6,linestyle=0,color='ffffff'XL,/device
     plots,[xstart+newxpix/2,xstart+newxpix/2+fix(newxpix*.05)],[bottom,bottom],thick=0.6,linestyle=0,color='ffffff'XL,/device
     plots,[xstart+newxpix-fix(newxpix*.05),xstart+newxpix],[top,top],thick=0.6,linestyle=0,color='ffffff'XL,/device
     plots,[xstart+newxpix-fix(newxpix*.05),xstart+newxpix],[bottom,bottom],thick=0.6,linestyle=0,color='ffffff'XL,/device
  endif

  if n_elements(reextract) ne 0 then begin    
     ;give option to save  
     ;use current spec1d name, with serendip replacing the target name
     splitname = strsplit(info.spec1dfile,'.',/extract)
     serendipname = splitname[0]+'.'+splitname[1]+'.'+splitname[2]+'.serendip.'+splitname[4]
     get_serendip_radec, event, info, serendipra=serendipra, serendipdec=serendipdec
     print,'Serendip at RA/DEC: '+strcompress(string(serendipra),/remove_all)+', '+strcompress(string(serendipdec),/remove_all)
     if serendipra ne -99 then begin
        show_serendip_posradec, event, info, serendipra, serendipdec
     endif
     serendipfile = dialog_pickfile(/WRITE, title='Save reextracted spec1d fits file', file=serendipname)
     if serendipfile ne '' then begin
       mwrfits,new_spec1d_box, serendipfile
       mwrfits,new_spec1d_opt, serendipfile
       print, 'Reextracted spec1d saved.'
       splname = strsplit(serendipfile,'/',/extract)
       info.reextracted = 1
       ;update ra, dec
       info.targetRA = serendipra
       info.targetDEC = serendipdec
       ;save everything
       info.spec1dfile = splname[n_elements(splname)-1]
       specpro_set_state, event, info
       create_reextracted_infofile, event, info, serendipra, serendipdec
     endif
  end

end ;make2Dplot.pro

;-----------------------------------------------------------------------
pro photplot,file,redshift, event, info, SED=sed

   specpro_get_state, event, info
   fmt='A,F,F,F,F'   
   readcol,file,F=fmt,filter,l,dl,m,dm,/SILENT

   goodphot = where(dm GT 0 and dm LT 1,goodcnt)

   ;lambda min and max in angstroms
   lmin=0.14*10000L
   lmax=10*10000L

   ;adjustments to min/max data for limts
   dmin=0.5
   dmax=2
   
   if goodcnt GT 0 then begin
   
     ;get flux min/max, this is backwards because of magnitudes
     plotmin = min(m(goodphot)+dm(goodphot))
     plotmin = round(plotmin) - dmin

     plotmax = max(m(goodphot)+dm(goodphot))
     goodlim = where( m GT 50 and dm gt 1,nlim)

     ;check if our photometry limits are higher than the best measured flux
     if nlim GT 0 then begin
       limmax = max(dm(goodlim))
       if limmax GT plotmax THEN plotmax = limmax
     endif

     plotmax = round(plotmax) + dmax
   
     ;plot data 
     ;Change to allow plotting past 1e5 angstroms
     if max(l(goodphot)) gt 1e5 then begin
        plot,l(goodphot),m(goodphot),xrange=[lmin,max(l(goodphot))+1e4],yrange=[plotmax,plotmin],psym=1,XTITLE='Wavelength (A)',YTITLE='AB Magnitude',/XLOG
        oploterror,l(goodphot),m(goodphot),dl(goodphot)/2.0,dm(goodphot),psym=1
     endif else begin
        plot,l(goodphot),m(goodphot),xrange=[lmin,lmax],yrange=[plotmax,plotmin],psym=1,XTITLE='Wavelength (A)',YTITLE='AB Magnitude',/XLOG
        oploterror,l(goodphot),m(goodphot),dl(goodphot)/2.0,dm(goodphot),psym=1
     endelse

     ;plot limits if they exist
     if nlim GT 0 then begin
       plotsym,1
       oplot,l(goodlim),dm(goodlim),psym=8
     endif

     ;plot the breaks
     breaks=[912.0,1215.6,4000.0,16000];
     breakname=['912','1216','4000','1.6um']
     for I=0,n_elements(breaks)-1 do begin
      breaklam=(breaks[I]*(redshift+1.0))
      name='!6' + breakname[I] + '!X'
      oplot,[breaklam,breaklam],[(plotmax-dmax),(plotmin+dmin)],color='ffa300'XL, linestyle=1,THICK=5
      xyouts,breaklam,(plotmax-dmax/3.0),name,orientation=90,charsize=1.3
     endfor
     ;plot the spec limits
     if info.lambdamin ne 0 then begin
       lims = [info.lambdamin, info.lambdamax]
       limname=['minlam','maxlam']
       for i=0,1 do begin
         limlam=lims[i]
         name='!6' + limname[i] + '!X'
         oplot,[limlam,limlam],[(plotmax-dmax),(plotmin+dmin)],color='00ff00'XL, linestyle=1,THICK=1
         ;xyouts,limlam,(plotmax+dmax),name,orientation=90,charsize=1.3
       endfor
     endif

     ;plot the SED if its there
     if keyword_set(sed) then begin

        rlaws=[4]
        if info.spectempidx GT 0 then begin
           types = ['NONE','SB','PGAL','PGAL','SGAL','SGAL','SB','QSO','PGAL','SGAL','SB','SB','QSO','QSO','STAR','STAR','STAR','STAR','STAR','STAR']
           fit = fit_sed(sed,1,0.02,rlaws,redshift,l,dl,m,dm,TYPE=types[info.spectempidx])
       endif else begin
           fit = fit_sed(sed,1,0.02,rlaws,redshift,l,dl,m,dm)
       endelse
        oplot,fit.wave,-2.5*alog10(fit.flux)+8.9,linestyle=0,color='0000FF'XL,thick=1
        string=fit.name  + '!C' + 'e(B-V)=' + strmid(strtrim(string(fit.ebv),2),0,4)
        xyouts,(lmin+(lmax-lmin)*2.0/5.0),(plotmax-dmax/3.0),string,orientation=0,charsize=1.3

     endif

  endif else begin
    print,'No valid photometry found!'
    erase
  endelse

end ;photplot


;-----------------------------------------------------------------------
;The specpro.pro main program. Must be last to compile.
;-----------------------------------------------------------------------
pro specpro, slit, small=small, space=space, basic=basic, histogram=histogram, OH=OH
  
;Create the GUI and initialize info structure to hold information accessible by all routines.
  
;Slit is an optional input. If not set, specpro waits for user to
;select a slit.

if keyword_set(histogram) then begin
   hist=1
endif else begin
   hist=0
endelse
     
if keyword_set(small) ne 1 and keyword_set(basic) ne 1 then begin

;Create the top level base
tlb = widget_base(column=2,title='SpecPro')

baseleft = widget_base(tlb,column=1,frame=1)
baseright = widget_base(tlb,column=1)

;------------------------
;Create left side widgets
leftsize = 25
blank0 = widget_label(baseleft,value = '')

slitbase = widget_base(baseleft,row=1)
enter_slit_label = widget_label(slitbase,value = '   Enter slit number: ')
enter_slit = widget_text(slitbase, value = '0', event_pro='select_slit', xsize = 10, /editable)
blank1 = widget_label(baseleft, value = '')
baseleft7 = widget_base(baseleft,row=1)
previous_slit_button = widget_button(baseleft7,value='Previous',event_pro='previous_slit',xsize=110, /align_center)
next_slit_button = widget_button(baseleft7,value='Next',event_pro='next_slit', xsize=110, /align_center)

;Create text widget to enter redshift guesses
zbase = widget_base(baseleft,row=1)
zlabel = widget_label(zbase, value = 'Enter redshift guess: ')
zin = widget_text(zbase, value = '0.0', xsize = 10, event_pro='redshift_update', /editable)

;Create buttons to increment/decrement redshift 
baseleft8 = widget_base(baseleft,row=1)
zdowncoarse_button = widget_button(baseleft8,value=' << ', event_pro='decreasez_coarse',xsize=53,/align_center)
zdown_button = widget_button(baseleft8,value=' < ',event_pro='decreasez',xsize=53, /align_center)
zup_button = widget_button(baseleft8,value=' > ',event_pro='increasez', xsize=53, /align_center)
zupcoarse_button = widget_button(baseleft8,value=' >> ', event_pro='increasez_coarse',xsize=53,/align_center)

;Create new base widget for stamps
stampbase = widget_base(baseleft,column=2)
stampbaseleft = widget_base(stampbase,column=1)
stampbaseright = widget_base(stampbase,colum=1)

;Create stamp image widgets
u_label = widget_label(stampbaseleft,value='            ')
u_stamp = widget_draw(stampbaseleft,/align_center,event_pro='stamp_button_event',/button_events)
b_label = widget_label(stampbaseleft,value='            ')
b_stamp = widget_draw(stampbaseleft,/align_center,event_pro='stamp_button_event',/button_events)
g_label = widget_label(stampbaseleft,value='            ')
g_stamp = widget_draw(stampbaseleft,/align_center,event_pro='stamp_button_event',/button_events)
v_label = widget_label(stampbaseleft,value='            ')
v_stamp = widget_draw(stampbaseleft,/align_center,event_pro='stamp_button_event',/button_events)
r_label = widget_label(stampbaseleft,value='            ')
r_stamp = widget_draw(stampbaseleft,/align_center,event_pro='stamp_button_event',/button_events)
i_label = widget_label(stampbaseleft,value='            ')
i_stamp = widget_draw(stampbaseleft,/align_center,event_pro='stamp_button_event',/button_events)
z_label = widget_label(stampbaseright,value='           ')
z_stamp = widget_draw(stampbaseright,/align_center,event_pro='stamp_button_event',/button_events)
j_label = widget_label(stampbaseright,value='           ')
j_stamp = widget_draw(stampbaseright,/align_center,event_pro='stamp_button_event',/button_events)
h_label = widget_label(stampbaseright,value='           ')
h_stamp = widget_draw(stampbaseright,/align_center,event_pro='stamp_button_event',/button_events)
k_label = widget_label(stampbaseright,value='           ')
k_stamp = widget_draw(stampbaseright,/align_center,event_pro='stamp_button_event',/button_events)

;Create optional stamp window

optionstamplist = make_array(100,/string,value='             ')
optional_label = widget_label(stampbaseright, value = 'Other: ', xsize = 100)
optional = widget_draw(stampbaseright,/align_center,event_pro='stamp_button_event',/button_events)
optionstamplabel = widget_label(stampbaseright, value='Other options:', /align_center)
stampoptionlist = widget_list(stampbaseright, value=optionstamplist, /align_center, ysize=6, event_pro='set_option_stamp')

;Contrast slider for stamp plot
;contrast2D_label = widget_label(baseleft, value='  Increase contrast:', al)
contraststamp_slider = widget_slider(baseleft,value=0,min=0,max=100,xsize=150,title='Stamp contrast',event_pro='contrast_stamp_update')

;-------------------------
;Create right side widgets
;Spec1d plot
spec1d_label = widget_label(baseright,value='Spec1D')
spec1d = widget_draw(baseright,xsize=1100,ysize=320,frame=1,event_pro='zoom_spec1D_event',/button_events)

baseright2 = widget_base(baseright,row=1,frame=1)

;Make rebin list
rebinsizes =['  1','  2','  4','  8',' 16',' 32',' 64','128']
rebin_list = widget_droplist(baseright2, value=rebinsizes, title='Bin:', event_pro='rebin_update') 

;Make rebin list
smoothsizes =['  1','  3','  5','  7','  9',' 11',' 13',' 15']
smooth_list = widget_droplist(baseright2, value=smoothsizes, title='Smooth:',  $
                    event_pro='smooth_update')    

;Make button for spec templates
spectemplist = ['None', 'VVDS LBG', 'VVDS Elliptical','VVDS S0', $ 
                'VVDS Early Spiral', 'VVDS Spiral','VVDS Starburst', 'SDSS Quasar', $
                'Red Galaxy','Green Galaxy','Blue Galaxy','LBG Shapley','SDSS LoBAL',  $
                'SDSS HiBAL', 'A0 star', 'F0 star', $
                'G0 star', 'K0 star', 'M0 star', 'M6 star']

spectemp_list = widget_droplist(baseright2,value=spectemplist,title='Template:',event_pro='spec_template_update')

;make button for template scaling adjustment
spectempscalelist = ['1x','2x','3x','4x','5x']
spectempscale_list = widget_droplist(baseright2,value=spectempscalelist,title='Template scale:', event_pro='update_template_scale')

;list of solutions
corrlist = ['None']
autoz_list = widget_droplist(baseright2,value=corrlist,title='Auto-z solution:',event_pro= 'autoz_solution_update')
 
;Add a centering button. Does a redshift auto determination near current value
blankspace = widget_label(baseright2, value='')
recenter_button = widget_button(baseright2, value='Auto center', xsize=95, event_pro='recenter_update') 

blankspace = widget_label(baseright2, value=' ')
;Add a button to pull up another 1D file for currently selected slit
alt1d_button = widget_button(baseright2, value='Alt 1D file', xsize=95, event_pro='pick_alt1d_file') 

blankspace = widget_label(baseright2, value='')
newbasenonex = widget_base(baseright2,row=1,/nonexclusive) 
autoupdate_button = widget_button(newbasenonex, value='Auto-z OFF', event_pro='auto_update_set')

;Create another base for line templates
newbase = widget_base(baseright,row=1,frame=1)
buttonbase = widget_base(newbase, row=1, /nonexclusive)
emission_button = widget_button(buttonbase, value='Emission', event_pro='line_template_update')
highz_button = widget_button(buttonbase, value='High-z', event_pro='line_template_update')
seyfert_button = widget_button(buttonbase, value='Seyfert emission', event_pro='line_template_update')
qsoabs_button = widget_button(buttonbase, value='QSO absorption', event_pro='line_template_update')
qsoem_button = widget_button(buttonbase, value='QSO emission', event_pro='line_template_update')
prominent_qsoem_button = widget_button(buttonbase, value = 'QSO emission (subset)', event_pro='line_template_update')
elliptical_button = widget_button(buttonbase, value='Elliptical', event_pro='line_template_update')
prominent_elliptical_button = widget_button(buttonbase, value='Elliptical (subset)', event_pro='line_template_update')
;sky_base = widget_base(baseright2,row=1, /nonexclusive)
;blankspace = widget_label(baseright2, value='   ')
sky_button_on = widget_button(buttonbase, value='Show error', xsize=80, event_pro='show_sky_update')
blankspace = widget_label(newbase, value='  ')
spec1d_reset = widget_button(newbase, value='Reset zoom', xsize=100,event_pro='reset_zoom', /align_center)

;Spec2d plot
spec2d_label = widget_label(baseright,value='Spec2D')
spec2d = widget_draw(baseright,xsize=1100,ysize=150,frame=1,event_pro='spec2d_button_event',/button_events)

;Create new bases for the bottom right
baseright00 = widget_base(baseright,row=1,frame=1)

;Drop list allowing user to view either calibrated or uncalibrated 1D
;objectlist = ['Calibrated 1D','Uncalibrated 1D']
;object_list = widget_droplist(baseright00,value=objectlist,title='Showing:', event_pro='calibrated_update')
;Instead of the calibrated/uncalibrated, have the QUIT button here

;Create QUIT button
quit_buffer = widget_label(baseright00, value = '')
quit_button = widget_button(baseright00,value='QUIT',event_pro='quit_application',/align_center,xsize=90)

;Buttons for 2D spec
baseright3 = widget_base(baseright00,row=1, /nonexclusive)
baseright4 = widget_base(baseright,row=1,frame=1)

;Contrast slider for 2D plot
contrast2D_label = widget_label(baseright00, value=' Spec2D Contrast:')
contrast2D_slider = widget_slider(baseright00,value=0,min=0,max=100,xsize=150,event_pro='contrast_2d_update')

;Make sigma list, to control scaling of 2D contrast
sigma_buffer = widget_label(baseright00, value = ' ')
sigmas =string([10,9,8,7,6,5,4,3,2,1])
sigma_list = widget_droplist(baseright00, value=sigmas, title='Max sigma:',  $
                    event_pro='sigma_update') 

;Create a button to "like" the current z
blankspace = widget_label(baseright00, value = '')  
like_button = widget_button(baseright00,value='Like z', event_pro = 'like_redshift', xsize=80)

;Create a button for no z found
blankspace = widget_label(baseright00, value = ' ')  
noz_button = widget_button(baseright00,value='No z', event_pro = 'no_redshift', xsize=80)

blankspace = widget_label(baseright00, value='  ')
;Add a button to pull up another 2D file for currently selected slit
alt2d_button = widget_button(baseright00, value='Alt 2D file', xsize=95, event_pro='pick_alt2d_file')

blankspace = widget_label(baseright00, value='')

;select2Dfile_button = widget_button(baseright00, value='Select 2D File', xsize=120, event_pro='select_2d_file')
baseforsed = widget_base(baseright00,row=1)
baseforsed2 = widget_base(baseforsed, row=1,/nonexclusive)

sedfit_button = widget_button(baseforsed2,value='SED fit ON',event_pro='fit_sed_update',/align_right)

;features_button = widget_button(baseright3,value='Show telluric features',event_pro='display_features',/align_center)
extract_button = widget_button(baseright3,value='Show extraction',event_pro='display_extract_pos',/align_center)


;Create even more bases for bottom right
baseright5 = widget_base(baseright4,column=1)
baseright6 = widget_base(baseright4,column=1)

;get the spec1d mask name
fileslit0 = file_search('*spec1d.*.???.*.fits')
if n_elements(fileslit0) ne 1 then fileslit0 = fileslit0[0]
if fileslit0 ne '' then begin
   if n_elements(fileslit0) gt 1 then fileslit0 = fileslit0[0] ;; in case there is a serendipitous source
   strslit0 = strsplit(fileslit0[0], '.', /extract)
   spec1dmaskname = strslit0[1]
endif else begin
   spec1dmaskname = 'Unknown'
endelse

;get the stamp mask name
fileslit0 = file_search('stamps.*.???.*.fits')
if n_elements(fileslit0) gt 1 then fileslit0 = fileslit0[0] ;; in case there is a serendipitous source
if fileslit0 ne '' then begin
  strslit0 = strsplit(fileslit0, '.', /extract)
  stampmaskname = strslit0[1]
endif else begin
  stampmaskname = spec1dmaskname
endelse

;Display information about source
blank2 = widget_label(baseright5, value = '')
mask_name = widget_label(baseright5, value='Mask name:       '+stampmaskname, $
                        /align_left, frame=1, xsize=310)
target_name = widget_label(baseright5,value='Source Name:',/align_left,frame=1, xsize=310)
target_id = widget_label(baseright5,value='Source ID:',/align_left, frame=1, xsize=310)
target_ra = widget_label(baseright5,value='RA: ',/align_left, frame=1, xsize=310)
target_dec = widget_label(baseright5,value='DEC: ',/align_left, frame=1, xsize=310)
target_zphot = widget_label(baseright5,value='zphot: ',/align_left, frame=1, xsize=310)
target_zpdf = widget_label(baseright5,value='zpdf: ',/align_left, frame=1, xsize=310)
target_zpdf_low = widget_label(baseright5,value='zpdf_low: ',/align_left, frame=1, xsize=310)
target_zpdf_up = widget_label(baseright5,value='zpdf_up: ',/align_left, frame=1, xsize=310)

;Create output fields
;outputlabel = widget_label(baseright5,value='Output',/align_center,xsize=170)
outputbase1 = widget_base(baseright5,row=1)
zlabel2 = widget_label(outputbase1, value = 'Redshift:')
zoutput = widget_text(outputbase1, value = '', xsize = 8, /editable)
zlabel3 = widget_label(outputbase1, value = '      Confidence:')
zconfidence =  widget_text(outputbase1, value = '', xsize = 8, /editable)
outputbase2 = widget_base(baseright5,row=1)
noteslabel = widget_label(outputbase2,value = 'Notes:   ')
notes = widget_text(outputbase2, value='', xsize=8, /editable)
initlabel = widget_label(outputbase2, value = '        Initials:')
initials = widget_text(outputbase2, value='', xsize = 8, /editable, /align_center)

outputbase3 = widget_base(baseright5,row=1) 
buffer = widget_label(outputbase3, value = '   ')
output_zinfo_button = widget_button(outputbase3, value='Save Redshift', event_pro='save_output_to_file',                                          /align_center, xsize=120)
;Button for saving 1D to ascii
buffer = widget_label(outputbase3, value = '   ')
oned2ascii = widget_button(outputbase3,value='Save Spec1D',xsize=120, event_pro='save_1d_to_ascii')

;SED
sed_label = widget_label(baseright6,value='SED', /align_center)
sed = widget_draw(baseright6, xsize=750,ysize=300,frame=1,/align_center, $
                  event_pro='phot_button_event',/button_events)

widget_control, tlb, /realize

;this array will hold line template file names, if they are selected
linetemplates = ['0','0','0','0','0','0','0','0'] 

;store current 1D min, max xy vals
xrange = [0.,0.]
yrange = [0.,0.]

;set the default output file name
outfilename = stampmaskname+'_zinfo.dat'

;see if we're in space
if keyword_set(space) then begin
   space = 1
endif else begin
   space = 0
endelse

;display OH lines or not
if keyword_set(OH) then begin
   OH = 1
endif else begin
   OH = 0
endelse

sed_templates = read_seds(getenv('SPVIEW')+'/','sed.list',500,10e4,256,2048)

info = {  $
          tlb_id: tlb,  $ 
          spec1d_id:spec1d,  $ 
          spec2d_id:spec2d,  $
          sed_id:sed,  $
          target_name_id:target_name,  $
          target_ra_id:target_ra,  $
          target_dec_id:target_dec,  $
          target_zphot_id:target_zphot,  $
          target_zpdf_id:target_zpdf,  $
          target_zpdf_low_id:target_zpdf_low,  $
          target_zpdf_up_id:target_zpdf_up,  $
          target_id:target_id,  $
          targetRA:0.0,  $
          targetDEC:0.0,  $
          u_stamp_id:u_stamp,  $
          b_stamp_id:b_stamp,  $
          g_stamp_id:g_stamp,  $
          v_stamp_id:v_stamp,  $
          r_stamp_id:r_stamp,  $
          i_stamp_id:i_stamp,  $
          z_stamp_id:z_stamp,  $
          j_stamp_id:j_stamp,  $ 
          h_stamp_id:h_stamp,  $
          k_stamp_id:k_stamp,  $
          u_label_id:u_label,  $
          b_label_id:b_label,  $
          g_label_id:g_label,  $
          v_label_id:v_label,  $
          r_label_id:r_label,  $
          i_label_id:i_label,  $
          z_label_id:z_label,  $
          j_label_id:j_label,  $ 
          h_label_id:h_label,  $
          k_label_id:k_label,  $
          optional_stamp_id:optional,  $
          stampoptionlist_id:stampoptionlist,  $
          enter_slit_id:enter_slit,  $
          zoutput_id:zoutput,  $
          zconfidence_id:zconfidence,  $
          initials_id:initials,  $
          notes_id:notes,  $
          spectemp_list_id:spectemp_list,  $
          autoz_list_id:autoz_list,  $
          redshift:0.0,  $ 
          spec1dfile:'',  $
          manualselect1d:0,  $
          stampfile:'',  $
          spec2dfile:'',  $
          manualselect2d:0,  $
          photfile:'',  $
          infofile:'',  $
          slitnumber:0,  $
          optionstamp:'', $
          optionstamplist:optionstamplist,  $
          optional_label:optional_label, $
          zin:zin, $
          linetemplates:linetemplates,  $
          rebin:1,  $
          smooth:1,  $
          rebin_list:rebin_list,  $
          smooth_list:smooth_list,  $
          showfeatures:1,  $
          showspecpos:0,  $
          zoom:0,  $
          x1Dpress:0, $
          x1Drelease:0,  $
          y1Dpress:0,  $
          y1Drelease:0,  $
          x2Dpress:0, $
          x2Drelease:0,  $
          y2Dpress:0,  $
          y2Drelease:0,  $
          reextractpix:0,  $
          reextractflag:0,  $
          xrange:xrange,  $
          yrange:yrange,  $
          contrast2D:0,  $
          stampcontrast:0,  $
          extractpos:0.0,  $
          extractwidth:0.0, $
          showsky:0,  $
          drawslit:0,  $  ;controls whether slit is drawn on stamp images
          spec2Dptr:ptr_new(),  $
          spec1Dptr:ptr_new(),  $ 
          stampsptr:ptr_new(),  $
          spectemplateptr:ptr_new(),  $
          corner1_ra:0.0,  $
          corner2_ra:0.0,  $
          corner3_ra:0.0,  $
          corner4_ra:0.0,  $
          corner1_dec:0.0,  $
          corner2_dec:0.0,  $      
          corner3_dec:0.0,  $
          corner4_dec:0.0,  $
          templatescale:1.0,  $
          spec1dmaskname:spec1dmaskname,  $
          stampmaskname:stampmaskname,  $
          usingsmallversion:0,  $
          usingbasicversion:0,  $
          outputdatasaved:0,  $
          outputfilesaved:0,  $
          outfilename:outfilename,  $
          upsigma:10,  $ ;controls the contrast scaling on 2D plot
          spectempidx:0,  $
          missingstampfile:0,  $
          missing1dfile:0,  $
          missing2dfile:0,  $
          missingphotfile:0,  $
          missinginfofile:0,  $
          corr_number:0,  $ ;determines which solution (1-6) to show from auto-z
          corr_results:fltarr(3,6),  $ ;holds z, z_err, and chi2 from auto-z
          autozflag:0,  $ ;tells whether autoz has been computed
          lambdamin:0,  $
          lambdamax:0,   $  ;used in photplot
          usephotzforautoz:0,  $
          slitPA:90.0,  $
          slitwid:1.0,  $
          sed_templates:sed_templates,  $
          sed_fit_flag:0,  $
          auto_update_flag:1,  $
          RAclicked:0.0,  $
          DECclicked:0.0,  $
          reextracted:0,  $
          space:space,  $
          hist:hist,  $
          OH:OH  $
       }


;If user has given a slit number, call select_slit to initialize
if n_elements(slit) ne 0 then begin
   info.slitnumber = slit
   initevent = {top:info.tlb_id}
   select_slit,initevent,info,init=1
endif else begin
   infoptr = ptr_new(info)
   ;Store pointer to info structure in top level base
   widget_control, tlb, set_uvalue=infoptr
endelse

;Start managing events
xmanager,'specpro', tlb, /no_block

;----------------------------------------------------------------------------------------------
endif else if keyword_set(small) ne 0 and keyword_set(basic) eq 0 then begin 
;Using small version

;Create the top level base
tlb = widget_base(column=2,title='SpecPro')

baseleft = widget_base(tlb,column=1,frame=1)
baseright = widget_base(tlb,column=1)

;------------------------
;Create left side widgets
leftsize = 25
blank0 = widget_label(baseleft,value = '')

slitbase = widget_base(baseleft,row=1)
enter_slit_label = widget_label(slitbase,value = '   Enter slit number: ')
enter_slit = widget_text(slitbase, value = '0', event_pro='select_slit', xsize = 10, /editable)
blank1 = widget_label(baseleft, value = '')
baseleft7 = widget_base(baseleft,row=1)
previous_slit_button = widget_button(baseleft7,value='Previous',event_pro='previous_slit',xsize=110, /align_center)
next_slit_button = widget_button(baseleft7,value='Next',event_pro='next_slit', xsize=110, /align_center)

;Create text widget to enter redshift guesses
zbase = widget_base(baseleft,row=1)
zlabel = widget_label(zbase, value = 'Enter redshift guess: ')
zin = widget_text(zbase, value = '0.0', xsize = 10, event_pro='redshift_update', /editable)

;Create buttons to increment/decrement redshift 
baseleft8 = widget_base(baseleft,row=1)
zdowncoarse_button = widget_button(baseleft8,value=' << ', event_pro='decreasez_coarse',xsize=53,/align_center)
zdown_button = widget_button(baseleft8,value=' < ',event_pro='decreasez',xsize=53, /align_center)
zup_button = widget_button(baseleft8,value=' > ',event_pro='increasez', xsize=53, /align_center)
zupcoarse_button = widget_button(baseleft8,value=' >> ', event_pro='increasez_coarse',xsize=53,/align_center)

;Create new base widget for stamps
stampbase = widget_base(baseleft,column=2)
stampbaseleft = widget_base(stampbase,column=1)
stampbaseright = widget_base(stampbase,colum=1)

;Create stamp image widgets
u_label = widget_label(stampbaseleft,value='            ')
u_stamp = widget_draw(stampbaseleft,/align_center,event_pro='stamp_button_event',/button_events)
b_label = widget_label(stampbaseleft,value='            ')
b_stamp = widget_draw(stampbaseleft,/align_center,event_pro='stamp_button_event',/button_events)
v_label = widget_label(stampbaseleft,value='            ')
v_stamp = widget_draw(stampbaseleft,/align_center,event_pro='stamp_button_event',/button_events)
r_label = widget_label(stampbaseright,value='           ')
r_stamp = widget_draw(stampbaseright,/align_center,event_pro='stamp_button_event',/button_events)
i_label = widget_label(stampbaseright,value='           ')
i_stamp = widget_draw(stampbaseright,/align_center,event_pro='stamp_button_event',/button_events)


optionstamplist = make_array(100,/string,value='             ')

optional_label = widget_label(stampbaseright, value = 'Other: ', xsize = 100)
optional = widget_draw(stampbaseright,/align_center,event_pro='stamp_button_event',/button_events)
optionstamplabel = widget_label(baseleft, value='Other options:', /align_center)
stampoptionlist = widget_list(baseleft, value=optionstamplist, /align_center, ysize=3, event_pro='set_option_stamp')

;Contrast slider for stamp plot
contraststamp_slider = widget_slider(baseleft,value=0,min=0,max=100,xsize=150,title='Stamp contrast',event_pro='contrast_stamp_update')

;-------------------------
;Create right side widgets
;Spec1d plot
;spec1d_label = widget_label(baseright,value='Spec1D')
spec1d = widget_draw(baseright,xsize=1024,ysize=220,frame=1,event_pro='zoom_spec1D_event',/button_events)

baseright2 = widget_base(baseright,row=1,frame=1)

;Make rebin list
rebinsizes =['  1','  2','  4','  8',' 16',' 32',' 64','128']
rebin_list = widget_droplist(baseright2, value=rebinsizes, title='Bin:', event_pro='rebin_update') 

;Make rebin list
smoothsizes =['  1','  3','  5','  7','  9',' 11',' 13',' 15']
smooth_list = widget_droplist(baseright2, value=smoothsizes, title='Smooth:',  $
                    event_pro='smooth_update')    

;Make button for spec templates
spectemplist = ['None', 'VVDS LBG', 'VVDS Elliptical','VVDS S0', $ 
                'VVDS Early Spiral', 'VVDS Spiral','VVDS Starburst', 'SDSS Quasar', $
                'Red Galaxy','Green Galaxy','Blue Galaxy','LBG Shapley','SDSS LoBAL',  $
                'SDSS HiBAL', 'A0 star', 'F0 star', $
                'G0 star', 'K0 star', 'M0 star', 'M6 star']

spectemp_list = widget_droplist(baseright2,value=spectemplist,title='Template:',event_pro='spec_template_update')

;list of solutions
corrlist = ['None']
autoz_list = widget_droplist(baseright2,value=corrlist,title='Auto-z solution:',event_pro= 'autoz_solution_update')

;make button for template scaling adjustment
spectempscalelist = ['1x','2x','3x','4x','5x']
spectempscale_list = widget_droplist(baseright2,value=spectempscalelist,title='Temp scale:', event_pro='update_template_scale')

;Add a centering button. Does a redshift auto determination near current value
blankspace = widget_label(baseright2, value='')
recenter_button = widget_button(baseright2, value='Auto center', xsize=95, event_pro='recenter_update') 

blankspace = widget_label(baseright2, value='')
;Add a button to pull up another 1D file for currently selected slit
alt1d_button = widget_button(baseright2, value='Alt 1D', xsize=95, event_pro='pick_alt1d_file') 

blankspace = widget_label(baseright2, value='')
newbasenonex = widget_base(baseright2,row=1,/nonexclusive) 
autoupdate_button = widget_button(newbasenonex, value='Auto-z OFF', event_pro='auto_update_set')

;Create another base for line templates
newbase = widget_base(baseright,row=1,frame=1)
buttonbase = widget_base(newbase, row=1, /nonexclusive)
emission_button = widget_button(buttonbase, value='Emission', event_pro='line_template_update')
highz_button = widget_button(buttonbase, value='High-z', event_pro='line_template_update')
seyfert_button = widget_button(buttonbase, value='Seyfert emission', event_pro='line_template_update')
qsoabs_button = widget_button(buttonbase, value='QSO absorption', event_pro='line_template_update')
qsoem_button = widget_button(buttonbase, value='QSO emission', event_pro='line_template_update')
prominent_qsoem_button = widget_button(buttonbase, value='QSO emission (subset)', event_pro='line_template_update')
elliptical_button = widget_button(buttonbase, value='Elliptical', event_pro='line_template_update')
prominent_elliptical_button = widget_button(buttonbase, value='Elliptical (subset)', event_pro='line_template_update')
;sky_base = widget_base(baseright2,row=1, /nonexclusive)
sky_button_on = widget_button(buttonbase, value='Error', xsize=50, event_pro='show_sky_update')
blankspace = widget_label(newbase, value='  ')
spec1d_reset = widget_button(newbase, value='Reset zoom', xsize=90,event_pro='reset_zoom', /align_center) 

;Spec2d plot
;spec2d_label = widget_label(baseright,value='Spec2D')
spec2d = widget_draw(baseright,xsize=1024,ysize=100,frame=1,event_pro='spec2d_button_event',/button_events)

;Create new bases for the bottom right
baseright00 = widget_base(baseright,row=1,frame=1)

;Create QUIT button
quit_buffer = widget_label(baseright00, value = '')
quit_button = widget_button(baseright00,value='QUIT',event_pro='quit_application',/align_center,xsize=90)

;Buttons for 2D spec
baseright3 = widget_base(baseright00,row=1, /nonexclusive)
baseright4 = widget_base(baseright,row=1,frame=1)

;Contrast slider for 2D plot
contrast2D_label = widget_label(baseright00, value=' Spec2D Contrast:')
contrast2D_slider = widget_slider(baseright00,value=0,min=0,max=100,xsize=150,event_pro='contrast_2d_update')

;Make sigma list, to control scaling of 2D contrast
sigma_buffer = widget_label(baseright00, value = '')
sigmas =string([10,9,8,7,6,5,4,3,2,1])
sigma_list = widget_droplist(baseright00, value=sigmas, title='Max sigma:',  $
                    event_pro='sigma_update') 

;Create a button to "like" the current z
blankspace = widget_label(baseright00, value = '')  
like_button = widget_button(baseright00,value='Like z', event_pro = 'like_redshift', xsize=80)

;Create a button for no z found
blankspace = widget_label(baseright00, value = '')  
noz_button = widget_button(baseright00,value='No z', event_pro = 'no_redshift', xsize=80)

blankspace = widget_label(baseright00, value='')
;Add a button to pull up another 2D file for currently selected slit
alt2d_button = widget_button(baseright00, value='Alt 2D', xsize=95, event_pro='pick_alt2d_file')

baseforsed = widget_base(baseright00,row=1)
baseforsed2 = widget_base(baseforsed, row=1,/nonexclusive)

sedfit_button = widget_button(baseforsed2,value='SED fit ON',event_pro='fit_sed_update',/align_right)

;features_button = widget_button(baseright3,value='Show telluric features',event_pro='display_features',/align_center)
extract_button = widget_button(baseright3,value='Show extraction',event_pro='display_extract_pos',/align_center)

;Create even more bases for bottom right
baseright5 = widget_base(baseright4,column=1)
baseright6 = widget_base(baseright5,row=1)

;get the spec1d mask name
fileslit0 = file_search('*spec1d.*.???.*.fits')
if n_elements(fileslit0) ne 1 then fileslit0 = fileslit0[0]
if fileslit0 ne '' then begin
   if n_elements(fileslit0) gt 1 then fileslit0 = fileslit0[0] ;; in case there is a serendipitous source
   strslit0 = strsplit(fileslit0[0], '.', /extract)
   spec1dmaskname = strslit0[1]
endif else begin
   spec1dmaskname = 'Unknown'
endelse

;get the stamp mask name
;spawn, 'ls stamps.*.fits', fileslit0
fileslit0 = file_search('stamps.*.???.*.fits')
if n_elements(fileslit0) gt 1 then fileslit0 = fileslit0[0] ;; in case there is a serendipitous source
if fileslit0 ne '' then begin
  strslit0 = strsplit(fileslit0, '.', /extract)
  stampmaskname = strslit0[1]
endif else begin
  stampmaskname = spec1dmaskname
endelse

;Display information about source
blank2 = widget_label(baseright5, value = '')
mask_name = widget_label(baseleft, value='Mask name:       '+stampmaskname, $
                        /align_left, xsize=220)
target_name = widget_label(baseright6,value='Source Name:',/align_left, xsize=200)
target_id = widget_label(baseright6,value='Source ID:',/align_left, xsize=200)
baseright7 = widget_base(baseright5,row=1)
target_ra = widget_label(baseright7,value='RA: ',/align_left, xsize=200)
target_dec = widget_label(baseright7,value='DEC: ',/align_left, xsize=200)
baseright8 = widget_base(baseright5,row=1)
target_zphot = widget_label(baseright8,value='zphot: ',/align_left, xsize=200)
target_zpdf = widget_label(baseright8,value='zpdf: ',/align_left, xsize=200)
baseright9=widget_base(baseright5,row=1)
target_zpdf_low = widget_label(baseright9,value='zpdf_low: ',/align_left, xsize=200)
target_zpdf_up = widget_label(baseright9,value='zpdf_up: ',/align_left, xsize=200)

;Create output fields
;outputlabel = widget_label(baseright5,value='Output',/align_center,xsize=170)
outputbase1 = widget_base(baseright5,row=1)
zlabel2 = widget_label(outputbase1, value = 'Redshift:')
zoutput = widget_text(outputbase1, value = '', xsize = 16, /editable)
zlabel3 = widget_label(outputbase1, value = '      Confidence:')
zconfidence =  widget_text(outputbase1, value = '', xsize = 16, /editable)
outputbase2 = widget_base(baseright5,row=1)
noteslabel = widget_label(outputbase2,value = 'Notes:   ')
notes = widget_text(outputbase2, value='', xsize=16, /editable)
initlabel = widget_label(outputbase2, value = '        Initials:')
initials = widget_text(outputbase2, value='', xsize = 16, /editable, /align_center)

outputbase3 = widget_base(baseright5,row=1) 
buffer = widget_label(outputbase3, value = '          ')
output_zinfo_button = widget_button(outputbase3, value='Save Redshift', event_pro='save_output_to_file',                                          /align_center, xsize=130)
;Button for saving 1D to ascii
buffer = widget_label(outputbase3, value = '     ')
oned2ascii = widget_button(outputbase3,value='Save Spec1D',xsize=130, event_pro='save_1d_to_ascii')

;SED
sedbase = widget_base(baseright4,column=1)
;sed_label = widget_label(sedbase,value='SED', /align_center)
sed = widget_draw(sedbase, xsize=600,ysize=200,frame=1,/align_center, $
                  event_pro='save_image_to_file',/button_events)

widget_control, tlb, /realize

;this array will hold line template file names, if they are selected
linetemplates = ['0','0','0','0','0','0','0','0'] 

;store current 1D min, max xy vals
xrange = [0.,0.]
yrange = [0.,0.]

;set the default output file name
outfilename = stampmaskname+'_zinfo.dat'

;see if we're in space
if keyword_set(space) then begin
   space = 1
endif else begin
   space = 0
endelse

;display OH lines or not
if keyword_set(OH) then begin
   OH = 1
endif else begin
   OH = 0
endelse

sed_templates = read_seds(getenv('SPVIEW')+'/','sed.list',500,10e4,256,2048)

info = {  $
          tlb_id: tlb,  $ 
          spec1d_id:spec1d,  $ 
          spec2d_id:spec2d,  $
          sed_id:sed,  $
          target_name_id:target_name,  $
          target_ra_id:target_ra,  $
          target_dec_id:target_dec,  $
          target_zphot_id:target_zphot,  $
          target_zpdf_id:target_zpdf,  $
          target_zpdf_low_id:target_zpdf_low,  $
          target_zpdf_up_id:target_zpdf_up,  $
          target_id:target_id,  $
          targetRA:0.0,  $
          targetDEC:0.0,  $
          u_stamp_id:u_stamp,  $
          b_stamp_id:b_stamp,  $
          v_stamp_id:v_stamp,  $
          r_stamp_id:r_stamp,  $
          i_stamp_id:i_stamp,  $
          u_label_id:u_label,  $
          b_label_id:b_label,  $
          v_label_id:v_label,  $
          r_label_id:r_label,  $
          i_label_id:i_label,  $
          optional_stamp_id:optional,  $
          stampoptionlist_id:stampoptionlist,  $
          enter_slit_id:enter_slit,  $
          zoutput_id:zoutput,  $
          zconfidence_id:zconfidence,  $
          initials_id:initials,  $
          notes_id:notes,  $
          spectemp_list_id:spectemp_list,  $
          autoz_list_id:autoz_list,  $
          redshift:0.0,  $ 
          spec1dfile:'',  $
          manualselect1d:0,  $
          stampfile:'',  $
          spec2dfile:'',  $
          manualselect2d:0,  $
          photfile:'',  $
          infofile:'',  $
          slitnumber:0,  $
          optionstamp:'', $
          optionstamplist:optionstamplist,  $
          optional_label:optional_label, $
          zin:zin, $
          linetemplates:linetemplates,  $
          rebin:1,  $
          smooth:1,  $
          rebin_list:rebin_list,  $
          smooth_list:smooth_list,  $
          showfeatures:1,  $
          showspecpos:0,  $
          zoom:0,  $
          x1Dpress:0, $
          x1Drelease:0,  $
          y1Dpress:0,  $
          y1Drelease:0,  $
          x2Dpress:0, $
          x2Drelease:0,  $
          y2Dpress:0,  $
          y2Drelease:0,  $
          reextractpix:0,  $
          reextractflag:0,  $
          xrange:xrange,  $
          yrange:yrange,  $
          contrast2D:0,  $
          stampcontrast:0,  $
          extractpos:0.0,  $
          extractwidth:0.0, $
          showsky:0,  $
          drawslit:0,  $  ;controls whether slit is drawn on stamp images
          spec2Dptr:ptr_new(),  $
          spec1Dptr:ptr_new(),  $ 
          stampsptr:ptr_new(),  $
          spectemplateptr:ptr_new(),  $
          corner1_ra:0.0,  $
          corner2_ra:0.0,  $
          corner3_ra:0.0,  $
          corner4_ra:0.0,  $
          corner1_dec:0.0,  $
          corner2_dec:0.0,  $      
          corner3_dec:0.0,  $
          corner4_dec:0.0,  $
          templatescale:1.0,  $
          spec1dmaskname:spec1dmaskname,  $
          stampmaskname:stampmaskname,  $
          usingsmallversion:1,  $
          usingbasicversion:0,  $
          outputdatasaved:0,  $
          outputfilesaved:0,  $
          outfilename:outfilename,  $
          upsigma:10,  $ ;controls the contrast scaling on 2D plot
          spectempidx:0,  $
          missingstampfile:0,  $
          missing1dfile:0,  $
          missing2dfile:0,  $
          missingphotfile:0,  $
          missinginfofile:0,  $
          corr_number:0,  $ ;determines which solution (1-6) to display from auto-z
          corr_results:fltarr(3,6),  $ ;holds z, z_err, and chi2 from auto-z
          autozflag:0,  $ ;tells whether autoz has been computed
          lambdamin:0,  $
          lambdamax:0,   $  ;used in photplot
          usephotzforautoz:0,  $
          slitPA:90.0,  $
          slitwid:1.0,  $
          sed_templates:sed_templates,  $
          sed_fit_flag:0,  $
          auto_update_flag:1,  $
          RAclicked:0.0,  $
          DECclicked:0.0,  $
          reextracted:0,  $
          space:space,  $
          hist:hist,  $
          OH:OH  $
       }


;If user has given a slit number, call select_slit to initialize
if n_elements(slit) ne 0 then begin
   info.slitnumber = slit
   initevent = {top:info.tlb_id}
   select_slit,initevent,info,init=1
endif else begin
   infoptr = ptr_new(info)
   ;Store point to info structure in top level base
   widget_control, tlb, set_uvalue=infoptr
endelse

;Start managing events
xmanager,'specpro', tlb, /no_block

;------------------------------------------------------------------------------------------------
endif else if (keyword_set(basic) ne 0 and keyword_set(small) eq 0) then begin 
;Using basic version

 ;Create the top level base
 tlb = widget_base(column=2,title='SpecPro')

 baseleft = widget_base(tlb,column=1,frame=1)
 baseright = widget_base(tlb,column=1)

 ;------------------------
 ;Create left side widgets
 leftsize = 25
 blank0 = widget_label(baseleft,value = '')

 slitbase = widget_base(baseleft,row=1)
 enter_slit_label = widget_label(slitbase,value = '   Enter slit number: ')
 enter_slit = widget_text(slitbase, value = '0', event_pro='select_slit', xsize = 10, /editable)
 blank1 = widget_label(baseleft, value = '')
 baseleft7 = widget_base(baseleft,row=1)
 previous_slit_button = widget_button(baseleft7,value='Previous',event_pro='previous_slit',xsize=110, /align_center)
 next_slit_button = widget_button(baseleft7,value='Next',event_pro='next_slit', xsize=110, /align_center)

 ;Create text widget to enter redshift guesses
 zbase = widget_base(baseleft,row=1)
 zlabel = widget_label(zbase, value = 'Enter redshift guess: ')
 zin = widget_text(zbase, value = '0.0', xsize = 10, event_pro='redshift_update', /editable)

 ;Create buttons to increment/decrement redshift 
 baseleft8 = widget_base(baseleft,row=1)
 zdowncoarse_button = widget_button(baseleft8,value=' << ', event_pro='decreasez_coarse',xsize=53,/align_center)
 zdown_button = widget_button(baseleft8,value=' < ',event_pro='decreasez',xsize=53, /align_center)
 zup_button = widget_button(baseleft8,value=' > ',event_pro='increasez', xsize=53, /align_center)
 zupcoarse_button = widget_button(baseleft8,value=' >> ', event_pro='increasez_coarse',xsize=53,/align_center)

 ;get the spec1d mask name
 fileslit0 = file_search('*spec1d.*.???.*.fits')
 if n_elements(fileslit0) ne 1 then fileslit0 = fileslit0[0]
 if fileslit0 ne '' then begin
    if n_elements(fileslit0) gt 1 then fileslit0 = fileslit0[0] ;; in case there is a serendipitous source
    strslit0 = strsplit(fileslit0, '.', /extract)
    spec1dmaskname = strslit0[1]
 endif else begin
    spec1dmaskname = 'Unknown'
 endelse

 ;get the stamp mask name
 fileslit0 = file_search('stamps.*.???.*.fits')
 if n_elements(fileslit0) gt 1 then fileslit0 = fileslit0[0] ;; in case there is a serendipitous source
 if fileslit0 ne '' then begin
   strslit0 = strsplit(fileslit0, '.', /extract)
   stampmaskname = strslit0[1]
 endif else begin
   stampmaskname = spec1dmaskname
 endelse

 ;Display information about source
 blank2 = widget_label(baseleft, value = '')
 mask_name = widget_label(baseleft, value='Mask name:       '+stampmaskname, $
                         /align_left, frame=1, xsize=220)
 target_name = widget_label(baseleft,value='Source Name:',/align_left,frame=1, xsize=220)
 target_id = widget_label(baseleft,value='Source ID:',/align_left, frame=1, xsize=220)
 target_ra = widget_label(baseleft,value='RA: ',/align_left, frame=1, xsize=220)
 target_dec = widget_label(baseleft,value='DEC: ',/align_left, frame=1, xsize=220)
 target_zphot = widget_label(baseleft,value='zphot: ',/align_left, frame=1, xsize=220)
 target_zpdf = widget_label(baseleft,value='zpdf: ',/align_left, frame=1, xsize=220)
 target_zpdf_low = widget_label(baseleft,value='zpdf_low: ',/align_left, frame=1, xsize=220)
 target_zpdf_up = widget_label(baseleft,value='zpdf_up: ',/align_left, frame=1, xsize=220)

 ;Create output fields
 zlabel2 = widget_label(baseleft, value = 'Redshift:')
 zoutput = widget_text(baseleft, value = '', xsize = 20, /editable, /align_center)
 zlabel3 = widget_label(baseleft, value = 'Confidence:')
 zconfidence =  widget_text(baseleft, value = '', xsize = 20, /editable, /align_center)
 outputbase2 = widget_base(baseleft,row=1)
 initlabel = widget_label(baseleft, value = 'Initials:')
 initials = widget_text(baseleft, value='', xsize = 20, /editable, /align_center)
 noteslabel = widget_label(baseleft,value = 'Notes:   ')
 notes = widget_text(baseleft, value='', xsize=20, /editable, /align_center)
  
 outputbase = widget_base(baseleft, row = 1)
 output_zinfo_button = widget_button(outputbase, value='Save Redshift', event_pro='save_output_to_file',/align_center, xsize=110)
 ;Button for saving 1D to ascii
 oned2ascii = widget_button(outputbase,value='Save Spec1D', event_pro='save_1d_to_ascii',/align_center, xsize=110)
 ;-------------------------
 ;Create right side widgets
 ;Spec1d plot
 spec1d_label = widget_label(baseright,value='Spec1D')
 spec1d = widget_draw(baseright,xsize=1100,ysize=320,frame=1,event_pro='zoom_spec1D_event',/button_events)

 baseright2 = widget_base(baseright,row=1,frame=1)

 ;Make rebin list
 rebinsizes =['  1','  2','  4','  8',' 16',' 32',' 64','128']
 rebin_list = widget_droplist(baseright2, value=rebinsizes, title='Bin:', event_pro='rebin_update') 

 ;Make rebin list
 smoothsizes =['  1','  3','  5','  7','  9',' 11',' 13',' 15']
 smooth_list = widget_droplist(baseright2, value=smoothsizes, title='Smooth:',  $
                     event_pro='smooth_update')    
 ;Make button for spec templates
 spectemplist = ['None', 'VVDS LBG', 'VVDS Elliptical','VVDS S0', $ 
                 'VVDS Early Spiral', 'VVDS Spiral','VVDS Starburst', 'SDSS Quasar', $
                 'Red Galaxy','Green Galaxy','Blue Galaxy','LBG Shapley','SDSS LoBAL',  $
                 'SDSS HiBAL', 'A0 star', 'F0 star', $
                 'G0 star', 'K0 star', 'M0 star', 'M6 star']
 
 spectemp_list = widget_droplist(baseright2,value=spectemplist,title='Template:',event_pro='spec_template_update')
 
 ;list of solutions
 corrlist = ['None']
 autoz_list = widget_droplist(baseright2,value=corrlist,title='Auto-z solution:',event_pro= 'autoz_solution_update')

 ;make button for template scaling adjustment
 spectempscalelist = ['1x','2x','3x','4x','5x']
 spectempscale_list = widget_droplist(baseright2,value=spectempscalelist,title='Temp scale:', event_pro='update_template_scale')

 ;Add a centering button. Does a redshift auto determination near current value
 blankspace = widget_label(baseright2, value='')
 recenter_button = widget_button(baseright2, value='Auto center', xsize=95, event_pro='recenter_update') 

 blankspace = widget_label(baseright2, value='')
 ;Add a button to pull up another 1D file for currently selected slit
 alt1d_button = widget_button(baseright2, value='Alt 1D', xsize=95, event_pro='pick_alt1d_file') 

 blankspace = widget_label(baseright2, value='')
 newbasenonex = widget_base(baseright2,row=1,/nonexclusive) 
 autoupdate_button = widget_button(newbasenonex, value='Auto-z OFF', event_pro='auto_update_set')
 
 ;Create another base for line templates
 newbase = widget_base(baseright,row=1,frame=1)
 buttonbase = widget_base(newbase, row=1, /nonexclusive)
 emission_button = widget_button(buttonbase, value='Emission', event_pro='line_template_update')
 highz_button = widget_button(buttonbase, value='High-z', event_pro='line_template_update')
 seyfert_button = widget_button(buttonbase, value='Seyfert emission', event_pro='line_template_update')
 qsoabs_button = widget_button(buttonbase, value='QSO absorption', event_pro='line_template_update')
 qsoem_button = widget_button(buttonbase, value='QSO emission', event_pro='line_template_update')
 prominent_qsoem_button = widget_button(buttonbase, value = 'QSO emission (subset)', event_pro='line_template_update')
 elliptical_button = widget_button(buttonbase, value='Elliptical', event_pro='line_template_update')
 prominent_elliptical_button = widget_button(buttonbase, value='Elliptical (subset)', event_pro='line_template_update')
 ;sky_base = widget_base(baseright2,row=1, /nonexclusive)
 ;blankspace = widget_label(baseright2, value='   ')
 sky_button_on = widget_button(buttonbase, value='Show error', xsize=80, event_pro='show_sky_update')
 blankspace = widget_label(newbase, value=' ')
 spec1d_reset = widget_button(newbase, value='Reset zoom', xsize=100,event_pro='reset_zoom', /align_center)

 ;Spec2d plot
 spec2d_label = widget_label(baseright,value='Spec2D')
 spec2d = widget_draw(baseright,xsize=1100,ysize=150,frame=1,event_pro='spec2d_button_event',/button_events)

 ;Create new bases for the bottom right
 baseright00 = widget_base(baseright,row=1,frame=1)

 ;Drop list allowing user to view either calibrated or uncalibrated 1D
 ;objectlist = ['Calibrated 1D','Uncalibrated 1D']
 ;object_list = widget_droplist(baseright00,value=objectlist,title='Showing:', event_pro='calibrated_update')
 ;Instead of the calibrated/uncalibrated, have the QUIT button here

 ;Buttons for 2D spec
 baseright3 = widget_base(baseright00,row=1, /nonexclusive)
 baseright4 = widget_base(baseright,row=1,frame=1)

 ;Create QUIT button
 quit_buffer = widget_label(baseright00, value = '')
 quit_button = widget_button(baseright00,value='QUIT',event_pro='quit_application',/align_center,xsize=90)
 
 ;Buttons for 2D spec
 baseright3 = widget_base(baseright00,row=1, /nonexclusive)
 baseright4 = widget_base(baseright,row=1,frame=1)

 ;Contrast slider for 2D plot
 contrast2D_label = widget_label(baseright00, value=' Spec2D Contrast:')
 contrast2D_slider = widget_slider(baseright00,value=0,min=0,max=100,xsize=150,event_pro='contrast_2d_update')

 ;Make sigma list, to control scaling of 2D contrast
 sigma_buffer = widget_label(baseright00, value = '')
 sigmas =string([10,9,8,7,6,5,4,3,2,1])
 sigma_list = widget_droplist(baseright00, value=sigmas, title='Max sigma:',  $
                     event_pro='sigma_update') 

 ;Create a button to "like" the current z
 blankspace = widget_label(baseright00, value = '')  
 like_button = widget_button(baseright00,value='Like z', event_pro = 'like_redshift', xsize=80)

 ;Create a button for no z found
 blankspace = widget_label(baseright00, value = '')  
 noz_button = widget_button(baseright00,value='No z', event_pro = 'no_redshift', xsize=80)

 blankspace = widget_label(baseright00, value='')
 ;Add a button to pull up another 2D file for currently selected slit
 alt2d_button = widget_button(baseright00, value='Alt 2D', xsize=95, event_pro='pick_alt2d_file')

 extract_button = widget_button(baseright3,value='Show extraction',event_pro='display_extract_pos',/align_center)

 widget_control, tlb, /realize

 ;this array will hold line template file names, if they are selected
 linetemplates = ['0','0','0','0','0','0','0','0'] 

 ;store current 1D min, max xy vals
 xrange = [0.,0.]
 yrange = [0.,0.]

 ;see if we're in space
 if keyword_set(space) then begin
    space = 1
 endif else begin
    space = 0
 endelse

;display OH lines or not
if keyword_set(OH) then begin
   OH = 1
endif else begin
   OH = 0
endelse

 ;set the default output file name
 outfilename = stampmaskname+'_zinfo.dat'

 info = {  $
           tlb_id: tlb,  $ 
           spec1d_id:spec1d,  $ 
           spec2d_id:spec2d,  $
           target_name_id:target_name,  $
           target_ra_id:target_ra,  $
           target_dec_id:target_dec,  $
           target_zphot_id:target_zphot,  $
           target_zpdf_id:target_zpdf,  $
           target_zpdf_low_id:target_zpdf_low,  $
           target_zpdf_up_id:target_zpdf_up,  $
           target_id:target_id,  $
           targetRA:0.0,  $
           targetDEC:0.0,  $
           enter_slit_id:enter_slit,  $
           zoutput_id:zoutput,  $
           zconfidence_id:zconfidence,  $
           initials_id:initials,  $
           notes_id:notes,  $
           spectemp_list_id:spectemp_list,  $
           autoz_list_id:autoz_list,  $
           redshift:0.0,  $ 
           spec1dfile:'',  $
           manualselect1d:0,  $
           stampfile:'',  $
           spec2dfile:'',  $
           manualselect2d:0,  $
           photfile:'',  $
           infofile:'',  $
           slitnumber:0,  $
           zin:zin, $
           linetemplates:linetemplates,  $
           rebin:1,  $
           smooth:1,  $
           rebin_list:rebin_list,  $
           smooth_list:smooth_list,  $
           showspecpos:0,  $
           zoom:0,  $
           x1Dpress:0, $
           x1Drelease:0,  $
           y1Dpress:0,  $
           y1Drelease:0,  $
           x2Dpress:0, $
           x2Drelease:0,  $
           y2Dpress:0,  $
           y2Drelease:0,  $
           reextractpix:0,  $
           reextractflag:0,  $
           xrange:xrange,  $
           yrange:yrange,  $
           contrast2D:0,  $
           stampcontrast:0,  $
           extractpos:0.0,  $
           extractwidth:0.0, $
           showsky:0,  $
           drawslit:0,  $  ;controls whether slit is drawn on stamp images
           spec2Dptr:ptr_new(),  $
           spec1Dptr:ptr_new(),  $ 
           stampsptr:ptr_new(),  $
           spectemplateptr:ptr_new(),  $
           corner1_ra:0.0,  $
           corner2_ra:0.0,  $
           corner3_ra:0.0,  $
           corner4_ra:0.0,  $
           corner1_dec:0.0,  $
           corner2_dec:0.0,  $      
           corner3_dec:0.0,  $
           corner4_dec:0.0,  $
           templatescale:1.0,  $
           spec1dmaskname:spec1dmaskname,  $
           stampmaskname:stampmaskname,  $
           usingsmallversion:0,  $
           usingbasicversion:1,  $
           outputdatasaved:0,  $
           outputfilesaved:0,  $
           outfilename:outfilename,  $
           upsigma:10,  $ ;controls the contrast scaling on 2D plot
           spectempidx:0,  $
           missingstampfile:0,  $
           missing1dfile:0,  $
           missing2dfile:0,  $
           missingphotfile:0,  $
           missinginfofile:0,  $
           corr_number:0,  $ ;determines which solution (1-6) to show from auto-z
           corr_results:fltarr(3,6),  $ ;holds z, z_err, and chi2 from auto-z
           autozflag:0,  $ ;tells whether autoz has been computed
           lambdamin:0,  $
           lambdamax:0,   $  ;used in photplot
           usephotzforautoz:0,  $
           auto_update_flag:1,  $
           reextracted:0,  $
           RAclicked:0,  $
           DECclicked:0,  $
           space:space,  $
           hist:hist,  $
           OH:OH  $
        }


 ;If user has given a slit number, call select_slit to initialize
 if n_elements(slit) ne 0 then begin
    info.slitnumber = slit
    initevent = {top:info.tlb_id}
    select_slit,initevent,info,init=1
 endif else begin
    infoptr = ptr_new(info)
    ;Store pointer to info structure in top level base
    widget_control, tlb, set_uvalue=infoptr
 endelse

 ;Start managing events
 xmanager,'specpro', tlb, /no_block

endif else begin ;using basic version
   print, 'Select basic or small version, but not both.'
   return
endelse


end ;specpro.pro
