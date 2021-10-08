/* $Id: segyhdr.h,v 1.4 1998/03/11 14:48:51 seistool Exp seistool $ */
#ifndef SEGYHDR_H
#define SEGYHDR_H

/* Definition of TAG for trace_headers */

#define SATAG_TRACEHEADER	33002
#define TRCBYTES		240

/*
* declarations for:
*      typedef struct {} segyhdr - th trace identification header
*
*  8-March -1996 changed to bitfields for Cray and better portability
*/

typedef struct {	/* segy - trace identification header */

	signed tracl   :32;	/* trace sequence number within line */

	signed tracr   :32;	/* trace sequence number within reel */

	signed fldr    :32;	/* field record number */

	signed tracf   :32;	/* trace number within field record */

	signed ep      :32;	/* energy source point number */

	signed cdp     :32;	/* CDP ensemble number */

	signed cdpt    :32;	/* trace number within CDP ensemble */

	signed trid    :16;	/* trace identification code:
			1 = seismic data
			2 = dead
			3 = dummy
			4 = time break
			5 = uphole
			6 = sweep
			7 = timing
			8 = water break
			9---, N = optional use (N = 32,767)

			Following are CWP id flags:

			 9 = autocorrelation

			10 = Fourier transformed - no packing
			     xr[0],xi[0], ..., xr[N-1],xi[N-1]

			11 = Fourier transformed - unpacked Nyquist
			     xr[0],xi[0],...,xr[N/2],xi[N/2]

			12 = Fourier transformed - packed Nyquist
	 		     even N:
			     xr[0],xr[N/2],xr[1],xi[1], ...,
				xr[N/2 -1],xi[N/2 -1]
				(note the exceptional second entry)
			     odd N:
			     xr[0],xr[(N-1)/2],xr[1],xi[1], ...,
				xr[(N-1)/2 -1],xi[(N-1)/2 -1],xi[(N-1)/2]
				(note the exceptional second & last entries)

			13 = Complex signal in the time domain
			     xr[0],xi[0], ..., xr[N-1],xi[N-1]

			14 = Fourier transformed - amplitude/phase
			     a[0],p[0], ..., a[N-1],p[N-1]

			15 = Complex time signal - amplitude/phase
			     a[0],p[0], ..., a[N-1],p[N-1]

			16 = Real part of complex trace from 0 to Nyquist

			17 = Imag part of complex trace from 0 to Nyquist

			18 = Amplitude of complex trace from 0 to Nyquist

			19 = Phase of complex trace from 0 to Nyquist

			21 = Wavenumber time domain (k-t)

			22 = Wavenumber frequency (k-omega)

			23 = Envelope of the complex time trace

			24 = Phase of the complex time trace

			25 = Frequency of the complex time trace

			30 = Depth-Range (z-x) traces

			101 = Seismic data packed to bytes (by supack1)
			
			102 = Seismic data packed to 2 bytes (by supack2)
			*/

	signed nvs    :16;   /* number of vertically summed traces (see vscode
			   in bhed structure) */

	signed nhs    :16;   /* number of horizontally summed traces (see vscode
			   in bhed structure) */

	signed duse   :16;   /* data use:
				1 = production
				2 = test */

	signed offset :32; /* distance from source point to receiver
			   group (negative if opposite to direction
			   in which the line was shot) */

	signed gelev  :32; /* receiver group elevation from sea level
			   (above sea level is positive) */

	signed selev  :32; /* source elevation from sea level
			   (above sea level is positive) */

	signed sdepth :32; /* source depth (positive) */

	signed gdel   :32; /* datum elevation at receiver group */

	signed sdel   :32; /* datum elevation at source */

	signed swdep  :32; /* water depth at source */

	signed gwdep  :32; /* water depth at receiver group */

	signed scalel :16; /* scale factor for previous 7 entries
			   with value plus or minus 10 to the
			   power 0, 1, 2, 3, or 4 (if positive,
			   multiply, if negative divide) */

	signed scalco :16; /* scale factor for next 4 entries
			   with value plus or minus 10 to the
			   power 0, 1, 2, 3, or 4 (if positive,
			   multiply, if negative divide) */

	signed  sx    :32;   /* X source coordinate */

	signed  sy    :32;   /* Y source coordinate */

	signed  gx    :32;   /* X group coordinate */

	signed  gy    :32;   /* Y group coordinate */

	signed counit :16;   /* coordinate units code:
				for previous four entries
				1 = length (meters or feet)
				2 = seconds of arc (in this case, the
				X values are longitude and the Y values
				are latitude, a positive value designates
				the number of seconds east of Greenwich
				or north of the equator */

	signed wevel  :16;	/* weathering velocity */

	signed swevel :16;	/* subweathering velocity */

	signed sut    :16;	/* uphole time at source */

	signed gut    :16;	/* uphole time at receiver group */

	signed sstat  :16;	/* source static correction */

	signed gstat  :16;	/* group static correction */

	signed tstat  :16;	/* total static applied */

	signed laga   :16; /* lag time A, time in ms between end of 240-
			   byte trace identification header and time
			   break, positive if time break occurs after
			   end of header, time break is defined as
			   the initiation pulse which maybe recorded
			   on an auxiliary trace or as otherwise
			   specified by the recording system */

	signed lagb   :16; /* lag time B, time in ms between the time break
			   and the initiation time of the energy source,
			   may be positive or negative */

	signed delrt  :16; /* delay recording time, time in ms between
			   initiation time of energy source and time
			   when recording of data samples begins
			   (for deep water work if recording does not
			   start at zero time) */

	signed muts   :16; /* mute time--start */

	signed mute   :16; /* mute time--end */

	unsigned ns   :16; /* number of samples in this trace */

	unsigned dt   :16; /* sample interval; in micro-seconds */

	signed gain   :16; /* gain type of field instruments code:
				1 = fixed
				2 = binary
				3 = floating point
				4 ---- N = optional use */

	signed igc    :16; /* instrument gain constant */

	signed igi    :16; /* instrument early or initial gain */

	signed corr   :16; /* correlated:
				1 = no
				2 = yes */

	signed sfs    :16; /* sweep frequency at start */

	signed sfe    :16; /* sweep frequency at end */

	signed slen   :16; /* sweep length in ms */

	signed styp   :16; /* sweep type code:
				1 = linear
				2 = cos-squared
				3 = other */

	signed stas   :16; /* sweep trace length at start in ms */

	signed stae   :16; /* sweep trace length at end in ms */

	signed tatyp  :16; /* taper type: 1=linear, 2=cos^2, 3=other */

	signed afilf  :16; /* alias filter frequency if used */

	signed afils  :16; /* alias filter slope */

	signed nofilf :16; /* notch filter frequency if used */

	signed nofils :16; /* notch filter slope */

	signed lcf    :16; /* low cut frequency if used */

	signed hcf    :16; /* high cut frequncy if used */

	signed lcs    :16; /* low cut slope */

	signed hcs    :16; /* high cut slope */

	signed year   :16; /* year data recorded */

	signed day    :16; /* day of year */

	signed hour   :16; /* hour of day (24 hour clock) */

	signed minute :16; /* minute of hour */

	signed sec    :16; /* second of minute */

	signed timbas :16; /* time basis code:
				1 = local
				2 = GMT
				3 = other */

	signed trwf   :16; /* trace weighting factor, defined as 1/2^N
			   volts for the least sigificant bit */

	signed grnors :16; /* geophone group number of roll switch
			   position one */

	signed grnofr :16; /* geophone group number of trace one within
			   original field record */

	signed grnlof :16; /* geophone group number of last trace within
			   original field record */

	signed gaps   :16;  /* gap size (total number of groups dropped) */

	signed otrav  :16;  /* overtravel taper code:
				1 = down (or behind)
				2 = up (or ahead) */

	/* local assignments */
#ifdef NEWSU
        unsigned machid :32; /* double word alignment needed on the Cray */
#endif

	float d1;	/* sample spacing for non-seismic data */

	float f1;	/* first sample location for non-seismic data */

	float d2;	/* sample spacing between traces */

	float f2;	/* first trace location */

	float ungpow;	/* negative of power used for dynamic
			   range compression */

	float unscale;	/* reciprocal of scaling factor to normalize
			   range */
	signed ntr   :32;   /* number of traces */

	signed mark  :16;   /* mark selected traces */

#ifdef NEWSU
        signed unass :16;   /* unassigned values */

#ifndef FLOAT64
	signed unass1 :16;   /* unassigned values */
        signed unass2 :16;   /* unassigned values */
        signed unass3 :16;   /* unassigned values */
        signed unass4 :16;   /* unassigned values */
        signed unass5 :16;   /* unassigned values */
        signed unass6 :16;   /* unassigned values */
        signed unass7 :16;   /* unassigned values */
        signed unass8 :16;   /* unassigned values */
        signed unass9 :16;   /* unassigned values */
        signed unass10 :16;   /* unassigned values */
        signed unass11 :16;   /* unassigned values */
        signed unass12 :16;   /* unassigned values */
#endif
#else
	short unass[15];
#endif

} segyhdr;

/* The following refer to the trid field in segy.h		*/
/* CHARPACK represents byte packed seismic data from supack1	*/
#define		CHARPACK	101
/* SHORTPACK represents 2 byte packed seismic data from supack2	*/
#define		SHORTPACK	102

/* TREAL represents real time traces 				*/
#define		TREAL		1
/* TDEAD represents dead time traces 				*/
#define		TDEAD		2
/* TDUMMY represents dummy time traces 				*/
#define		TDUMMY		3
/* TBREAK represents time break traces 				*/
#define		TBREAK		4
/* UPHOLE represents uphole traces 				*/
#define		UPHOLE		5
/* SWEEP represents sweep traces 				*/
#define		SWEEP		6
/* TIMING represents timing traces 				*/
#define		TIMING		7
/* WBREAK represents timing traces 				*/
#define		WBREAK		8

/* TCMPLX represents complex time traces 			*/
#define		TCMPLX		13
/* TAMPH represents time domain data in amplitude/phase form	*/
#define		TAMPH		15
/* FPACK represents packed frequency domain data 		*/
#define		FPACK		12
/* FUNPACKNYQ represents complex frequency domain data 		*/
#define		FUNPACKNYQ	11
/* FCMPLX represents complex frequency domain data 		*/
#define		FCMPLX		10
/* FAMPH represents freq domain data in amplitude/phase form	*/
#define		FAMPH		14
/* REALPART represents the real part of a trace to Nyquist	*/
#define		REALPART	16
/* IMAGPART represents the real part of a trace to Nyquist	*/
#define		IMAGPART	17
/* AMPLITUDE represents the amplitude of a trace to Nyquist	*/
#define		AMPLITUDE	18
/* PHASE represents the phase of a trace to Nyquist		*/
#define		PHASE		19
/* KT represents wavenumber-time domain data 			*/
#define		KT		21
/* KOMEGA represents wavenumber-frequency domain data		*/
#define		KOMEGA		22
/* ENVELOPE represents the envelope of the complex time trace	*/
#define		ENVELOPE	23
/* INSTPHASE represents the phase of the complex time trace	*/
#define		INSTPHASE	24
/* INSTFREQ represents the frequency of the complex time trace	*/
#define		INSTFREQ	25
/* DEPTH represents traces in depth-range (z-x)			*/
#define		TRID_DEPTH	30
/* 3C data...  v,h1,h2=(11,12,13)+32 so a bitmask will convert  */
/* between conventions */
/* TVERT represents the vertical component */
#define     TVERT          43
/* TVHOZ1 represents the horizontal-1 component */
#define     THORZ1         44
/* TVHOZ2 represents the horizontal-2 component */
#define     THORZ2         45
/* TRADIAL represents the radial component */
#define     TRADIAL        46
/* TTRANS represents the transverse component */
#define     TTRANS         47

#endif /* SEGYHDR_H */
