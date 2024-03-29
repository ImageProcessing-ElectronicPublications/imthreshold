//  Zlib license
//
// ImThreshold header library.
//
//  Copyright (C) 2017:
//  zvezdochiot <zvezdochiot@user.sourceforge.net>
//  Homepage: https://sourceforge.net/projects/imthreshold/

#include <inttypes.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <memory.h>

#ifndef IMTHRESHOLD_H
#define IMTHRESHOLD_H
#define MAXVAL 256

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

#ifndef PI
#define PI  ((float)3.14159265358979323846264338327950288419716939937510)
#endif
#ifndef RADGRD
#define RADGRD ((float)57.295779513082320876798154814105)
#endif
#ifndef ABS
#define ABS(a)    ((a) < 0 ? (-(a)) : (a))
#endif
#ifndef MIN
#define MIN(a,b)  ((a) < (b) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a,b)  ((a) > (b) ? (a) : (b))
#endif
#ifndef TRIM
#define TRIM(x,a,b) (MIN(FA_MAX((x),(a)),(b)))
#endif

typedef uint8_t BYTE;
typedef uint16_t WORD;

typedef struct
{
    BYTE c[3];
    WORD s;
}
IMTpixel;

typedef struct
{
    float c[3];
    unsigned long int n;
}
IMTcluster;

typedef struct
{
    unsigned int width;
    unsigned int height;
    unsigned int bpp;
    unsigned int min;
    unsigned int max;
    float mid;
    float mean;
    float std;
    float wb;
}
IMTinfo;

enum Scaler
{
    SCALER_NEAREST,
    SCALER_BICONT,
    SCALER_BICUBIC,
    SCALER_BILINE,
    SCALER_BIAKIMA,
    SCALER_CROP,
    SCALER_GSAMPLE
};

////////////////////////////////////////////////////////////////////////////////

BYTE ByteClamp(int);
WORD Byte3Clamp(int);
unsigned IndexClamp(int, unsigned);
IMTpixel IMTset (BYTE, BYTE, BYTE);
IMTpixel IMTcalcS (IMTpixel);
IMTpixel IMTccorS (IMTpixel);
bool IMTequal (IMTpixel, IMTpixel);
IMTpixel IMTcopy (IMTpixel);
IMTpixel IMTnear (IMTpixel, IMTpixel, IMTpixel);
IMTpixel** IMTalloc (unsigned, unsigned);
void IMTfree (IMTpixel**, unsigned);
BYTE** BWalloc (unsigned, unsigned);
void BWfree (BYTE**, unsigned);
WORD** TLalloc (unsigned, unsigned);
void TLfree (WORD**, unsigned);
IMTpixel IMTdiffS (IMTpixel);
IMTpixel IMTrefilter1p (IMTpixel, IMTpixel);
IMTpixel IMTinterpolation (IMTpixel**, unsigned, unsigned, float, float);
BYTE IMTmax (IMTpixel**, unsigned, unsigned);
BYTE IMTmin (IMTpixel**, unsigned, unsigned);
float IMTmean (IMTpixel**, unsigned, unsigned);
float IMTdev (IMTpixel**, float, unsigned, unsigned);
float IMTwb (IMTpixel**, float, unsigned, unsigned);
int IMTwbauto (IMTpixel**, unsigned, unsigned, int);
void IMTHist (IMTpixel**, unsigned long long*, unsigned, unsigned, unsigned, unsigned, unsigned);
int IMTHistBiMod (unsigned long long*, unsigned, float);
IMTpixel IMTRGBtoRYB4 (IMTpixel, int);
IMTpixel IMTRGBtoYCbCr (IMTpixel, int);
IMTpixel IMTRGBtoHSV (IMTpixel, int);
void IMTFilterRGBtoRYB4 (IMTpixel**, unsigned, unsigned, int);
void IMTFilterRGBtoYCbCr (IMTpixel**, unsigned, unsigned, int);
void IMTFilterRGBtoHSV (IMTpixel**, unsigned, unsigned, int);
char* IMTFilterRGBtoCSP (IMTpixel**, unsigned, unsigned, char*, int);
void IMTFilterSCCor (IMTpixel**, unsigned, unsigned);
void IMTFilterSMirror (IMTpixel**, unsigned, unsigned);
int IMTFilterSNorm (IMTpixel**, unsigned, unsigned);
void IMTFilterSCScale (IMTpixel**, unsigned, unsigned, float, IMTpixel);
void IMTFilterSEdge (IMTpixel**, IMTpixel**, unsigned, unsigned);
void IMTFilterSNegate (IMTpixel**, unsigned, unsigned);
void IMTFilterSMean (IMTpixel**, IMTpixel**, unsigned, unsigned, float);
void IMTFilterSSelect (IMTpixel**, unsigned, unsigned, unsigned);
void IMTFilterInvertBW (BYTE**, unsigned, unsigned);
unsigned IMTdist (IMTpixel, IMTpixel);
float IMTdist3c2p (IMTpixel, float*);
IMTpixel IMTmeanIc (IMTpixel**, unsigned, unsigned, unsigned, unsigned);
IMTpixel IMTmaxIc (IMTpixel**, unsigned, unsigned, unsigned, unsigned);
IMTpixel IMTminIc (IMTpixel**, unsigned, unsigned, unsigned, unsigned);
IMTpixel IMTaverageIc (IMTpixel**, IMTpixel, unsigned, unsigned, unsigned, unsigned, float);
void IMTFilterDespeck2 (BYTE**, unsigned, unsigned, unsigned);
void IMTFilterDHatch (BYTE**, unsigned, unsigned, unsigned);
unsigned IMTFilterDMag2 (BYTE**, unsigned, unsigned, unsigned);
void IMTFilterDNeuro2 (BYTE**, unsigned, unsigned, unsigned, float, unsigned);
float IMTFilterDphist (BYTE**, unsigned, unsigned, unsigned);
void IMTFilterDMinMax (BYTE**, unsigned, unsigned, unsigned);
void IMTFilterDSmearing (BYTE**, unsigned, unsigned, unsigned);
void IMTFilterDWGrid (BYTE**, unsigned int, unsigned int, unsigned int);
float IMTFilterDeNoiseDiff1p (IMTpixel**, unsigned, unsigned, float);
float IMTFilterDeNoiseDiff (IMTpixel**, unsigned, unsigned, unsigned, float);
void IMTBlurMask (IMTpixel**, BYTE**, unsigned, unsigned, int);
void IMTReduceBW (BYTE**, BYTE**, unsigned, unsigned, unsigned, unsigned, unsigned, BYTE, BYTE);
IMTpixel IMTFilterIllumCorr (IMTpixel**, IMTpixel**, IMTpixel**, unsigned, unsigned);
float IMTFilterFindSkew (BYTE**, unsigned, unsigned);
void IMTFilterRotate (IMTpixel**, IMTpixel**, unsigned, unsigned, float);
IMTinfo IMTFilterInfo (IMTpixel**, unsigned, unsigned, unsigned);
void IMTFilterInvert (IMTpixel**, unsigned, unsigned);
float IMTFilterIRange (IMTpixel**, unsigned, unsigned, int, float);
void IMTFilterCopy (IMTpixel**, IMTpixel**, unsigned, unsigned);
void IMTFilterAdSmooth (IMTpixel**, IMTpixel**, unsigned, unsigned, float);
void IMTFilterAutoLevel (IMTpixel**, IMTpixel**, unsigned int, unsigned int, unsigned int);
void IMTFilterAutoWhite (IMTpixel**, IMTpixel**, unsigned, unsigned, float);
void IMTFilterInpaint (IMTpixel**, BYTE**, IMTpixel**, unsigned, unsigned, unsigned);
void IMTFilterSeparate (IMTpixel**, BYTE**, IMTpixel**, unsigned, unsigned, unsigned);
void IMTFilterSeparateBGFGL (IMTpixel**, BYTE**, IMTpixel**, IMTpixel**, unsigned, unsigned, unsigned, unsigned, unsigned, float);
void IMTFilterSeparateBGFGC (IMTpixel**, BYTE**, IMTpixel**, IMTpixel**, unsigned, unsigned, unsigned, unsigned, unsigned, float);
void IMTFilterSeparateDelta (IMTpixel**, BYTE**, IMTpixel**, unsigned, unsigned, int, float);
void IMTFilterEqualize (IMTpixel**, unsigned, unsigned);
IMTpixel IMTFilterGreyWorld (IMTpixel**, unsigned, unsigned);
IMTpixel IMTFilterGreyNorm (IMTpixel**, unsigned, unsigned);
float IMTFilterLevelMean (IMTpixel**, IMTpixel**, unsigned, unsigned, int, float, float, int, int);
float IMTFilterLevelSigma (IMTpixel**, IMTpixel**, unsigned, unsigned, float, float);
float IMTFilterLevelSize (IMTpixel**, IMTpixel**, unsigned, unsigned, int, int);
void IMTFilterMathAverage (IMTpixel**, IMTpixel**, unsigned, unsigned, int);
void IMTFilterMathBlur (IMTpixel**, IMTpixel**, unsigned, unsigned, unsigned, int);
void IMTFilterMathDistance (IMTpixel**, IMTpixel**, unsigned, unsigned, int);
void IMTFilterMathDivide (IMTpixel**, IMTpixel**, unsigned, unsigned, int);
void IMTFilterMathGeometric (IMTpixel**, IMTpixel**, unsigned, unsigned, int);
void IMTFilterMathHarmonic (IMTpixel**, IMTpixel**, unsigned, unsigned, int);
void IMTFilterMathMax (IMTpixel**, IMTpixel**, unsigned, unsigned, int);
void IMTFilterMathMin (IMTpixel**, IMTpixel**, unsigned, unsigned, int);
void IMTFilterMathMinus (IMTpixel**, IMTpixel**, unsigned, unsigned, int);
void IMTFilterMathMirror (IMTpixel**, IMTpixel**, unsigned, unsigned, int);
void IMTFilterMathMultiply (IMTpixel**, IMTpixel**, unsigned, unsigned, int);
void IMTFilterMathNorm (IMTpixel**, IMTpixel**, unsigned, unsigned, int);
void IMTFilterMathOverlay (IMTpixel**, IMTpixel**, unsigned int, unsigned int, int);
void IMTFilterMathPlus (IMTpixel**, IMTpixel**, unsigned, unsigned, int);
void IMTFilterMathScreen (IMTpixel**, IMTpixel**, unsigned int, unsigned int, int);
float IMTFilterMathSharpenBadMetric (IMTpixel**, IMTpixel**, unsigned, unsigned);
void IMTFilterMathSqr (IMTpixel**, unsigned, unsigned);
void IMTFilterMathThreshold (IMTpixel**, IMTpixel**, BYTE**, unsigned, unsigned, int);
float IMTFilterMirror (IMTpixel**, IMTpixel**, unsigned, unsigned);
float IMTFilterMirrorPart (IMTpixel**, IMTpixel**, unsigned, unsigned, float);
void IMTFilterMirrorMean (IMTpixel**, unsigned, unsigned);
void IMTFilterMonoColor (IMTpixel**, unsigned, unsigned);
void IMTFilterMorph (IMTpixel**, IMTpixel**, unsigned, unsigned, int, bool);
IMTpixel IMTmeanIcM (IMTpixel**, bool**, unsigned, unsigned, unsigned, unsigned, unsigned, unsigned);
IMTpixel IMTmeanPtIcM (IMTpixel**, IMTpixel, float*, bool**, unsigned, unsigned, unsigned, unsigned, unsigned, unsigned);
float IMTwbIcM (IMTpixel**, IMTpixel, bool**, unsigned, unsigned, unsigned, unsigned, unsigned, unsigned);
IMTpixel IMTmeanThIcM (IMTpixel**, IMTpixel, bool**, unsigned, unsigned, unsigned, unsigned, unsigned, unsigned, int);
IMTpixel IMTmeanRadIcM (IMTpixel**, IMTpixel, float**, bool**, unsigned, unsigned, unsigned, unsigned, unsigned, unsigned);
IMTpixel IMTmeanBlIcM (IMTpixel**, IMTpixel, float**, float*, bool**, unsigned, unsigned, unsigned, unsigned, unsigned, unsigned);
IMTpixel IMTmeanMinIcM (IMTpixel**, IMTpixel, bool**, unsigned, unsigned, unsigned, unsigned, unsigned, unsigned);
IMTpixel IMTmeanMaxIcM (IMTpixel**, IMTpixel, bool**, unsigned, unsigned, unsigned, unsigned, unsigned, unsigned);
unsigned IMTFilterKMeans (IMTpixel**, unsigned, unsigned, unsigned, unsigned);
void IMTFilterPeron (IMTpixel**, IMTpixel**, unsigned, unsigned, float, float);
float IMTFilterPosterize (IMTpixel**, unsigned, unsigned, unsigned);
float IMTFilterQuant (IMTpixel**, unsigned, unsigned, unsigned);
void IMTFilterPMean (IMTpixel**, IMTpixel**, unsigned, unsigned, float, int, bool);
int IMTFilterRetinex (IMTpixel**, IMTpixel**, unsigned, unsigned, int, float);
float IMTFilterRS (IMTpixel**, unsigned, unsigned);
void IMTSelGaussInitMatrix (float, float*, int);
void IMTFilterSelGauss (IMTpixel**, IMTpixel**, unsigned, unsigned, float, int);
float IMTFilterShrink (IMTpixel**, unsigned, unsigned, int);
void IMTFilterSobel (IMTpixel**, IMTpixel**, unsigned, unsigned);
void IMTFilterUnRipple (IMTpixel**, IMTpixel**, unsigned, unsigned, int, float);
void IMTGaussLineMatrix (float*, float);
void IMTFilterGaussBlur (IMTpixel**, IMTpixel**, unsigned, unsigned, float);
int IMTFilterClusterBiModC (IMTpixel**, IMTpixel**, unsigned, unsigned, unsigned);
float IMTFilterClusterBWC (IMTpixel**, IMTpixel**, unsigned, unsigned, unsigned);
void IMTFilterUnsharpMask (IMTpixel**, IMTpixel**, unsigned, unsigned, float, int);
float IMTFilterNoiseVariance (IMTpixel**, unsigned, unsigned, int);
void IMTFilterWiener (IMTpixel**, IMTpixel**, unsigned, unsigned, int, float);
int IMTFilterWhiteFill (IMTpixel**, unsigned, unsigned);
void IMTFilterReverse (IMTpixel**, IMTpixel**, unsigned int, unsigned int);
IMTpixel IMTInterpolateBiCubic (IMTpixel**, int, int, float, float);
IMTpixel IMTInterpolateBiLine (IMTpixel**, int, int, float, float);
float IMTInterpolateAkima (float*, float);
IMTpixel IMTInterpolateBiAkima (IMTpixel**, int, int, float, float);
void IMTFilterSBicub (IMTpixel**, IMTpixel**, unsigned, unsigned, unsigned, unsigned);
void IMTFilterSBicont (IMTpixel**, IMTpixel**, unsigned, unsigned, unsigned, unsigned);
void IMTFilterSBilin (IMTpixel**, IMTpixel**, unsigned, unsigned, unsigned, unsigned);
void IMTFilterSBiakima (IMTpixel**, IMTpixel**, unsigned, unsigned, unsigned, unsigned);
void IMTFilterSCrop (IMTpixel**, IMTpixel**, unsigned, unsigned, unsigned, unsigned);
void IMTFilterSize (IMTpixel**, IMTpixel**, int, unsigned, unsigned, unsigned, unsigned);
int IMTFilterSBWMag2 (BYTE**, BYTE**, unsigned, unsigned, unsigned, unsigned);
int IMTFilterSBWReduce2 (BYTE**, BYTE**, unsigned, unsigned, unsigned, unsigned);
void IMTFilterSGsample (IMTpixel**, IMTpixel**, unsigned, unsigned, unsigned, unsigned);
void IMTFilterSHRIS (IMTpixel**, IMTpixel**, unsigned, unsigned, int);
void IMTFilterSGSampleUp (IMTpixel**, IMTpixel**, unsigned, unsigned, int);
void IMTFilterSReduce (IMTpixel**, IMTpixel**, unsigned, unsigned, unsigned);
void IMTFilterSNearest (IMTpixel**, IMTpixel**, unsigned, unsigned, unsigned, unsigned);
void IMTFilterSScaleNX (IMTpixel**, IMTpixel**, unsigned int, unsigned int, int);
void IMTFilterSFRP (IMTpixel** p_im, IMTpixel** d_im, unsigned height, unsigned width, int smode);
int IMTFilterThreshold (IMTpixel**, BYTE**, unsigned, unsigned, int);
int IMTFilterThresholdLayer (IMTpixel**, WORD**, BYTE**, unsigned, unsigned);
void IMTFilterTLayerToImg (WORD**, IMTpixel**, unsigned, unsigned);
int IMTFilterTAbutaleb (IMTpixel**, BYTE**, unsigned, unsigned, int);
int IMTFilterTBernsenLayer(IMTpixel**, WORD**, unsigned, unsigned, int, unsigned);
int IMTFilterTBernsen (IMTpixel**, BYTE**, unsigned, unsigned, int, unsigned);
int IMTFilterTBHTValue (IMTpixel**, unsigned, unsigned);
int IMTFilterTBHT (IMTpixel**, BYTE**, unsigned, unsigned);
int IMTFilterTBiModValueBound (IMTpixel**, unsigned, unsigned, int, int);
int IMTFilterTBiModValueIcP (IMTpixel**, unsigned, unsigned, unsigned, unsigned, float);
int IMTFilterTBiModValueIc (IMTpixel**, unsigned, unsigned, unsigned, unsigned);
int IMTFilterTBiModValueP (IMTpixel**, unsigned, unsigned, float);
int IMTFilterTBiModValue (IMTpixel**, unsigned, unsigned);
int IMTFilterTBiMod (IMTpixel**, BYTE**, unsigned, unsigned, int);
int IMTFilterTBiModP (IMTpixel**, IMTpixel**, unsigned, unsigned, unsigned);
int IMTFilterTBiModC (IMTpixel**, BYTE**, unsigned, unsigned, int);
int IMTFilterTBiModLayer (IMTpixel**, WORD**, unsigned, unsigned, int, float, int);
int IMTFilterTBiModRegion (IMTpixel**, BYTE**, unsigned, unsigned, int, float, int);
int IMTFilterTBlurDiv (IMTpixel**, BYTE**, unsigned int, unsigned int, int, float, float);
int IMTFilterTColor (IMTpixel**, BYTE**, unsigned, unsigned, int);
int IMTFilterTChistianLayer (IMTpixel**, WORD**, unsigned, unsigned, int, float, int, int, float);
int IMTFilterTChistian (IMTpixel**, BYTE**, unsigned, unsigned, int, float, int, int, float);
int IMTFilterTDalg (IMTpixel**, BYTE**, unsigned, unsigned, int, int, int, int);
int IMTFilterTDither (IMTpixel**, BYTE**, unsigned, unsigned);
int IMTFilterTDithH (IMTpixel**, BYTE**, unsigned, unsigned, int, int, unsigned);
int IMTFilterTDithO (IMTpixel**, BYTE**, unsigned, unsigned, int, int);
int IMTFilterTDithBayer (IMTpixel**, BYTE**, unsigned, unsigned, int);
int IMTFilterTDithDots (IMTpixel**, BYTE**, unsigned, unsigned, int);
int IMTFilterTDjVuL (IMTpixel**, BYTE**, IMTpixel**, IMTpixel**, unsigned, unsigned, unsigned, unsigned, unsigned, int, float, float, unsigned);
int IMTFilterTEdge (IMTpixel**, BYTE**, unsigned int, unsigned int, int, float);
int IMTFilterTEdgeDiv (IMTpixel**, BYTE**, unsigned int, unsigned int, int, float, float, float);
int IMTFilterTEdgePlus (IMTpixel**, BYTE**, unsigned int, unsigned int, int, float, float);
int IMTFilterTEntValue (IMTpixel**, unsigned, unsigned);
int IMTFilterTEnt (IMTpixel**, BYTE**, unsigned, unsigned);
int IMTFilterTEqBrightValue (IMTpixel**, unsigned, unsigned);
int IMTFilterTEqBright (IMTpixel**, BYTE**, unsigned, unsigned);
int IMTFilterGatosBG (IMTpixel**, BYTE**, IMTpixel**, unsigned, unsigned, int);
int IMTFilterTGatos (IMTpixel**, BYTE**, IMTpixel**, BYTE**, unsigned, unsigned, float, float, float);
int IMTFilterTGradValue (IMTpixel**, unsigned, unsigned);
int IMTFilterTGrad (IMTpixel**, BYTE**, unsigned, unsigned);
int IMTFilterTGravureLayer (IMTpixel**, WORD**, unsigned, unsigned, int, float, int, int, float);
int IMTFilterTGravure (IMTpixel**, BYTE**, unsigned, unsigned, int, float, int, int, float);
int IMTFilterTHalftone2 (IMTpixel**, BYTE**, unsigned, unsigned);
int IMTFilterTJanniValue (IMTpixel**, unsigned, unsigned);
int IMTFilterTJanni (IMTpixel**, BYTE**, unsigned, unsigned);
int IMTFilterTKMeans (IMTpixel**, BYTE**, unsigned, unsigned, unsigned, unsigned);
int IMTFilterTMscaleLayer (IMTpixel**, WORD**, unsigned, unsigned, unsigned, float, float, float);
int IMTFilterTMscale (IMTpixel**, BYTE**, unsigned, unsigned, unsigned, float, float, float);
int IMTFilterTNiblackLayer (IMTpixel**, WORD**, unsigned, unsigned, int, float, int, int, float);
int IMTFilterTNiblack (IMTpixel**, BYTE**, unsigned, unsigned, int, float, int, int, float);
int IMTFilterTOtsuValue (IMTpixel**, unsigned, unsigned);
int IMTFilterTOtsu (IMTpixel**, BYTE**, unsigned, unsigned);
int IMTFilterTQuadModValue (IMTpixel**, unsigned, unsigned);
int IMTFilterTQuadMod (IMTpixel**, BYTE**, unsigned, unsigned, int);
int IMTFilterTRotValue (IMTpixel**, unsigned, unsigned, bool);
int IMTFilterTRot (IMTpixel**, BYTE**, unsigned, unsigned, bool);
int IMTFilterTSauvolaLayer (IMTpixel**, WORD**, unsigned, unsigned, int, float, int, int, int, float);
int IMTFilterTSauvola (IMTpixel**, BYTE**, unsigned, unsigned, int, float, int, int, int, float);
int IMTFilterTSizeLayer (IMTpixel**, WORD**, unsigned, unsigned, int, int, int, int);
int IMTFilterTSize (IMTpixel**, BYTE**, unsigned, unsigned, int, int, int, int);
int IMTFilterTText (IMTpixel**, BYTE**, unsigned, unsigned, unsigned, unsigned);
int IMTFilterTTsaiValue (IMTpixel**, unsigned, unsigned, int);
int IMTFilterTTsai (IMTpixel**, BYTE**, unsigned, unsigned, int);
int IMTFilterTWhiteRohrer (IMTpixel**, BYTE**, unsigned, unsigned, int, int, int, int, int, int);

#ifdef __cplusplus
} // extern "C"
#endif // _\cplusplus

#endif // IMTHRESHOLD_H
