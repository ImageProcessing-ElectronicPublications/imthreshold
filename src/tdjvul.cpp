//    Zlib license
//
// DjVuL thresholds an image (Multi-scale binarization).
//
//    Copyright (C) 2017:
//    zvezdochiot    <zvezdochiot@user.sourceforge.net>

#include <unistd.h>
#include <FreeImage.h>
#include "imthresholdfreeimage.h"

////////////////////////////////////////////////////////////////////////////////

void ImthresholdIMTFilterDjVuLTitle()
{
    printf("ImThreshold.\n");
    printf("BookScanLib Project: http://djvu-soft.narod.ru/\n\n");
    printf("DjVuL thresholds an image (Multi-scale binarization).\n");
    printf("Homepage: https://sourceforge.net/projects/imthreshold/.\n\n");
}

void ImthresholdIMTFilterDjVuLUsage()
{
    printf("Usage : imthreshold-tdjvul [options] <input_file> <output_file>(BW) [fg_file] [bg_file] [fg_mask] [bg_mask]\n\n");
    printf("options:\n");
    printf("          -a N.N  anisotropic (float, optional, default = 0.0)\n");
    printf("          -b N    base block size (int, optional, default = 3)\n");
    printf("          -c N    clean fg, bg blur radius (int, optional, default = 0)\n");
    printf("          -d N    despeckle aperture size (int, optional, default = 0)\n");
    printf("          -f N    foreground divide (int, optional, default = 2)\n");
    printf("          -i      invert station (bool, optional, default = false)\n");
    printf("          -l N    level (int, optional, default = 0 [auto])\n");
    printf("          -m str  name metod:\n");
    printf("                    'djvuc'\n");
    printf("                    'djvul (default)'\n");
    printf("          -o N.N  overlay (float, optional, default = 0.5)\n");
    printf("          -p N    posterize fg (int, optional, default = 0)\n");
    printf("          -q str  colorspace:\n");
    printf("                    'rgb' (default)\n");
    printf("                    'ryb4'\n");
    printf("                    'ycbcr'\n");
    printf("                    'hsv'\n");
    printf("          -w N    w/b mode (int, optional, default = 0 [auto], >0-white, <0-black)\n");
    printf("          -h      this help\n");
}

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
    // call this ONLY when linking with FreeImage as a static library
#ifdef FREEIMAGE_LIB
    FreeImage_Initialise();
#endif // FREEIMAGE_LIB

    int opt;
    float anisotropic = 0.0;
    float doverlay = 0.5;
    int bgs = 3;
    int fgs = 2;
    int level = 0;
    int wbmode = 0;
    int fclean = 0;
    int fdespeckle = 0;
    unsigned fposter = 0;
    bool finvs = false;
    bool fhelp = false;
    char *csp, *cspn;
    char *namefilter;
    namefilter = (char*)"djvul";
    csp = (char*)"rgb";
    while ((opt = getopt(argc, argv, ":a:b:c:d:f:il:m:o:p:q:w:h")) != -1)
    {
        switch(opt)
        {
        case 'a':
            anisotropic = atof(optarg);
            break;
        case 'b':
            bgs = atof(optarg);
            break;
        case 'c':
            fclean = atof(optarg);
            break;
        case 'd':
            fdespeckle = atof(optarg);
            break;
        case 'f':
            fgs = atof(optarg);
            break;
        case 'i':
            finvs = true;
            break;
        case 'l':
            level = atof(optarg);
            break;
        case 'm':
            namefilter = optarg;
            break;
        case 'o':
            doverlay = atof(optarg);
            break;
        case 'p':
            fposter = atof(optarg);
            break;
        case 'q':
            csp = optarg;
            break;
        case 'w':
            wbmode = atof(optarg);
            break;
        case 'h':
            fhelp = true;
            break;
        case ':':
            printf("option needs a value\n");
            break;
        case '?':
            printf("unknown option: %c\n", optopt);
            break;
        }
    }

    ImthresholdIMTFilterDjVuLTitle();

    if(optind + 2 > argc || fhelp || bgs < 1 || fgs < 1 || level < 0)
    {
        ImthresholdIMTFilterDjVuLUsage();
        return 0;
    }
    const char *src_filename = argv[optind];
    const char *output_filename = argv[optind + 1];
    const char *fg_filename;
    if(optind + 2 < argc)
    {
        fg_filename = argv[optind + 2];
    }
    const char *bg_filename;
    if(optind + 3 < argc)
    {
        bg_filename = argv[optind + 3];
    }
    const char *fg_maskname;
    if(optind + 4 < argc)
    {
        fg_maskname = argv[optind + 4];
    }
    const char *bg_maskname;
    if(optind + 5 < argc)
    {
        bg_maskname = argv[optind + 5];
    }

    FreeImage_SetOutputMessage(FreeImageErrorHandler);

    if (fclean < 0)
    {
        fclean = -fclean;
    }
    printf("Input= %s\n", src_filename);
    FIBITMAP *dib = ImthresholdGenericLoader(src_filename, 0);
    if (dib)
    {
        if (FreeImage_GetImageType(dib) == FIT_BITMAP)
        {
            FIBITMAP* dst_dib;
            FIBITMAP* fg_dib;
            FIBITMAP* bg_dib;
            FIBITMAP* fgm_dib;
            FIBITMAP* bgm_dib;
            if (FreeImage_GetBPP(dib) == 1)
            {
                dst_dib = ImthresholdFilterNone(dib);
            }
            else
            {
                printf("Level= %d\n", level);
                printf("Anisotropic= %f\n", anisotropic);
                printf("Overlay= %f\n", doverlay);
                unsigned width = FreeImage_GetWidth(dib);
                unsigned height = FreeImage_GetHeight(dib);
                unsigned widthbg = (width + bgs - 1) / bgs;
                unsigned heightbg = (height + bgs - 1) / bgs;
                unsigned widthfg = (widthbg + fgs - 1) / fgs;
                unsigned heightfg = (heightbg + fgs - 1) / fgs;

                IMTpixel** p_im = IMTalloc(height, width);
                BYTE** m_im = BWalloc(height, width);
                IMTpixel** fg_im = IMTalloc(heightfg, widthfg);
                IMTpixel** bg_im = IMTalloc(heightbg, widthbg);
                BYTE** fgm_im = BWalloc(heightfg, widthfg);
                BYTE** bgm_im = BWalloc(heightbg, widthbg);

                ImthresholdGetData(dib, p_im);
                FreeImage_Unload(dib);

                cspn = IMTFilterRGBtoCSP(p_im, height, width, csp, 1);
                printf("ColorSpace= %s\n", cspn);

                if (fposter != 0)
                {
                    printf("Posterize= %d\n", fposter);
                }
                if (strcmp(namefilter, "djvuc") == 0)
                {
                    printf("Filter= %s\n", namefilter);
                    wbmode = IMTwbauto (p_im, height, width, wbmode);
                    (void) IMTFilterTBiMod (p_im, m_im, height, width, 0);
                    if (wbmode < 0)
                    {
                        IMTFilterInvertBW(m_im, height, width);
                    }
                    IMTFilterSeparateBGFGC (p_im, m_im, fg_im, bg_im, height, width, bgs, fgs, level, doverlay);
                }
                else
                {
                    printf("Filter= djvul\n");
                    wbmode = IMTFilterTDjVuL(p_im, m_im, fg_im, bg_im, height, width, bgs, fgs, level, wbmode, anisotropic, doverlay, fposter);
                }
                if (wbmode > 0)
                {
                    printf("Mode= white\n");
                }
                else
                {
                    printf("Mode= black\n");
                }
                IMTfree(p_im, height);

                cspn = IMTFilterRGBtoCSP(fg_im, heightfg, widthfg, csp, -1);
                cspn = IMTFilterRGBtoCSP(bg_im, heightbg, widthbg, csp, -1);
                printf("ColorSpace= %s\n", cspn);

                if (fdespeckle > 0)
                {
                    printf("Despeckle= %d\n", fdespeckle);
                    IMTFilterDespeck2(m_im, height, width, fdespeckle);
                }
                IMTReduceBW(m_im, fgm_im, height, width, heightfg, widthfg, bgs * fgs, 0, 255);
                IMTReduceBW(m_im, bgm_im, height, width, heightbg, widthbg, bgs, 255, 0);
                dst_dib = FreeImage_Allocate(width, height, 1);
                if (finvs)
                {
                    printf("InvertStation= true\n");
                    IMTFilterInvertBW (m_im, height, width);
                }
                ImthresholdSetDataBW(dst_dib, m_im);
                BWfree(m_im, height);
                if (fclean > 0)
                {
                    printf("Blur= %d\n", fclean);
                    if (fposter == 0)
                    {
                        IMTBlurMask(fg_im, fgm_im, heightfg, widthfg, fclean);
                    }
                    IMTBlurMask(bg_im, bgm_im, heightbg, widthbg, fclean);
                }
                if (finvs)
                {
                    fg_dib = FreeImage_Allocate(widthbg, heightbg, 24);
                    ImthresholdSetData(fg_dib, bg_im);
                    bg_dib = FreeImage_Allocate(widthfg, heightfg, 24);
                    ImthresholdSetData(bg_dib, fg_im);
                }
                else
                {
                    fg_dib = FreeImage_Allocate(widthfg, heightfg, 24);
                    ImthresholdSetData(fg_dib, fg_im);
                    bg_dib = FreeImage_Allocate(widthbg, heightbg, 24);
                    ImthresholdSetData(bg_dib, bg_im);
                }
                IMTfree(fg_im, heightfg);
                IMTfree(bg_im, heightbg);
                if (finvs)
                {
                    fgm_dib = FreeImage_Allocate(widthbg, heightbg, 1);
                    ImthresholdSetDataBW(fgm_dib, bgm_im);
                    bgm_dib = FreeImage_Allocate(widthfg, heightfg, 1);
                    ImthresholdSetDataBW(bgm_dib, fgm_im);
                }
                else
                {
                    fgm_dib = FreeImage_Allocate(widthfg, heightfg, 1);
                    ImthresholdSetDataBW(fgm_dib, fgm_im);
                    bgm_dib = FreeImage_Allocate(widthbg, heightbg, 1);
                    ImthresholdSetDataBW(bgm_dib, bgm_im);
                }
                BWfree(fgm_im, heightfg);
                BWfree(bgm_im, heightbg);

            }
            if (dst_dib)
            {
                FREE_IMAGE_FORMAT out_fif = FreeImage_GetFIFFromFilename(output_filename);
                if(out_fif != FIF_UNKNOWN)
                {
                    FreeImage_Save(out_fif, dst_dib, output_filename, 0);
                    printf("Output= %s\n", output_filename);
                }
                FreeImage_Unload(dst_dib);
            }
            if (fg_dib)
            {
                if(optind + 2 < argc)
                {
                    FREE_IMAGE_FORMAT out_fif = FreeImage_GetFIFFromFilename(fg_filename);
                    if(out_fif != FIF_UNKNOWN)
                    {
                        FreeImage_Save(out_fif, fg_dib, fg_filename, 0);
                        printf("Output= %s\n", fg_filename);
                    }
                }
                FreeImage_Unload(fg_dib);
            }
            if (bg_dib)
            {
                if(optind + 3 < argc)
                {
                    FREE_IMAGE_FORMAT out_fif = FreeImage_GetFIFFromFilename(bg_filename);
                    if(out_fif != FIF_UNKNOWN)
                    {
                        FreeImage_Save(out_fif, bg_dib, bg_filename, 0);
                        printf("Output= %s\n", bg_filename);
                    }
                }
                FreeImage_Unload(bg_dib);
            }
            if (fgm_dib)
            {
                if(optind + 4 < argc)
                {
                    FREE_IMAGE_FORMAT out_fif = FreeImage_GetFIFFromFilename(fg_maskname);
                    if(out_fif != FIF_UNKNOWN)
                    {
                        FreeImage_Save(out_fif, fgm_dib, fg_maskname, 0);
                        printf("Output= %s\n", fg_maskname);
                    }
                }
                FreeImage_Unload(fgm_dib);
            }
            if (bgm_dib)
            {
                if(optind + 5 < argc)
                {
                    FREE_IMAGE_FORMAT out_fif = FreeImage_GetFIFFromFilename(bg_maskname);
                    if(out_fif != FIF_UNKNOWN)
                    {
                        FreeImage_Save(out_fif, bgm_dib, bg_maskname, 0);
                        printf("Output= %s\n", bg_maskname);
                    }
                }
                FreeImage_Unload(bgm_dib);
            }
            printf("\n");
        }
        else
        {
            printf("%s\n", "Unsupported format type.");
            FreeImage_Unload(dib);
        }
    }

    // call this ONLY when linking with FreeImage as a static library
#ifdef FREEIMAGE_LIB
    FreeImage_DeInitialise();
#endif // FREEIMAGE_LIB

    return 0;
}
