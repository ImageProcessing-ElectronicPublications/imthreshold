//    Zlib license
//
// BGFG separate based Multi-scale separate.
//
//    Copyright (C) 2017:
//    zvezdochiot    <zvezdochiot@user.sourceforge.net>

#include <unistd.h>
#include <FreeImage.h>
#include "imthresholdfreeimage.h"

////////////////////////////////////////////////////////////////////////////////

void ImthresholdFilterMathTitle()
{
    printf("ImThreshold.\n");
    printf("BookScanLib Project: http://djvu-soft.narod.ru/\n\n");
    printf("Math opertion of two image.\n");
    printf("Homepage: https://sourceforge.net/projects/imthreshold/.\n\n");
}

void ImthresholdFilterMathUsage()
{
    printf("Usage : imthreshold-math [options] <input_file> <input_file> <output_file>\n\n");
    printf("options:\n");
    printf("          -m str  name metod:\n");
    printf("                    'average'\n");
    printf("                    'distance'\n");
    printf("                    'divide'\n");
    printf("                    'sbmetric'\n");
    printf("                    'geometric'\n");
    printf("                    'harmonic'\n");
    printf("                    'max'\n");
    printf("                    'min'\n");
    printf("                    'minus'\n");
    printf("                    'mirror'\n");
    printf("                    'multiply'\n");
    printf("                    'norm'\n");
    printf("                    'plus'\n");
    printf("                    'threshold (default)'\n");
    printf("          -d N    delta (int, optional, default = 0)\n");
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
    int delta = 0;
    bool fhelp = false;
    char *namefilter;
    namefilter="threshold";
    while ((opt = getopt(argc, argv, ":m:d:h")) != -1)
    {
        switch(opt)
        {
            case 'm':
                namefilter = optarg;
                break;
            case 'd':
                delta = atof(optarg);
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

    ImthresholdFilterMathTitle ();

    if(optind + 3 > argc || fhelp)
    {
        ImthresholdFilterMathUsage ();
        return 0;
    }
    const char *src_filename = argv[optind];
    const char *math_filename = argv[optind + 1];
    const char *dst_filename = argv[optind + 2];

    FreeImage_SetOutputMessage (FreeImageErrorHandler);

    printf ("Input= %s\n", src_filename);
    FIBITMAP *dib = ImthresholdGenericLoader (src_filename, 0);
    if (dib)
    {
        if (FreeImage_GetImageType(dib) == FIT_BITMAP)
        {
            unsigned width = FreeImage_GetWidth (dib);
            unsigned height = FreeImage_GetHeight (dib);

            IMTpixel** p_im = IMTalloc(height, width);

            ImthresholdGetData(dib, p_im);
            FreeImage_Unload(dib);

            printf("Input= %s\n", math_filename);
            FIBITMAP *dib = ImthresholdGenericLoader(math_filename, 0);
            if (dib)
            {
                if (FreeImage_GetImageType(dib) == FIT_BITMAP)
                {
                    unsigned widthm = FreeImage_GetWidth (dib);
                    unsigned heightm = FreeImage_GetHeight (dib);

                    if ((height == heightm) && (width == widthm))
                    {
                        FIBITMAP* dst_dib;

                        IMTpixel** m_im = IMTalloc(height, width);

                        ImthresholdGetData(dib, m_im);
                        
                        if (strcmp(namefilter, "average") == 0)
                        {
                            printf("Filter= %s\n", namefilter);
                            printf("Delta= %d\n", delta);

                            IMTFilterMathAverage (p_im, m_im, height, width, delta);
                            dst_dib = FreeImage_Allocate (width, height, 24);
                            ImthresholdSetData (dst_dib, p_im);
                        }
                        else if (strcmp(namefilter, "distance") == 0)
                        {
                            printf("Filter= %s\n", namefilter);
                            printf("Delta= %d\n", delta);

                            IMTFilterMathDistance (p_im, m_im, height, width, delta);
                            dst_dib = FreeImage_Allocate (width, height, 24);
                            ImthresholdSetData (dst_dib, p_im);
                        }
                        else if (strcmp(namefilter, "divide") == 0)
                        {
                            printf("Filter= %s\n", namefilter);
                            printf("Delta= %d\n", delta);

                            IMTFilterMathDivide (p_im, m_im, height, width, delta);
                            dst_dib = FreeImage_Allocate (width, height, 24);
                            ImthresholdSetData (dst_dib, p_im);
                        }
                        else if (strcmp(namefilter, "geometric") == 0)
                        {
                            printf("Filter= %s\n", namefilter);
                            printf("Delta= %d\n", delta);

                            IMTFilterMathGeometric (p_im, m_im, height, width, delta);
                            dst_dib = FreeImage_Allocate (width, height, 24);
                            ImthresholdSetData (dst_dib, p_im);
                        }
                        else if (strcmp(namefilter, "harmonic") == 0)
                        {
                            printf("Filter= %s\n", namefilter);
                            printf("Delta= %d\n", delta);

                            IMTFilterMathHarmonic (p_im, m_im, height, width, delta);
                            dst_dib = FreeImage_Allocate (width, height, 24);
                            ImthresholdSetData (dst_dib, p_im);
                        }
                        else if (strcmp(namefilter, "max") == 0)
                        {
                            printf("Filter= %s\n", namefilter);
                            printf("Delta= %d\n", delta);

                            IMTFilterMathMax (p_im, m_im, height, width, delta);
                            dst_dib = FreeImage_Allocate (width, height, 24);
                            ImthresholdSetData (dst_dib, p_im);
                        }
                        else if (strcmp(namefilter, "min") == 0)
                        {
                            printf("Filter= %s\n", namefilter);
                            printf("Delta= %d\n", delta);

                            IMTFilterMathMin (p_im, m_im, height, width, delta);
                            dst_dib = FreeImage_Allocate (width, height, 24);
                            ImthresholdSetData (dst_dib, p_im);
                        }
                        else if (strcmp(namefilter, "minus") == 0)
                        {
                            printf("Filter= %s\n", namefilter);
                            printf("Delta= %d\n", delta);

                            IMTFilterMathMinus (p_im, m_im, height, width, delta);
                            dst_dib = FreeImage_Allocate (width, height, 24);
                            ImthresholdSetData (dst_dib, p_im);
                        }
                        else if (strcmp(namefilter, "mirror") == 0)
                        {
                            printf("Filter= %s\n", namefilter);
                            printf("Delta= %d\n", delta);

                            IMTFilterMathMirror (p_im, m_im, height, width, delta);
                            dst_dib = FreeImage_Allocate (width, height, 24);
                            ImthresholdSetData (dst_dib, p_im);
                        }
                        else if (strcmp(namefilter, "multiply") == 0)
                        {
                            printf("Filter= %s\n", namefilter);
                            printf("Delta= %d\n", delta);

                            IMTFilterMathMultiply (p_im, m_im, height, width, delta);
                            dst_dib = FreeImage_Allocate (width, height, 24);
                            ImthresholdSetData (dst_dib, p_im);
                        }
                        else if (strcmp(namefilter, "norm") == 0)
                        {
                            printf("Filter= %s\n", namefilter);
                            printf("Delta= %d\n", delta);

                            IMTFilterMathNorm (p_im, m_im, height, width, delta);
                            dst_dib = FreeImage_Allocate (width, height, 24);
                            ImthresholdSetData (dst_dib, p_im);
                        }
                        else if (strcmp(namefilter, "plus") == 0)
                        {
                            printf("Filter= %s\n", namefilter);
                            printf("Delta= %d\n", delta);

                            IMTFilterMathPlus (p_im, m_im, height, width, delta);
                            dst_dib = FreeImage_Allocate (width, height, 24);
                            ImthresholdSetData (dst_dib, p_im);
                        }
                        else if (strcmp(namefilter, "sbmetric") == 0)
                        {
                            float emetric;
                            printf("Filter= %s\n", namefilter);

                            emetric = IMTFilterMathSharpenBadMetric (p_im, m_im, height, width);
                            printf("SharpenBad= %f\n", emetric);
                            dst_dib = FreeImage_Allocate (width, height, 24);
                            ImthresholdSetData (dst_dib, p_im);
                        }
                        else
                        {
                            printf("Filter= threshold\n");
                            printf("Delta= %d\n", delta);

                            BYTE** d_im = BWalloc(height, width);

                            IMTFilterMathThreshold (p_im, m_im, d_im, height, width, delta);
                            dst_dib = FreeImage_Allocate (width, height, 1);
                            ImthresholdSetDataBW (dst_dib, d_im);
                            
                            BWfree(d_im, height);
                        }
                        FreeImage_Unload(dib);
                        IMTfree(m_im, height);

                        if (dst_dib)
                        {
                            FREE_IMAGE_FORMAT out_fif = FreeImage_GetFIFFromFilename (dst_filename);
                            if(out_fif != FIF_UNKNOWN)
                            {
                                FreeImage_Save (out_fif, dst_dib, dst_filename, 0);
                                printf("Output= %s\n", dst_filename);
                            }
                            FreeImage_Unload(dst_dib);
                        }
                    } else {
                        printf("%s\n", "Size math uncorect.");
                        FreeImage_Unload (dib);
                    }
                } else {
                    printf("%s\n", "Unsupported color mode.");
                    FreeImage_Unload (dib);
                }
            }
            IMTfree(p_im, height);

            printf("\n");
        } else {
            printf("%s\n", "Unsupported format type.");
            FreeImage_Unload (dib);
        }
    }

    // call this ONLY when linking with FreeImage as a static library
#ifdef FREEIMAGE_LIB
    FreeImage_DeInitialise ();
#endif // FREEIMAGE_LIB

    return 0;
}
