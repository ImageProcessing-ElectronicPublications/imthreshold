//    Zlib license
//
// Half Reverse Interpolate Scale image (HRIS).
//
//    Copyright (C) 2017:
//    zvezdochiot    <zvezdochiot@user.sourceforge.net>

#include <unistd.h>
#include <FreeImage.h>
#include "imthresholdfreeimage.h"

////////////////////////////////////////////////////////////////////////////////

void ImthresholdFilterSHRISTitle()
{
    printf("ImThreshold.\n");
    printf("BookScanLib Project: http://djvu-soft.narod.ru/\n\n");
    printf("Half Reverse Interpolate Scale image (HRIS).\n");
    printf("Homepage: https://sourceforge.net/projects/imthreshold/.\n\n");
}

void ImthresholdFilterSHRISUsage()
{
    printf("Usage : imthreshold-shris [options] <input_file> <output_file>\n\n");
    printf("options:\n");
    printf("          -f str  name filter:\n");
    printf("                    'hris' (default)\n");
    printf("                    'gsample'\n");
    printf("                    'scalenx'\n");
    printf("                    'iris'\n");
    printf("                    'frp' (long!)\n");
    printf("          -m      mode (int, optional, default = 2, {2,3})\n");
    printf("          -q str  colorspace:\n");
    printf("                    'rgb' (default)\n");
    printf("                    'ryb4'\n");
    printf("                    'ycbcr'\n");
    printf("                    'hsv'\n");
    printf("          -r      reduce scale (bool, optional)\n");
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
    unsigned smode = 2;
    bool reduce = false;
    bool fhelp = false;
    char *namefilter, *csp, *cspn;
    unsigned width, height, width2, height2;
    int hr, wr, hr2, wr2;
    IMTpixel **p_im, **d_im;
    IMTpixel **s_im, **m_im, **x_im;
    FIBITMAP *dib, *dst_dib;
    FREE_IMAGE_FORMAT out_fif;
    namefilter = (char*)"hris";
    csp = (char*)"rgb";
    while ((opt = getopt(argc, argv, ":f:m:q:rh")) != -1)
    {
        switch(opt)
        {
        case 'f':
            namefilter = optarg;
            break;
        case 'm':
            smode = atof(optarg);
            break;
        case 'q':
            csp = optarg;
            break;
        case 'r':
            reduce = true;
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

    ImthresholdFilterSHRISTitle();

    if(optind + 2 > argc || fhelp)
    {
        ImthresholdFilterSHRISUsage();
        return 0;
    }
    const char *src_filename = argv[optind];
    const char *output_filename = argv[optind + 1];

    if (smode < 2)
    {
        smode = 2;
    }
    if ((strcmp(namefilter, "gsample") > 0) && (strcmp(namefilter, "frp") > 0) && !(reduce))
    {
        if (smode > 3)
        {
            smode = 3;
        }
    }

    FreeImage_SetOutputMessage(FreeImageErrorHandler);

    printf("Input= %s\n", src_filename);
    dib = ImthresholdGenericLoader(src_filename, 0);
    if (dib)
    {
        if (FreeImage_GetImageType(dib) == FIT_BITMAP)
        {
            width = FreeImage_GetWidth(dib);
            height = FreeImage_GetHeight(dib);
            printf("Mode= %d\n", smode);
            if (reduce > 0)
            {
                width2 = (width + smode - 1) / smode;
                height2 = (height + smode - 1) / smode;
            }
            else
            {
                width2 = width * smode;
                height2 = height * smode;
            }
            printf("Width= %d\n", width2);
            printf("Height= %d\n", height2);
            p_im = IMTalloc(height, width);
            d_im = IMTalloc(height2, width2);

            ImthresholdGetData(dib, p_im);
            FreeImage_Unload(dib);

            cspn = IMTFilterRGBtoCSP(p_im, height, width, csp, 1);
            printf("ColorSpace= %s\n", cspn);

            if (reduce)
            {
                printf("Scale= Reduce.\n");
                IMTFilterSReduce(p_im, d_im, height, width, smode);
            }
            else
            {
                if (strcmp(namefilter, "iris") == 0)
                {
                    printf("Scale= Up %s\n", namefilter);
                    wr = (width + smode - 1) / smode;
                    hr = (height + smode - 1) / smode;
                    wr2 = wr * smode;
                    hr2 = hr * smode;
                    s_im = IMTalloc(hr, wr);
                    m_im = IMTalloc(hr2, wr2);
                    x_im = IMTalloc(height2, width2);
                    IMTFilterSReduce(p_im, s_im, height, width, smode);
                    IMTFilterSHRIS(s_im, m_im, hr, wr, smode);
                    IMTfree(s_im, hr);
                    IMTFilterMathMinus (m_im, p_im, height, width, 127);
                    IMTFilterInvert(m_im, hr2, wr2);
                    IMTFilterSHRIS(m_im, x_im, height, width, smode);
                    IMTfree(m_im, hr2);
                    IMTFilterSHRIS(p_im, d_im, height, width, smode);
                    IMTFilterMathPlus (x_im, d_im, height2, width2, -127);
                    IMTFilterMathAverage (x_im, d_im, height2, width2, 0);
                    IMTFilterMathAverage (d_im, x_im, height2, width2, 0);
                    IMTfree(x_im, height2);
                }
                else if (strcmp(namefilter, "gsample") == 0)
                {
                    printf("Scale= Up %s\n", namefilter);
                    IMTFilterSGSampleUp(p_im, d_im, height, width, smode);
                }
                else if (strcmp(namefilter, "scalenx") == 0)
                {
                    printf("Scale= Up %s\n", namefilter);
                    IMTFilterSScaleNX(p_im, d_im, height, width, smode);
                }
                else if (strcmp(namefilter, "frp") == 0)
                {
                    printf("Scale= Up %s\n", namefilter);
                    IMTFilterSFRP(p_im, d_im, height, width, smode);
                }
                else
                {
                    printf("Scale= Up hris\n");
                    IMTFilterSHRIS(p_im, d_im, height, width, smode);
                }
            }
            IMTfree(p_im, height);

            cspn = IMTFilterRGBtoCSP(d_im, width2, height2, csp, -1);
            printf("ColorSpace= %s\n", cspn);

            dst_dib = FreeImage_Allocate(width2, height2, 24);
            ImthresholdSetData(dst_dib, d_im);
            IMTfree(d_im, height2);

            if (dst_dib)
            {
                out_fif = FreeImage_GetFIFFromFilename(output_filename);
                if(out_fif != FIF_UNKNOWN)
                {
                    FreeImage_Save(out_fif, dst_dib, output_filename, 0);
                    printf("Output= %s\n\n", output_filename);
                }
                FreeImage_Unload(dst_dib);
            }
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

