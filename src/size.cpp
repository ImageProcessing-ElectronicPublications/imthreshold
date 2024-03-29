//    Zlib license
//
// Resize image.
//
//    Copyright (C) 2017:
//    zvezdochiot    <zvezdochiot@user.sourceforge.net>

#include <unistd.h>
#include <FreeImage.h>
#include "imthresholdfreeimage.h"

////////////////////////////////////////////////////////////////////////////////

void ImthresholdFilterSizeTitle()
{
    printf("ImThreshold.\n");
    printf("BookScanLib Project: http://djvu-soft.narod.ru/\n\n");
    printf("Resize image.\n");
    printf("Homepage: https://sourceforge.net/projects/imthreshold/.\n\n");
}

void ImthresholdFilterSizeUsage()
{
    printf("Usage : imthreshold-size [options] <input_image> <output_image>\n\n");
    printf("options:\n");
    printf("          -f str  name filter:\n");
    printf("                    'bicubic'\n");
    printf("                    'biline'\n");
    printf("                    'bicont'\n");
    printf("                    'biakima'\n");
    printf("                    'crop'\n");
    printf("                    'gsample' (default)\n");
    printf("                    'nearest'\n");
    printf("          -p N.N  part prefilter RIS (float, optional, default = 1.0)\n");
    printf("          -q str  colorspace:\n");
    printf("                    'rgb' (default)\n");
    printf("                    'ryb4'\n");
    printf("                    'ycbcr'\n");
    printf("                    'hsv'\n");
    printf("          -r N.N  ratio (float, optional, default = 1.0)\n");
    printf("          -w N    new width (int, optional, default = [auto])\n");
    printf("          -z N    new height (int, optional, default = [auto])\n");
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
    float ratio = 1.0f, ppart = 1.0f;
    int newh = 0, neww = 0;
    bool fhelp = false;
    char *namefilter, *csp, *cspn;
    float ims = 0.0f, imratio;
    unsigned width, height, new_width, new_height, hr, wr;
    int imscaler;
    IMTpixel **p_im, **d_im, **c_im;
    IMTpixel **s_im, **m_im, **x_im;
    FIBITMAP *dib, *dst_dib;
    FREE_IMAGE_FORMAT out_fif;
    namefilter = (char*)"gsample";
    csp = (char*)"rgb";
    while ((opt = getopt(argc, argv, ":f:p:q:r:w:z:h")) != -1)
    {
        switch(opt)
        {
        case 'f':
            namefilter = optarg;
            break;
        case 'p':
            ppart = atof(optarg);
            break;
        case 'q':
            csp = optarg;
            break;
        case 'r':
            ratio = atof(optarg);
            break;
        case 'w':
            neww = atof(optarg);
            break;
        case 'z':
            newh = atof(optarg);
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

    ImthresholdFilterSizeTitle();

    if(optind + 2 > argc || fhelp || (ratio <= 0 && (newh <= 0 &&  neww <= 0)))
    {
        ImthresholdFilterSizeUsage();
        return 0;
    }
    const char *src_filename = argv[optind];
    const char *output_filename = argv[optind + 1];

    FreeImage_SetOutputMessage(FreeImageErrorHandler);

    printf("Input= %s\n", src_filename);
    dib = ImthresholdGenericLoader(src_filename, 0);
    if (dib)
    {
        if (FreeImage_GetImageType(dib) == FIT_BITMAP)
        {
            width = FreeImage_GetWidth(dib);
            height = FreeImage_GetHeight(dib);
            if (newh <= 0 && neww <= 0)
            {
                new_width = width * ratio;
                new_height = height * ratio;
            }
            else if (newh > 0 && neww > 0)
            {
                new_width = neww;
                new_height = newh;
            }
            else if (newh > 0)
            {
                new_width = (width * newh) / height;
                new_height = newh;
            }
            else
            {
                new_width = neww;
                new_height = (height * neww) / width;
            }
            p_im = IMTalloc(height, width);
            d_im = IMTalloc(new_height, new_width);

            printf("Width= %d\n", new_width);
            printf("Height= %d\n", new_height);

            ImthresholdGetData(dib, p_im);
            FreeImage_Unload(dib);

            cspn = IMTFilterRGBtoCSP(p_im, height, width, csp, 1);
            printf("ColorSpace= %s\n", cspn);

            if (strcmp(namefilter, "bicubic") == 0)
            {
                printf("Filter= %s\n", namefilter);
                imscaler = SCALER_BICUBIC;
            }
            else if (strcmp(namefilter, "biline") == 0)
            {
                printf("Filter= %s\n", namefilter);
                imscaler = SCALER_BILINE;
            }
            else if (strcmp(namefilter, "bicont") == 0)
            {
                printf("Filter= %s\n", namefilter);
                imscaler = SCALER_BICONT;
            }
            else if (strcmp(namefilter, "biakima") == 0)
            {
                printf("Filter= %s\n", namefilter);
                imscaler = SCALER_BIAKIMA;
            }
            else if (strcmp(namefilter, "crop") == 0)
            {
                printf("Filter= %s\n", namefilter);
                imscaler = SCALER_CROP;
            }
            else if (strcmp(namefilter, "nearest") == 0)
            {
                printf("Filter= %s\n", namefilter);
                imscaler = SCALER_NEAREST;
            }
            else
            {
                printf("Filter= gsample\n");
                imscaler = SCALER_GSAMPLE;
            }
            printf("Prefilter= %f\n", ppart);
            IMTFilterSize(p_im, d_im, imscaler, height, width, new_height, new_width);
            if (ppart > 0.0f)
            {
                c_im = IMTalloc(height, width);
                IMTFilterSize(d_im, c_im, imscaler, new_height, new_width, height, width);
                ims = IMTFilterMirrorPart(p_im, c_im, height, width, ppart);
                IMTfree(c_im, height);
                printf("RIS= %f\n", ims);
                IMTFilterSize(p_im, d_im, imscaler, height, width, new_height, new_width);
            }
            else if (ppart < 0.0f)
            {
                imratio = new_height;
                imratio /= height;
                imratio *= new_width;
                imratio /= width;
                if (imratio > 0.0f)
                {
                    imratio = sqrt(imratio);
                    wr = (width + imratio - 1) / imratio;
                    hr = (height + imratio - 1) / imratio;
                    s_im = IMTalloc(hr, wr);
                    m_im = IMTalloc(height, width);
                    x_im = IMTalloc(new_height, new_width);
                    IMTFilterSize(p_im, s_im, imscaler, height, width, hr, wr);
                    IMTFilterSize(s_im, m_im, imscaler, hr, wr, height, width);
                    IMTfree(s_im, hr);
                    IMTFilterMathMinus (m_im, p_im, height, width, 127);
                    ims = IMTFilterIRange (m_im, height, width, 127, 0.25f * ppart);
                    printf("IRIS= %f\n", ims);
                    IMTFilterSize(m_im, x_im, imscaler, height, width, new_height, new_width);
                    IMTfree(m_im, height);
                    IMTFilterSize(m_im, x_im, imscaler, height, width, new_height, new_width);
                    IMTFilterMathPlus (d_im, x_im, new_height, new_width, -127);
                }
            }
            IMTfree(p_im, height);

            cspn = IMTFilterRGBtoCSP(d_im, new_height, new_width, csp, -1);
            printf("ColorSpace= %s\n", cspn);

            dst_dib = FreeImage_Allocate(new_width, new_height, 24);
            ImthresholdSetData(dst_dib, d_im);
            IMTfree(d_im, new_height);

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
