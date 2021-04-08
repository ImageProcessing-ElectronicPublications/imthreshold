//    Zlib license
//
// Deskew image.
//
//    Copyright (C) 2017:
//    zvezdochiot    <zvezdochiot@user.sourceforge.net>

#include <unistd.h>
#include <FreeImage.h>
#include "imthresholdfreeimage.h"

////////////////////////////////////////////////////////////////////////////////

void ImthresholdFilterDeskewTitle()
{
    printf("ImThreshold.\n");
    printf("BookScanLib Project: http://djvu-soft.narod.ru/\n\n");
    printf("The Deskew algorithm adapted for the FreeImage\n");
    printf("Homepage: https://sourceforge.net/projects/imthreshold/.\n\n");
}

////////////////////////////////////////////////////////////////////////////////

void ImthresholdFilterDeskewUsage()
{
    printf("Usage : imthreshold-deskew [options] <input_image> <output_image>\n\n");
    printf("options:\n");
    printf("          -f str  name filter:\n");
    printf("                    'bimod' (default)\n");
    printf("                    'grad'\n");
    printf("                    'niblack'\n");
    printf("                    'sauvola'\n");
    printf("          -h      this help\n");
}

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
    // call this ONLY when linking with FreeImage as a static library
#ifdef FREEIMAGE_LIB
    FreeImage_Initialise();
#endif // FREEIMAGE_LIB

    char *namefilter = "bimod";
    int opt;
    bool fhelp = false;
    while ((opt = getopt(argc, argv, ":f:h")) != -1)
    {
        switch(opt)
        {
            case 'f':
                namefilter = optarg;
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

    ImthresholdFilterDeskewTitle();

    if(optind + 2 > argc || fhelp)
    {
        ImthresholdFilterDeskewUsage();
        return 0;
    }
    const char *src_filename = argv[optind];
    const char *output_filename = argv[optind + 1];

    FreeImage_SetOutputMessage(FreeImageErrorHandler);

    printf("Input= %s\n", src_filename);
    FIBITMAP *dib = ImthresholdGenericLoader(src_filename, 0);
    if (dib)
    {
        if (FreeImage_GetImageType(dib) == FIT_BITMAP)
        {
            FIBITMAP* dst_dib;
            unsigned width = FreeImage_GetWidth(dib);
            unsigned height = FreeImage_GetHeight(dib);
            float fskew;

            IMTpixel** p_im = IMTalloc(height, width);
            BYTE** s_im = BWalloc(height, width);
            IMTpixel** d_im = IMTalloc(height, width);

            ImthresholdGetData(dib, p_im);
            FreeImage_Unload(dib);
            if (strcmp(namefilter, "grad") == 0)
            {
                printf("Filter= %s\n", namefilter);
                (void)IMTFilterTGrad(p_im, s_im, height, width);
            }
            else if (strcmp(namefilter, "niblack") == 0)
            {
                printf("Filter= %s\n", namefilter);
                (void)IMTFilterTNiblack(p_im, s_im, height, width, 7, 0.2, 0, 255, -5.0f);
            }
            else if (strcmp(namefilter, "sauvola") == 0)
            {
                printf("Filter= %s\n", namefilter);
                (void)IMTFilterTSauvola(p_im, s_im, height, width, 7, 0.2, 128, 0, 255, -5.0f);
            }
            else
            {
                printf("Filter= %s\n", "bimod");
                (void)IMTFilterTBiMod(p_im, s_im, height, width, 0);
            }
            fskew = IMTFilterFindSkew(s_im, height, width);
            BWfree(s_im, height);
            printf("Angle= %f\n", fskew);
            if (fskew > 0.0f || fskew < 0.0f)
                IMTFilterRotate(p_im, d_im, height, width, fskew);
            else
                IMTFilterCopy(p_im, d_im, height, width);
            IMTfree(p_im, height);

            dst_dib = FreeImage_Allocate(width, height, 24);
            ImthresholdSetData(dst_dib, d_im);
            IMTfree(d_im, height);

            if (dst_dib)
            {
                FREE_IMAGE_FORMAT out_fif = FreeImage_GetFIFFromFilename(output_filename);
                if(out_fif != FIF_UNKNOWN)
                {
                    FreeImage_Save(out_fif, dst_dib, output_filename, 0);
                    printf("Output= %s\n\n", output_filename);
                }
                FreeImage_Unload(dst_dib);
            }
        } else {
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

