//    Zlib license
//
// Text thresholding image.
//
//    Copyright (C) 2017:
//    zvezdochiot    <zvezdochiot@user.sourceforge.net>

#include <unistd.h>
#include <FreeImage.h>
#include "imthresholdfreeimage.h"

////////////////////////////////////////////////////////////////////////////////

void ImthresholdFilterTTextTitle()
{
    printf("ImThreshold.\n");
    printf("BookScanLib Project: http://djvu-soft.narod.ru/\n\n");
    printf("Text thresholding image.\n");
    printf("Homepage: https://sourceforge.net/projects/imthreshold/.\n\n");
}

void ImthresholdFilterTTextUsage()
{
    printf("Usage : imthreshold-ttext [options] <input_file> <output_file>(BW) [textg_file]\n\n");
    printf("options:\n");
    printf("          -c N    contour amplitude (int, optional, default = 5)\n");
    printf("          -q str  colorspace:\n");
    printf("                    'rgb' (default)\n");
    printf("                    'ryb4'\n");
    printf("                    'ycbcr'\n");
    printf("                    'hsv'\n");
    printf("          -r N    radius (int, optional, default = 5)\n");
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
    int contour = 5;
    int radius = 5;
    bool fhelp = false;
    int threshold = 0;
    char *csp, *cspn;
    csp = (char*)"rgb";
    while ((opt = getopt(argc, argv, ":c:q:r:h")) != -1)
    {
        switch(opt)
        {
        case 'c':
            contour = atof(optarg);
            break;
        case 'q':
            csp = optarg;
            break;
        case 'r':
            radius = atof(optarg);
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

    ImthresholdFilterTTextTitle();

    if(optind + 2 > argc || contour < 1 || radius < 1 || fhelp)
    {
        ImthresholdFilterTTextUsage();
        return 0;
    }
    const char *src_filename = argv[optind];
    const char *output_filename = argv[optind + 1];
    const char *txtg_filename;
    if(optind + 2 < argc)
    {
        txtg_filename = argv[optind + 2];
    }

    FreeImage_SetOutputMessage(FreeImageErrorHandler);

    printf("Input= %s\n", src_filename);
    FIBITMAP *dib = ImthresholdGenericLoader(src_filename, 0);
    if (dib)
    {
        if (FreeImage_GetImageType(dib) == FIT_BITMAP)
        {
            FIBITMAP* dst_dib;
            FIBITMAP* txt_dib;
            unsigned width = FreeImage_GetWidth(dib);
            unsigned height = FreeImage_GetHeight(dib);
            unsigned y, x;
            IMTpixel** p_im = IMTalloc(height, width);
            BYTE** d_im = BWalloc(height, width);

            printf("Contour= %d\n", contour);
            printf("Radius= %d\n", radius);

            ImthresholdGetData(dib, p_im);
            FreeImage_Unload(dib);

            cspn = IMTFilterRGBtoCSP(p_im, height, width, csp, 1);
            printf("ColorSpace= %s\n", cspn);

            threshold = IMTFilterTText(p_im, d_im, height, width, contour, radius);
            printf("Threshold= %d\n", threshold / 3);
            dst_dib = FreeImage_Allocate(width, height, 1);
            ImthresholdSetDataBW(dst_dib, d_im);
            for (y = 0; y < height; y++)
            {
                for (x = 0; x < width; x++)
                {
                    if (d_im[y][x] == 0)
                    {
                        p_im[y][x] = IMTset(255, 255, 255);
                    }
                }
            }
            BWfree(d_im, height);
            txt_dib = FreeImage_Allocate(width, height, 24);
            ImthresholdSetData(txt_dib, p_im);
            IMTfree(p_im, height);

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
            if (txt_dib)
            {
                if(optind + 2 < argc)
                {
                    FREE_IMAGE_FORMAT out_fif = FreeImage_GetFIFFromFilename(txtg_filename);
                    if(out_fif != FIF_UNKNOWN)
                    {
                        FreeImage_Save(out_fif, txt_dib, txtg_filename, 0);
                        printf("Output= %s\n", txtg_filename);
                    }
                }
                FreeImage_Unload(txt_dib);
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
