//    Zlib license
//
// Thresholds an image according to Gatos et al.'s method.
//
//    Copyright (C) 2017:
//    zvezdochiot    <zvezdochiot@user.sourceforge.net>

#include <unistd.h>
#include <FreeImage.h>
#include "imthresholdfreeimage.h"

////////////////////////////////////////////////////////////////////////////////

void ImthresholdFilterTGatosTitle()
{
    printf("ImThreshold.\n");
    printf("BookScanLib Project: http://djvu-soft.narod.ru/\n\n");
    printf("Thresholds an image according to Gatos et al.'s method.\n");
    printf("Homepage: https://sourceforge.net/projects/imthreshold/.\n\n");
}

void ImthresholdFilterTGatosUsage()
{
    printf("Usage : imthreshold-tgatos [options] <input_file> <output_file>(BW) [bg_file] [niblack_file](BW)\n\n");
    printf("options:\n");
    printf("          -m str  mode {niblack, sauvola, chistian, bimod, dalg, default = niblack)\n");
    printf("          -r N    radius (int, optional, default = 7)\n");
    printf("          -s N.N  sensitivity (float, optional, default = -0.2)\n");
    printf("          -f N    dynamic range (int, optional, default = 128)\n");
    printf("          -l N    lower bound (int, optional, default = 20)\n");
    printf("          -u N    upper bound (int, optional, default = 150)\n");
    printf("          -d N.N  delta (float, optional, default = 0.0)\n");
    printf("          -q N.N  q (float, optional, default = 0.6)\n");
    printf("          -1 N.N  p1 (float, optional, default = 0.5)\n");
    printf("          -2 N.N  p2 (float, optional, default = 0.8)\n");
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
    int fmode = 0;
    int radius = 7;
    float sensitivity = -0.2;
    int dynamic_range = 128;
    int lower_bound = 20;
    int upper_bound = 150;
    float delta = 0.0;
    float q = 0.6;
    float p1 = 0.5;
    float p2 = 0.8;
    bool fhelp = false;
    int threshold = 0;
    while ((opt = getopt(argc, argv, ":m:r:s:f:l:u:d:q:1:2:h")) != -1)
    {
        switch(opt)
        {
        case 'm':
            if (strcmp(optarg, "niblack") == 0)
            {
                fmode = 0;
            }
            if (strcmp(optarg, "sauvola") == 0)
            {
                fmode = 1;
            }
            if (strcmp(optarg, "chistian") == 0)
            {
                fmode = 2;
            }
            if (strcmp(optarg, "bimod") == 0)
            {
                fmode = 3;
            }
            if (strcmp(optarg, "dalg") == 0)
            {
                fmode = 4;
            }
            break;
        case 'r':
            radius = atof(optarg);
            break;
        case 's':
            sensitivity = atof(optarg);
            break;
        case 'f':
            dynamic_range = atof(optarg);
            break;
        case 'l':
            lower_bound = atof(optarg);
            break;
        case 'u':
            upper_bound = atof(optarg);
            break;
        case 'd':
            delta = atof(optarg);
            break;
        case 'q':
            q = atof(optarg);
            break;
        case '1':
            p1 = atof(optarg);
            break;
        case '2':
            p2 = atof(optarg);
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

    ImthresholdFilterTGatosTitle();

    if(optind + 2 > argc || fhelp || radius <= 0)
    {
        ImthresholdFilterTGatosUsage();
        return 0;
    }
    const char *src_filename = argv[optind];
    const char *output_filename = argv[optind + 1];
    const char *bg_filename;
    if(optind + 2 < argc)
    {
        bg_filename = argv[optind + 2];
    }
    const char *bin_filename;
    if(optind + 3 < argc)
    {
        bin_filename = argv[optind + 3];
    }

    FreeImage_SetOutputMessage(FreeImageErrorHandler);

    printf("Input= %s\n", src_filename);
    if (fmode == 0)
    {
        printf("Mode= Niblack\n");
    }
    if (fmode == 1)
    {
        printf("Mode= Sauvola\n");
    }
    if (fmode == 2)
    {
        printf("Mode= Chistian\n");
    }
    if (fmode == 3)
    {
        printf("Mode= BiMod\n");
    }
    if (fmode == 4)
    {
        printf("Mode= D-alg\n");
    }
    FIBITMAP *dib = ImthresholdGenericLoader(src_filename, 0);
    if (dib)
    {
        if (FreeImage_GetImageType(dib) == FIT_BITMAP)
        {
            FIBITMAP* bin_dib;
            FIBITMAP* bg_dib;
            FIBITMAP* dst_dib;
            unsigned width = FreeImage_GetWidth(dib);
            unsigned height = FreeImage_GetHeight(dib);

            IMTpixel** p_im = IMTalloc(height, width);
            BYTE** d_im = BWalloc(height, width);
            IMTpixel** bg_im = IMTalloc(height, width);
            BYTE** g_im = BWalloc(height, width);

            switch(fmode)
            {
            case 1:
                printf("Dynamic= %d\n", dynamic_range);
            case 0:
            case 2:
                printf("Sensitivity= %f\n", sensitivity);
                printf("Lower= %d\n", lower_bound);
                printf("Upper= %d\n", upper_bound);
            case 4:
                printf("Radius= %d\n", radius);
            case 3:
                printf("Delta= %f\n", delta);
            }

            ImthresholdGetData(dib, p_im);
            FreeImage_Unload(dib);
            switch(fmode)
            {
            case 0:
                threshold = IMTFilterTNiblack(p_im, d_im, height, width, radius, sensitivity, lower_bound, upper_bound, delta);
                break;
            case 1:
                threshold = IMTFilterTSauvola(p_im, d_im, height, width, radius, sensitivity, dynamic_range, lower_bound, upper_bound, delta);
                break;
            case 2:
                threshold = IMTFilterTChistian(p_im, d_im, height, width, radius, sensitivity, lower_bound, upper_bound, delta);
                break;
            case 3:
                threshold = IMTFilterTBiMod(p_im, d_im, height, width, delta);
                break;
            case 4:
                threshold = IMTFilterTDalg(p_im, d_im, height, width, radius, lower_bound, upper_bound, delta);
                break;
            }
            threshold = IMTFilterGatosBG(p_im, d_im, bg_im, height, width, radius);
            printf("q= %f\n", q);
            printf("p1= %f\n", p1);
            printf("p2= %f\n", p2);
            threshold = IMTFilterTGatos(p_im, d_im, bg_im, g_im, height, width, q, p1, p2);
            printf("Threshold= %d\n", threshold / 3);
            IMTfree(p_im, height);
            bin_dib = FreeImage_Allocate(width, height, 1);
            ImthresholdSetDataBW(bin_dib, d_im);
            BWfree(d_im, height);
            bg_dib = FreeImage_Allocate(width, height, 24);
            ImthresholdSetData(bg_dib, bg_im);
            IMTfree(bg_im, height);
            dst_dib = FreeImage_Allocate(width, height, 1);
            ImthresholdSetDataBW(dst_dib, g_im);
            BWfree(g_im, height);

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
            if (bg_dib)
            {
                if(optind + 2 < argc)
                {
                    FREE_IMAGE_FORMAT out_fif = FreeImage_GetFIFFromFilename(bg_filename);
                    if(out_fif != FIF_UNKNOWN)
                    {
                        FreeImage_Save(out_fif, bg_dib, bg_filename, 0);
                        printf("Output= %s\n", bg_filename);
                    }
                }
                FreeImage_Unload(bg_dib);
                printf("\n");
            }
            if (bin_dib)
            {
                if(optind + 3 < argc)
                {
                    FREE_IMAGE_FORMAT out_fif = FreeImage_GetFIFFromFilename(bin_filename);
                    if(out_fif != FIF_UNKNOWN)
                    {
                        FreeImage_Save(out_fif, bin_dib, bin_filename, 0);
                        printf("Output= %s\n", bin_filename);
                    }
                }
                FreeImage_Unload(bin_dib);
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
