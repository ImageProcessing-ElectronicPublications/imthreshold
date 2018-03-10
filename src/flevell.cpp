//    Zlib license
//
// LevelL filter image.
//
//    Copyright (C) 2017:
//    zvezdochiot    <zvezdochiot@user.sourceforge.net>

#include <unistd.h>
#include <FreeImage.h>
#include "imthresholdfreeimage.h"

////////////////////////////////////////////////////////////////////////////////

void ImthresholdFilterLevelLTitle()
{
    printf("ImThreshold.\n");
    printf("BookScanLib Project: http://djvu-soft.narod.ru/\n\n");
    printf("LevelL filter image.\n");
    printf("Homepage: https://sourceforge.net/projects/imthreshold/.\n\n");
}

void ImthresholdFilterLevelLUsage()
{
    printf("Usage : imthreshold-flevell [options] <input_file> <output_file>\n\n");
    printf("options:\n");
    printf("          -l N    level (int, optional, default = 10)\n");
    printf("          -n N    break level (int, optional, default = 0)\n");
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
    int level = 10, num = 0;
    bool fhelp = false;
    double iml = 0;
    while ((opt = getopt(argc, argv, ":l:n:h")) != -1)
    {
        switch(opt)
        {
            case 'l':
                level = atof(optarg);
                break;
            case 'n':
                num = atof(optarg);
                break;
            case 'h':
                fhelp = 1;
                break;
            case ':':
                printf("option needs a value\n");
                break;
            case '?':
                printf("unknown option: %c\n", optopt);
                break;
        }
    }

    ImthresholdFilterLevelLTitle();

    if(optind + 2 > argc || fhelp)
    {
        ImthresholdFilterLevelLUsage();
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
            unsigned y;

            IMTpixel** p_im;
            p_im = (IMTpixel**)malloc(height * sizeof(IMTpixel*));
            for (y = 0; y < height; y++) {p_im[y] = (IMTpixel*)malloc(width * sizeof(IMTpixel));}
            IMTpixel** d_im;
            d_im = (IMTpixel**)malloc(height * sizeof(IMTpixel*));
            for (y = 0; y < height; y++) {d_im[y] = (IMTpixel*)malloc(width * sizeof(IMTpixel));}

            printf("Level= %d\n", level);
            printf("Num= %d\n", num);

            ImthresholdGetData(dib, p_im);
            FreeImage_Unload(dib);
            iml = IMTFilterLevelL(p_im, d_im, height, width, level, num);
            printf("Mean= %f\n", iml);
            for (y = 0; y < height; y++){free(p_im[y]);}
            free(p_im);
            dst_dib = FreeImage_Allocate(width, height, 24);
            ImthresholdSetData(dst_dib, d_im);
            for (y = 0; y < height; y++){free(d_im[y]);}
            free(d_im);

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
