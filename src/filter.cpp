//  Zlib license
//
// Filter image.
//
//  Copyright (C) 2017:
//  zvezdochiot <zvezdochiot@user.sourceforge.net>

#include <unistd.h>
#include <FreeImage.h>
#include "imthresholdfreeimage.h"

////////////////////////////////////////////////////////////////////////////////

void ImthresholdFilterTitle()
{
    printf("ImThreshold.\n");
    printf("BookScanLib Project: http://djvu-soft.narod.ru/\n\n");
    printf("Filter image.\n");
    printf("Homepage: https://sourceforge.net/projects/imthreshold/.\n\n");
}

////////////////////////////////////////////////////////////////////////////////

void ImthresholdFilterUsage()
{
    printf("Usage : imthreshold-filter [options] <input_image> <output_image>\n\n");
    printf("options:\n");
    printf("          -a N.N  amount (float, optional, default = 0.5)\n");
    printf("          -c N.N  contour factor (float, optional, default = -1[auto])\n");
    printf("          -d N    max delta (int, optional, default = 50)\n");
    printf("          -f str  name filter:\n");
    printf("                    'adsmooth'\n");
    printf("                    'autolevel'\n");
    printf("                    'autowhite'\n");
    printf("                    'bimod'\n");
    printf("                    'blur'\n");
    printf("                    'bwc'\n");
    printf("                    'cluster'\n");
    printf("                    'denoise'\n");
    printf("                    'greynorm'\n");
    printf("                    'greyworld'\n");
    printf("                    'illumc'\n");
    printf("                    'invert'\n");
    printf("                    'kmeans'\n");
    printf("                    'levelmean'\n");
    printf("                    'levelsigma'\n");
    printf("                    'levelsize'\n");
    printf("                    'mirror'\n");
    printf("                    'monocolor'\n");
    printf("                    'mdilate'\n");
    printf("                    'merose'\n");
    printf("                    'none' (default)\n");
    printf("                    'peron'\n");
    printf("                    'posterize'\n");
    printf("                    'quant'\n");
    printf("                    'retinex'\n");
    printf("                    'rs'\n");
    printf("                    'selgauss'\n");
    printf("                    'shrink'\n");
    printf("                    'unripple'\n");
    printf("                    'unsharp'\n");
    printf("                    'whitefill'\n");
    printf("                    'wiener'\n");
    printf("          -i N    num divide (int, optional, default = 1)\n");
    printf("          -k N    K num (int, optional, default = 2)\n");
    printf("          -l N    lower bound (int, optional, default = 32)\n");
    printf("          -n N.N  noise size (float, optional, default = -1.0[auto])\n");
    printf("          -p N    posterize divide factor (int, optional, default = 16)\n");
    printf("          -q str  colorspace:\n");
    printf("                    'rgb' (default)\n");
    printf("                    'ryb4'\n");
    printf("                    'ycbcr'\n");
    printf("                    'hsv'\n");
    printf("          -r N.N  radius (float, optional, default = 3.0)\n");
    printf("          -s N.N  sigma (float, optional, default = 4.0)\n");
    printf("          -t N.N  threshold (float, optional, default = 0.0)\n");
    printf("          -u N    upper bound (int, optional, default = 223)\n");
    printf("          -h      this help\n");
}

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
    // call this ONLY when linking with FreeImage as a static library
#ifdef FREEIMAGE_LIB
    FreeImage_Initialise();
#endif // FREEIMAGE_LIB

    unsigned width, height;
    int opt;
    float amount = 0.5f;
    float contour = -1.0f;
    int maxdelta = 50;
    unsigned ndiv = 1;
    int knum = 2;
    int lower_bound = 32;
    float noise = -1.0f;
    int posterdiv = 16;
    float radius = 5.0f;
    float sigma = 4.0f;
    float threshold = 0.0f;
    int upper_bound = 223;
    bool frf = false;
    bool fhelp = false;
    char *namefilter, *csp, *cspn;
    namefilter="none";
    csp = "rgb";
    while ((opt = getopt(argc, argv, ":a:c:d:f:i:k:l:n:p:q:r:s:t:u:h")) != -1)
    {
        switch(opt)
        {
        case 'a':
            amount = atof(optarg);
            break;
        case 'c':
            contour = atof(optarg);
            break;
        case 'd':
            maxdelta = atof(optarg);
            break;
        case 'f':
            namefilter = optarg;
            break;
        case 'i':
            ndiv = atof(optarg);
            break;
        case 'k':
            knum = atof(optarg);
            break;
        case 'l':
            lower_bound = atof(optarg);
            break;
        case 'n':
            noise = atof(optarg);
            break;
        case 'p':
            posterdiv = atof(optarg);
            break;
        case 'q':
            csp = optarg;
            break;
        case 'r':
            radius = atof(optarg);
            break;
        case 's':
            sigma = atof(optarg);
            break;
        case 't':
            threshold = atof(optarg);
            break;
        case 'u':
            upper_bound = atof(optarg);
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

    ImthresholdFilterTitle();

    if (optind + 2 > argc || fhelp)
    {
        ImthresholdFilterUsage();
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
            unsigned radiusint;
            float imsh = 0.0f, ks = 0.0f, kwb = 0.5f, kdenoises;
            int thres = 0, imdm, kimi = 0, niter;
            IMTpixel dim;
            FIBITMAP* dst_dib;
            width = FreeImage_GetWidth(dib);
            height = FreeImage_GetHeight(dib);
            IMTpixel** p_im = IMTalloc(height, width);
            IMTpixel** d_im = IMTalloc(height, width);
            ImthresholdGetData(dib, p_im);
            FreeImage_Unload(dib);

            if (radius < 0.0f)
            {
                frf = true;
                radius = -radius;
            }
            radiusint = (unsigned)(radius + 0.5f);
            radiusint = (radiusint > 0) ? radiusint : 1;
            posterdiv = (posterdiv < 2) ? 2 : posterdiv;
            knum = (knum < 2) ? 2 : knum;

            cspn = IMTFilterRGBtoCSP(p_im, height, width, csp, 1);
            printf("ColorSpace= %s\n", cspn);

            if (strcmp(namefilter, "adsmooth") == 0)
            {
                printf("Filter= %s\n", namefilter);
                printf("Radius= %f\n", radius);
                IMTFilterAdSmooth(p_im, d_im, height, width, radius);
            }
            else if (strcmp(namefilter, "autolevel") == 0)
            {
                printf("Filter= %s\n", namefilter);
                printf("Radius= %d\n", radiusint);
                IMTFilterAutoLevel (p_im, d_im, height, width, radiusint);
            }
            else if (strcmp(namefilter, "autowhite") == 0)
            {
                printf("Filter= %s\n", namefilter);
                printf("Radius= %f\n", radius);
                IMTFilterAutoWhite(p_im, d_im, height, width, radius);
            }
            else if (strcmp(namefilter, "bimod") == 0)
            {
                printf("Filter= %s\n", namefilter);
                printf("Count= %d\n", posterdiv);
                thres = IMTFilterTBiModP(p_im, d_im, height, width, posterdiv);
                thres /= 3;
                printf("Threshold= %d\n", thres);
            }
            else if (strcmp(namefilter, "bwc") == 0)
            {
                printf("Filter= %s\n", namefilter);
                printf("Radius= %d\n", radiusint);
                kwb = IMTFilterClusterBWC(p_im, d_im, height, width, radiusint);
                printf("W/b= %f\n", kwb);
            }
            else if (strcmp(namefilter, "blur") == 0)
            {
                printf("Filter= %s\n", namefilter);
                radius = (radius > 0.0f) ? radius : 1.0f;
                printf("Radius= %f\n", radius);
                IMTFilterGaussBlur(p_im, d_im, height, width, radius);
            }
            else if (strcmp(namefilter, "cluster") == 0)
            {
                printf("Filter= %s\n", namefilter);
                printf("K= %d\n", knum);
                imdm = IMTFilterClusterBiModC (p_im, d_im, height, width, knum);
                printf("Dmax= %d\n", imdm);
            }
            else if (strcmp(namefilter, "denoise") == 0)
            {
                printf("Filter= %s\n", namefilter);
                printf("Radius= %d\n", radiusint);
                printf("KUse= %f\n", threshold);
                IMTFilterCopy (p_im, d_im, height, width);
                kdenoises = IMTFilterDeNoiseDiff (d_im, height, width, radiusint, threshold);
                printf("KDeNoise= %f\n", kdenoises);
            }
            else if (strcmp(namefilter, "greyworld") == 0)
            {
                printf("Filter= %s\n", namefilter);
                IMTFilterCopy (p_im, d_im, height, width);
                dim = IMTFilterGreyWorld(d_im, height, width);
                printf("Normalize= %d,%d,%d\n", dim.c[0], dim.c[1], dim.c[2]);
            }
            else if (strcmp(namefilter, "greynorm") == 0)
            {
                printf("Filter= %s\n", namefilter);
                IMTFilterCopy (p_im, d_im, height, width);
                dim = IMTFilterGreyNorm(d_im, height, width);
                printf("Normalize= %d,%d,%d: %d\n", dim.c[0], dim.c[1], dim.c[2], dim.s);
            }
            else if (strcmp(namefilter, "illumc") == 0)
            {
                printf("Filter= %s\n", namefilter);
                IMTpixel** b_im = IMTalloc(height, width);

                radius = (radius > 0.0f) ? radius : 1.0f;
                printf("Radius= %f\n", radius);
                printf("Threshold= %f\n", threshold);
                IMTFilterGaussBlur(p_im, b_im, height, width, radius);
                dim = IMTFilterIllumCorr(p_im, b_im, d_im, height, width);
                printf("Dominante= %d,%d,%d\n", dim.c[0], dim.c[1], dim.c[2]);
                IMTfree(b_im, height);
            }
            else if (strcmp(namefilter, "invert") == 0)
            {
                printf("Filter= %s\n", namefilter);
                IMTFilterCopy (p_im, d_im, height, width);
                IMTFilterInvert(d_im, height, width);
            }
            else if (strcmp(namefilter, "kmeans") == 0)
            {
                printf("Filter= %s\n", namefilter);
                printf("K= %d\n", knum);
                IMTFilterCopy (p_im, d_im, height, width);
                threshold = IMTFilterKMeans(d_im, height, width, knum, threshold);
                printf("Iters= %d\n", (int)threshold);
            }
            else if (strcmp(namefilter, "levelmean") == 0)
            {
                printf("Filter= %s\n", namefilter);
                printf("Radius= %d\n", radiusint);
                printf("Threshold= %f\n", threshold);
                printf("Lower= %d\n", lower_bound);
                printf("Upper= %d\n", upper_bound);
                contour = IMTFilterLevelMean(p_im, d_im, height, width, radiusint, contour, threshold, lower_bound, upper_bound);
                printf("Contour= %f\n", contour);
            }
            else if (strcmp(namefilter, "levelsigma") == 0)
            {
                printf("Filter= %s\n", namefilter);
                printf("Average= %f\n", amount);
                printf("Part= %f\n", threshold);
                ks = IMTFilterLevelSigma(p_im, d_im, height, width, amount, threshold);
                printf("Multy= %f\n", ks);
            }
            else if (strcmp(namefilter, "levelsize") == 0)
            {
                printf("Filter= %s\n", namefilter);
                printf("Radius= %d\n", radiusint);
                printf("Delta= %f\n", threshold);
                ks = IMTFilterLevelSize(p_im, d_im, height, width, radiusint, threshold);
                printf("Distance= %f\n", ks);
            }
            else if (strcmp(namefilter, "mirror") == 0)
            {
                printf("Filter= %s\n", namefilter);
                IMTFilterCopy (p_im, d_im, height, width);
                IMTFilterMirrorMean(d_im, height, width);
            }
            else if (strcmp(namefilter, "monocolor") == 0)
            {
                printf("Filter= %s\n", namefilter);
                IMTFilterCopy (p_im, d_im, height, width);
                IMTFilterMonoColor(d_im, height, width);
            }
            else if (strcmp(namefilter, "mdilate") == 0)
            {
                printf("Filter= %s\n", namefilter);
                printf("Radius= %f\n", radius);
                IMTFilterMorph(p_im, d_im, height, width, radius, true);
            }
            else if (strcmp(namefilter, "merose") == 0)
            {
                printf("Filter= %s\n", namefilter);
                printf("Radius= %f\n", radius);
                IMTFilterMorph(p_im, d_im, height, width, radius, false);
            }
            else if (strcmp(namefilter, "peron") == 0)
            {
                printf("Filter= %s\n", namefilter);
                printf("Radius= %f\n", radius);
                noise = (noise > 0.0f) ? noise : (IMTFilterNoiseVariance (p_im, height, width, radiusint) / 2);
                printf("Noise= %f\n", noise);
                IMTFilterPeron(p_im, d_im, height, width, radius, noise);
            }
            else if (strcmp(namefilter, "posterize") == 0)
            {
                printf("Filter= %s\n", namefilter);
                printf("Posterize= %d\n", posterdiv);
                IMTFilterCopy (p_im, d_im, height, width);
                imsh = IMTFilterPosterize(d_im, height, width, posterdiv);
                printf("Delta= %f\n", imsh);
            }
            else if (strcmp(namefilter, "quant") == 0)
            {
                printf("Filter= %s\n", namefilter);
                printf("Quant= %d\n", posterdiv);
                IMTFilterCopy (p_im, d_im, height, width);
                imsh = IMTFilterQuant(d_im, height, width, posterdiv);
                printf("Delta= %f\n", imsh);
            }
            else if (strcmp(namefilter, "retinex") == 0)
            {
                printf("Filter= %s\n", namefilter);
                radiusint = (radiusint > 0) ? radiusint : 1;
                printf("Radius= %d\n", radiusint);
                printf("Sigma= %f\n", sigma);
                kimi = IMTFilterRetinex(p_im, d_im, height, width, radiusint, sigma);
                printf("K= %d\n", kimi);
            }
            else if (strcmp(namefilter, "rs") == 0)
            {
                printf("Filter= %s\n", namefilter);
                IMTFilterCopy (p_im, d_im, height, width);
                imsh = IMTFilterRS(d_im, height, width);
                printf("Stdev= %f\n", imsh);
            }
            else if (strcmp(namefilter, "selgauss") == 0)
            {
                printf("Filter= %s\n", namefilter);
                radius = (radius > 0.0f) ? radius : 1.0f;
                printf("Radius= %f\n", radius);
                printf("MaxDelta= %d\n", maxdelta);
                IMTFilterSelGauss(p_im, d_im, height, width, radius, maxdelta);
            }
            else if (strcmp(namefilter, "shrink") == 0)
            {
                printf("Filter= %s\n", namefilter);
                threshold = (threshold <= 0) ? 4 : threshold;
                printf("Threshold= %f\n", threshold);
                IMTFilterCopy (p_im, d_im, height, width);
                imsh = IMTFilterShrink(d_im, height, width, threshold);
                printf("Shrink= %f\n", imsh);
            }
            else if (strcmp(namefilter, "unripple") == 0)
            {
                printf("Filter= %s\n", namefilter);
                radiusint = (radiusint > 0) ? radiusint : 1;
                printf("Radius= %d\n", radiusint);
                printf("MaxDelta= %d\n", maxdelta);
                IMTFilterUnRipple(p_im, d_im, height, width, radiusint, maxdelta);
            }
            else if (strcmp(namefilter, "unsharp") == 0)
            {
                printf("Filter= %s\n", namefilter);
                IMTpixel** b_im = IMTalloc(height, width);

                radius = (radius > 0.0f) ? radius : 1.0f;
                printf("Radius= %f\n", radius);
                printf("Amount= %f\n", amount);
                printf("Threshold= %f\n", threshold);
                IMTFilterCopy (p_im, d_im, height, width);
                IMTFilterGaussBlur(d_im, b_im, height, width, radius);
                IMTFilterUnsharpMask(d_im, b_im, height, width, amount, threshold);
                IMTfree(b_im, height);
            }
            else if (strcmp(namefilter, "whitefill") == 0)
            {
                printf("Filter= %s\n", namefilter);
                IMTFilterCopy (p_im, d_im, height, width);
                niter = IMTFilterWhiteFill(d_im, height, width);
                printf("Iteration= %d\n", niter);
            }
            else if (strcmp(namefilter, "wiener") == 0)
            {
                printf("Filter= %s\n", namefilter);
                printf("Radius= %d\n", radiusint);
                noise = (noise > 0.0f) ? noise : IMTFilterNoiseVariance (p_im, height, width, radiusint);
                printf("Noise= %f\n", noise);
                IMTFilterWiener(p_im, d_im, height, width, radiusint, noise);
            }
            else
            {
                printf("Filter= none\n");
                IMTFilterCopy (p_im, d_im, height, width);
            }

            if (frf)
            {
                printf("Postprocess = Reverse\n");
                IMTFilterReverse (p_im, d_im, height, width);
            }
            IMTfree(p_im, height);

            cspn = IMTFilterRGBtoCSP(d_im, height, width, csp, -1);
            printf("ColorSpace= %s\n", cspn);

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

