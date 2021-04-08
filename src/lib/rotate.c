//  Zlib license
//
// ImThreshold library.
//
//  Copyright (C) 2017:
//  zvezdochiot <zvezdochiot@user.sourceforge.net>
//  Homepage: https://sourceforge.net/projects/imthreshold/

#include "../imthreshold.h"

/////////////////////////////////////////////////////////////////////////////////////////////

unsigned int PageTools_Next_Pow2(unsigned int n)
{
    unsigned int retval = 1;
    while(retval < n)
        retval <<= 1;
    return retval;
}

// ----------------------------------------------------------

void PageTools_Radon(BYTE** p_im, unsigned h, unsigned w, int sign, unsigned int sharpness[])
{
    unsigned char *bitcount = (unsigned char*)malloc(sizeof(unsigned char) * 256);
    unsigned short int *p1, *p2; // Stored columnwise
    unsigned int w2 = PageTools_Next_Pow2(w);
    unsigned int s = h * w2;
    p1 = (unsigned short int*)malloc(sizeof(unsigned short int) * s);
    p2 = (unsigned short int*)malloc(sizeof(unsigned short int) * s);
    // Fill in the first table
    memset(p1, 0, sizeof(unsigned short int) * s);
    unsigned int ir, ic, i, j, cnt;

    for(i = 0; i < 256; i++)
    {
        j = i, cnt = 0;
        do cnt += j & 1; while(j >>= 1);
        bitcount[i] = cnt;
    }

    for(ir = 0; ir < h; ir++)
    {
        unsigned char val;

        for(ic = 0; ic < w; ic++)
        {
            if(sign > 0)
                val = p_im[ir][w - 1 - ic];
            else
                val = p_im[ir][ic];
            p1[h * ic + ir] = bitcount[val];
        }
    }

    // Iterate
    unsigned short int *x1 = p1, *x2 = p2;
    unsigned int step = 1;
    while(step < w2)
    {
        for(i = 0; i < w2; i += 2 * step)
        {
            for(j = 0; j < step; j++)
            {
                // Columns-sources:
                unsigned short int *s1 = x1 + h * (i + j), *s2 = x1 + h * (i + j + step);

                // Columns-targets:
                unsigned short int *t1 = x2 + h * (i + 2 * j), *t2 = x2 + h * (i + 2 * j + 1);
                unsigned int m;
                for(m = 0; m < h; m++)
                {
                    t1[m] = s1[m];
                    t2[m] = s1[m];
                    if(m + j < h)
                        t1[m] += s2[m+j];
                    if(m + j + 1 < h)
                        t2[m] += s2[m + j + 1];
                }
            }
        }
        // Swap the tables:
        unsigned short int *aux = x1;
        x1 = x2;
        x2 = aux;
        // Increase the step:
        step += step;
    }
    // Now, compute the sum of squared finite differences:
    for(ic = 0; ic < w2; ic++)
    {
        unsigned int acc = 0;
        unsigned short int *col = x1 + h * ic;
        for(ir = 0; ir + 1 < h; ir++)
        {
            int diff = (int)(col[ir])-(int)(col[ir + 1]);
            acc += diff * diff;
        }
        sharpness[w2 - 1 + sign * ic] = acc;

    }
    free(p1);
    free(p2);
    free(bitcount);
}

// ----------------------------------------------------------

float IMTFilterFindSkew(BYTE** p_im, unsigned h, unsigned w)
{
    unsigned int w2 = PageTools_Next_Pow2(w);

    unsigned int ssize = 2 * w2 - 1; // Size of sharpness table
    unsigned int *sharpness = malloc(sizeof(unsigned int) * ssize);
    PageTools_Radon(p_im, h, w, 1, sharpness);
    PageTools_Radon(p_im, h, w, -1, sharpness);
    unsigned int i, imax = 0, vmax=0;
    double sum = 0.;
    for(i = 0; i < ssize; i++)
    {
        unsigned int s = sharpness[i];

        if(s > vmax)
        {
            imax = i;
            vmax = s;
        }
        sum += s;
    }

    if (vmax <= 3 * sum / h)
        return 0; // Heuristics !!!
    int iskew = (int)imax - (int)w2 + 1;
    free(sharpness);
    return atan((float)iskew / w2) * RADGRD;
}

////////////////////////////////////////////////////////////////////////////////

void IMTFilterRotate (IMTpixel** p_im, IMTpixel** d_im, unsigned height, unsigned width, float angle)
{
    unsigned y, x;
    float yt, xt, yr, xr, ktc, kts;

    angle /= RADGRD;
    kts = sin(angle);
    ktc = cos(angle);
    for (y = 0; y < height; y++)
    {
        yt = (float)y;
        yt -= height / 2;
        for (x = 0; x < width; x++)
        {
            xt = (float)x;
            xt -= (float)width / 2;
            yr = ktc * yt - kts * xt;
            yr += (float)height / 2;
            xr = ktc * xt + kts * yt;
            xr += (float)width / 2;
            if (yr >= 0.0 && yr < height && xr >= 0.0 && xr < width)
            {
                d_im[y][x] = IMTinterpolation(p_im, height, width, yr, xr);
            } else {
                d_im[y][x] = IMTset(255,255,255);
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
