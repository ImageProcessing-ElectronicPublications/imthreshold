//  Zlib license
//
// ImThreshold library.
//
//  Copyright (C) 2017:
//  zvezdochiot <zvezdochiot@user.sourceforge.net>
//  Homepage: https://sourceforge.net/projects/imthreshold/

#include "../imthreshold.h"

////////////////////////////////////////////////////////////////////////////////

void IMTFilterRGBtoRYB4 (IMTpixel** p_im, unsigned height, unsigned width, int direct)
{
    unsigned y, x;
    IMTpixel pim;

    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            pim = p_im[y][x];
            pim = IMTRGBtoRYB4 (pim, direct);
            p_im[y][x] = pim;
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void IMTFilterRGBtoYCbCr (IMTpixel** p_im, unsigned height, unsigned width, int direct)
{
    unsigned y, x;
    IMTpixel pim;

    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            pim = p_im[y][x];
            pim = IMTRGBtoYCbCr (pim, direct);
            p_im[y][x] = pim;
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void IMTFilterRGBtoHSV (IMTpixel** p_im, unsigned height, unsigned width, int direct)
{
    unsigned y, x;
    IMTpixel pim;

    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            pim = p_im[y][x];
            pim = IMTRGBtoHSV (pim, direct);
            p_im[y][x] = pim;
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

char* IMTFilterRGBtoCSP (IMTpixel** p_im, unsigned height, unsigned width, char* csp, int direct)
{
    char *colorspace = "RGB";
    if (strcmp(csp, "ryb4") == 0)
    {
        if (direct >= 0) colorspace = "RYB4";
        IMTFilterRGBtoRYB4(p_im, height, width, direct);
    }
    else if (strcmp(csp, "ycbcr") == 0)
    {
        if (direct >= 0) colorspace = "YCbCr";
        IMTFilterRGBtoYCbCr(p_im, height, width, direct);
    }
    else if (strcmp(csp, "hsv") == 0)
    {
        if (direct >= 0) colorspace = "HSV";
        IMTFilterRGBtoHSV(p_im, height, width, direct);
    }
    return colorspace;
}

////////////////////////////////////////////////////////////////////////////////

void IMTFilterSCCor (IMTpixel** p_im, unsigned height, unsigned width)
{
    unsigned y, x;
    IMTpixel pim;

    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            pim = p_im[y][x];
            pim = IMTccorS(pim);
            p_im[y][x] = pim;
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void IMTFilterSMirror (IMTpixel** p_im, unsigned height, unsigned width)
{
    unsigned y, x, tim, pim;

    tim = (unsigned)IMTFilterTBiModValueIc (p_im, 0, 0, height, width);
    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            pim = (unsigned)p_im[y][x].s;
            if (pim > tim)
            {
                pim = 765 - pim;
            }
            pim *= 2;
            p_im[y][x].s = Byte3Clamp(pim);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterSNorm (IMTpixel** p_im, unsigned height, unsigned width)
{
    unsigned y, x, im, immin, immax, imd, threshold = 0;

    immin = 765;
    immax = 0;
    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            im = (unsigned)p_im[y][x].s;
            if (im < immin)
            {
                immin = im;
            }
            if (im > immax)
            {
                immax = im;
            }
        }
    }
    imd = immax - immin;
    if (imd > 0)
    {
        for (y = 0; y < height; y++)
        {
            for (x = 0; x < width; x++)
            {
                im = (unsigned)p_im[y][x].s;
                im -= immin;
                im *= 765;
                im /= imd;
                p_im[y][x].s = Byte3Clamp(im);
            }
        }
    }
    threshold = immax + immin;
    threshold /= 6;

    return (int)threshold;
}

////////////////////////////////////////////////////////////////////////////////

void IMTFilterSCScale (IMTpixel** p_im, unsigned height, unsigned width, float sensitivity, IMTpixel center)
{
    unsigned y, x, d;
    int im;
    float imx;

    sensitivity += 1.0f;
    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            for (d = 0; d < 3; d++)
            {
                im = p_im[y][x].c[d];
                imx = im - center.c[d];
                imx *= sensitivity;
                im = (int)(imx + 0.5f);
                im += center.c[d];
                p_im[y][x].c[d] = ByteClamp(im);
            }
            p_im[y][x] = IMTcalcS(p_im[y][x]);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void IMTFilterSEdge (IMTpixel** p_im, IMTpixel** b_im, unsigned height, unsigned width)
{
    unsigned y, x;
    int ims, imsd;
    IMTpixel pim, bim;

    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            pim = p_im[y][x];
            bim = b_im[y][x];
            ims = (int)pim.s;
            ims -= (int)bim.s;
            imsd = (int)IMTdist(pim, bim);
            if (ims < 0)
            {
                ims = 384 - imsd;
            }
            else
            {
                ims = 384 + imsd;
            }
            p_im[y][x].s = Byte3Clamp(ims);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void IMTFilterSNegate (IMTpixel** p_im, unsigned height, unsigned width)
{
    unsigned y, x;
    int ims;
    IMTpixel pim;

    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            pim = p_im[y][x];
            ims = 765 - (int)pim.s;
            p_im[y][x].s = Byte3Clamp(ims);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void IMTFilterSMean (IMTpixel** p_im, IMTpixel** b_im, unsigned height, unsigned width, float sensitivity)
{
    unsigned y, x;
    float ims;
    IMTpixel pim, bim;

    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            pim = p_im[y][x];
            bim = b_im[y][x];
            ims = (float)pim.s * sensitivity;
            ims += (float)bim.s * (1.0 - sensitivity);
            ims *= 0.5f;
            p_im[y][x].s = Byte3Clamp((int)(ims + 0.5f));
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void IMTFilterSSelect (IMTpixel** p_im, unsigned height, unsigned width, unsigned comp)
{
    unsigned y, x;
    comp = (comp > 2) ? 2 : comp;

    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            p_im[y][x].s = p_im[y][x].c[comp];
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

IMTpixel IMTFilterIllumCorr (IMTpixel** p_im, IMTpixel** b_im, IMTpixel** d_im, unsigned height, unsigned width)
{
    unsigned y, x, d, im, imb;
    float res;
    int i, k;
    IMTpixel mean; // dominant color
    BYTE val;
    long max, palette[768] = {0};
    int max_pos[3];

    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            for (d = 0; d < 3; d++)
            {
                im = (unsigned)b_im[y][x].c[d];
                im += (d * 256);
                palette[im]++;
            }
        }
    }

    // finding the dominant color
    k = 0;
    for (d = 0; d < 3; d++)
    {
        max = 0;
        for (i = 0; i < 256; i++)
        {
            if (palette[k] > max)
            {
                max = palette[k];
                max_pos[d] = i;
            }
            k++;
        }
    }
    mean.s = 0;
    for (d = 0; d < 3; d++)
    {
        val = ByteClamp(max_pos[d]); // dominant color
        mean.c[d] = val;
        mean.s += val;
    }

    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            for (d = 0; d < 3; d++)
            {
                im = (unsigned)p_im[y][x].c[d];
                im++;
                imb = (unsigned)b_im[y][x].c[d];
                imb++;
                res = (float)im;
                res /= imb;
                res *= mean.c[d];
                res += im;
                res /= 2;
                val = ByteClamp((int)(res + 0.5));
                d_im[y][x].c[d] = val;
            }
            d_im[y][x] = IMTcalcS(d_im[y][x]);
        }
    }
    return mean;
}

////////////////////////////////////////////////////////////////////////////////

void IMTFilterInvert (IMTpixel** p_im, unsigned height, unsigned width)
{
    unsigned x, y, d;
    BYTE val;
    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            for (d = 0; d < 3; d++)
            {
                val = ByteClamp(255 - p_im[y][x].c[d]);
                p_im[y][x].c[d] = val;
            }
            p_im[y][x] = IMTcalcS(p_im[y][x]);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

float IMTFilterIRange (IMTpixel** p_im, unsigned height, unsigned width, int delta, float mult)
{
    unsigned x, y, d;
    float imx, ims = 0.0f;
    BYTE val;
    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            for (d = 0; d < 3; d++)
            {
                imx = p_im[y][x].c[d];
                imx -= delta;
                imx *= mult;
                ims += ((imx < 0.0f) ? -imx : imx);
                imx += delta;
                val = ByteClamp((int)(imx + 0.5f));
                p_im[y][x].c[d] = val;
            }
            p_im[y][x] = IMTcalcS(p_im[y][x]);
        }
    }
    ims /= 3;
    ims /= width;
    ims /= height;

    return ims;
}

////////////////////////////////////////////////////////////////////////////////

void IMTFilterCopy (IMTpixel** p_im, IMTpixel** d_im, unsigned height, unsigned width)
{
    unsigned y, x;

    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            d_im[y][x] = p_im[y][x];
        }
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////

void IMTFilterEqualize (IMTpixel** p_im, unsigned height, unsigned width)
{
    unsigned int y, x, d, i, val;
    unsigned long int histogram[256] = {0}, szi;

    szi = (height * width) >> 8;
    for (d = 0; d < 3; d++)
    {
        for (i = 0; i < 256; i++)
        {
            histogram[i] = 0;
        }
        for (y = 0; y < height; y++)
        {
            for (x = 0; x < width; x++)
            {
                val = p_im[y][x].c[d];
                histogram[val]++;
            }
        }
        for (i = 1; i < 256; i++)
        {
            histogram[i] += histogram[i - 1];
        }
        for (i = 0; i < 256; i++)
        {
            histogram[i] += (szi >> 1);
            histogram[i] /= szi;
            histogram[i] = (histogram[i] < 255) ? histogram[i] : 255;
        }
        for (y = 0; y < height; y++)
        {
            for (x = 0; x < width; x++)
            {
                val = p_im[y][x].c[d];
                val = histogram[val];
                p_im[y][x].c[d] = val;
            }
        }
    }
    for ( y = 0; y < height; y++ )
    {
        for ( x = 0; x < width; x++ )
        {
            p_im[y][x] = IMTcalcS (p_im[y][x]);
        }
    }

}

/////////////////////////////////////////////////////////////////////////////////////////////

IMTpixel IMTFilterGreyNorm (IMTpixel** p_im, unsigned height, unsigned width)
{
    unsigned y, x, d, n;
    int im;
    float di, dist, sd, mins = 768.0, maxs = 0.0, means, ds, ming, maxg, meang, dg;
    IMTpixel imd;
    float pc[3] = {0.227, 0.453, 0.320}, mc[3], sc[3], ac[3];
    BYTE grey;

    n = height * width;
    ming = mins;
    maxg = maxs;
    for (d = 0; d < 3; d++)
    {
        mc[d] = 0;
        sc[d] = 0;
    }
    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            for (d = 0; d < 3; d++)
            {
                im = (int)p_im[y][x].c[d];
                mc[d] += (float)im;
            }
            im = (int)p_im[y][x].s;
            di = (float)im;
            if (di < mins)
            {
                mins = di;
            }
            if (di > maxs)
            {
                maxs = di;
            }
        }
    }
    mins /= 3;
    maxs /= 3;
    means = (maxs + mins) / 2;
    ds = (maxs > mins) ? (maxs - mins) : ((means > 127.5) ? (255.0 - means) : means);
    for (d = 0; d < 3; d++)
    {
        mc[d] /= n;
    }
    sd = 0;
    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            dist = 0;
            for (d = 0; d < 3; d++)
            {
                im = (int)p_im[y][x].c[d];
                di = (float)im;
                di -= mc[d];
                di *= pc[d];
                if (di < 0)
                {
                    di = -di;
                }
                sc[d] += di;
                di *= di;
                dist += di;
            }
            dist = sqrt(dist);
            sd += dist;
        }
    }
    sd /= n;
    for (d = 0; d < 3; d++)
    {
        sc[d] /= n;
    }
    if (sd < 0.0 || sd > 0.0)
    {
        for (d = 0; d < 3; d++)
        {
            ac[d] = sc[d] / sd;
        }
    }
    else
    {
        for (d = 0; d < 3; d++)
        {
            ac[d] = 1.0;
        }
    }
    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            dist = 0;
            for (d = 0; d < 3; d++)
            {
                im = (int)p_im[y][x].c[d];
                di = (float)im;
                di -= mc[d];
                di *= pc[d];
                di *= ac[d];
                dist += di;
            }
            if (dist < ming)
            {
                ming = dist;
            }
            if (dist > maxg)
            {
                maxg = dist;
            }
        }
    }
    meang = (maxg + ming) / 2;
    dg = (maxg > ming) ? (maxg - ming) : 1.0;
    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            dist = 0;
            for (d = 0; d < 3; d++)
            {
                im = (int)p_im[y][x].c[d];
                di = (float)im;
                di -= mc[d];
                di *= pc[d];
                di *= ac[d];
                dist += di;
            }
            dist -= meang;
            dist /= dg;
            dist *= ds;
            dist += means;
            grey = ByteClamp((int)(dist + 0.5));
            for (d = 0; d < 3; d++)
            {
                p_im[y][x].c[d] = grey;
            }
        }
    }
    for ( y = 0; y < height; y++ )
    {
        for ( x = 0; x < width; x++ )
        {
            p_im[y][x] = IMTcalcS (p_im[y][x]);
        }
    }

    for (d = 0; d < 3; d++)
    {
        imd.c[d] = ByteClamp((int)(ac[d] * 255.0 + 0.5));
    }
    imd.s = Byte3Clamp((int)means);

    return imd;
}

////////////////////////////////////////////////////////////////////////////////

IMTpixel IMTFilterGreyWorld (IMTpixel** p_im, unsigned height, unsigned width)
{
    unsigned y, x, d;
    int im;
    float si, st, kt, vt;
    IMTpixel imd;

    si = 0;
    for ( y = 0; y < height; y++ )
    {
        for ( x = 0; x < width; x++ )
        {
            im = (int)p_im[y][x].s;
            si += (float)im;
        }
    }
    si /= 3;
    for (d = 0; d < 3; d++)
    {
        st=0;
        for ( y = 0; y < height; y++ )
        {
            for ( x = 0; x < width; x++ )
            {
                im = (int)p_im[y][x].c[d];
                st += (float)im;
            }
        }
        kt = (st < 0.0 || st > 0.0) ? (si / st) : 1.0;
        imd.c[d] = ByteClamp((int)(255 * kt / 2));
        for ( y = 0; y < height; y++ )
        {
            for ( x = 0; x < width; x++ )
            {
                im = (int)p_im[y][x].c[d];
                vt = (float)(im);
                vt *= kt;
                p_im[y][x].c[d] = ByteClamp((int)(vt + 0.5));
            }
        }
    }
    for ( y = 0; y < height; y++ )
    {
        for ( x = 0; x < width; x++ )
        {
            p_im[y][x] = IMTcalcS (p_im[y][x]);
        }
    }

    return imd;
}

////////////////////////////////////////////////////////////////////////////////

float IMTFilterLevelMean (IMTpixel** p_im, IMTpixel** d_im, unsigned height, unsigned width, int radius, float contour, float thres, int lower_bound, int upper_bound)
{
    unsigned r, x, y, x0, y0, x2, y2, xf, yf, d, n;
    int im, mean, t;
    float imx, val, km;
    int d_bound;

    if (upper_bound < lower_bound)
    {
        upper_bound += lower_bound;
        lower_bound = upper_bound - lower_bound;
        upper_bound -= lower_bound;
    }
    if (contour < 0)
    {
        if ((upper_bound - lower_bound + 1) > 0)
        {
            contour = 256.0 / (upper_bound - lower_bound + 1);
        }
        else
        {
            contour = 1.0;
        }
    }
    if (thres < 0)
    {
        thres = -thres;
    }

    d_bound = upper_bound - lower_bound;
    r = (unsigned)((radius < 0) ? -radius : radius);

    for (y = 0; y < height; y++)
    {
        y0 = ((y > r) ? (y - r) : 0);
        y2 = (((y + r + 1) < height) ? (y + r + 1) : height);
        for (x = 0; x < width; x++)
        {
            x0 = ((x > r) ? (x - r) : 0);
            x2 = (((x + r + 1) < width) ? (x + r + 1) : width);
            mean = 0;
            n = 0;
            for (yf = y0; yf< y2; yf++)
            {
                for (xf = x0; xf < x2; xf++)
                {
                    mean += (int)p_im[yf][xf].s;
                    n++;
                }
            }
            mean /= n;
            im = (int)p_im[y][x].s;
            km = 1.0;
            if (im > 0)
            {
                t = im - mean;
                t = (t > -thres && t < thres) ? 0 : t;
                t *= contour;
                imx = (float)mean / 3.0;
                imx -= lower_bound;
                imx *= 255.0;
                imx /= d_bound;
                if (imx < 0.0)
                {
                    imx = 0.0;
                }
                if (imx > 255.0)
                {
                    imx = 255.0;
                }
                imx *= 3.0;
                imx += t;
                km = (float)imx / im;
            }
            for (d = 0; d < 3; d++)
            {
                val = (float)p_im[y][x].c[d];
                val *= km;
                d_im[y][x].c[d] = ByteClamp((int)(val + 0.5));
            }
            d_im[y][x] = IMTcalcS (d_im[y][x]);
        }
    }

    return contour;
}

////////////////////////////////////////////////////////////////////////////////

float IMTFilterLevelSigma (IMTpixel** p_im, IMTpixel** d_im, unsigned height, unsigned width, float thres, float fpart)
{
    unsigned y, x, d;
    int im;
    float imx, tx, apart, k, ks;

    if (thres < 0.0)
    {
        thres = -thres;
    }
    if (thres > 1.0)
    {
        thres = 1.0;
    }
    if (fpart < -1.0)
    {
        fpart = -1.0;
    }
    if (fpart > 1.0)
    {
        fpart = 1.0;
    }
    apart = fpart;
    if (apart < 0.0)
    {
        apart = -apart;
    }

    ks = 0.0;
    for ( y = 0; y < height; y++ )
    {
        for ( x = 0; x < width; x++ )
        {
            im = (int)p_im[y][x].s;
            imx = (float)im;
            imx /= 765.0;
            if (fpart > 0.0)
            {
                if (imx < thres)
                {
                    tx = thres;
                    imx = tx - sqrt(tx * tx - imx * imx);
                }
                else
                {
                    tx = 1.0 - thres;
                    imx = 1.0 - imx;
                    imx = 1.0 - tx + sqrt(tx * tx - imx * imx);
                }
            }
            else
            {
                if (imx < thres)
                {
                    tx = thres;
                    imx -= thres;
                    imx = sqrt(tx * tx - imx * imx);
                }
                else
                {
                    tx = 1.0 - thres;
                    imx -= thres;
                    imx = 1.0 - sqrt(tx * tx - imx * imx);
                }
            }
            imx *= 765.0;
            k = imx + 1.0;
            k /= (im + 1.0);
            k -= 1.0;
            k *= apart;
            k += 1.0;
            ks += k;
            for (d = 0; d < 3; d++)
            {
                im = (int)p_im[y][x].c[d];
                imx = (float)im;
                imx *= k;
                d_im[y][x].c[d] = ByteClamp((int)(imx + 0.5));
            }
            d_im[y][x] = IMTcalcS (d_im[y][x]);
        }
    }
    ks /= (float)height;
    ks /= (float)width;

    return ks;
}

///////////////////////////////////////////////////////////////////////////////

float IMTFilterLevelSize (IMTpixel** p_im, IMTpixel** d_im, unsigned height, unsigned width, int radius, int delta)
{
    unsigned x, y, d, sn = 0;
    IMTpixel mim;
    float st = 0;
    int c0, c;
    BYTE val;

    if (radius < 0)
    {
        radius = -radius;
    }

    IMTpixel** s_im = IMTalloc(radius, radius);
    IMTFilterSGsample(p_im, s_im, height, width, radius, radius);
    IMTFilterSBicub(s_im, d_im, radius, radius, height, width);
    IMTfree(s_im, radius);
    mim = IMTmeanIc (p_im, 0, 0, height, width);

    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            for (d = 0; d < 3; d++)
            {
                c0 = (int)p_im[y][x].c[d];
                c = c0;
                c += (int)mim.c[d];
                c -= (int)d_im[y][x].c[d];
                c += delta;
                val = ByteClamp(c);
                c = (int)val;
                c = (c > c0) ? (c - c0) : (c0 -c);
                st += c;
                sn++;
                d_im[y][x].c[d] = val;
            }
        }
    }
    sn = (sn < 0.0 || sn > 0.0) ? sn : 1.0;
    st /= sn;

    return st;
}

////////////////////////////////////////////////////////////////////////////////

float IMTFilterMirror (IMTpixel** p_im, IMTpixel** m_im, unsigned height, unsigned width)
{
    unsigned y, x, d;
    int im, imm, imd;
    float ims = 0.0;

    for ( y = 0; y < height; y++ )
    {
        for ( x = 0; x < width; x++ )
        {
            for (d = 0; d < 3; d++)
            {
                im = (int)p_im[y][x].c[d];
                imm = (int)m_im[y][x].c[d];
                imd = im - imm;
                im += imd;
                ims += (float)imd;
                p_im[y][x].c[d] = ByteClamp(im);
            }
            p_im[y][x] = IMTcalcS (p_im[y][x]);
        }
    }
    ims /= 3.0;
    ims /= (float)width;
    ims /= (float)height;

    return ims;
}

///////////////////////////////////////////////////////////////////////////////

float IMTFilterMirrorPart (IMTpixel** p_im, IMTpixel** m_im, unsigned height, unsigned width, float part)
{
    unsigned y, x, d;
    int im, imm, imd;
    float imx, ims = 0.0;

    for ( y = 0; y < height; y++ )
    {
        for ( x = 0; x < width; x++ )
        {
            for (d = 0; d < 3; d++)
            {
                im = (int)p_im[y][x].c[d];
                imm = (int)m_im[y][x].c[d];
                imd = im - imm;
                imx = (float)imd;
                imx *= part;
                imd = (int)(imx + 0.5);
                im += imd;
                ims += (float)((imd < 0) ? -imd : imd);
                p_im[y][x].c[d] = ByteClamp(im);
            }
            p_im[y][x] = IMTcalcS (p_im[y][x]);
        }
    }
    ims /= 3;
    ims /= width;
    ims /= height;

    return ims;
}

///////////////////////////////////////////////////////////////////////////////

void IMTFilterMirrorMean (IMTpixel** p_im, unsigned height, unsigned width)
{
    unsigned y, x, d;
    int im, imm;

    for ( y = 0; y < height; y++ )
    {
        for ( x = 0; x < width; x++ )
        {
            for (d = 0; d < 3; d++)
            {
                im = (int)p_im[y][x].c[d];
                imm = (im < 128) ? (im + 127) : (im - 128);
                p_im[y][x].c[d] = ByteClamp(imm);
            }
            p_im[y][x] = IMTcalcS (p_im[y][x]);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

unsigned IMTFilterKMeans (IMTpixel** IMTim, unsigned height, unsigned width, unsigned ncluster, unsigned iters)
{
    int y, x, ym, xm, d, itr, i, h, l, dist, distmin, cluster, changes;
    float fcluster;
    IMTpixel p, pm;
    IMTcluster *means = (IMTcluster*)malloc(ncluster * sizeof(IMTcluster));

    for (i = 0; i < ncluster; i++)
    {
        h = 255 * i / ncluster;
        p.c[0] = h;
        p.c[1] = 255;
        p.c[2] = 255;
        p = IMTcalcS (p);
        p = IMTRGBtoHSV (p, -1);
        distmin = 256 * 3;
        ym = 0;
        xm = 0;
        for (y = 0; y < height; y++)
        {
            for (x = 0; x < width; x++)
            {
                dist = IMTdist(p, IMTim[y][x]);
                if (dist < distmin)
                {
                    distmin = dist;
                    ym = y;
                    xm = x;
                }

            }
        }
        for (d = 0; d < 3; d++)
        {
            means[i].c[d] = (float)IMTim[ym][xm].c[d];
        }
        means[i].n = 1;
    }

    fcluster = 1.0f / sqrt((float)ncluster);
    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            IMTim[y][x].s = 0;
        }
    }
    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            distmin = 256 * 3;
            l = 0;
            for (i = 0; i < ncluster; i++)
            {
                for (d = 0; d < 3; d++)
                {
                    pm.c[d] = (unsigned)(means[i].c[d] + 0.5f);
                }
                dist = IMTdist(pm, IMTim[y][x]);
                if (dist < distmin)
                {
                    distmin = dist;
                    l = i;
                }
            }
            IMTim[y][x].s = l;
        }
    }

    for (itr = 0; itr < iters; itr++)
    {
        for (i = 0; i < ncluster; i++)
        {
            for (d = 0; d < 3; d++)
            {
                means[i].c[d] = 0.0f;
            }
            means[i].n = 0;
        }
        for (y = 0; y < height; y++)
        {
            for (x = 0; x < width; x++)
            {
                cluster = IMTim[y][x].s;
                for (d = 0; d < 3; d++)
                {
                    means[cluster].c[d] += (float)IMTim[y][x].c[d];
                }
                means[cluster].n++;
            }
        }
        changes = 0;
        for (i = 0; i < ncluster; i++)
        {
            if (means[i].n > 0)
            {
                for (d = 0; d < 3; d++)
                {
                    means[i].c[d] /= (float)means[i].n;
                }
            }
            else
            {
                h = 255 * i / ncluster;
                p.c[0] = h;
                p.c[1] = 255;
                p.c[2] = 255;
                p = IMTcalcS (p);
                p = IMTRGBtoHSV (p, -1);
                distmin = 256 * 3;
                l = 0;
                for (i = 0; i < ncluster; i++)
                {
                    for (d = 0; d < 3; d++)
                    {
                        pm.c[d] = (unsigned)(means[i].c[d] + 0.5f);
                    }
                    dist = IMTdist(p, pm);
                    if (dist < distmin)
                    {
                        distmin = dist;
                        l = i;
                    }
                }
                for (d = 0; d < 3; d++)
                {
                    p.c[d] *= fcluster;
                    p.c[d] += (unsigned)(means[l].c[d] * (1.0f - fcluster) + 0.5f);
                }
                distmin = 256 * 3;
                ym = 0;
                xm = 0;
                for (y = 0; y < height; y++)
                {
                    for (x = 0; x < width; x++)
                    {
                        dist = IMTdist(p, IMTim[y][x]);
                        if (dist < distmin)
                        {
                            distmin = dist;
                            ym = y;
                            xm = x;
                        }

                    }
                }
                for (d = 0; d < 3; d++)
                {
                    means[i].c[d] = (float)IMTim[ym][xm].c[d];
                }
                means[i].n = 1;
                changes++;
            }
        }
        for (y = 0; y < height; y++)
        {
            for (x = 0; x < width; x++)
            {
                distmin = 256 * 3;
                cluster = 0;
                for (i = 0; i < ncluster; i++)
                {
                    for (d = 0; d < 3; d++)
                    {
                        pm.c[d] = (unsigned)(means[i].c[d] + 0.5f);
                    }
                    dist = IMTdist(pm, IMTim[y][x]);
                    if (dist < distmin)
                    {
                        distmin = dist;
                        cluster = i;
                    }
                }
                if (cluster != IMTim[y][x].s)
                {
                    changes++;
                    IMTim[y][x].s = cluster;
                }
            }
        }
        if (changes == 0)
        {
            break;
        }
    }
    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            cluster = IMTim[y][x].s;
            for (d = 0; d < 3; d++)
            {
                IMTim[y][x].c[d] = ByteClamp((int)(means[cluster].c[d] + 0.5f));
            }
            IMTim[y][x] = IMTcalcS (IMTim[y][x]);
        }
    }
    free(means);

    return iters;
}

////////////////////////////////////////////////////////////////////////////////

float IMTFilterPosterize (IMTpixel** p_im, unsigned height, unsigned width, unsigned thres)
{
    unsigned y, x, d, n;
    unsigned imin, imax, imst, imc;
    float xmin, xmax, xd, ims = 0.0, ime, sumc = 0.0;
    IMTpixel im, immin, immax;

    if (thres < 1)
    {
        thres = 1;
    }
    if (thres > 256)
    {
        thres = 256;
    }

    if (thres > 1)
    {
        n=0;
        for (d = 0; d < 3; d++)
        {
            immin.c[d] = 255;
            immax.c[d] = 0;
        }
        for (y = 0; y < height; y++)
        {
            for (x = 0; x < width; x++)
            {
                im = p_im[y][x];
                for (d = 0; d < 3; d++)
                {
                    if (im.c[d] < immin.c[d])
                    {
                        immin.c[d] = im.c[d];
                    }
                    if (im.c[d] > immax.c[d])
                    {
                        immax.c[d] = im.c[d];
                    }
                }
            }
        }
        for (d = 0; d < 3; d++)
        {
            imin = (unsigned)immin.c[d];
            xmin = (float)imin;
            xd = (float)(immax.c[d] - immin.c[d] + 1);
            xd /= 256.0;
            xd *= thres;
            imst = (unsigned)immax.c[d];
            imst++;
            while (imin < imst)
            {
                xmax = xmin + xd;
                imax = (unsigned)(xmax + 0.5);
                n = 0;
                sumc = 0.0;
                for (y = 0; y < height; y++)
                {
                    for (x = 0; x < width; x++)
                    {
                        imc = (unsigned)p_im[y][x].c[d];
                        if ((imc >= imin) && (imc < imax))
                        {
                            sumc += (float)imc;
                            n++;
                        }
                    }
                }
                if (n > 0)
                {
                    sumc /= n;
                    sumc += 0.5;
                    for (y = 0; y < height; y++)
                    {
                        for (x = 0; x < width; x++)
                        {
                            imc = (unsigned)p_im[y][x].c[d];
                            if ((imc >= imin) && (imc < imax))
                            {
                                ime = (float)imc;
                                ime -= sumc;
                                ime = (ime < 0.0) ? -ime : ime;
                                ims += ime;
                                p_im[y][x].c[d] = ByteClamp((int)(sumc + 0.5));
                            }
                        }
                    }
                }
                xmin += xd;
                imin = (unsigned)(xmin + 0.5);
            }
        }
        ims /= 3;
        ims /= width;
        ims /= height;
    }
    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            p_im[y][x] = IMTcalcS (p_im[y][x]);
        }
    }

    return ims;
}

////////////////////////////////////////////////////////////////////////////////

float IMTFilterQuant (IMTpixel** p_im, unsigned height, unsigned width, unsigned quant)
{
    unsigned y, x, d, k, l, n;
    unsigned im, imt, imd;
    unsigned histogram[256];
    unsigned histthres[256];
    unsigned histquant[256];
    float histstep, histsum, histsumk, imds = 0.0;

    quant = (quant > 0) ? quant : 1;
    n = height * width;
    n = (n > 0) ? n : 1;
    histstep = (float)n;
    histstep /= quant;
    for (d = 0; d < 3; d++)
    {
        for (k = 0; k < 256; k++)
        {
            histogram[k] = 0;
            histthres[k] = 0;
            histquant[k] = 0;
        }
        for (y = 0; y < height; y++)
        {
            for (x = 0; x < width; x++)
            {
                im = (unsigned)p_im[y][x].c[d];
                histogram[im]++;
            }
        }
        histsum = 0.0;
        for (k = 0; k < 256; k++)
        {
            histsum += histogram[k];
            l = (histsum > 0.0) ? (unsigned)((histsum - 1.0)/ histstep) : 0;
            histthres[k] = l;
        }
        for (l = 0; l < quant; l++)
        {
            histsum = 0;
            histsumk = 0;
            for (k = 0; k < 256; k++)
            {
                if (histthres[k] == l)
                {
                    histsum += histogram[k];
                    histsumk += (histogram[k] * k);
                }
            }
            histsumk = (histsum > 0) ? (histsumk / histsum) : histsumk;
            for (k = 0; k < 256; k++)
            {
                if (histthres[k] == l)
                {
                    histquant[k] = (unsigned)histsumk;
                }
            }
        }
        for (y = 0; y < height; y++)
        {
            for (x = 0; x < width; x++)
            {
                im = (unsigned)p_im[y][x].c[d];
                imt = histquant[im];
                imd = (im < imt) ? (imt - im) : (im - imt);
                p_im[y][x].c[d] = ByteClamp((int)imt);
                imds += imd;
            }
        }
    }
    imds /= n;
    imds /= 3;
    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            p_im[y][x] = IMTcalcS (p_im[y][x]);
        }
    }

    return imds;
}

////////////////////////////////////////////////////////////////////////////////

void IMTFilterMonoColor (IMTpixel** p_im, unsigned height, unsigned width)
{
    unsigned x, y, d, maxc;
    BYTE maxcv;
    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            maxc = 0;
            maxcv = p_im[y][x].c[0];
            for (d = 1; d < 3; d++)
            {
                if (p_im[y][x].c[d] > maxcv)
                {
                    maxc = d;
                    maxcv = p_im[y][x].c[d];
                }
            }
            for (d = 0; d < 3; d++)
            {
                p_im[y][x].c[d] = 0;
            }
            p_im[y][x].c[maxc] = maxcv;
            p_im[y][x].s = (WORD)maxcv;
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterWhiteFill (IMTpixel** p_im, unsigned height, unsigned width)
{
    unsigned y, x, yy, xx;
    int sw, niter = 0;
    WORD im, imf;

    sw = 1;
    while (sw > 0)
    {
        niter++;
        sw = 0;
        for ( x = 0; x < width; x++ )
        {
            for ( y = 0; y < height; y++ )
            {
                im = p_im[y][x].s;
                if (im != 765)
                {
                    if (x > 0)
                    {
                        imf = p_im[y][x - 1].s;
                        if (imf == 765)
                        {
                            p_im[y][x - 1] = p_im[y][x];
                        }
                    }
                    if (x < width - 1)
                    {
                        imf = p_im[y][x + 1].s;
                        if (imf == 765)
                        {
                            p_im[y][x + 1] = p_im[y][x];
                        }
                    }
                    if (y > 0)
                    {
                        imf = p_im[y - 1][x].s;
                        if (imf == 765)
                        {
                            p_im[y - 1][x] = p_im[y][x];
                        }
                    }
                    if (y < height - 1)
                    {
                        imf = p_im[y + 1][x].s;
                        if (imf == 765)
                        {
                            p_im[y + 1][x] = p_im[y][x];
                        }
                    }
                }
                else
                {
                    sw++;
                }
            }
        }
        sw = 0;
        for ( xx = width; xx > 0; xx-- )
        {
            x = xx - 1;
            for ( yy = height; yy > 0; yy-- )
            {
                y = yy - 1;
                im = p_im[y][x].s;
                if (im != 765)
                {
                    if (x > 0)
                    {
                        imf = p_im[y][x - 1].s;
                        if (imf == 765)
                        {
                            p_im[y][x - 1] = p_im[y][x];
                        }
                    }
                    if (x < width - 1)
                    {
                        imf = p_im[y][x + 1].s;
                        if (imf == 765)
                        {
                            p_im[y][x + 1] = p_im[y][x];
                        }
                    }
                    if (y > 0)
                    {
                        imf = p_im[y - 1][x].s;
                        if (imf == 765)
                        {
                            p_im[y - 1][x] = p_im[y][x];
                        }
                    }
                    if (y < height - 1)
                    {
                        imf = p_im[y + 1][x].s;
                        if (imf == 765)
                        {
                            p_im[y + 1][x] = p_im[y][x];
                        }
                    }
                }
                else
                {
                    sw++;
                }
            }
        }
    }
    return niter;
}

////////////////////////////////////////////////////////////////////////////////
