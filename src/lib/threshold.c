//  Zlib license
//
// ImThreshold library.
//
//  Copyright (C) 2017:
//  zvezdochiot <zvezdochiot@user.sourceforge.net>
//  Homepage: https://sourceforge.net/projects/imthreshold/

#include "../imthreshold.h"

////////////////////////////////////////////////////////////////////////////////

void IMTHist (IMTpixel** p_im, unsigned long long* histogram, unsigned y0, unsigned x0, unsigned y1, unsigned x1, unsigned histsize)
{
    unsigned y, x, im;

    for (im = 0; im < histsize; im++)
    {
        histogram[im] = 0;
    }
    for (y = y0; y < y1; y++)
    {
        for (x = x0; x < x1; x++)
        {
            im = (unsigned)p_im[y][x].s;
            histogram[im]++;
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

int IMTHistBiMod (unsigned long long* histogram, unsigned histsize, float part)
{
    unsigned k, T, Tn;
    unsigned long long im, iw, ib, Tw, Tb;
    int threshold = 0;

    part = (part < 0.0f || part > 1.0f) ? 0.5f : part;
    T = (unsigned)(part * (float)histsize + 0.5f);
    Tn = 0;
    while ( T != Tn )
    {
        Tn = T;
        Tb = Tw = ib = iw = 0;
        for (k = 0; k < T; k++)
        {
            im = histogram[k];
            Tb += (im * k);
            ib += im;
        }
        for (k = T; k < histsize; k++)
        {
            im = histogram[k];
            Tw += (im * k);
            iw += im;
        }
        Tb /= (ib > 1) ? ib : 1;
        Tw /= (iw > 1) ? iw : 1;
        if (iw == 0 && ib == 0)
        {
            T = Tn;
        }
        else if (iw == 0)
        {
            T = (unsigned)Tb;
        }
        else if (ib == 0)
        {
            T = (unsigned)Tw;
        }
        else
        {
            T = (unsigned)(part * (float)Tw + (1.0f - part) * (float)Tb + 0.5f);
        }
    }
    threshold = (int)T;

    return threshold;
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterThreshold (IMTpixel** p_im, BYTE** d_im, unsigned height, unsigned width, int threshold)
{
    unsigned y, x;
    BYTE val;
    int im;

    for (y = 0; y < height; y++ )
    {
        for (x = 0; x < width; x++)
        {
            im = p_im[y][x].s;
            val = (BYTE) ((im >= threshold) ? 255 : 0 );
            d_im[y][x] = val;
        }
    }

    return threshold;
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterThresholdLayer (IMTpixel** p_im, WORD** t_im, BYTE** d_im, unsigned height, unsigned width)
{
    unsigned y, x;
    BYTE val;
    int im, imt, threshold = 0;
    float imts;

    imts = 0;
    for (y = 0; y < height; y++ )
    {
        for (x = 0; x < width; x++)
        {
            im = p_im[y][x].s;
            imt = t_im[y][x];
            imts += imt;
            val = (BYTE) ((im >= imt) ? 255 : 0 );
            d_im[y][x] = val;
        }
    }
    imts /= height;
    imts /= width;
    threshold = (int)(imts + 0.5);

    return threshold;
}

////////////////////////////////////////////////////////////////////////////////

void IMTFilterTLayerToImg (WORD** t_im, IMTpixel** d_im, unsigned height, unsigned width)
{
    unsigned y, x;
    BYTE val;
    int im;

    for (y = 0; y < height; y++ )
    {
        for (x = 0; x < width; x++)
        {
            im = t_im[y][x];
            im += 1;
            im /= 3;
            val = ByteClamp(im);
            d_im[y][x] = IMTset(val, val, val);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterTAbutaleb (IMTpixel** p_im, BYTE** d_im, unsigned height, unsigned width, int radius)
{
    unsigned x, y, i, n, st = 0, sn = 0;
    int y1, x1, y2, x2, yf, xf;
    float imx, imm;
    int h = height;
    int w = width;
    unsigned a, b, s, t;
    BYTE val;
    int hw = 256; // square histogram width
    int hist_size = hw*hw*sizeof(float);
    float one_over_area = 1.0 / (width * height);
    float P_sum, p, H_sum, H_end, Phi, P, H;
    float Phi_max = 1;
    float Phi_max_sub = 1;
    float tiny = 1e-6;
    unsigned threshold = 0, avg_threshold = 0;

    float* histogram = (float*)malloc(hist_size);
    float* P_histogram = (float*)malloc(hist_size);
    float* H_histogram = (float*)malloc(hist_size);

    memset(histogram,0,hist_size);
    memset(P_histogram,0,hist_size);
    memset(H_histogram,0,hist_size);

    if (radius < 0)
    {
        radius = -radius;
    }

    for (y = 0; y < height; y++)
    {
        y1 = y;
        y1 -= radius;
        if (y1 < 0)
        {
            y1 = 0;
        }
        y2 = y + 1;
        y2 += radius;
        if (y2 > h)
        {
            y2 = h;
        }
        for (x = 0; x < width; x++)
        {
            x1 = x;
            x1 -= radius;
            if (x1 < 0)
            {
                x1 = 0;
            }
            x2 = x + 1;
            x2 += radius;
            if (x2 > w)
            {
                x2 = w;
            }
            imm = 0;
            n = 0;
            for (yf = y1; yf < y2; yf++)
            {
                for (xf = x1; xf < x2; xf++)
                {
                    imx = (float)p_im[yf][xf].s;
                    imm += imx;
                    n++;
                }
            }
            if (n > 0)
            {
                imm /= n;
            }
            imm /= 3;
            imx = (float)p_im[y][x].s;
            imx /= 3;
            a = (unsigned)imx;
            b = (unsigned)imm;
            i = (a * hw + b);
            histogram[i]++;
        }
    }

    for (b = 0; b < 256; b++)
    {
        for (a = 0; a < 256; a++)
        {
            i = (a * hw + b);
            histogram[i] = histogram[i] * one_over_area;
        }
    }
    P_sum = 0.0;
    for (s = 0; s < 256; s++)
    {
        i = s * hw;
        P_sum += histogram[i];
        P_histogram[i] = P_sum;
    }
    for (t = 1; t < 256; t++)
    {
        P_sum = 0.0;
        for (s = 0; s < 256; s++)
        {
            i = s * hw + t;
            P_sum += histogram[i];
            P_histogram[i] = P_histogram[i - 1] + P_sum;
        }
    }
    H_sum = 0.0;
    for (s = 0; s < 256; s++)
    {
        i = s * hw;
        p = histogram[i];
        if (p != 0)
        {
            H_sum -= (p * log(p));
        }
        H_histogram[i] = H_sum;
    }
    for (t = 1; t < 256; t++)
    {
        H_sum = 0.0;
        for (s = 0; s < 256; ++s)
        {
            i = s * hw + t;
            p = histogram[i];
            if (p != 0)
            {
                H_sum -= (p * log(p));
            }
            H_histogram[i] = H_histogram[i - 1] + H_sum;
        }
    }

    while (1 + Phi_max_sub > 1)
    {
        Phi_max = Phi_max_sub;
        Phi_max_sub = Phi_max / 10;
    }
    H_end = H_histogram[255*hw+255];
    threshold = 0;
    avg_threshold = 0;

    for (s = 0; s < 256; s++)
    {
        for (t = 0; t < 256; t++)
        {
            i = s * hw + t;
            P = P_histogram[i];
            H = H_histogram[i];
            if ((P > tiny) && ((1.0 - P) > tiny))
            {
                Phi = log(P * (1.0 - P)) + H / P + (H_end - H) / (1.0 - P);
                if (Phi > Phi_max)
                {
                    Phi_max = Phi;
                    threshold = s;
                    avg_threshold = t;
                }
            }
        }
    }

    for (y = 0; y < height; y++)
    {
        y1 = y;
        y1 -= radius;
        if (y1 < 0)
        {
            y1 = 0;
        }
        y2 = y + 1;
        y2 += radius;
        if (y2 > h)
        {
            y2 = h;
        }
        for (x = 0; x < width; x++)
        {
            x1 = x;
            x1 -= radius;
            if (x1 < 0)
            {
                x1 = 0;
            }
            x2 = x + 1;
            x2 += radius;
            if (x2 > w)
            {
                x2 = w;
            }
            imm = 0;
            n = 0;
            for (yf = y1; yf < y2; yf++)
            {
                for (xf = x1; xf < x2; xf++)
                {
                    imx = (float)p_im[yf][xf].s;
                    imm += imx;
                    n++;
                }
            }
            if (n > 0)
            {
                imm /= n;
            }
            imm /= 3;
            imx = (float)p_im[y][x].s;
            imx /= 3;
            if (imx <= (float)threshold && imm <= (float)avg_threshold)
            {
                val = 0;
            }
            else
            {
                val = 255;
                if (imx > (float)threshold)
                {
                    st += threshold;
                    sn++;
                }
                if (imm > (float)avg_threshold)
                {
                    st += avg_threshold;
                    sn++;
                }
            }
            d_im[y][x] = val;
        }
    }
    if (sn > 0)
    {
        threshold = st * 3 / sn;
    }

    free(histogram);
    free(P_histogram);
    free(H_histogram);

    return threshold;
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterTBernsenLayer(IMTpixel** p_im, WORD** t_im, unsigned height, unsigned width, int radius, unsigned contrast_limit)
{
    int x, y, i, j, xt, yt;
    unsigned im, t;
    unsigned c, minimum, maximum;
    WORD val;
    int h = height;
    int w = width;
    int threshold = 0, st = 0, sn = 0;

    if (radius < 0)
    {
        radius = -radius;
    }

    for (y = 0; y < h; y++)
    {
        for (x = 0; x < w; x++)
        {
            minimum = 768;
            maximum = 0;
            for (i = -radius; i <= radius; i++)
            {
                yt = y + i;
                if (yt >= 0 && yt < h)
                {
                    for (j = -radius; j <= radius; j++)
                    {
                        xt = x + j;
                        if (xt >= 0 && xt < w)
                        {
                            im = p_im[yt][xt].s;
                            minimum = MIN(minimum, im);
                            maximum = MAX(maximum, im);
                        }
                    }
                }
            }
            c = maximum - minimum;
            if (c < contrast_limit * 3)
            {
                val = p_im[y][x].s;
            }
            else
            {
                t = (maximum + minimum) / 2;
                st += t;
                sn++;
                val = Byte3Clamp((int)t);
            }
            t_im[y][x] = val;
        }
    }
    if (sn > 0)
    {
        threshold =  st / sn;
    }

    return threshold;
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterTBernsen(IMTpixel** p_im, BYTE** d_im, unsigned height, unsigned width, int radius, unsigned contrast_limit)
{
    int threshold = 0;
    WORD** t_im = TLalloc(height, width);

    threshold = IMTFilterTBernsenLayer(p_im, t_im, height, width, radius, contrast_limit);
    threshold = IMTFilterThresholdLayer (p_im, t_im, d_im, height, width);

    TLfree(t_im, height);

    return threshold;
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterTBHTValue (IMTpixel** p_im, unsigned height, unsigned width)
{
    int Tmax = 768;
    float wl = 0, wr = 0, wd = 0;
    int i, il = 0, ir = Tmax  - 1, threshold = Tmax / 2;
    unsigned long long histogram[768] = {0};

    IMTHist (p_im, histogram, 0, 0, height, width, Tmax);
    wd = 0;
    for (i = 0; i < Tmax; i++)
    {
        wd += (float)histogram[i];
    }
    wd *= 2.0;
    for (i = il; i <= threshold; i++)
    {
        wl += (float)histogram[i];
    }
    for (i = threshold + 1; i <= ir; i++)
    {
        wr += (float)histogram[i];
    }

    while (il <= ir)
    {
        threshold = (il + ir) / 2;
        if (wr > wl)
        {
            wr -= (float)histogram[ir];
            ir--;
            wr += (float)histogram[threshold] * 0.5;
            wl -= (float)histogram[threshold] * 0.5;
        }
        else
        {
            wl -= (float)histogram[il];
            il++;
            wl += (float)histogram[threshold + 1] * 0.5;
            wr -= (float)histogram[threshold + 1] * 0.5;
        }
    }

    return threshold;
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterTBHT (IMTpixel** p_im, BYTE** d_im, unsigned height, unsigned width)
{
    int threshold = 384;

    threshold = IMTFilterTBHTValue (p_im, height, width);
    threshold = IMTFilterThreshold (p_im, d_im, height, width, threshold);

    return threshold;
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterTBiModValueBound (IMTpixel** p_im, unsigned height, unsigned width, int lower_bound, int upper_bound)
{
    unsigned y, x;
    int im;
    float T, Tw, Tb, Tn, iw, ib;
    int threshold = 0;

    if (upper_bound < lower_bound)
    {
        upper_bound += lower_bound;
        lower_bound = upper_bound - lower_bound;
        upper_bound -= lower_bound;
    }

    T = (lower_bound + upper_bound);
    T *= 0.5;
    Tn = 0.0;
    while ( T != Tn )
    {
        Tn = T;
        Tb = Tw = ib = iw = 0;
        for (y = 0; y < height; y++ )
        {
            for (x = 0; x < width; x++)
            {
                im = p_im[y][x].s;
                if ((im >= lower_bound) && (im <= upper_bound))
                {
                    if ( im > T)
                    {
                        Tw += im;
                        iw++;
                    }
                    else
                    {
                        Tb += im;
                        ib++;
                    }
                }
            }
        }
        if (iw == 0 && ib == 0)
        {
            T = Tn;
        }
        else if (iw == 0)
        {
            T = Tb/ib;
        }
        else if (ib == 0)
        {
            T = Tw/iw;
        }
        else
        {
            T = ((Tw/iw) + (Tb/ib)) * 0.5;
        }
    }
    threshold = (int)(T + 0.5);

    return threshold;
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterTBiModValueIcP (IMTpixel** p_im, unsigned y0, unsigned x0, unsigned y1, unsigned x1, float part)
{
    unsigned Tmax = 768;
    int threshold = 0;
    unsigned long long histogram[768] = {0};

    IMTHist (p_im, histogram, y0, x0, y1, x1, Tmax);
    threshold = IMTHistBiMod (histogram, Tmax, part);

    return threshold;
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterTBiModValueIc (IMTpixel** p_im, unsigned y0, unsigned x0, unsigned y1, unsigned x1)
{
    return IMTFilterTBiModValueIcP (p_im, y0, x0, y1, x1, 0.5);
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterTBiModValueP (IMTpixel** p_im, unsigned height, unsigned width, float part)
{
    return IMTFilterTBiModValueIcP (p_im, 0, 0, height, width, part);
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterTBiModValue (IMTpixel** p_im, unsigned height, unsigned width)
{
    return IMTFilterTBiModValueIc (p_im, 0, 0, height, width);
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterTBiMod (IMTpixel** p_im, BYTE** d_im, unsigned height, unsigned width, int delta)
{
    int threshold = 0;

    threshold = IMTFilterTBiModValue (p_im, height, width);
    threshold += delta;
    threshold = IMTFilterThreshold (p_im, d_im, height, width, threshold);

    return threshold;
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterTBiModP (IMTpixel** p_im, IMTpixel** d_im, unsigned height, unsigned width, unsigned count)
{
    int threshold = 0, t, Tmax = 256;
    unsigned i, y, x, d, im;
    float part;
    BYTE val;
    unsigned long long histogram[256];

    count = (count < 2) ? 2 : count;
    threshold = 0;
    for (d = 0; d < 3; d++)
    {
        IMTFilterSSelect (p_im, height, width, d);
        IMTHist (p_im, histogram, 0, 0, height, width, Tmax);
        for (i = 1; i < count; i++)
        {
            part = (float)i;
            part /= (float)count;
            val = ByteClamp((int)(255 * i / (count - 1)));
            t = IMTHistBiMod (histogram, Tmax, part);
            threshold += t;
            for (y = 0; y < height; y++ )
            {
                for (x = 0; x < width; x++)
                {
                    im = p_im[y][x].c[d];
                    if (im >= t)
                    {
                        d_im[y][x].c[d] = val;
                    }
                }
            }
        }
    }
    threshold /= (count - 1);
    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            d_im[y][x] = IMTcalcS(d_im[y][x]);
        }
    }

    return threshold;
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterTBiModLayer (IMTpixel** p_im, WORD** t_im, unsigned height, unsigned width, int radius, float sensitivity, int delta)
{
    unsigned int y, x, y0, x0, y1, x1, r;
    float TG, T, Ts;
    WORD val;
    int threshold = 0;

    radius = (radius < 0) ? -radius : radius;
    r = radius;
    TG = IMTFilterTBiModValueIc (p_im, 0, 0, height, width);
    if (r > 0)
    {
        Ts = 0;
        TG *= (1 - sensitivity);
        for (y = 0; y < height; y++ )
        {
            y0 = (y < r) ? 0 : y - r;
            y1 = y + r;
            y1 = (y1 < height) ? y1 : height;
            for (x = 0; x < width; x++ )
            {
                x0 = (x < r) ? 0 : (x - r);
                x1 = x + r;
                x1 = (x1 < width) ? x1 : width;
                T = IMTFilterTBiModValueIc (p_im, y0, x0, y1, x1);
                T *= sensitivity;
                T += TG;
                T += delta;
                Ts += T;
                val = Byte3Clamp((int)(T + 0.5));
                t_im[y][x] = val;
            }
        }
        Ts /= height;
        Ts /= width;
        threshold = (int)(Ts + 0.5);
    }
    else
    {
        threshold = TG;
        threshold += delta;
        val = Byte3Clamp((int)(threshold + 0.5));
        for (y = 0; y < height; y++)
        {
            for (x = 0; x < width; x++)
            {
                t_im[y][x] = val;
            }
        }
    }

    return threshold;
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterTBiModRegion (IMTpixel** p_im, BYTE** d_im, unsigned height, unsigned width, int radius, float sensitivity, int delta)
{
    int threshold = 0;
    WORD** t_im = TLalloc(height, width);

    threshold = IMTFilterTBiModLayer (p_im, t_im, height, width, radius, sensitivity, delta);
    threshold = IMTFilterThresholdLayer (p_im, t_im, d_im, height, width);

    TLfree(t_im, height);

    return threshold;
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterTBiModC (IMTpixel** p_im, BYTE** d_im, unsigned height, unsigned width, int delta)
{
    unsigned y, x, d, im, ims;
    unsigned thres[3];
    BYTE val;
    int threshold = 0;
    unsigned Tmax = 256;
    unsigned long long histogram[256] = {0};

    for (d = 0; d < 3; d++)
    {
        for (im = 0; im < Tmax; im++)
        {
            histogram[im] = 0;
        }
        for (y = 0; y < height; y++)
        {
            for (x = 0; x < width; x++)
            {
                im = (unsigned)p_im[y][x].c[d];
                histogram[im]++;
            }
        }
        thres[d] = IMTHistBiMod (histogram, Tmax, 0.5f) + delta;
    }
    threshold = 0;
    for (d = 0; d < 3; d++)
    {
        threshold += thres[d];
    }
    for (y = 0; y < height; y++ )
    {
        for (x = 0; x < width; x++)
        {
            ims = 0;
            for (d = 0; d < 3; d++)
            {
                im = p_im[y][x].c[d];
                ims += ((im >= thres[d]) ? 255 : 0);
            }
            val = (BYTE) ((ims >= 384) ? 255 : 0 );
            d_im[y][x] = val;
        }
    }

    return threshold;
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterTBlurDiv (IMTpixel** p_im, BYTE** d_im, unsigned int height, unsigned int width, int radius, float sensitivity, float delta)
{
    /* DEPRECATED */
    return IMTFilterTEdgeDiv (p_im, d_im, height, width, radius, 0.0f, sensitivity, delta);
}

///////////////////////////////////////////////////////////////////////////////

int IMTFilterTColor (IMTpixel** p_im, BYTE** d_im, unsigned height, unsigned width, int delta)
{
    unsigned y, x;
    int threshold;

    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            p_im[y][x] = IMTdiffS (p_im[y][x]);
        }
    }
    threshold = IMTFilterTBiModValue (p_im, height, width);
    threshold += delta;
    threshold = IMTFilterThreshold (p_im, d_im, height, width, threshold);

    return threshold;
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterTChistianLayer (IMTpixel** p_im, WORD** t_im, unsigned height, unsigned width, int radius, float sensitivity, int lower_bound, int upper_bound, float delta)
{
    unsigned int x, y, n, y1, x1, y2, x2, yf, xf;
    float imx, imMg, imVg, imm, imv, t, st = 0, sn = 0;
    int threshold = 0;
    WORD val;

    radius = (radius < 0) ? -radius : radius;
    if (upper_bound < lower_bound)
    {
        upper_bound += lower_bound;
        lower_bound = upper_bound - lower_bound;
        upper_bound -= lower_bound;
    }
    lower_bound *= 3;
    upper_bound *= 3;

    imMg = 765.0f;
    imVg = 0.0f;
    for (y = 0; y < height; y++)
    {
        y1 = (y < radius) ? 0 : (y - radius);
        y2 = y + radius + 1;
        y2 = (y2 < height) ? y2 : height;
        for (x = 0; x < width; x++)
        {
            x1 = (x < radius) ? 0 : (x - radius);
            x2 = x + radius + 1;
            x2 = (x2 < width) ? x2 : width;
            imm = 0;
            imv = 0;
            n = 0;
            for (yf = y1; yf < y2; yf++)
            {
                for (xf = x1; xf < x2; xf++)
                {
                    imx = (float)p_im[yf][xf].s;
                    imm += imx;
                    imv += (imx * imx);
                    n++;
                }
            }
            if (n == 0)
            {
                n = 1;
            }
            imm /= n;
            imv /= n;
            imv -= (imm * imm);
            if (imv < 0)
            {
                imv = -imv;
            }
            imv = sqrt(imv);
            if (imm < imMg)
            {
                imMg = imm;
            }
            if (imv > imVg)
            {
                imVg = imv;
            }
        }
    }
    if (imVg == 0)
    {
        imVg = 1;
    }
    for (y = 0; y < height; y++)
    {
        y1 = (y < radius) ? 0 : (y - radius);
        y2 = y + radius + 1;
        y2 = (y2 < height) ? y2 : height;
        for (x = 0; x < width; x++)
        {
            x1 = (x < radius) ? 0 : (x - radius);
            x2 = x + radius + 1;
            x2 = (x2 < width) ? x2 : width;
            imm = 0;
            imv = 0;
            n = 0;
            for (yf = y1; yf < y2; yf++)
            {
                for (xf = x1; xf < x2; xf++)
                {
                    imx = (float)p_im[yf][xf].s;
                    imm += imx;
                    imv += (imx * imx);
                    n++;
                }
            }
            if (n == 0)
            {
                n = 1;
            }
            imm /= n;
            imv /= n;
            imv -= (imm * imm);
            if (imv < 0)
            {
                imv = -imv;
            }
            imv = sqrt(imv);
            imx = (float)p_im[y][x].s;
            imx /= 3;
            if (imx < lower_bound)
            {
                t = lower_bound;
            }
            else if (imx > upper_bound)
            {
                t = upper_bound;
            }
            else
            {
                t = (imm + sensitivity * (imMg - imm + (imv / imVg) * (imm - imMg)) + delta);
            }
            val = Byte3Clamp((int)(t + 0.5));
            st += t;
            sn++;
            t_im[y][x] = val;
        }
    }
    sn = (sn > 0) ? sn : 1;
    threshold = (st / sn);

    return threshold;
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterTChistian (IMTpixel** p_im, BYTE** d_im, unsigned height, unsigned width, int radius, float sensitivity, int lower_bound, int upper_bound, float delta)
{
    int threshold = 0;
    WORD** t_im = TLalloc(height, width);

    threshold = IMTFilterTChistianLayer (p_im, t_im, height, width, radius, sensitivity, lower_bound, upper_bound, delta);
    threshold = IMTFilterThresholdLayer (p_im, t_im, d_im, height, width);

    TLfree(t_im, height);

    return threshold;
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterTDalg (IMTpixel** p_im, BYTE** d_im, unsigned height, unsigned width, int region_size, int lower_bound, int upper_bound, int delta)
{
    unsigned x, y, i, j, Tmax = 768;
    unsigned whg, wwn, iy0, ix0, iyn, ixn, tx, tt;
    unsigned wwidth = region_size;
    float sw, swt, dsr, dst;
    BYTE val;
    unsigned long long histogram[768] = {0};
    int threshold = 0;

    if (upper_bound < lower_bound)
    {
        upper_bound += lower_bound;
        lower_bound = upper_bound - lower_bound;
        upper_bound -= lower_bound;
    }
    lower_bound *= 3;
    upper_bound *= 3;

    whg = (height + wwidth - 1) / wwidth;
    wwn = (width + wwidth - 1) / wwidth;
    for (y = 0; y < whg; y++)
    {
        iy0 = y * wwidth;
        iyn = iy0 + wwidth;
        if (iyn > height)
        {
            iyn = height;
        }
        for (x = 0; x < wwn; x++)
        {
            ix0 = x * wwidth;
            ixn = ix0 + wwidth;
            if (ixn > width)
            {
                ixn = width;
            }
            sw = 0;
            for (j = iy0; j < iyn; j++)
            {
                for (i = ix0; i < ixn; i++)
                {
                    tx = p_im[j][i].s;
                    sw += tx;
                }
            }
            tt = Tmax - 1;
            swt = 0;
            if (sw > 0)
            {
                IMTHist (p_im, histogram, iy0, ix0, iyn, ixn, Tmax);
                while ( swt < sw && tt > 0)
                {
                    dsr = sw - swt;
                    swt += (histogram[tt] * (Tmax - 1));
                    dst = swt - sw;
                    tt--;
                }
                if (dst > dsr)
                {
                    tt++;
                }
            }
            tt += delta;
            for (j = iy0; j < iyn; j++)
            {
                for (i = ix0; i < ixn; i++)
                {
                    tx = p_im[j][i].s;
                    if (tx < lower_bound)
                    {
                        val = (BYTE)0;
                    }
                    else if (tx > upper_bound)
                    {
                        val = (BYTE)255;
                    }
                    else
                    {
                        val = (BYTE)((tx > tt) ? 255 : 0);
                    }
                    d_im[j][i] = val;
                }
            }
            threshold += tt;
        }
    }
    threshold /= whg;
    threshold /= wwn;

    return threshold;
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterTDither (IMTpixel** p_im, BYTE** d_im, unsigned height, unsigned width)
{
    unsigned x, y, dx, dy, i, j, l, m, n;
    unsigned blw, blh, s = 0, sn = 0;
    int threshold = 383;

    int dithermn[16] = {0};
    int dithermdy[16] = {0};
    int dithermdx[16] = {0};
    int erm[16] = {0};
    int val0, val;
    int erv, w;
    BYTE vali;
    static int ditherm[] =
    {
        14, 13, 1, 2,
        4, 6, 11, 9,
        0, 3, 15, 12,
        10, 8, 5, 7
    };


    blw = (width + 3)/4;
    blh = (height + 3)/4;
    for ( y = 0; y < 16; y++ )
    {
        x = ditherm[y];
        dithermn[x] = y;
        dy = y / 4;
        dx = y - dy * 4;
        dithermdy[x] = dy;
        dithermdx[x] = dx;
    }

    for ( i = 0; i < blh; i++ )
    {
        for ( j = 0; j < blw; j++)
        {
            for ( l = 0; l < 16; l++ )
            {
                erm[l] = 0;
            }
            for ( l = 0; l < 16; l++ )
            {
                m = dithermn[l];
                dy = dithermdy[l];
                y = i * 4 + dy;
                if ( y < height)
                {
                    dx = dithermdx[l];
                    x = j * 4 + dx;
                    if ( x < width)
                    {
                        val0 = p_im[y][x].s;
                        val0 += erm[l];
                        val = (( val0 >= threshold ) ? 255 : 0 );
                        erv = val0 - val * 3;
                        vali = (BYTE) (val);
                        if (val > 0)
                        {
                            s++;
                        }
                        sn++;
                        d_im[y][x] = vali;
                        w = 0;
                        if ( dx > 0 )
                        {
                            n = ditherm[m - 1];
                            if (n > l)
                            {
                                w++;
                                w++;
                            }
                        }
                        if (dx < 3 && x < width)
                        {
                            n = ditherm[m + 1];
                            if (n > l)
                            {
                                w++;
                                w++;
                            }
                        }
                        if (dy > 0)
                        {
                            n = ditherm[m - 4];
                            if (n > l)
                            {
                                w++;
                                w++;
                            }
                            if (dx > 0)
                            {
                                n = ditherm[m - 5];
                                if (n > l)
                                {
                                    w++;
                                }
                            }
                            if (dx < 3 && x < width)
                            {
                                n = ditherm[m - 3];
                                if (n > l)
                                {
                                    w++;
                                }
                            }
                        }
                        if (dy < 3 && y < height)
                        {
                            n = ditherm[m + 4];
                            if (n > l)
                            {
                                w++;
                                w++;
                            }
                            if (dx > 0)
                            {
                                n = ditherm[m + 3];
                                if (n > l)
                                {
                                    w++;
                                }
                            }
                            if (dx < 3 && x < width)
                            {
                                n = ditherm[m + 5];
                                if (n > l)
                                {
                                    w++;
                                }
                            }
                        }
                        if (w > 0)
                        {
                            erv /= w;
                            if (dx > 0)
                            {
                                n = ditherm[m - 1];
                                if (n > l)
                                {
                                    erm[n] += (2 * erv);
                                }
                            }
                            if (dx < 3 && x < width)
                            {
                                n = ditherm[m + 1];
                                if (n > l)
                                {
                                    erm[n] += (2 * erv);
                                }
                            }
                            if (dy > 0)
                            {
                                n = ditherm[m - 4];
                                if (n > l)
                                {
                                    erm[n] += (2 * erv);
                                }
                                if (dx > 0)
                                {
                                    n = ditherm[m - 5];
                                    if (n > l)
                                    {
                                        erm[n] += erv;
                                    }
                                }
                                if (dx < 3 && x < width)
                                {
                                    n = ditherm[m - 3];
                                    if (n > l)
                                    {
                                        erm[n] += erv;
                                    }
                                }
                            }
                            if (dy < 3 && y < height)
                            {
                                n = ditherm[m + 4];
                                if (n > l)
                                {
                                    erm[n] += (2 * erv);
                                }
                                if (dx > 0)
                                {
                                    n = ditherm[m + 3];
                                    if (n > l)
                                    {
                                        erm[n] += erv;
                                    }
                                }
                                if (dx < 3 && x < width)
                                {
                                    n = ditherm[m + 5];
                                    if (n > l)
                                    {
                                        erm[n] += erv;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    if (sn > 0)
    {
        threshold = 765 * s / sn;
    }

    return threshold;
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterTDithH (IMTpixel** p_im, BYTE** d_im, unsigned height, unsigned width, int kpg, int delta, unsigned wwidth)
{
    unsigned x, y, i, j, k, l, lmin = 0;
    unsigned whg, wwn, iy0, ix0, iy, ix, tt;
    unsigned herrmin, ww, kw;
    BYTE val;
    int tx, imm, threshold = 0;
    // Knuth D.E. dither matrix
    int hdith[4][4] =
    {
        {  1,  5, 10, 14 },
        {  3,  7,  8, 12 },
        { 13,  9,  6,  2 },
        { 15, 11,  4,  0 }
    };
    int tdith[4][4] =
    {
        { 1, 5, 7, 0 },
        { 3, 6, 2, 0 },
        { 8, 4, 0, 0 },
        { 0, 0, 0, 0 }
    };
    int qdith[4][4] =
    {
        { 1, 2, 0, 0 },
        { 3, 0, 0, 0 },
        { 0, 0, 0, 0 },
        { 0, 0, 0, 0 }
    };
    int dith[4][4];
    int hdithy[2][17], hdithx[2][17], herr, herrp, herrg;
    wwidth = (wwidth < 2) ? 2 : wwidth;
    wwidth = (wwidth > 4) ? 4 : wwidth;
    ww = wwidth * wwidth;
    kw = 768 / ww;
    if (wwidth == 2)
    {
        for (y = 0; y < wwidth; y++)
        {
            for (x = 0; x < wwidth; x++)
            {
                dith[y][x] = qdith[y][x];
            }
        }
    }
    else if (wwidth == 3)
    {
        for (y = 0; y < wwidth; y++)
        {
            for (x = 0; x < wwidth; x++)
            {
                dith[y][x] = tdith[y][x];
            }
        }
    }
    else
    {
        for (y = 0; y < wwidth; y++)
        {
            for (x = 0; x < wwidth; x++)
            {
                dith[y][x] = hdith[y][x];
            }
        }
    }
    for (y = 0; y < wwidth; y++)
    {
        for (x = 0; x < wwidth; x++)
        {
            l = dith[y][x] + 1;
            hdithy[0][l] = y;
            hdithx[0][l] = x;
            hdithy[1][l] = y;
            hdithx[1][l] = wwidth - x - 1;
        }
    }
    whg = (height + wwidth - 1) / wwidth;
    wwn = (width + wwidth - 1) / wwidth;
    for (y = 0; y < whg; y++)
    {
        iy0 = y * wwidth;
        for (x = 0; x < wwn; x++)
        {
            ix0 = x * wwidth;
            k = (y + x) % 2;
            imm = 0;
            for (l = 1; l < (ww + 1); l++)
            {
                j = hdithy[k][l];
                i = hdithx[k][l];
                iy = iy0 + j;
                ix = ix0 + i;
                tx = (iy < height && ix < width) ? p_im[iy][ix].s : 765;
                tx += delta;
                imm += tx;
            }
            imm /= ww;
            herrp = herrg = 0;
            for (l = 1; l < (ww + 1); l++)
            {
                j = hdithy[k][l];
                i = hdithx[k][l];
                iy = iy0 + j;
                ix = ix0 + i;
                tx = (iy < height && ix < width) ? p_im[iy][ix].s : 765;
                tx += delta;
                herrg += (765 - tx);
                tx -= imm;
                tx *= kpg;
                tx += imm;
                tx = (int)Byte3Clamp(tx);
                herrp += (765 - tx);
            }
            herrmin = herrp + herrg;
            lmin = 0;
            for (l = 1; l < (ww + 1); l++)
            {
                j = hdithy[k][l];
                i = hdithx[k][l];
                iy = iy0 + j;
                ix = ix0 + i;
                tx = (iy < height && ix < width) ? p_im[iy][ix].s : 765;
                tx += delta;
                tx -= imm;
                tx *= kpg;
                tx += imm;
                tx = (int)Byte3Clamp(tx);
                herrp += (tx + tx - 765);
                herrg -= 765;
                herr = (herrp < 0) ? (-herrp) : herrp;
                herr += (herrg < 0) ? (-herrg) : herrg;
                if (herr < herrmin)
                {
                    herrmin = herr;
                    lmin = l;
                }
            }
            for (l = 1; l < (ww + 1); l++)
            {
                j = hdithy[k][l];
                i = hdithx[k][l];
                iy = iy0 + j;
                ix = ix0 + i;
                if (iy < height && ix < width)
                {
                    val = (BYTE) ( ( l > lmin ) ? 255 : 0 );
                    d_im[iy][ix] = val;
                }
            }
            tt = kw * (ww - lmin);
            threshold += tt;
        }
    }
    threshold /= whg;
    threshold /= wwn;

    return threshold;
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterTDithO (IMTpixel** p_im, BYTE** d_im, unsigned height, unsigned width, int kpg, int delta)
{
    unsigned x, y, i, j, l, lmin = 0;
    unsigned whg, wwn, iy0, ix0, iy, ix, tt;
    unsigned wwidth = 8, herrmin, ww = wwidth * wwidth, kw = 768 / ww;
    BYTE val;
    int tx, imm, threshold = 0;
    // Knuth D.E. dither matrix
    int odith[8][8] =
    {
        { 35, 48, 40, 32, 28, 15, 23, 31 },
        { 43, 59, 56, 52, 20,  4,  7, 11 },
        { 51, 62, 60, 44, 12,  1,  3, 19 },
        { 38, 46, 54, 36, 25, 17,  9, 27 },
        { 29, 14, 22, 30, 34, 49, 41, 33 },
        { 21,  5,  6, 10, 42, 58, 57, 53 },
        { 13,  0,  2, 18, 50, 63, 61, 45 },
        { 24, 16,  8, 26, 39, 47, 55, 37 }
    };
    int odithy[65], odithx[65], herr, herrp, herrg;
    for (y = 0; y < wwidth; y++)
    {
        for (x = 0; x < wwidth; x++)
        {
            l = odith[y][x] + 1;
            odithy[l] = y;
            odithx[l] = x;
        }
    }
    whg = (height + wwidth - 1) / wwidth;
    wwn = (width + wwidth - 1) / wwidth;
    for (y = 0; y < whg; y++)
    {
        iy0 = y * wwidth;
        for (x = 0; x < wwn; x++)
        {
            ix0 = x * wwidth;
            imm = 0;
            for (l = 1; l < (ww + 1); l++)
            {
                j = odithy[l];
                i = odithx[l];
                iy = iy0 + j;
                ix = ix0 + i;
                tx = (iy < height && ix < width) ? p_im[iy][ix].s : 765;
                tx += delta;
                imm += tx;
            }
            imm /= ww;
            herrp = herrg = 0;
            for (l = 1; l < (ww + 1); l++)
            {
                j = odithy[l];
                i = odithx[l];
                iy = iy0 + j;
                ix = ix0 + i;
                tx = (iy < height && ix < width) ? p_im[iy][ix].s : 765;
                tx += delta;
                herrg += (765 - tx);
                tx -= imm;
                tx *= kpg;
                tx += imm;
                tx = (int)Byte3Clamp(tx);
                herrp += (765 - tx);
            }
            herrmin = herrp + herrg;
            lmin = 0;
            for (l = 1; l < (ww + 1); l++)
            {
                j = odithy[l];
                i = odithx[l];
                iy = iy0 + j;
                ix = ix0 + i;
                tx = (iy < height && ix < width) ? p_im[iy][ix].s : 765;
                tx += delta;
                tx -= imm;
                tx *= kpg;
                tx += imm;
                tx = (int)Byte3Clamp(tx);
                herrp += (tx + tx - 765);
                herrg -= 765;
                herr = (herrp < 0) ? (-herrp) : herrp;
                herr += (herrg < 0) ? (-herrg) : herrg;
                if (herr < herrmin)
                {
                    herrmin = herr;
                    lmin = l;
                }
            }
            for (l = 1; l < (ww + 1); l++)
            {
                j = odithy[l];
                i = odithx[l];
                iy = iy0 + j;
                ix = ix0 + i;
                if (iy < height && ix < width)
                {
                    val = (BYTE) ( ( l > lmin ) ? 255 : 0 );
                    d_im[iy][ix] = val;
                }
            }
            tt = kw * (ww - lmin);
            threshold += tt;
        }
    }
    threshold /= whg;
    threshold /= wwn;

    return threshold;
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterTDithBayer (IMTpixel** p_im, BYTE** d_im, unsigned height, unsigned width, int delta)
{
    unsigned int y, x, yb, xb;
    unsigned int wwidth = 8;
    int thres, pix;
    int threshold = 383;
    // Bayer dither matrix
    int bdith[8][8] =
    {
        { 0, 32,  8, 40,  2, 34, 10, 42},
        {48, 16, 56, 24, 50, 18, 58, 26},
        {12, 44,  4, 36, 14, 46,  6, 38},
        {60, 28, 52, 20, 62, 30, 54, 22},
        { 3, 35, 11, 43,  1, 33,  9, 41},
        {51, 19, 59, 27, 49, 17, 57, 25},
        {15, 47,  7, 39, 13, 45,  5, 37},
        {63, 31, 55, 23, 61, 29, 53, 21}
    };

    for (y = 0; y < wwidth; y++)
    {
        for (x = 0; x < wwidth; x++)
        {
            bdith[y][x] *= 12;
        }
    }
    for (y = 0; y < height; y++)
    {
        yb = y % wwidth;
        for (x = 0; x < width; x++)
        {
            xb = x % wwidth;
            thres = bdith[yb][xb];
            pix = p_im[y][x].s;
            pix += delta;
            d_im[y][x] = (BYTE)((pix < thres) ? 0 : 255);
        }
    }

    return threshold;
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterTDithDots (IMTpixel** p_im, BYTE** d_im, unsigned height, unsigned width, int delta)
{
    unsigned int y, x, yb, xb;
    unsigned int wwidth = 8;
    int thres, pix;
    int threshold = 383;
    // Dots dither matrix
    int ddith[8][8] =
    {
        {13,  9,  5, 12, 18, 22, 26, 19},
        { 6,  1,  0,  8, 25, 30, 31, 23},
        {10,  2,  3,  4, 21, 29, 28, 27},
        {14,  7, 11, 15, 17, 24, 20, 16},
        {18, 22, 26, 19, 13,  9,  5, 12},
        {25, 30, 31, 23,  6,  1,  0,  8},
        {21, 29, 28, 27, 10,  2,  3,  4},
        {17, 24, 20, 16, 14,  7, 11, 15}
    };

    for (y = 0; y < wwidth; y++)
    {
        for (x = 0; x < wwidth; x++)
        {
            ddith[y][x] *= 24;
        }
    }
    for (y = 0; y < height; y++)
    {
        yb = y % wwidth;
        for (x = 0; x < width; x++)
        {
            xb = x % wwidth;
            thres = ddith[yb][xb];
            pix = p_im[y][x].s;
            pix += delta;
            d_im[y][x] = (BYTE)((pix < thres) ? 0 : 255);
        }
    }

    return threshold;
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterTDjVuL (IMTpixel** p_im, BYTE** m_im, IMTpixel** fg_im, IMTpixel** bg_im, unsigned height, unsigned width, unsigned bgs, unsigned fgs, unsigned level, int wbmode, float anisotropic, float doverlay, unsigned fposter)
{
    unsigned widthbg = (width + bgs - 1) / bgs;
    unsigned heightbg = (height + bgs - 1) / bgs;
    unsigned widthfg = (widthbg + fgs - 1) / fgs;
    unsigned heightfg = (heightbg + fgs - 1) / fgs;
    unsigned whcp, y, x, d, l, i, j, y0, x0, y1, x1, y0b, x0b, y1b, x1b, yb, xb, yf, xf, blsz;
    float immean, imwb;
    BYTE fgbase, bgbase;
    unsigned cnth, cntw;
    IMTpixel pim, gim, tim, fgim, bgim;
    float fgdist, bgdist, kover, fgpart, bgpart, lpart, fgk = 1.0;
    unsigned maskbl, maskover, bgsover, fgsum[3], bgsum[3], fgnum, bgnum;

    IMTpixel** fgt_im = IMTalloc(heightbg, widthbg);

    whcp = height;
    whcp += width;
    whcp /= 2;
    blsz = 1;
    if (level == 0)
    {
        while (bgs * blsz < whcp)
        {
            level++;
            blsz *= 2;
        }
    }
    else
    {
        for (l = 0; l < level; l++)
        {
            blsz *= 2;
        }
    }
    blsz /= 2;
    immean = IMTmean(p_im, height, width);
    imwb = IMTwb(p_im, immean, height, width);
    if (anisotropic == 0)
    {
        fgk = sqrt(1.5 - imwb);
    }
    if (doverlay < 0)
    {
        doverlay = 0;
    }
    kover = doverlay + 1.0;
    if (wbmode == 0)
    {
        if (imwb < 0.5)
        {
            wbmode = -1;
        }
        else
        {
            wbmode = 1;
        }
    }
    if (wbmode < 0)
    {
        fgbase = 255;
        bgbase = 0;
    }
    else
    {
        fgbase = 0;
        bgbase = 255;
    }
    for (y = 0; y < heightbg; y++)
    {
        for (x = 0; x < widthbg; x++)
        {
            fgt_im[y][x] = IMTset(fgbase, fgbase, fgbase);
            bg_im[y][x] = IMTset(bgbase, bgbase, bgbase);
        }
    }
    for (l = 0; l < level; l++)
    {
        cnth = (heightbg + blsz - 1) / blsz;
        cntw = (widthbg + blsz - 1) / blsz;
        maskbl = bgs * blsz;
        maskover = (kover * maskbl);
        bgsover = (kover * blsz);
        lpart = (float)(level - l) / (float)level;
        for (i = 0; i < cnth; i++)
        {
            y0 = i * maskbl;
            y1 = (((y0 + maskover) < height) ? (y0 + maskover) : height);
            y0b = i * blsz;
            y1b = (((y0b + bgsover) < heightbg) ? (y0b + bgsover) : heightbg);
            for (j = 0; j < cntw; j++)
            {
                x0 = j * maskbl;
                x1 = (((x0 + maskover) < width) ? (x0 + maskover) : width);
                x0b = j * blsz;
                x1b = (((x0b + bgsover) < widthbg) ? (x0b + bgsover) : widthbg);

                gim = IMTmeanIc(p_im, y0, x0, y1, x1);

                fgim = IMTmeanIc(fgt_im, y0b, x0b, y1b, x1b);
                bgim = IMTmeanIc(bg_im, y0b, x0b, y1b, x1b);

                fgdist = IMTdist(gim, fgim);
                bgdist = IMTdist(gim, bgim);

                fgk = (fgdist + bgdist);
                if (fgk > 0)
                {
                    fgk = (bgdist - fgdist) / fgk;
                    fgk *= anisotropic;
                    fgk = exp(fgk);
                }
                else
                {
                    fgk = 1.0;
                }
                for (d = 0; d < 3; d++)
                {
                    fgsum[d] = 0;
                    bgsum[d] = 0;
                }
                fgnum = 0;
                bgnum = 0;
                for (y = y0; y < y1; y++)
                {
                    for (x = x0; x < x1; x++)
                    {
                        pim = p_im[y][x];
                        //tim = IMTrefilter1p (pim, gim);
                        tim = pim;
                        fgdist = IMTdist(tim, fgim);
                        bgdist = IMTdist(tim, bgim);
                        if (fgdist * fgk < bgdist)
                        {
                            for (d = 0; d < 3; d++)
                            {
                                fgsum[d] += (int)pim.c[d];
                            }
                            fgnum++;
                        }
                        else
                        {
                            for (d = 0; d < 3; d++)
                            {
                                bgsum[d] += (int)pim.c[d];
                            }
                            bgnum++;
                        }
                    }
                }
                if (fgnum > 0)
                {
                    for (d = 0; d < 3; d++)
                    {
                        fgsum[d] /= fgnum;
                        fgim.c[d] = fgsum[d];
                    }
                }
                if (bgnum > 0)
                {
                    for (d = 0; d < 3; d++)
                    {
                        bgsum[d] /= bgnum;
                        bgim.c[d] = bgsum[d];
                    }
                }
                pim = IMTmeanIc(p_im, y0, x0, y1, x1);

                fgdist = IMTdist(pim, fgim);
                bgdist = IMTdist(pim, bgim);
                fgpart = 1.0;
                bgpart = 1.0;
                if ((fgdist + bgdist)> 0)
                {
                    fgpart += (2.0 * fgdist / (fgdist + bgdist));
                    bgpart += (2.0 * bgdist / (fgdist + bgdist));
                }
                fgpart *= lpart;
                bgpart *= lpart;
                fgim = IMTaverageIc(fgt_im, fgim, y0b, x0b, y1b, x1b, fgpart);
                bgim = IMTaverageIc(bg_im, bgim, y0b, x0b, y1b, x1b, bgpart);
            }
        }
        blsz /= 2;
    }
    IMTFilterSReduce (fgt_im, fg_im, heightbg, widthbg, fgs);
    if (fposter != 0)
    {
        (void)IMTFilterPosterize(fg_im, heightfg, widthfg, fposter);
    }
    for (y = 0; y < height; y++)
    {
        yb = y / bgs;
        yf = yb / fgs;
        for (x = 0; x < width; x++)
        {
            xb = x /bgs;
            xf = xb / fgs;
            pim = p_im[y][x];
            fgim = fg_im[yf][xf];
            bgim = bg_im[yb][xb];
            bgdist = IMTdist(pim, bgim);
            fgdist = IMTdist(pim, fgim);
            m_im[y][x] = (BYTE) ((fgdist < bgdist) ? 0 : 255);
        }
    }

    IMTfree(fgt_im, heightbg);

    return wbmode;
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterTEdge (IMTpixel** p_im, BYTE** d_im, unsigned int height, unsigned int width, int radius, float delta)
{
    int threshold = 0;
    IMTpixel** b_im = IMTalloc(height, width);

    IMTFilterGaussBlur (p_im, b_im, height, width, radius);
    IMTFilterSEdge (p_im, b_im, height, width);

    IMTfree(b_im, height);

    threshold = IMTFilterTBiMod (p_im, d_im, height, width, delta);

    return threshold;
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterTEdgeDiv (IMTpixel** p_im, BYTE** d_im, unsigned int height, unsigned int width, int radius, float sensitivity, float sensitivityblur, float delta)
{
    unsigned int y, x;
    float origin, blur, edge, edgeplus, edgenorm, retval;
    int threshold = 0;
    IMTpixel** b_im = IMTalloc(height, width);

    IMTFilterGaussBlur (p_im, b_im, height, width, radius);
    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            origin = p_im[y][x].s;
            blur = b_im[y][x].s;
            retval = origin;
            if (sensitivity > 0.0f)
            {
                edge = (origin + 1.0f) / (blur + 1.0f) - 0.5f;
                edgeplus = origin * edge;
                retval = sensitivity * edgeplus + (1.0f - sensitivity) * origin;
            }
            if (sensitivityblur > 0.0f)
            {
                edge = (blur + 1.0f) / (retval + 1.0f) - 0.5f;
                edgenorm = sensitivityblur * edge + (1.0f - sensitivityblur);
                retval = (edgenorm > 0.0f) ? (origin / edgenorm) : origin;
            }
            retval = (retval < 0.0f) ? 0.0f : (retval < 765.0f) ? retval : 765.0f;
            b_im[y][x].s = (WORD)(retval + 0.5f);
        }
    }
    threshold = IMTFilterTBiMod (b_im, d_im, height, width, delta);
    IMTfree(b_im, height);

    return threshold;
}

///////////////////////////////////////////////////////////////////////////////

int IMTFilterTEdgePlus (IMTpixel** p_im, BYTE** d_im, unsigned int height, unsigned int width, int radius, float sensitivity, float delta)
{
    /* DEPRECATED */
    return IMTFilterTEdgeDiv (p_im, d_im, height, width, radius, sensitivity, 0.0f, delta);
}

///////////////////////////////////////////////////////////////////////////////

int IMTFilterTEntValue (IMTpixel** p_im, unsigned height, unsigned width)
{
    unsigned i, t, cn = 768;
    float sum = 0, pTB, hhB, pTW, hhW, jMax, j, hn;
    float pT[768] = {0};
    float epsilon = 0;
    float hB[768] = {0};
    float hW[768] = {0};
    int threshold = 0;
    unsigned long long histogram[768] = {0};

    IMTHist (p_im, histogram, 0, 0, height, width, cn);
    sum = (float)(width * height);
    hn = 1.0f / sum;

    pT[0] = hn * histogram[0];
    for (i = 1; i < cn; i++)
    {
        pT[i] = pT[i - 1] + hn * (float)histogram[i];
    }

    for (t = 0; t < cn; t++)
    {
        pTB = pT[t];
        if (pTB > epsilon)
        {
            hhB = 0;
            for (i = 0; i <= t; i++)
            {
                if (histogram[i] > 0)
                {
                    hhB -= (hn * (float)histogram[i]) / pTB * log(hn * (float)histogram[i] / pTB);
                }
            }
            hB[t] = hhB;
        }
        else
        {
            hB[t] = 0;
        }

        pTW = 1 - pT[t];
        if (pTW > epsilon)
        {
            hhW = 0;
            for (i = t + 1; i < cn; ++i)
            {
                if (histogram[i] > 0)
                {
                    hhW -= hn * (float)histogram[i] / pTW * log(hn * (float)histogram[i] / pTW);
                }
            }
            hW[t] = hhW;
        }
        else
        {
            hW[t] = 0;
        }
    }

    // Find histogram index with maximum entropy
    jMax = hB[0] + hW[0];

    for (t = 1; t < cn; ++t)
    {
        j = hB[t] + hW[t];
        if (j > jMax)
        {
            jMax = j;
            threshold = t;
        }
    }

    return threshold;
}
////////////////////////////////////////////////////////////////////////////////

int IMTFilterTEnt (IMTpixel** p_im, BYTE** d_im, unsigned height, unsigned width)
{
    int threshold = 0;

    threshold = IMTFilterTEntValue (p_im, height, width);
    threshold = IMTFilterThreshold (p_im, d_im, height, width, threshold);

    return threshold;
}
////////////////////////////////////////////////////////////////////////////////

int IMTFilterTEqBrightValue (IMTpixel** p_im, unsigned height, unsigned width)
{
    unsigned i, tt, Tmax = 768;
    float imx, sw, swt, dsr = 0, dst = 0;
    unsigned long long histogram[768] = {0};
    int threshold = 0;

    IMTHist (p_im, histogram, 0, 0, height, width, Tmax);
    sw = 0;
    for (i = 0; i < Tmax; i++)
    {
        sw += (float)(histogram[i] * i);
    }
    tt = Tmax - 1;
    swt = 0;
    while ( swt < sw && tt > 0)
    {
        dsr = sw - swt;
        imx = histogram[tt];
        swt += (imx * (Tmax - 1));
        dst = swt - sw;
        tt--;
    }
    if (dst > dsr)
    {
        tt++;
        swt -= (imx * (Tmax - 1));
    }
    threshold = tt;

    return threshold;
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterTEqBright (IMTpixel** p_im, BYTE** d_im, unsigned height, unsigned width)
{
    int threshold = 0;

    threshold = IMTFilterTEqBrightValue (p_im, height, width);
    threshold = IMTFilterThreshold (p_im, d_im, height, width, threshold);

    return threshold;
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterGatosBG (IMTpixel** p_im, BYTE** d_im, IMTpixel** bg_im, unsigned height, unsigned width, int radius)
{
    int x, y, i, j, xt, yt;
    unsigned d;
    float imx;
    float bin_sum;
    float p_sum[3] = {0.0};
    BYTE val;
    int h = height;
    int w = width;
    int threshold = 0, st = 0, sn = 0;

    radius = (radius < 0) ? -radius : radius;
    for (y = 0; y < h; y++)
    {
        for (x = 0; x < w; x++)
        {
            if (d_im[y][x] != 0)
            {
                // white
                bg_im[y][x] = p_im[y][x];
            }
            else
            {
                for (d = 0; d < 3; d++)
                {
                    p_sum[d] = 0;
                }
                bin_sum = 0;
                for (i = -radius; i <= radius; i++)
                {
                    yt = y + i;
                    if (yt >= 0 && yt < h)
                    {
                        for (j = -radius; j <= radius; j++)
                        {
                            xt = x + j;
                            if (xt >= 0 && xt < w)
                            {
                                if (d_im[yt][xt] != 0) // white
                                {
                                    bin_sum++;
                                    for (d = 0; d < 3; d++)
                                    {
                                        imx = (float)p_im[yt][xt].c[d];
                                        p_sum[d] += imx;
                                    }
                                }
                            }
                        }
                    }
                }
                if (bin_sum > 0)
                {
                    for (d = 0; d < 3; d++)
                    {
                        imx = p_sum[d] / (float)bin_sum;
                        val = ByteClamp((int)(imx + 0.5));
                        bg_im[y][x].c[d] = val;
                    }
                    bg_im[y][x] = IMTcalcS (bg_im[y][x]);
                }
                else
                {
                    bg_im[y][x] = IMTset(255, 255, 255);
                    st++;
                }
            }
            sn++;
        }
    }

    if (sn > 0)
    {
        threshold = 768 * st / sn;
    }

    return threshold;
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterTGatos (IMTpixel** p_im, BYTE** d_im, IMTpixel** bg_im, BYTE** g_im, unsigned height, unsigned width, float q, float p1, float p2)
{
    unsigned x, y, t, s;
    float imx, imk;
    int sum = 0;
    unsigned delta_numerator = 0;
    unsigned delta_denominator = 0;
    int bin_sum = 0;
    float delta, b;
    BYTE val;
    int threshold = 0, st = 0, sn = 0;

    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            t = (unsigned)bg_im[y][x].s;
            s = (unsigned)p_im[y][x].s;
            if (t > s)
            {
                delta_numerator += (t - s);
            }
        }
    }

    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            if (d_im[y][x] == 0)
            {
                // black
                delta_denominator++;
            }
        }
    }

    delta = (float)delta_numerator;
    delta /= (float)delta_denominator;
    delta /= 3;

    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            if(d_im[y][x] != 0) // white
            {
                bin_sum++;
                sum += (float)bg_im[y][x].s;
            }
        }
    }
    sum /= 3;

    b = sum / bin_sum;

    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            imx = (float)p_im[y][x].s;
            imk = (float)bg_im[y][x].s;
            imx /= 3;
            imk /= 3;
            imx = imk -imx;
            imk = q * delta * (((1 - p2) / (1 + exp(((-4 * imk) / (b * (1 - p1))) + ((2 * (1 + p1)) / (1 - p1))))) + p2);
            val = (BYTE) ((imx > imk) ? 0 : 255);
            if (val > 0)
            {
                st++;
            }
            g_im[y][x] = val;
            sn++;
        }
    }
    if (sn > 0)
    {
        threshold = 765 * st / sn;
    }

    return threshold;
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterTGradValue (IMTpixel** p_im, unsigned height, unsigned width)
{
    int x, y, xp, xn, yp, yn;
    int im, imp, imn;
    float GX, GY, G, SG = 0, SGI = 0, SI = 0, SN = 0, S = 0;
    int threshold = 0;
    int h = (int)height;
    int w = (int)width;

    for (y = 0; y < h; y++)
    {
        yp = y -1;
        if (yp < 0)
        {
            yp = 0;
        }
        yn = y + 1;
        if (yn >= h)
        {
            yn = h - 1;
        }
        for (x = 0; x < w; x++)
        {
            im = (int)p_im[y][x].s;
            xp = x -1;
            if (xp < 0)
            {
                xp = 0;
            }
            xn = x + 1;
            if (xn >= w)
            {
                xn = w - 1;
            }
            imp = (int)p_im[y][xp].s;
            imn = (int)p_im[y][xn].s;
            GX = (float)(imn - imp);
            if (GX < 0)
            {
                GX = -GX;
            }
            imp = (int)p_im[yp][x].s;
            imn = (int)p_im[yn][x].s;
            GY = (float)(imn - imp);
            if (GY < 0)
            {
                GY = -GY;
            }
            G = sqrt(GX * GX + GY * GY);
            SG += G;
            SGI += ((float)im * G);
            SI += (float)im;
            SN++;
        }
    }

    if (SG==0)
    {
        S = SI / SN;
    }
    else
    {
        S = SGI / SG;
    }

    threshold = (int)(S);

    return threshold;
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterTGrad (IMTpixel** p_im, BYTE** d_im, unsigned height, unsigned width)
{
    int threshold = 0;

    threshold = IMTFilterTGradValue (p_im, height, width);
    threshold = IMTFilterThreshold (p_im, d_im, height, width, threshold);

    return threshold;
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterTGravureLayer (IMTpixel** p_im, WORD** t_im, unsigned height, unsigned width, int radius, float sensitivity, int lower_bound, int upper_bound, float delta)
{
    unsigned int x, y, n;
    unsigned int y1, x1, y2, x2, yf, xf;
    float imx, imMg, imVg, imm, imv, t, st = 0, sn = 0;
    int threshold = 0;
    WORD val;

    radius = (radius < 0) ? -radius : radius;
    if (upper_bound < lower_bound)
    {
        upper_bound += lower_bound;
        lower_bound = upper_bound - lower_bound;
        upper_bound -= lower_bound;
    }
    lower_bound *= 3;
    upper_bound *= 3;

    imMg = IMTFilterTBiModValue (p_im, height, width);
    imVg = 0;
    n = 0;
    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            imx = (float)p_im[y][x].s;
            imx -= imMg;
            imVg += (imx * imx);
            n++;
        }
    }
    if (n == 0)
    {
        n = 1;
    }
    imVg /= n;
    if (imVg < 0)
    {
        imVg = -imVg;
    }
    imVg = sqrt(imVg);
    if (imVg == 0)
    {
        imVg = 1;
    }
    for (y = 0; y < height; y++)
    {
        y1 = (y < radius) ? 0 : (y - radius);
        y2 = y + radius + 1;
        y2 = (y2 < height) ? y2 : height;
        for (x = 0; x < width; x++)
        {
            x1 = (x < radius) ? 0 : (x - radius);
            x2 = x + radius + 1;
            x2 = (x2 < width) ? x2 : width;
            imm = 0;
            imv = 0;
            n = 0;
            for (yf = y1; yf < y2; yf++)
            {
                for (xf = x1; xf < x2; xf++)
                {
                    imx = (float)p_im[yf][xf].s;
                    imm += imx;
                    imv += (imx * imx);
                    n++;
                }
            }
            if (n == 0)
            {
                n = 1;
            }
            imm /= n;
            imv /= n;
            imv -= (imm * imm);
            if (imv < 0)
            {
                imv = -imv;
            }
            imv = sqrt(imv);
            imx = (float)p_im[y][x].s;
            imx /= 3;
            if (imx < lower_bound)
            {
                t = lower_bound;
            }
            else if (imx > upper_bound)
            {
                t = upper_bound;
            }
            else
            {
                t = (imm + sensitivity * (imMg - imm + (imv / imVg) * (imm - imMg)) + delta);
            }
            val = Byte3Clamp((int)(t + 0.5));
            st += t;
            sn++;
            t_im[y][x] = val;
        }
    }
    if (sn == 0)
    {
        sn = 1;
    }
    threshold =  (st / sn);

    return threshold;
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterTGravure (IMTpixel** p_im, BYTE** d_im, unsigned height, unsigned width, int radius, float sensitivity, int lower_bound, int upper_bound, float delta)
{
    int threshold = 0;
    WORD** t_im = TLalloc(height, width);

    threshold = IMTFilterTGravureLayer (p_im, t_im, height, width, radius, sensitivity, lower_bound, upper_bound, delta);
    threshold = IMTFilterThresholdLayer (p_im, t_im, d_im, height, width);

    TLfree(t_im, height);

    return threshold;
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterTHalftone2 (IMTpixel** p_im, BYTE** d_im, unsigned height, unsigned width)
{
    unsigned y, x, im, n, s, threshold = 0;

    n = 0;
    s = 0;
    for (y = 0; y < height; y++ )
    {
        for (x = 0; x < width; x++)
        {
            im = p_im[y][x].s;
            n++;
            if (im < 153)
            {
                d_im[2 * y][2 * x] = 0;
                d_im[2 * y + 1][2 * x] = 0;
                d_im[2 * y][2 * x + 1] = 0;
                d_im[2 * y + 1][2 * x + 1] = 0;
            }
            else if (im < 306)
            {
                d_im[2 * y][2 * x] = 255;
                d_im[2 * y + 1][2 * x] = 0;
                d_im[2 * y][2 * x + 1] = 0;
                d_im[2 * y + 1][2 * x + 1] = 0;
                s++;
            }
            else if (im < 459)
            {
                d_im[2 * y][2 * x] = 255;
                d_im[2 * y + 1][2 * x] = 0;
                d_im[2 * y][2 * x + 1] = 0;
                d_im[2 * y + 1][2 * x + 1] = 255;
                s++;
                s++;
            }
            else if (im < 612)
            {
                d_im[2 * y][2 * x] = 255;
                d_im[2 * y][2 * x + 1] = 255;
                d_im[2 * y + 1][2 * x] = 0;
                d_im[2 * y + 1][2 * x + 1] = 255;
                s++;
                s++;
                s++;
            }
            else
            {
                d_im[2 * y][2 * x] = 255;
                d_im[2 * y][2 * x + 1] = 255;
                d_im[2 * y + 1][2 * x] = 255;
                d_im[2 * y + 1][2 * x + 1] = 255;
                s++;
                s++;
                s++;
                s++;
            }
        }
    }
    threshold = s * 64 / n;

    return threshold;
}


////////////////////////////////////////////////////////////////////////////////

int IMTFilterTJanniValue (IMTpixel** p_im, unsigned height, unsigned width)
{
    unsigned i, cn = 768;
    unsigned gmin=0, gmax=256, gmid;
    float spg = 0;
    unsigned long long histogram[768] = {0};
    int threshold = 0;

    IMTHist (p_im, histogram, 0, 0, height, width, cn);

    i = 0;
    while (histogram[i] == 0 && i < cn)
    {
        i++;
    }
    gmin = i;
    i = cn - 1;
    while (histogram[i] == 0 && i > 0)
    {
        i--;
    }
    gmax = i;

    gmid = (gmin + gmax) / 2;

    for (i = (gmid + 1); i <= gmax; i++)
    {
        spg += (float)histogram[i];
    }
    spg /= (float)(height * width);
    threshold = (int)(gmin + ((gmax - gmin) * spg));
    threshold = (threshold + gmid) / 2;

    return threshold;
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterTJanni (IMTpixel** p_im, BYTE** d_im, unsigned height, unsigned width)
{
    int threshold = 0;

    threshold = IMTFilterTJanniValue (p_im, height, width);
    threshold = IMTFilterThreshold (p_im, d_im, height, width, threshold);

    return threshold;
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterTKMeans (IMTpixel** p_im, BYTE** d_im, unsigned height, unsigned width, unsigned knum, unsigned iters)
{
    unsigned knum2 = knum / 2;
    int y, x, i, j;
    float kptb, kptn;

    unsigned d, knm, ct, t, st = 0, sn = 0;
    float kdmin, cstep;
    unsigned imm[3] = {0};
    BYTE val;
    int h = (int)height;
    int w = (int)width;
    int threshold = 0;
    int k = 0;
    int dl, dls = 1;
    float** kpt;
    float** ksum;
    float* kdist;

    cstep = 255.0 / (knum - 1);

    kpt = (float**)malloc(knum * sizeof(float*));
    for (y = 0; y < (int)knum; y++)
    {
        kpt[y] = (float*)malloc(4 * sizeof(float));
    }
    ksum = (float**)malloc(knum * sizeof(float*));
    for (y = 0; y < (int)knum; y++)
    {
        ksum[y] = (float*)malloc(4 * sizeof(float));
    }
    kdist = (float*)malloc(knum * sizeof(float));

    for (i = 0; i < (int)knum; i++)
    {
        ct = (unsigned)(cstep * i + 0.5);
        for (j = 0; j < 3; j++)
        {
            kpt[i][j] = ct;
        }
    }
    for (y = 0; y < h; y++)
    {
        for (x = 0; x < w; x++)
        {
            for (d = 0; d < 3; d++)
            {
                imm[d] += (unsigned)p_im[y][x].c[d];
            }
        }
    }
    for (d = 0; d < 3; d++)
    {
        imm[d] /= height;
        imm[d] /= width;
    }

    while (dls > 0 && k < (int)iters)
    {
        for (i = 0; i < (int)knum; i++)
        {
            for (j = 0; j < 4; j++)
            {
                ksum[i][j] = 0;
            }
        }
        for (y = 0; y < h; y++)
        {
            for (x = 0; x < w; x++)
            {
                for (i = 0; i < (int)knum; i++)
                {
                    kdist[i] = IMTdist3c2p(p_im[y][x], kpt[i]);
                }
                kdmin = 10.0;
                knm = 0;
                for (i = 0; i < (int)knum; i++)
                {
                    if (kdist[i] < kdmin)
                    {
                        kdmin = kdist[i];
                        knm = (unsigned)i;
                    }
                }
                for (j = 0; j < 3; j++)
                {
                    ksum[knm][j] += (float)p_im[y][x].c[j];
                }
                ksum[knm][3]++;
            }
        }
        dls = 0;
        for (i = 0; i < (int)knum; i++)
        {
            if (ksum[i][3] > 0)
            {
                for (j = 0; j < 3; j++)
                {
                    kptb = kpt[i][j];
                    kptn = (ksum[i][j] / ksum[i][3]);
                    kptn = (kptb + kptn) / 2;
                    dl = kptb - kptn;
                    dl *= dl;
                    dls += dl;
                    kpt[i][j] = kptn;
                }
            }
            else
            {
                for (j = 0; j < 3; j++)
                {
                    kptb = kpt[i][j];
                    kptn = (kptb * iters + imm[j]) / (iters + 1);
                    dl = kptb - kptn;
                    dl *= dl;
                    dls += dl;
                    kpt[i][j] = kptn;
                }
            }
        }
        k++;
    }

    for (j = 0; j < 3; j++)
    {
        ksum[0][j] = 0;
        for (i = 0; i < (int)knum; i++)
        {
            ksum[0][j] += kpt[i][j];
        }
        ksum[0][j] /= knum;
    }
    kptb = 0;
    for (j = 0; j < 3; j++)
    {
        kptb += (kpt[knum2][j] - ksum[0][j]);
    }
    for (i = 0; i < (int)knum; i++)
    {
        kpt[i][3] = 0;
        for (j = 0; j < 3; j++)
        {
            kpt[i][3] += kpt[i][j];
        }
    }
    for (y = 0; y < h; y++)
    {
        for (x = 0; x < w; x++)
        {
            for (i = 0; i < (int)knum; i++)
            {
                kdist[i] = IMTdist3c2p(p_im[y][x], kpt[i]);
            }
            kdmin = 10;
            knm = 0;
            for (i = 0; i < (int)knum; i++)
            {
                if (kdist[i] < kdmin)
                {
                    kdmin = kdist[i];
                    knm = (unsigned)i;
                }
            }
            if (knum - knum2 * 2 == 0)
            {
                val = (BYTE) ((knm >= knum2) ? 255 : 0);
                t = (unsigned)kpt[knum2][3];
                st += t;
                sn++;
            }
            else
            {
                if (knm == knum2)
                {
                    val = (BYTE) ((kptb >= 0) ? 255 : 0);
                }
                else
                {
                    val = (BYTE) ((knm >= knum2) ? 255 : 0);
                }
            }
            d_im[y][x] = val;
        }
    }
    if (sn > 0)
    {
        threshold = (int)(st / sn + 0.5);
    }

    free(kdist);
    for (y = 0; y < (int)knum; y++)
    {
        free(kpt[y]);
    }
    free(kpt);
    for (y = 0; y < (int)knum; y++)
    {
        free(ksum[y]);
    }
    free(ksum);

    return threshold;
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterTMscaleLayer (IMTpixel** p_im, WORD** t_im, unsigned height, unsigned width, unsigned radius, float sensitivity, float doverlay, float delta)
{
    unsigned whcp, y, x, l, i, j, y0, x0, y1, x1, blsz, rsz;
    float immean, st, kover, sensdiv, senspos, sensinv;
    unsigned pim, immin, immax, imt, cnth, cntw, level = 0;
    unsigned maskbl, maskover;
    WORD tim;
    int threshold = 0;

    whcp = height;
    whcp += width;
    whcp /= 2;
    blsz = 1;
    while (blsz < whcp)
    {
        level++;
        blsz *= 2;
    }
    blsz /= 2;
    rsz = 1;
    while (rsz < radius)
    {
        if (level > 1)
        {
            level--;
        }
        rsz *= 2;
    }
    immin = (unsigned)p_im[0][0].s;
    immax = immin;
    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            pim = (unsigned)p_im[y][x].s;
            if (pim < immin)
            {
                immin = pim;
            }
            if (pim > immax)
            {
                immax = pim;
            }
        }
    }
    immean = (float)(immax + immin);
    immean /= 2.0;
    immean += 0.5;
    tim = Byte3Clamp((int)immean);
    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            t_im[y][x] = tim;
        }
    }
    if (doverlay < 0)
    {
        doverlay = 0;
    }
    kover = doverlay + 1.0;
    if (sensitivity < 0)
    {
        sensitivity = -sensitivity;
        sensdiv = sensitivity;
        sensdiv += 1.0;
        sensinv = 1.0 / sensdiv;
        senspos = sensitivity / sensdiv;
    }
    else
    {
        sensdiv = sensitivity;
        sensdiv += 1.0;
        senspos = 1.0 / sensdiv;
        sensinv = sensitivity / sensdiv;
    }
    for (l = 0; l < level; l++)
    {
        cnth = (height + blsz - 1) / blsz;
        cntw = (width + blsz - 1) / blsz;
        maskbl = blsz;
        maskover = (unsigned)(kover * maskbl);
        for (i = 0; i < cnth; i++)
        {
            y0 = i * maskbl;
            y1 = y0 + maskover;
            if (y1 > height)
            {
                y1 = height;
            }
            for (j = 0; j < cntw; j++)
            {
                x0 = j * maskbl;
                x1 = x0 + maskover;
                if (x1 > width)
                {
                    x1 = width;
                }

                immin = (unsigned)p_im[y0][x0].s;
                immax = immin;
                for (y = y0; y < y1; y++)
                {
                    for (x = x0; x < x1; x++)
                    {
                        pim = (unsigned)p_im[y][x].s;
                        if (pim < immin)
                        {
                            immin = pim;
                        }
                        if (pim > immax)
                        {
                            immax = pim;
                        }
                    }
                }
                immean = (float)(immax + immin);
                immean /= 2.0;
                immean *= sensinv;
                for (y = y0; y < y1; y++)
                {
                    for (x = x0; x < x1; x++)
                    {
                        imt = (unsigned)t_im[y][x];
                        imt *= senspos;
                        imt += immean;
                        imt += 0.5;
                        t_im[y][x] = Byte3Clamp((int)imt);
                    }
                }
            }
        }
        blsz /= 2;
    }
    st = 0;
    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            imt = (unsigned)t_im[y][x];
            imt += delta;
            t_im[y][x] = Byte3Clamp((int)imt);
            st += imt;
        }
    }
    st /= height;
    st /= width;
    st += 0.5;
    threshold = (int)st;

    return threshold;
}

///////////////////////////////////////////////////////////////////////////////

int IMTFilterTMscale (IMTpixel** p_im, BYTE** d_im, unsigned height, unsigned width, unsigned radius, float sensitivity, float doverlay, float delta)
{
    int threshold = 0;
    WORD** t_im = TLalloc(height, width);

    threshold = IMTFilterTMscaleLayer (p_im, t_im, height, width, radius, sensitivity, doverlay, delta);
    threshold = IMTFilterThresholdLayer (p_im, t_im, d_im, height, width);

    TLfree(t_im, height);

    return threshold;
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterTNiblackLayer (IMTpixel** p_im, WORD** t_im, unsigned height, unsigned width, int radius, float sensitivity, int lower_bound, int upper_bound, float delta)
{
    unsigned int y, x, n, y1, x1, y2, x2, yf, xf;
    float imx, imm, imv, t, st = 0, sn = 0;
    int threshold = 0;
    WORD val;

    radius = (radius < 0) ? -radius : radius;
    if (upper_bound < lower_bound)
    {
        upper_bound += lower_bound;
        lower_bound = upper_bound - lower_bound;
        upper_bound -= lower_bound;
    }
    lower_bound *= 3;
    upper_bound *= 3;

    for (y = 0; y < height; y++)
    {
        y1 = (y < radius) ? 0 : (y - radius);
        y2 = y + radius + 1;
        y2 = (y2 < height) ? y2 : height;
        for (x = 0; x < width; x++)
        {
            x1 = (x < radius) ? 0 : (x - radius);
            x2 = x + radius + 1;
            x2 = (x2 < width) ? x2 : width;
            imm = 0;
            imv = 0;
            n = 0;
            for (yf = y1; yf < y2; yf++)
            {
                for (xf = x1; xf < x2; xf++)
                {
                    imx = (float)p_im[yf][xf].s;
                    imm += imx;
                    imv += (imx * imx);
                    n++;
                }
            }
            if (n == 0)
            {
                n = 1;
            }
            imm /= n;
            imv /= n;
            imv -= (imm * imm);
            if (imv < 0)
            {
                imv = -imv;
            }
            imv = sqrt(imv);
            imx = (float)p_im[y][x].s;
            if (imx < (float)lower_bound)
            {
                t = (float)lower_bound;
            }
            else if (imx > (float)upper_bound)
            {
                t = (float)upper_bound;
            }
            else
            {
                t = (imm + sensitivity * imv + delta);
            }
            val = Byte3Clamp((int)(t + 0.5));
            st += t;
            sn++;
            t_im[y][x] = val;
        }
    }
    sn = (sn < 0.0 || sn > 0.0) ? sn : 1.0;
    threshold =  (int)(st / sn + 0.5);

    return threshold;
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterTNiblack (IMTpixel** p_im, BYTE** d_im, unsigned height, unsigned width, int radius, float sensitivity, int lower_bound, int upper_bound, float delta)
{
    int threshold = 0;
    WORD** t_im = TLalloc(height, width);

    threshold = IMTFilterTNiblackLayer (p_im, t_im, height, width, radius, sensitivity, lower_bound, upper_bound, delta);
    threshold = IMTFilterThresholdLayer (p_im, t_im, d_im, height, width);

    TLfree(t_im, height);

    return threshold;
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterTOtsuValue (IMTpixel** p_im, unsigned height, unsigned width)
{
    unsigned i, Tmax = 768;
    float w = 0.0f, u = 0.0f, uT = 0.0f, hn;
    float work1, work2, work3 = 0.0f;
    int threshold = 0;
    unsigned long long histogram[768] = {0};

    IMTHist (p_im, histogram, 0, 0, height, width, Tmax);
    hn = (float)(width * height);
    for (i = 0; i < Tmax; i++)
    {
        uT += (float)(histogram[i] * (i + 1));
    }
    uT /= hn;

    for (i = 0; i < Tmax; i++)
    {
        w += (float)histogram[i];
        u += (float)((i + 1) * histogram[i]);
        work1 = (uT * w - u);
        work2 = (work1 * work1) / (w * (hn - w));

        if (work2 > work3)
        {
            work3 = work2;
            threshold = (int)i;
        }
    }

    threshold -= 1;

    if(threshold <= 0)
    {
        threshold = (int)(uT + 0.5);
    }

    return threshold;
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterTOtsu (IMTpixel** p_im, BYTE** d_im, unsigned height, unsigned width)
{
    int threshold = 0;

    threshold = IMTFilterTOtsuValue (p_im, height, width);
    threshold = IMTFilterThreshold (p_im, d_im, height, width, threshold);

    return threshold;
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterTQuadModValue (IMTpixel** p_im, unsigned height, unsigned width)
{
    int T, TG, Trw, Trb;
    int threshold = 0;

    TG = IMTFilterTBiModValue (p_im, height, width);
    Trw = IMTFilterTBiModValueBound (p_im, height, width, TG, 765);
    Trb = IMTFilterTBiModValueBound (p_im, height, width, 0, TG);
    T = IMTFilterTBiModValueBound (p_im, height, width, Trb, Trw);
    T += TG;
    T++;
    T /= 2;
    threshold = T;

    return threshold;
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterTQuadMod (IMTpixel** p_im, BYTE** d_im, unsigned height, unsigned width, int delta)
{
    int threshold = 0;

    threshold = IMTFilterTQuadModValue (p_im, height, width);
    threshold += delta;
    threshold = IMTFilterThreshold (p_im, d_im, height, width, threshold);

    return threshold;
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterTRotValue (IMTpixel** p_im, unsigned height, unsigned width, bool weight)
{
    unsigned Tmax = 768;
    unsigned i, il, ir;
    float wl, wr, kw, wi;
    int threshold = 0;
    unsigned long long histogram[768] = {0};

    IMTHist (p_im, histogram, 0, 0, height, width, Tmax);
    if (weight)
    {
        for (i = 0; i < Tmax; i++)
        {
            kw = (float)Tmax;
            wi = 2 *( kw - (float)i) / kw;
            histogram[i] *= wi;
        }
    }
    wl = wr = 0.0f;
    il = 0;
    ir = Tmax - 1;
    while (il < (ir - 1))
    {
        if (wl < wr)
        {
            il++;
            for (i = 0; i <= il; i++)
            {
                wl += (float)histogram[i];
            }
        }
        else
        {
            ir--;
            for (i = ir; i < Tmax; i++)
            {
                wr += (float)histogram[i];
            }
        }
    }

    threshold = il;

    return threshold;
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterTRot (IMTpixel** p_im, BYTE** d_im, unsigned height, unsigned width, bool weight)
{
    int threshold = 0;

    threshold = IMTFilterTRotValue (p_im, height, width, weight);
    threshold = IMTFilterThreshold (p_im, d_im, height, width, threshold);

    return threshold;
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterTSauvolaLayer (IMTpixel** p_im, WORD** t_im, unsigned height, unsigned width, int radius, float sensitivity, int dynamic_range, int lower_bound, int upper_bound, float delta)
{
    unsigned int x, y, n, y1, x1, y2, x2, yf, xf;
    float imx, imm, imv, ima, t, st = 0, sn = 0;
    int threshold = 0;
    WORD val;

    radius = (radius < 0) ? -radius : radius;
    if (upper_bound < lower_bound)
    {
        upper_bound += lower_bound;
        lower_bound = upper_bound - lower_bound;
        upper_bound -= lower_bound;
    }
    lower_bound *= 3;
    upper_bound *= 3;

    for (y = 0; y < height; y++)
    {
        y1 = (y < radius) ? 0 : (y - radius);
        y2 = y + radius + 1;
        y2 = (y2 < height) ? y2 : height;
        for (x = 0; x < width; x++)
        {
            x1 = (x < radius) ? 0 : (x - radius);
            x2 = x + radius + 1;
            x2 = (x2 < width) ? x2 : width;
            imm = 0;
            imv = 0;
            n = 0;
            for (yf = y1; yf < y2; yf++)
            {
                for (xf = x1; xf < x2; xf++)
                {
                    imx = (float)p_im[yf][xf].s;
                    imm += imx;
                    imv += (imx * imx);
                    n++;
                }
            }
            if (n == 0)
            {
                n = 1;
            }
            imm /= n;
            imv /= n;
            imv -= (imm * imm);
            imv = (imv < 0) ? -imv : imv;
            imv = sqrt(imv);
            ima = 1.0 - imv / (float)dynamic_range;
            imx = (float)p_im[y][x].s;
            if (imx < (float)lower_bound)
            {
                t = (float)lower_bound;
            }
            else if (imx > (float)upper_bound)
            {
                t = (float)upper_bound;
            }
            else
            {
                t = imm * (1.0 - sensitivity * ima) + delta;
            }
            val = Byte3Clamp((int)(t + 0.5));
            st += t;
            sn++;
            t_im[y][x] = val;
        }
    }
    sn = (sn < 0.0 || sn > 0.0) ? sn : 1.0;
    threshold = (int)(st / sn + 0.5);

    return threshold;
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterTSauvola (IMTpixel** p_im, BYTE** d_im, unsigned height, unsigned width, int radius, float sensitivity, int dynamic_range, int lower_bound, int upper_bound, float delta)
{
    int threshold = 0;
    WORD** t_im = TLalloc(height, width);

    threshold = IMTFilterTSauvolaLayer (p_im, t_im, height, width, radius, sensitivity, dynamic_range, lower_bound, upper_bound, delta);
    threshold = IMTFilterThresholdLayer (p_im, t_im, d_im, height, width);

    TLfree(t_im, height);

    return threshold;
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterTSizeLayer (IMTpixel** p_im, WORD** t_im, unsigned height, unsigned width, int radius, int lower_bound, int upper_bound, int delta)
{
    unsigned x, y, sn = 0;
    float st = 0;
    int im, t, threshold = 0;
    WORD val;

    radius = (radius < 0) ? -radius : radius;

    IMTpixel** s_im = IMTalloc(radius, radius);
    IMTpixel** b_im = IMTalloc(height, width);

    if (upper_bound < lower_bound)
    {
        upper_bound += lower_bound;
        lower_bound = upper_bound - lower_bound;
        upper_bound -= lower_bound;
    }
    lower_bound *= 3;
    upper_bound *= 3;

    IMTFilterSGsample(p_im, s_im, height, width, radius, radius);
    IMTFilterSBicub(s_im, b_im, radius, radius, height, width);
    IMTfree(s_im, radius);

    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            im = (int)p_im[y][x].s;
            if (im < lower_bound)
            {
                t = lower_bound;
            }
            else if (im > upper_bound)
            {
                t = upper_bound;
            }
            else
            {
                t = (int)b_im[y][x].s + delta;
            }
            val = Byte3Clamp((int)t);
            st += t;
            sn++;
            t_im[y][x] = val;
        }
    }
    sn = (sn < 0.0 || sn > 0.0) ? sn : 1.0;
    threshold = (int)(st / sn + 0.5);
    IMTfree(b_im, height);

    return threshold;
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterTSize (IMTpixel** p_im, BYTE** d_im, unsigned height, unsigned width, int radius, int lower_bound, int upper_bound, int delta)
{
    int threshold = 0;
    WORD** t_im = TLalloc(height, width);

    threshold = IMTFilterTSizeLayer (p_im, t_im, height, width, radius, lower_bound, upper_bound, delta);
    threshold = IMTFilterThresholdLayer (p_im, t_im, d_im, height, width);

    TLfree(t_im, height);

    return threshold;
}
////////////////////////////////////////////////////////////////////////////////

int IMTFilterTText (IMTpixel** p_im, BYTE** d_im, unsigned height, unsigned width, unsigned contour, unsigned radius)
{
    unsigned wr, ws, s, sn, n, z;
    int y, x, y0, x0, yk, xk, yf, xf, ys, xs, im, imf, ims, imd;
    BYTE val;
    int h = (int)height;
    int w = (int)width;
    int threshold = 0;

    WORD** b_im = TLalloc(height, width);

    radius = (radius < 0) ? -radius : radius;
    wr = 2 * radius;
    ws = wr * wr - 1;

    z = 0;
    for (y = 0; y < h; y++)
    {
        y0 = y - radius;
        yk = y + radius;
        if (y0 < 0)
        {
            y0 = 0;
        }
        if (yk > h)
        {
            yk = h;
        }
        for (x = 0; x < w; x++)
        {
            x0 = x - radius;
            xk = x + radius;
            if (x0 < 0)
            {
                x0 = 0;
            }
            if (xk > w)
            {
                xk = w;
            }
            n = 0;
            im = (int)p_im[y][x].s;
            for (yf = y0; yf < yk; yf++)
            {
                for (xf = x0; xf < xk; xf++)
                {
                    imf = (int)p_im[yf][xf].s;
                    sn = 0;
                    s = 0;
                    for (ys = yf - 1; ys < yf + 2; ys++)
                    {
                        if (ys >= 0 && ys < h)
                        {
                            for (xs = xf - 1; xs < xf + 2; xs++)
                            {
                                if (xs >= 0 && xs < w)
                                {
                                    if (xs != xf && ys != yf)
                                    {
                                        ims = (int)p_im[ys][xs].s;
                                        imd = imf - ims;
                                        if (imd < 0)
                                        {
                                            imd = -imd;
                                        }
                                        if (imd > (int)contour)
                                        {
                                            s++;
                                        }
                                        sn++;
                                    }
                                }
                            }
                        }
                    }
                    s *= 2;
                    if (sn > 0)
                    {
                        s /= sn;
                    }
                    n += s;
                }
            }
            if (z < n)
            {
                z = n;
            }
            b_im[y][x] = (WORD)n;
        }
    }
    z /= 2;
    threshold = (int)(384 * z / ws);
    for (y = 0; y < h; y++ )
    {
        for (x = 0; x < w; x++)
        {
            im = (int)b_im[y][x];
            val = (BYTE) ((im >= (int)z) ? 255 : 0);
            d_im[y][x] = val;
        }
    }

    TLfree(b_im, height);

    return threshold;
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterTTsaiValue (IMTpixel** p_im, unsigned height, unsigned width, int shift)
{

    unsigned cn = 768, i;
    int threshold = 0;
    unsigned long long histogram[768] = {0};
    float criterion = 0.0;
    float m1, m2, m3;
    float cd, c0, c1, z0, z1, pd, p0, p1;

    IMTHist (p_im, histogram, 0, 0, height, width, 768);

    m1 = m2 = m3 = 0.0;

    for (i = 0; i < cn; i++)
    {
        m1 += i * (float)histogram[i];
        m2 += i * i * (float)histogram[i];
        m3 += i * i * i * (float)histogram[i];
    }

    cd = m2 - m1 * m1;
    c0 = (-m2 * m2 + m1 * m3) / cd;
    c1 = (-m3 + m2 * m1) / cd;

    z0 = 0.5 * (-c1 - sqrt(c1 * c1 - 4.0 * c0));
    z1 = 0.5 * (-c1 + sqrt(c1 * c1 - 4.0 * c0));

    pd = z1 - z0;
    p0 = (z1 - m1) / pd;
    p1 = 1.0 - p0;

    for (threshold = 0; threshold < (int)cn; threshold++)
    {
        criterion += (float)histogram[threshold];
        if (criterion > p1) break;
    }

    if(threshold == 255)
    {
        threshold = 0;
    }

    threshold += shift * 3;

    return threshold;
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterTTsai (IMTpixel** p_im, BYTE** d_im, unsigned height, unsigned width, int shift)
{
    int threshold = 0;

    threshold = IMTFilterTTsaiValue (p_im, height, width, shift);
    threshold = IMTFilterThreshold (p_im, d_im, height, width, threshold);

    return threshold;
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterTWhiteRohrer(IMTpixel** p_im, BYTE** d_im, unsigned height, unsigned width, int x_lookahead, int y_lookahead, int bias_mode, int bias_factor, int f_factor, int g_factor)
{
    int x, y;
    int h = (int)height;
    int w = (int)width;
    BYTE val;

    int WR1_F_OFFSET = 255;
    int WR1_G_OFFSET = 255;
    int BIN_FOREGROUND = 0;
    int BIN_BACKGROUND = 255;
    int WR1_BIAS_CROSSOVER = 93;
    int WR1_BIAS = 20;
    float WR1_BLACK_BIAS_FACTOR = 0.0;
    float WR1_WHITE_BIAS_FACTOR = -0.25;
    int wr1_f_tab[512] =
    {
        -62,  -62,  -61,  -61,  -60,  -60,  -59,  -59,
            -58,  -58,  -57,  -57,  -56,  -56,  -54,  -54,
            -53,  -53,  -52,  -52,  -51,  -51,  -50,  -50,
            -49,  -49,  -48,  -48,  -47,  -47,  -46,  -46,
            -45,  -45,  -44,  -44,  -43,  -43,  -42,  -42,
            -41,  -41,  -41,  -41,  -40,  -40,  -39,  -39,
            -38,  -38,  -37,  -37,  -36,  -36,  -36,  -36,
            -35,  -35,  -34,  -34,  -33,  -33,  -33,  -33,
            -32,  -32,  -31,  -31,  -31,  -31,  -30,  -30,
            -29,  -29,  -29,  -29,  -28,  -28,  -27,  -27,
            -27,  -27,  -26,  -26,  -25,  -25,  -25,  -25,
            -24,  -24,  -24,  -24,  -23,  -23,  -23,  -23,
            -22,  -22,  -22,  -22,  -21,  -21,  -21,  -21,
            -20,  -20,  -20,  -20,  -19,  -19,  -19,  -19,
            -18,  -18,  -18,  -18,  -17,  -17,  -17,  -17,
            -16,  -16,  -16,  -16,  -16,  -16,  -15,  -15,
            -15,  -15,  -14,  -14,  -14,  -14,  -14,  -14,
            -13,  -13,  -13,  -13,  -13,  -13,  -12,  -12,
            -12,  -12,  -12,  -12,  -11,  -11,  -11,  -11,
            -11,  -11,  -10,  -10,  -10,  -10,  -10,  -10,
            -9,   -9,   -9,   -9,   -9,   -9,   -8,   -8,
            -8,   -8,   -8,   -8,   -8,   -8,   -7,   -7,
            -7,   -7,   -7,   -7,   -7,   -7,   -6,   -6,
            -6,   -6,   -6,   -6,   -6,   -6,   -5,   -5,
            -5,   -5,   -5,   -5,   -5,   -5,   -4,   -4,
            -3,   -3,   -2,   -2,   -2,   -2,   -2,   -2,
            -2,   -2,   -2,   -2,   -1,   -1,   -1,   -1,
            -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
            -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
            -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
            -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
            -1,   -1,   -1,   -1,   -1,   -1,    0,    0,
            1,    1,    2,    2,    2,    2,    2,    2,
            2,    2,    2,    2,    2,    2,    2,    2,
            2,    2,    2,    2,    2,    2,    2,    2,
            2,    2,    2,    2,    2,    2,    2,    2,
            2,    2,    2,    2,    2,    2,    2,    2,
            2,    2,    2,    2,    2,    2,    2,    2,
            2,    2,    2,    2,    2,    2,    2,    2,
            3,    3,    3,    3,    3,    3,    3,    3,
            3,    3,    3,    3,    3,    3,    3,    3,
            3,    3,    3,    3,    3,    3,    4,    4,
            4,    4,    4,    4,    4,    4,    4,    4,
            4,    4,    4,    4,    4,    4,    4,    4,
            4,    4,    4,    4,    4,    4,    4,    4,
            4,    4,    5,    5,    5,    5,    5,    5,
            5,    5,    5,    5,    5,    5,    5,    5,
            5,    5,    6,    6,    6,    6,    6,    6,
            6,    6,    6,    6,    6,    6,    6,    6,
            6,    6,    6,    6,    6,    6,    6,    6,
            6,    6,    6,    6,    6,    6,    6,    6,
            6,    6,    7,    7,    7,    7,    7,    7,
            7,    7,    7,    7,    7,    7,    7,    7,
            7,    7,    7,    7,    7,    7,    7,    7,
            7,    7,    7,    7,    8,    8,    8,    8,
            8,    8,    8,    8,    8,    8,    8,    8,
            8,    8,    8,    8,    8,    8,    8,    8,
            8,    8,    9,    9,    9,    9,    9,    9,
            9,    9,    9,    9,    9,    9,    9,    9,
            9,    9,    9,    9,    9,    9,    9,    9,
            9,    9,    9,    9,    9,    9,    9,    9,
            9,    9,   10,   10,   10,   10,   10,   10,
            10,   10,   10,   10,   10,   10,   10,   10,
            10,   10,   10,   10,   10,   10,   10,    0
        };
    int wr1_g_tab[512] =
    {
        -126, -126, -125, -125, -124, -124, -123, -123,
            -122, -122, -121, -121, -120, -120, -119, -119,
            -118, -118, -117, -117, -116, -116, -115, -115,
            -114, -114, -113, -113, -112, -112, -111, -111,
            -110, -110, -109, -109, -108, -108, -107, -107,
            -106, -106, -105, -105, -104, -104, -103, -103,
            -102, -102, -101, -101, -100, -100,  -99,  -99,
            -98,  -98,  -97,  -97,  -96,  -96,  -95,  -95,
            -94,  -94,  -93,  -93,  -92,  -92,  -91,  -91,
            -90,  -90,  -89,  -89,  -88,  -88,  -87,  -87,
            -86,  -86,  -85,  -85,  -84,  -84,  -83,  -83,
            -82,  -82,  -81,  -81,  -80,  -80,  -79,  -79,
            -78,  -78,  -77,  -77,  -76,  -76,  -75,  -75,
            -74,  -74,  -73,  -73,  -72,  -72,  -71,  -71,
            -70,  -70,  -69,  -69,  -68,  -68,  -67,  -67,
            -66,  -66,  -65,  -65,  -64,  -64,  -63,  -63,
            -61,  -61,  -59,  -59,  -57,  -57,  -54,  -54,
            -52,  -52,  -50,  -50,  -48,  -48,  -46,  -46,
            -44,  -44,  -42,  -42,  -41,  -41,  -39,  -39,
            -37,  -37,  -36,  -36,  -34,  -34,  -33,  -33,
            -31,  -31,  -30,  -30,  -29,  -29,  -27,  -27,
            -26,  -26,  -25,  -25,  -24,  -24,  -23,  -23,
            -22,  -22,  -21,  -21,  -20,  -20,  -19,  -19,
            -18,  -18,  -17,  -17,  -16,  -16,  -15,  -15,
            -14,  -14,  -14,  -14,  -13,  -13,  -12,  -12,
            -12,  -12,  -11,  -11,  -10,  -10,  -10,  -10,
            -9,   -9,   -8,   -8,   -8,   -8,   -7,   -7,
            -7,   -7,   -6,   -6,   -6,   -6,   -5,   -5,
            -5,   -5,   -4,   -4,   -2,   -2,   -2,   -2,
            -2,   -2,   -1,   -1,   -1,   -1,   -1,   -1,
            -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
            -1,   -1,   -1,   -1,   -1,   -1,    0,    0,
            1,    1,    1,    1,    1,    1,    1,    1,
            1,    1,    1,    1,    1,    1,    1,    1,
            1,    1,    2,    2,    2,    2,    2,    2,
            4,    4,    5,    5,    5,    5,    6,    6,
            6,    6,    7,    7,    7,    7,    8,    8,
            8,    8,    9,    9,   10,   10,   10,   10,
            11,   11,   12,   12,   12,   12,   13,   13,
            14,   14,   14,   14,   15,   15,   16,   16,
            17,   17,   18,   18,   19,   19,   20,   20,
            21,   21,   22,   22,   23,   23,   24,   24,
            25,   25,   26,   26,   27,   27,   29,   29,
            30,   30,   31,   31,   33,   33,   34,   34,
            36,   36,   37,   37,   39,   39,   41,   41,
            42,   42,   44,   44,   46,   46,   48,   48,
            50,   50,   52,   52,   54,   54,   57,   57,
            59,   59,   61,   61,   63,   63,   64,   64,
            65,   65,   66,   66,   67,   67,   68,   68,
            69,   69,   70,   70,   71,   71,   72,   72,
            73,   73,   74,   74,   75,   75,   76,   76,
            77,   77,   78,   78,   79,   79,   80,   80,
            81,   81,   82,   82,   83,   83,   84,   84,
            85,   85,   86,   86,   87,   87,   88,   88,
            89,   89,   90,   90,   91,   91,   92,   92,
            93,   93,   94,   94,   95,   95,   96,   96,
            97,   97,   98,   98,   99,   99,  100,  100,
            101,  101,  102,  102,  103,  103,  104,  104,
            105,  105,  106,  106,  107,  107,  108,  108,
            109,  109,  110,  110,  111,  111,  112,  112,
            113,  113,  114,  114,  115,  115,  116,  116,
            117,  117,  118,  118,  119,  119,  120,  120,
            121,  121,  122,  122,  123,  123,  124,  124,
            125,  125,  126,  126,  127,  127,  127,    0
        };

    int u, prevY, Y = 0;
    int f, g;
    int x_ahead, y_ahead;
    int t, tx, rx, tbias;
    int offset = WR1_BIAS;
    unsigned im, threshold = 0, sum, sqr_sum;
    float mu = 0.0, s_dev = 0.0;
    int nZ = (int)(2 * width + 1), n;
    int Z[nZ];
    int st = 0, sn = 0;

    x_lookahead = x_lookahead % width;

    if (bias_mode == 0)
    {
        sum = 0;
        for (y = 0; y < h; y++)
        {
            for (x = 0; x < w; x++)
            {
                im = (unsigned)p_im[y][x].s;
                im += 2;
                im /= 3;
                sum += im;
            }
        }
        mu = (float)sum / (width*height);

        sqr_sum = 0;
        for (y = 0; y < h; y++)
        {
            for (x = 0; x < w; x++)
            {
                im = (unsigned)p_im[y][x].s;
                im += 2;
                im /= 3;
                sqr_sum += (im * im);
            }
        }
        s_dev = sqrt((float)sqr_sum / (width*height) - mu * mu);
        offset = (int)(s_dev - 40);
    }
    else
    {
        offset = bias_mode;
    }

    for(n = 0; n < nZ; ++n)
    {
        Z[n] = 0;
    }
    Z[0] = prevY = (int)mu;

    for (y = 0; y< 1 + y_lookahead; y++)
    {
        if (y < y_lookahead)
        {
            t = (int)width;
        }
        else
        {
            t = x_lookahead;
        }
        for (x=0; x < t; x++)
        {
            im = (unsigned)p_im[y][x].s;
            im += 2;
            im /= 3;
            u = (int)im;
            f = -wr1_f_tab[WR1_F_OFFSET - (u - prevY)];
            Y = prevY + f;
            if (y == 1)
            {
                Z[x] = (int)mu;
            }
            else
            {
                g = -wr1_g_tab[WR1_G_OFFSET - (Y - Z[x])];
                Z[x] = Z[x] + g;
            }
        }

    }
    x_ahead = 1 + x_lookahead;
    y_ahead = 1 + y_lookahead;

    for (y = 0; y < h; y++)
    {
        for (x = 0; x < w; x++)
        {
            threshold = (unsigned)((bias_factor < 0) ? -bias_factor : bias_factor);
            tx = 256 - Z[x_ahead];
            tbias = -offset;
            if (tx < WR1_BIAS_CROSSOVER)
            {
                rx = tx - tbias - (int)(WR1_BLACK_BIAS_FACTOR * (WR1_BIAS_CROSSOVER - tx));
            }
            else if (tx >= WR1_BIAS_CROSSOVER)
            {
                rx = tx + tbias + (int)(WR1_WHITE_BIAS_FACTOR * (tx - WR1_BIAS_CROSSOVER));
            }
            else
            {
                rx = tx;
            }
            if (rx < BIN_FOREGROUND)
            {
                rx = BIN_FOREGROUND;
            }
            if (rx > BIN_BACKGROUND)
            {
                rx = BIN_BACKGROUND;
            }
            rx = 256 - rx;
            threshold *= rx;
            threshold /= 100;
            st += threshold;
            sn++;
            im = (unsigned)p_im[y][x].s;
            im += 2;
            im /= 3;
            val = (BYTE) ((im < threshold) ? 0 : 255);
            d_im[y][x] = val;
            x_ahead++;
            if (x_ahead >= w)
            {
                x_ahead = 1;
                y_ahead++;
            }
            if (y_ahead < h)
            {
                prevY = Y;
                im = (unsigned)p_im[y_ahead][x_ahead].s;
                im += 2;
                im /= 3;
                f = -wr1_f_tab[WR1_F_OFFSET - (im - prevY)];
                Y = prevY + f_factor * f / 100;
                g = -wr1_g_tab[WR1_G_OFFSET - (Y - Z[x_ahead])];
                Z[x_ahead] = Z[x_ahead] + g_factor * g / 100;
            }
            else
            {
                Z[x_ahead] = Z[x_ahead - 1];
            }
        }
    }
    if (sn > 0)
    {
        threshold = (unsigned)Byte3Clamp((int)(3 * st / sn));
    }

    return (int)threshold;
}

////////////////////////////////////////////////////////////////////////////////
