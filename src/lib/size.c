//  Zlib license
//
// ImThreshold library.
//
//  Copyright (C) 2017:
//  zvezdochiot <zvezdochiot@user.sourceforge.net>
//  Homepage: https://sourceforge.net/projects/imthreshold/

#include "../imthreshold.h"

////////////////////////////////////////////////////////////////////////////////

IMTpixel IMTInterpolateBiCubic (IMTpixel** p_im, int height, int width, float y, float x)
{
    int i, d, xi, yi, xf, yf;
    float d0, d2, d3, a0, a1, a2, a3;
    float dx, dy;
    float Cc, C[4];
    IMTpixel dim;

    yi = (int)y;
    yi = (yi < 0) ? 0 : (yi < height) ? yi : (height - 1);
    xi = (int)x;
    xi = (xi < 0) ? 0 : (xi < width) ? xi : (width - 1);
    dy = y - yi;
    dx = x - xi;
    for(d = 0; d < 3; d++)
    {
        for(i = -1; i < 3; i++)
        {
            yf = (int)y + i;
            yf = (yf < 0) ? 0 : (yf < height) ? yf : (height - 1);
            xf = (int)x;
            xf = (xf < 0) ? 0 : (xf < width) ? xf : (width - 1);
            a0 = p_im[yf][xf].c[d];
            xf = (int)x - 1;
            xf = (xf < 0) ? 0 : (xf < width) ? xf : (width - 1);
            d0 = p_im[yf][xf].c[d];
            d0 -= a0;
            xf = (int)x + 1;
            xf = (xf < 0) ? 0 : (xf < width) ? xf : (width - 1);
            d2 = p_im[yf][xf].c[d];
            d2 -= a0;
            xf = (int)x + 2;
            xf = (xf < 0) ? 0 : (xf < width) ? xf : (width - 1);
            d3 = p_im[yf][xf].c[d];
            d3 -= a0;
            a1 =  -1.0f / 3.0f * d0 + d2 - 1.0f / 6.0f * d3;
            a2 = 1.0f / 2.0f * d0 + 1.0f / 2.0f * d2;
            a3 = -1.0f / 6.0 * d0 - 1.0f / 2.0f * d2 + 1.0f / 6.0f * d3;
            C[i + 1] = a0 + (a1 + (a2 + a3 * dx) * dx) * dx;
        }
        d0 = C[0] - C[1];
        d2 = C[2] - C[1];
        d3 = C[3] - C[1];
        a0 = C[1];
        a1 =  -1.0f / 3.0f * d0 + d2 - 1.0f / 6.0f * d3;
        a2 = 1.0f / 2.0f * d0 + 1.0f / 2.0f * d2;
        a3 = -1.0f / 6.0f * d0 - 1.0f / 2.0f * d2 + 1.0f / 6.0f * d3;
        Cc = a0 + (a1 + (a2 + a3 * dy) * dy) * dy;
        dim.c[d] = ByteClamp((int)(Cc + 0.5f));
    }
    dim = IMTcalcS (dim);

    return dim;
}

////////////////////////////////////////////////////////////////////////////////

IMTpixel IMTInterpolateBiLine (IMTpixel** p_im, int height, int width, float y, float x)
{
    int d, xi, yi, xf, yf;
    float dx1, dy1, dx2, dy2;
    IMTpixel dim;

    yi = (int)y;
    yi = (yi < 0) ? 0 : (yi < height) ? yi : (height - 1);
    xi = (int)x;
    xi = (xi < 0) ? 0 : (xi < width) ? xi : (width - 1);
    dy1 = y - yi;
    dx1 = x - xi;
    dy2 = 1.0 - dy1;
    dx2 = 1.0 - dx1;
    yf = (int)y + 1;
    yf = (yf < 0) ? 0 : (yf < height) ? yf : (height - 1);
    xf = (int)x + 1;
    xf = (xf < 0) ? 0 : (xf < width) ? xf : (width - 1);
    for (d = 0; d < 3; d++)
    {
        dim.c[d] = ByteClamp((int) (dy2 * (dx2 * p_im[yi][xi].c[d] + dx1 * p_im[yi][xf].c[d]) + dy1 * (dx2 * p_im[yf][xi].c[d] + dx1 * p_im[yf][xf].c[d]) + 0.5));
    }
    dim = IMTcalcS (dim);

    return dim;
}

////////////////////////////////////////////////////////////////////////////////

float IMTInterpolateAkima (float* y, float dx)
{
    int i;
    float a, b, val;
    float m[5], t[2];

    for(i = 0; i < 5; i++)
    {
        m[i] = y[i + 1] - y[i];
    }
    for(i = 0; i < 2; i++)
    {
        a = (m[i + 2] > m[i + 3]) ? (m[i + 2] - m[i + 3]) : (m[i + 3] - m[i + 2]);
        b = (m[i] > m[i + 1]) ? (m[i] - m[i + 1]) : (m[i + 1] - m[i]);
        t[i] = ((a + b) > 0) ? ((a * m[i + 1] + b * m[i + 2]) / (a + b)) : (0.5f * (m[i + 1] + m[i + 2]));
    }
    val = y[2] + (t[0] + ((m[2] + m[2] + m[2] - t[0] - t[0] - t[1]) + (t[0] + t[1] - m[2] - m[2]) * dx) * dx) * dx;

    return val;
}

////////////////////////////////////////////////////////////////////////////////

IMTpixel IMTInterpolateBiAkima (IMTpixel** p_im, int height, int width, float y, float x)
{
    int i, j, d, xi, yi, xf, yf, xfa, yfa;
    float z, zz[6], zzz[6];
    IMTpixel dim;

    yi = (int)y;
    yi = (yi < 0) ? 0 : (yi < height) ? yi : (height - 1);
    xi = (int)x;
    xi = (xi < 0) ? 0 : (xi < width) ? xi : (width - 1);
    for(d = 0; d < 3; d++)
    {
        for(i = -2; i < 4; i++)
        {
            yf = (int)y + i;
            yfa = (yf < 0) ? 0 : (yf < height) ? yf : (height - 1);
            for(j = -2; j < 4; j++)
            {
                xf = (int)x + j;
                xfa = (xf < 0) ? 0 : (xf < width) ? xf : (width - 1);
                zz[j + 2] = p_im[yfa][xfa].c[d];
            }
            zzz[i + 2] = IMTInterpolateAkima(zz, (x - xi));
        }
        z = IMTInterpolateAkima(zzz, (y - yi));
        dim.c[d] = ByteClamp((int)(z + 0.5f));
    }
    dim = IMTcalcS (dim);

    return dim;
}

////////////////////////////////////////////////////////////////////////////////

void IMTFilterSBicub (IMTpixel** p_im, IMTpixel** d_im, unsigned height, unsigned width, unsigned new_height, unsigned new_width)
{
    unsigned y, x;
    float xFactor = (float)width / new_width;
    float yFactor = (float)height / new_height;
    float ox, oy;

    for (y = 0; y < new_height; y++ )
    {
        oy  = ((float)y + 0.5f) * yFactor - 0.5f;
        for (x = 0; x < new_width; x++ )
        {
            ox  = ((float)x  + 0.5f) * xFactor - 0.5f;
            d_im[y][x] = IMTInterpolateBiCubic (p_im, height, width, oy, ox);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void IMTFilterSBilin (IMTpixel** p_im, IMTpixel** d_im, unsigned height, unsigned width, unsigned new_height, unsigned new_width)
{
    unsigned y, x;
    float xFactor = (float)width / new_width;
    float yFactor = (float)height / new_height;
    float ox, oy;

    for (y = 0; y < new_height; y++)
    {
        oy  = ((float)y + 0.5f) * yFactor - 0.5f;
        for (x = 0; x < new_width; x++ )
        {
            ox  = ((float)x  + 0.5f) * xFactor - 0.5f;
            d_im[y][x] = IMTInterpolateBiLine (p_im, height, width, oy, ox);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void IMTFilterSBicont (IMTpixel** p_im, IMTpixel** d_im, unsigned height, unsigned width, unsigned new_height, unsigned new_width)
{

    unsigned d;
    float xFactor = (float)width / new_width;
    float yFactor = (float)height / new_height;
    float ox, oy, dx, dy;
    int y, x, i, j, k, ox1, oy1, ox2, oy2;
    int h = new_height, w = new_width;
    int ymax = height - 1, xmax = width - 1;
    IMTpixel mim, **pim;
    float sdist, pdist[16];

    pim = IMTalloc (4, 4);
    for (y = 0; y < h; y++ )
    {
        oy  = ((float)y + 0.5f) * yFactor - 0.5f;
        oy1 = (int)oy;
        dy  = oy - (float)oy1;
        for (x = 0; x < w; x++ )
        {
            ox  = ((float)x  + 0.5f) * xFactor - 0.5f;
            ox1 = (int)ox;
            dx  = ox - (float)ox1;
            mim = IMTInterpolateBiCubic (p_im, height, width, oy, ox);
            sdist = 0;
            k = 0;
            for (i = 0; i < 4; i++)
            {
                oy2 = oy1 + i - 1;
                if ( oy2 < 0 )
                {
                    oy2 = 0;
                }
                if ( oy2 > ymax )
                {
                    oy2 = ymax;
                }
                for (j = 0; j < 4; j++)
                {
                    ox2 = ox1 + j - 1;
                    if ( ox2 < 0.0f )
                    {
                        ox2 = 0.0f;
                    }
                    if ( ox2 > xmax )
                    {
                        ox2 = xmax;
                    }
                    pim[i][j] = p_im[oy2][ox2];
                    pdist[k] = IMTdist(mim, pim[i][j]);
                    sdist += pdist[k];
                    k++;
                }
            }
            sdist /= 16.0f;
            if (sdist == 0.0f)
            {
                for (k = 0; k < 16; k++)
                {
                    pdist[k] = 1.0f;
                }
            }
            else
            {
                for (k = 0; k < 16; k++)
                {
                    pdist[k] += sdist;
                    pdist[k] = sdist / pdist[k];
                }
                sdist = 0;
                for (k = 0; k < 16; k++)
                {
                    sdist += pdist[k];
                }
                sdist /= 16.0f;
                if (sdist > 0.0f)
                {
                    for (k = 0; k < 16; k++)
                    {
                        pdist[k] /= sdist;
                    }
                }
            }
            k = 0;
            for (i = 0; i < 4; i++)
            {
                for (j = 0; j < 4; j++)
                {
                    for (d = 0; d < 3; d++)
                    {
                        pim[i][j].c[d] = ByteClamp((int) ((float)pim[i][j].c[d] * pdist[k] + (float)mim.c[d] * (1.0f - pdist[k]) + 0.5f));
                    }
                    k++;
                }
            }
            d_im[y][x] = IMTInterpolateBiCubic (pim, 4, 4, (1.0f + dy), (1.0f + dx));
        }
    }
    IMTfree (pim, 4);
}

////////////////////////////////////////////////////////////////////////////////

void IMTFilterSBiakima (IMTpixel** p_im, IMTpixel** d_im, unsigned height, unsigned width, unsigned new_height, unsigned new_width)
{
    unsigned y, x;
    float xFactor = (float)width / new_width;
    float yFactor = (float)height / new_height;
    float ox, oy;

    for (y = 0; y < new_height; y++ )
    {
        oy  = ((float)y + 0.5f) * yFactor - 0.5f;
        for (x = 0; x < new_width; x++ )
        {
            ox  = ((float)x  + 0.5f) * xFactor - 0.5f;
            d_im[y][x] = IMTInterpolateBiAkima (p_im, height, width, oy, ox);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void IMTFilterSCrop (IMTpixel** p_im, IMTpixel** d_im, unsigned height, unsigned width, unsigned new_height, unsigned new_width)
{
    unsigned int y, x;
    int yo, xo, dy, dx;

    dy = ((int) height - (int) new_height) >> 1;
    dx = ((int) width - (int) new_width) >> 1;
    for (y = 0; y < new_height; y++)
    {
        yo  = (int) y + dy;
        yo = (yo < 0) ? 0 : ((yo < (int) height) ? yo : (height -1));
        for (x = 0; x < new_width; x++ )
        {
            xo  = (int) x + dx;
            xo = (xo < 0) ? 0 : ((xo < (int) width) ? xo : (width -1));
            d_im[y][x] = p_im[yo][xo];
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void IMTFilterSize (IMTpixel** p_im, IMTpixel** d_im, int scaler, unsigned height, unsigned width, unsigned new_height, unsigned new_width)
{
    switch(scaler)
    {
    case SCALER_NEAREST:
        IMTFilterSNearest(p_im, d_im, height, width, new_height, new_width);
        break;
    case SCALER_BICONT:
        IMTFilterSBicont(p_im, d_im, height, width, new_height, new_width);
        break;
    case SCALER_BICUBIC:
        IMTFilterSBicub(p_im, d_im, height, width, new_height, new_width);
        break;
    case SCALER_BILINE:
        IMTFilterSBilin(p_im, d_im, height, width, new_height, new_width);
        break;
    case SCALER_BIAKIMA:
        IMTFilterSBiakima(p_im, d_im, height, width, new_height, new_width);
        break;
    case SCALER_CROP:
        IMTFilterSCrop(p_im, d_im, height, width, new_height, new_width);
        break;
    case SCALER_GSAMPLE:
        IMTFilterSGsample(p_im, d_im, height, width, new_height, new_width);
        break;
    }
}

////////////////////////////////////////////////////////////////////////////////

void IMTFilterSScaleNX (IMTpixel** p_im, IMTpixel** d_im, unsigned int height, unsigned int width, int smode)
{
    unsigned int y, x, y0, x0, y1, x1, y2, x2;
    IMTpixel pA, pB, pC, pD, pE, pF, pG, pH, pI;
    IMTpixel pBH, pDF;
    IMTpixel pEA, pEC, pEG, pEI;
    IMTpixel rA, rB, rC, rD, rE, rF, rG, rH, rI;

    smode = (smode < 2) ? 2 : (smode > 3) ? 3 : smode;

    if (smode == 2)
    {
        for (y = 0; y < height; y++)
        {
            y0 = (y > 0) ? (y - 1) : 0;
            y1 = y + 1;
            y1 = (y1 < height) ? y1 : (height - 1);
            y2 = y + y;
            for (x = 0; x < width; x++)
            {
                x0 = (x > 0) ? (x - 1) : 0;
                x1 = x + 1;
                x1 = (x1 < width) ? x1 : (width - 1);
                x2 = x + x;

                pB = p_im[y0][x];
                pD = p_im[y][x0];
                pE = p_im[y][x];
                pF = p_im[y][x1];
                pH = p_im[y1][x];

                rA = IMTcopy(pE);
                rB = IMTcopy(pE);
                rC = IMTcopy(pE);
                rD = IMTcopy(pE);
                pBH = IMTnear(pE, pB, pH);
                pDF = IMTnear(pE, pD, pF);
                if (IMTequal(pBH, pE) && IMTequal(pDF, pE))
                {
                    rA = IMTnear(pE, pB, pD);
                    rB = IMTnear(pE, pB, pF);
                    rC = IMTnear(pE, pH, pD);
                    rD = IMTnear(pE, pH, pF);
                }
                d_im[y2 + 0][x2 + 0] = rA;
                d_im[y2 + 0][x2 + 1] = rB;
                d_im[y2 + 1][x2 + 0] = rC;
                d_im[y2 + 1][x2 + 1] = rD;
            }
        }
    }
    else
    {
        for (y = 0; y < height; y++)
        {
            y0 = (y > 0) ? (y - 1) : 0;
            y1 = y + 1;
            y1 = (y1 < height) ? y1 : (height - 1);
            y2 = y + y + y;
            for (x = 0; x < width; x++)
            {
                x0 = (x > 0) ? (x - 1) : 0;
                x1 = x + 1;
                x1 = (x1 < width) ? x1 : (width - 1);
                x2 = x + x + x;

                pA = p_im[y0][x0];
                pB = p_im[y0][x];
                pC = p_im[y0][x1];
                pD = p_im[y][x0];
                pE = p_im[y][x];
                pF = p_im[y][x1];
                pG = p_im[y1][x0];
                pH = p_im[y1][x];
                pI = p_im[y1][x1];

                rA = IMTcopy(pE);
                rB = IMTcopy(pE);
                rC = IMTcopy(pE);
                rD = IMTcopy(pE);
                rE = IMTcopy(pE);
                rF = IMTcopy(pE);
                rG = IMTcopy(pE);
                rH = IMTcopy(pE);
                rI = IMTcopy(pE);
                pBH = IMTnear(pE, pB, pH);
                pDF = IMTnear(pE, pD, pF);
                if (IMTequal(pBH, pE) && IMTequal(pDF, pE))
                {
                    rA = IMTnear(pE, pB, pD);
                    rC = IMTnear(pE, pB, pF);
                    rG = IMTnear(pE, pD, pH);
                    rI = IMTnear(pE, pF, pH);

                    pEA = IMTnear(rA, pE, pA);
                    pEC = IMTnear(rC, pE, pC);
                    pEG = IMTnear(rG, pE, pG);
                    pEI = IMTnear(rI, pE, pI);

                    rB = IMTnear(pE, pEA, pEC);
                    rD = IMTnear(pE, pEA, pEG);
                    rF = IMTnear(pE, pEC, pEI);
                    rH = IMTnear(pE, pEG, pEI);
                }


                d_im[y2 + 0][x2 + 0] = rA;
                d_im[y2 + 0][x2 + 1] = rB;
                d_im[y2 + 0][x2 + 2] = rC;
                d_im[y2 + 1][x2 + 0] = rD;
                d_im[y2 + 1][x2 + 1] = rE;
                d_im[y2 + 1][x2 + 2] = rF;
                d_im[y2 + 2][x2 + 0] = rG;
                d_im[y2 + 2][x2 + 1] = rH;
                d_im[y2 + 2][x2 + 2] = rI;
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterSBWMag2 (BYTE** d_im, BYTE** r_im, unsigned height, unsigned width, unsigned height2, unsigned width2)
{
    unsigned y, x, y1, x1, y2, x2, n, s, threshold = 0;
    int im, er, ym, xm, yp, xp;
    int h = (int)height;
    int w = (int)width;

    n = 0;
    s = 0;
    er = 0;
    for (y = 0; y < height; y++ )
    {
        ym = (int)y;
        ym--;
        if (ym < 0)
        {
            ym = 0;
        }
        yp = (int)y;
        yp++;
        if (yp >= h)
        {
            yp = h - 1;
        }
        y1 = y * 2;
        if (y1 >= height2)
        {
            y1 = height2 -1;
        }
        y2 = y1 + 1;
        if (y2 >= height2)
        {
            y2 = height2 -1;
        }
        for (x = 0; x < width; x++)
        {
            xm = (int)x;
            xm--;
            if (xm < 0)
            {
                xm = 0;
            }
            xp = (int)x;
            xp++;
            if (xp >= w)
            {
                xp = w - 1;
            }
            x1 = 2 * x;
            if (x1 >= width2)
            {
                x1 = width2 -1;
            }
            x2 = x1 + 1;
            if (x2 >= width2)
            {
                x2 = width2 -1;
            }
            im = (int)d_im[y][x];
            er = (16 * im);
            im = (int)d_im[ym][x];
            er -= (2 * im);
            im = (int)d_im[yp][x];
            er -= (2 * im);
            im = (int)d_im[y][xm];
            er -= (2 * im);
            im = (int)d_im[y][xp];
            er -= (2 * im);
            im = (int)d_im[ym][xm];
            er -= im;
            im = (int)d_im[yp][xm];
            er -= im;
            im = (int)d_im[ym][xp];
            er -= im;
            im = (int)d_im[yp][xp];
            er -= im;
            er /= 12;

            im = (int)d_im[y][x];
            im += (int)d_im[ym][x];
            im += (int)d_im[y][xm];
            im += (int)d_im[ym][xm];
            im += er;
            im /= 4;
            if (im == 127)
            {
                if ((int)d_im[y][x] > 127)
                {
                    r_im[y1][x1] = (BYTE)255;
                    s++;
                }
                else
                {
                    r_im[y1][x1] = (BYTE)0;
                }
            }
            else if (im < 127)
            {
                r_im[y1][x1] = (BYTE)0;
            }
            else
            {
                r_im[y1][x1] = (BYTE)255;
                s++;
            }
            n++;
            im = (int)d_im[y][x];
            im += (int)d_im[yp][x];
            im += (int)d_im[y][xm];
            im += (int)d_im[yp][xm];
            im += er;
            im /= 4;
            if (im == 127)
            {
                if ((int)d_im[y][x] > 127)
                {
                    r_im[y2][x1] = (BYTE)255;
                    s++;
                }
                else
                {
                    r_im[y2][x1] = (BYTE)0;
                }
            }
            else if (im < 127)
            {
                r_im[y2][x1] = (BYTE)0;
            }
            else
            {
                r_im[y2][x1] = (BYTE)255;
                s++;
            }
            n++;
            im = (int)d_im[y][x];
            im += (int)d_im[ym][x];
            im += (int)d_im[y][xp];
            im += (int)d_im[ym][xp];
            im += er;
            im /= 4;
            if (im == 127)
            {
                if ((int)d_im[y][x] > 127)
                {
                    r_im[y1][x2] = (BYTE)255;
                    s++;
                }
                else
                {
                    r_im[y1][x2] = (BYTE)0;
                }
            }
            else if (im < 127)
            {
                r_im[y1][x2] = (BYTE)0;
            }
            else
            {
                r_im[y1][x2] = (BYTE)255;
                s++;
            }
            n++;
            im = (int)d_im[y][x];
            im += (int)d_im[yp][x];
            im += (int)d_im[y][xp];
            im += (int)d_im[yp][xp];
            im += er;
            im /= 4;
            if (im == 127)
            {
                if ((int)d_im[y][x] > 127)
                {
                    r_im[y2][x2] = (BYTE)255;
                    s++;
                }
                else
                {
                    r_im[y2][x2] = (BYTE)0;
                }
            }
            else if (im < 127)
            {
                r_im[y2][x2] = (BYTE)0;
            }
            else
            {
                r_im[y2][x2] = (BYTE)255;
                s++;
            }
            n++;
        }
    }
    threshold = s * 256 / n;

    return threshold;
}

////////////////////////////////////////////////////////////////////////////////

int IMTFilterSBWReduce2 (BYTE** d_im, BYTE** r_im, unsigned height, unsigned width, unsigned height2, unsigned width2)
{
    unsigned y, x, y1, x1, y2, x2, im, n, s, threshold = 0;
    int er, era;

    n = 0;
    s = 0;
    er = 0;
    for (y = 0; y < height2; y++ )
    {
        y1 = 2 * y;
        if (y1 >= height)
        {
            y1 = height -1;
        }
        y2 = y1 + 1;
        if (y2 >= height)
        {
            y2 = height -1;
        }
        for (x = 0; x < width2; x++)
        {
            x1 = 2 * x;
            if (x1 >= width)
            {
                x1 = width -1;
            }
            x2 = x1 + 1;
            if (x2 >= width)
            {
                x2 = width -1;
            }
            im = d_im[y1][x1];
            im += d_im[y2][x1];
            im += d_im[y1][x2];
            im += d_im[y2][x2];
            im /= 255;
            n++;
            if (im == 2)
            {
                if ((im + er) < 3)
                {
                    r_im[y][x] = 0;
                    s++;
                    era = (int)im;
                    er = era;
                }
                else
                {
                    r_im[y][x] = 255;
                    s++;
                    era = (int)im;
                    er = -era;
                }
            }
            else if (im < 2)
            {
                r_im[y][x] = 0;
                era = (int)im;
                er = era;
            }
            else
            {
                r_im[y][x] = 255;
                s++;
                era = (int)im;
                er = -era;
            }
        }
    }
    threshold = s * 256 / n;

    return threshold;
}

////////////////////////////////////////////////////////////////////////////////

void IMTFilterSGsample (IMTpixel** p_im, IMTpixel** d_im, unsigned height, unsigned width, unsigned new_height, unsigned new_width)
{
    unsigned d;
    float xFactor = (float)width / new_width;
    float yFactor = (float)height / new_height;
    float y0, x0, yn, xn, y0t, x0t, ynt, xnt, yt, xt, dyt, dxt, fyt, fxt, s, ss, sim, gy, gx, gc;
    int y, x, y0i, x0i, yni, xni, i, j, yti, xti, ypr, xpr, ynx, xnx;
    int h = (int)height;
    int w = (int)width;
    int h2 = (int)new_height;
    int w2 = (int)new_width;

    for (y = 0; y < h2; y++)
    {
        y0 = (float)y;
        y0 *= yFactor;
        y0i = (int)y0;
        yn = y0 + yFactor;
        yni = (int)yn;
        if ((float)yni < yn)
        {
            yni++;
        }
        for (x = 0; x < w2; x++)
        {
            x0 = (float)x;
            x0 *= xFactor;
            x0i = (int)x0;
            xn = x0 + xFactor;
            xni = (int)xn;
            if ((float)xni < xn)
            {
                xni++;
            }
            for (d = 0; d < 3; d++)
            {
                sim = 0;
                ss = 0;
                for (i = y0i; i < yni; i++)
                {
                    y0t = (float)i;
                    if (y0t < y0)
                    {
                        y0t = y0;
                    }
                    ynt = (float)(i + 1);
                    if (ynt > yn)
                    {
                        ynt = yn;
                    }
                    yt = (y0t + ynt) / 2;
                    yti = (int)yt;
                    if (yti >= h)
                    {
                        yti = h - 1;
                    }
                    dyt = ynt - y0t;
                    fyt = yt - yti - 0.5;
                    ypr = yti - 1;
                    if (ypr < 0)
                    {
                        ypr = 0;
                    }
                    ynx = yti + 1;
                    if (ynx >= h)
                    {
                        ynx = h - 1;
                    }
                    for (j = x0i; j < xni; j++)
                    {
                        x0t = (float)j;
                        if (x0t < x0)
                        {
                            x0t = x0;
                        }
                        xnt = (float)(j + 1);
                        if (xnt > xn)
                        {
                            xnt = xn;
                        }
                        xt = (x0t + xnt) / 2;
                        xti = (int)xt;
                        if (xti >= w)
                        {
                            xti = w - 1;
                        }
                        dxt = xnt - x0t;
                        fxt = xt - xti - 0.5;
                        xpr = xti - 1;
                        if (xpr < 0)
                        {
                            xpr = 0;
                        }
                        xnx = xti + 1;
                        if (xnx >= w)
                        {
                            xnx = w - 1;
                        }
                        gy = (float)p_im[ynx][xti].c[d];
                        gy += (float)p_im[ynx][xti].c[d];
                        gy += (float)p_im[ynx][xnx].c[d];
                        gy += (float)p_im[ynx][xpr].c[d];
                        gy -= (float)p_im[ypr][xti].c[d];
                        gy -= (float)p_im[ypr][xti].c[d];
                        gy -= (float)p_im[ypr][xnx].c[d];
                        gy -= (float)p_im[ypr][xpr].c[d];
                        gy /= 8.0;
                        gx = (float)p_im[yti][xnx].c[d];
                        gx += (float)p_im[yti][xnx].c[d];
                        gx += (float)p_im[ynx][xnx].c[d];
                        gx += (float)p_im[ypr][xnx].c[d];
                        gx -= (float)p_im[yti][xpr].c[d];
                        gx -= (float)p_im[yti][xpr].c[d];
                        gx -= (float)p_im[ynx][xpr].c[d];
                        gx -= (float)p_im[ypr][xpr].c[d];
                        gx /= 8.0;
                        gc = (float)p_im[yti][xti].c[d];
                        gc += fyt * gy;
                        gc += fxt * gx;
                        s = dxt * dyt;
                        ss += s;
                        sim += (gc * s);
                    }
                }
                if (ss < 0.0 || ss > 0.0)
                {
                    sim /= ss;
                }
                d_im[y][x].c[d] = ByteClamp((int)(sim + 0.5));
            }
            d_im[y][x] = IMTcalcS (d_im[y][x]);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void IMTFilterSHRIS (IMTpixel** p_im, IMTpixel** d_im, unsigned height, unsigned width, int smode)
{
    unsigned d;
    int y, x, y2, x2, yp1, xp1, yp2, xp2, yn1, xn1, yn2, xn2;
    int h = (int)height;
    int w = (int)width;
    int h2 = h * smode;
    int w2 = w * smode;

    float imx;
    float b11, b12, b13, b21, b22, b23, b31, b32, b33;
    float r11, r12, r13, r21, r22, r23, r31, r32, r33;

    float k0[2];
    float k1[3];
    float k2[3];

    smode = (smode < 2) ? 2 : (smode > 3) ? 3 : smode;

    if (smode == 2)
    {
        k0[0] = 0.750;
        k0[1] = 0.250;
        k1[0] = k0[0] * k0[0];
        k1[1] = k0[0] * k0[1];
        k1[2] = k0[1] * k0[1];
        k2[0] = k1[0];
        k2[1] = k1[1] / 2;
        k2[2] = k1[2] / 4;
    }
    else
    {
        k0[0] = 2.0 / 3.0;
        k0[1] = 1.0 / 3.0;
        k1[0] = k0[0] * k0[0];
        k1[1] = k0[0] * k0[1];
        k1[2] = k0[1] * k0[1];
        k2[0] = (1 + 4 * k0[0] + 4 * k1[0]) / 9;
        k2[1] = (k0[1] + 2 * k1[1]) / 9;
        k2[2] = k1[2] / 9;
    }

    for (y = 0; y < h; y++)
    {
        yp1 = y - 1;
        if (yp1 < 0)
        {
            yp1 = 0;
        }
        yp2 = y - 2;
        if (yp2 < 0)
        {
            yp2 = 0;
        }
        yn1 = y + 1;
        if (yn1 > h - 1)
        {
            yn1 = h - 1;
        }
        yn2 = y + 2;
        if (yn2 > h - 1)
        {
            yn2 = h - 1;
        }
        for (x = 0; x < w; x++)
        {
            xp1 = x - 1;
            if (xp1 < 0)
            {
                xp1 = 0;
            }
            xp2 = x - 2;
            if (xp2 < 0)
            {
                xp2 = 0;
            }
            xn1 = x + 1;
            if (xn1 > w - 1)
            {
                xn1 = w - 1;
            }
            xn2 = x + 2;
            if (xn2 > w - 1)
            {
                xn2 = w - 1;
            }
            for (d = 0; d < 3; d++)
            {
                imx = (float)p_im[yp2][xp2].c[d];
                b11 = -(k2[2] * imx);
                imx = (float)p_im[yp2][xp1].c[d];
                b11 -= (k2[1] * imx);
                b12 = -(k2[2] * imx);
                imx = (float)p_im[yp2][x].c[d];
                b11 -= (k2[2] * imx);
                b12 -= (k2[1] * imx);
                b13 = -(k2[2] * imx);
                imx = (float)p_im[yp2][xn1].c[d];
                b12 -= (k2[2] * imx);
                b13 -= (k2[1] * imx);
                imx = (float)p_im[yp2][xn2].c[d];
                b13 -= (k2[2] * imx);
                imx = (float)p_im[yp1][xp2].c[d];
                b11 -= (k2[1] * imx);
                b21 = -(k2[2] * imx);
                imx = (float)p_im[yp1][xp1].c[d];
                b11 += ((2.0 - k2[0]) * imx);
                b12 -= (k2[1] * imx);
                b21 -= (k2[1] * imx);
                b22 = -(k2[2] * imx);
                imx = (float)p_im[yp1][x].c[d];
                b11 -= (k2[1] * imx);
                b12 += ((2.0 - k2[0]) * imx);
                b13 -= (k2[1] * imx);
                b21 -= (k2[2] * imx);
                b22 -= (k2[1] * imx);
                b23 = -(k2[2] * imx);
                imx = (float)p_im[yp1][xn1].c[d];
                b12 -= (k2[1] * imx);
                b13 += ((2.0 - k2[0]) * imx);
                b22 -= (k2[2] * imx);
                b23 -= (k2[1] * imx);
                imx = (float)p_im[yp1][xn2].c[d];
                b13 -= (k2[1] * imx);
                b23 -= (k2[2] * imx);
                imx = (float)p_im[y][xp2].c[d];
                b11 -= (k2[2] * imx);
                b21 -= (k2[1] * imx);
                b31 = -(k2[2] * imx);
                imx = (float)p_im[y][xp1].c[d];
                b11 -= (k2[1] * imx);
                b12 -= (k2[2] * imx);
                b21 += ((2.0 - k2[0]) * imx);
                b22 -= (k2[1] * imx);
                b31 -= (k2[1] * imx);
                b32 = -(k2[2] * imx);
                imx = (float)p_im[y][x].c[d];
                b11 -= (k2[2] * imx);
                b12 -= (k2[1] * imx);
                b13 -= (k2[2] * imx);
                b21 -= (k2[1] * imx);
                b22 += ((2.0 - k2[0]) * imx);
                b23 -= (k2[1] * imx);
                b31 -= (k2[2] * imx);
                b32 -= (k2[1] * imx);
                b33 = -(k2[2] * imx);
                imx = (float)p_im[y][xn1].c[d];
                b12 -= (k2[2] * imx);
                b13 -= (k2[1] * imx);
                b22 -= (k2[1] * imx);
                b23 += ((2.0 - k2[0]) * imx);
                b32 -= (k2[2] * imx);
                b33 -= (k2[1] * imx);
                imx = (float)p_im[y][xn2].c[d];
                b13 -= (k2[2] * imx);
                b23 -= (k2[1] * imx);
                b33 -= (k2[2] * imx);
                imx = (float)p_im[yn1][xp2].c[d];
                b21 -= (k2[2] * imx);
                b31 -= (k2[1] * imx);
                imx = (float)p_im[yn1][xp1].c[d];
                b21 -= (k2[1] * imx);
                b22 -= (k2[2] * imx);
                b31 += ((2.0 - k2[0]) * imx);
                b32 -= (k2[1] * imx);
                imx = (float)p_im[yn1][x].c[d];
                b21 -= (k2[2] * imx);
                b22 -= (k2[1] * imx);
                b23 -= (k2[2] * imx);
                b31 -= (k2[1] * imx);
                b32 += ((2.0 - k2[0]) * imx);
                b33 -= (k2[1] * imx);
                imx = (float)p_im[yn1][xn1].c[d];
                b22 -= (k2[2] * imx);
                b23 -= (k2[1] * imx);
                b32 -= (k2[1] * imx);
                b33 += ((2.0 - k2[0]) * imx);
                imx = (float)p_im[yn1][xn2].c[d];
                b23 -= (k2[2] * imx);
                b33 -= (k2[1] * imx);
                imx = (float)p_im[yn2][xp2].c[d];
                b31 -= (k2[2] * imx);
                imx = (float)p_im[yn2][xp1].c[d];
                b31 -= (k2[1] * imx);
                b32 -= (k2[2] * imx);
                imx = (float)p_im[yn2][x].c[d];
                b31 -= (k2[2] * imx);
                b32 -= (k2[1] * imx);
                b33 -= (k2[2] * imx);
                imx = (float)p_im[yn2][xn1].c[d];
                b32 -= (k2[2] * imx);
                b33 -= (k2[1] * imx);
                imx = (float)p_im[yn2][xn2].c[d];
                b33 -= (k2[2] * imx);

                if (smode == 2)
                {
                    r11 = (k1[0] * b22 + k1[1] * (b12 + b21) + k1[2] * b11);
                    r12 = (k1[0] * b22 + k1[1] * (b12 + b23) + k1[2] * b13);
                    r21 = (k1[0] * b22 + k1[1] * (b32 + b21) + k1[2] * b31);
                    r22 = (k1[0] * b22 + k1[1] * (b32 + b23) + k1[2] * b33);

                    y2 = y * 2;
                    x2 = x * 2;

                    d_im[y2][x2].c[d] = ByteClamp((int)(r11 + 0.5));
                    d_im[y2][x2 + 1].c[d] = ByteClamp((int)(r12 + 0.5));
                    d_im[y2 + 1][x2].c[d] = ByteClamp((int)(r21 + 0.5));
                    d_im[y2 + 1][x2 + 1].c[d] = ByteClamp((int)(r22 + 0.5));
                }
                else
                {
                    r11 = (k1[0] * b22 + k1[1] * (b12 + b21) + k1[2] * b11);
                    r12 = (k0[0] * b22 + k0[1] * b12);
                    r13 = (k1[0] * b22 + k1[1] * (b12 + b23) + k1[2] * b13);
                    r21 = (k0[0] * b22 + k0[1] * b21);
                    r22 = (float)p_im[y][x].c[d];
                    r23 = (k0[0] * b22 + k0[1] * b23);
                    r31 = (k1[0] * b22 + k1[1] * (b32 + b21) + k1[2] * b31);
                    r32 = (k0[0] * b22 + k0[1] * b32);
                    r33 = (k1[0] * b22 + k1[1] * (b32 + b23) + k1[2] * b33);

                    y2 = y * 3;
                    x2 = x * 3;

                    d_im[y2][x2].c[d] = ByteClamp((int)(r11 + 0.5));
                    d_im[y2][x2 + 1].c[d] = ByteClamp((int)(r12 + 0.5));
                    d_im[y2][x2 + 2].c[d] = ByteClamp((int)(r13 + 0.5));
                    d_im[y2 + 1][x2].c[d] = ByteClamp((int)(r21 + 0.5));
                    d_im[y2 + 1][x2 + 1].c[d] = ByteClamp((int)(r22 + 0.5));
                    d_im[y2 + 1][x2 + 2].c[d] = ByteClamp((int)(r23 + 0.5));
                    d_im[y2 + 2][x2].c[d] = ByteClamp((int)(r31 + 0.5));
                    d_im[y2 + 2][x2 + 1].c[d] = ByteClamp((int)(r32 + 0.5));
                    d_im[y2 + 2][x2 + 2].c[d] = ByteClamp((int)(r33 + 0.5));
                }
            }
        }
    }
    for (y = 0; y < h2; y++)
    {
        for (x = 0; x < w2; x++)
        {
            d_im[y][x] = IMTcalcS (d_im[y][x]);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void IMTFilterSGSampleUp (IMTpixel** p_im, IMTpixel** d_im, unsigned height, unsigned width, int smode)
{
    unsigned y, x, yf, xf, yr, xr, d, k, l, bpp = 3;
    int dx, dy, ir, isz;
    IMTpixel wm, wr;
    float imx, k0, kx, ky, kz, dd;

    if (smode < 2)
    {
        smode = 2;
    }
    ir = 2 * smode;
    isz = ir * ir;
    dd = 1.0 / (float)isz;

    for (y = 0; y < height; y++)
    {
        yr = y * smode;
        for (x = 0; x < width; x++)
        {
            xr = x * smode;
            for (d = 0; d < bpp; d++)
            {
                imx = (float)p_im[y][x].c[d];
                imx *= 36.0;
                yf = (y > 0) ? (y - 1) : y;
                xf = (x > 0) ? (x - 1) : x;
                imx -= (float)p_im[yf][xf].c[d];
                imx -= (float)p_im[yf][x].c[d];
                imx -= (float)p_im[yf][x].c[d];
                xf = (x < width - 1) ? (x + 1) : x;
                imx -= (float)p_im[yf][xf].c[d];
                imx -= (float)p_im[y][xf].c[d];
                imx -= (float)p_im[y][xf].c[d];
                xf = (x > 0) ? (x - 1) : x;
                imx -= (float)p_im[y][xf].c[d];
                imx -= (float)p_im[y][xf].c[d];
                yf = (y < height - 1) ? (y + 1) : y;
                imx -= (float)p_im[yf][xf].c[d];
                imx -= (float)p_im[yf][x].c[d];
                imx -= (float)p_im[yf][x].c[d];
                xf = (x < width - 1) ? (x + 1) : x;
                imx -= (float)p_im[yf][xf].c[d];
                imx /= 24.0;
                wm.c[d] = ByteClamp((int)(imx + 0.5));
            }
            for (k =  0; k < smode; k++)
            {
                dy = (int)(k + k - smode + 1);
                yf = (dy < 0) ? ((y > 0) ? (y - 1) : y) : ((y < height - 1) ? (y + 1) : y);
                for (l =  0; l < smode; l++)
                {
                    dx = (int)(l + l - smode + 1);
                    xf = (dx < 0) ? ((x > 0) ? (x - 1) : x) : ((x < width - 1) ? (x + 1) : x);
                    k0 = (float)isz;
                    kx = (float)ABS(dx * ir);
                    ky = (float)ABS(dy * ir);
                    kz = (float)ABS(dx * dy);
                    kx -= kz;
                    ky -= kz;
                    k0 -= kz;
                    k0 -= kx;
                    k0 -= ky;
                    for (d = 0; d < bpp; d++)
                    {
                        imx = (k0 * (float)wm.c[d]);
                        imx += (ky * (float)p_im[yf][x].c[d]);
                        imx += (kx * (float)p_im[y][xf].c[d]);
                        imx += (kz * (float)p_im[yf][xf].c[d]);
                        imx *= dd;
                        wr.c[d] = ByteClamp((int)(imx + 0.5));
                    }
                    d_im[yr + k][xr + l] = wr;
                }
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void IMTFilterSReduce (IMTpixel** p_im, IMTpixel** d_im, unsigned height, unsigned width, unsigned reduce)
{
    unsigned y, x, y0, x0, y1, x1;
    if (reduce > 1)
    {
        unsigned widthr = (width + reduce - 1) / reduce;
        unsigned heightr = (height + reduce - 1) / reduce;
        for (y = 0; y < heightr; y++)
        {
            y0 = y * reduce;
            y1 = y0 + reduce;
            if (y1 > height)
            {
                y1 = height;
            }
            for (x = 0; x < widthr; x++)
            {
                x0 = x * reduce;
                x1 = x0 + reduce;
                if (x1 > width)
                {
                    x1 = width;
                }
                d_im[y][x] = IMTmeanIc(p_im, y0, x0, y1, x1);
            }
        }
    }
    else
    {
        for (y = 0; y < height; y++)
        {
            for (x = 0; x < width; x++)
            {
                d_im[y][x] = p_im[y][x];
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void IMTFilterSNearest (IMTpixel** p_im, IMTpixel** d_im, unsigned height, unsigned width, unsigned new_height, unsigned new_width)
{
    float xFactor = (float) width / new_width;
    float yFactor = (float) height / new_height;
    int y, x, ox, oy;
    int h = new_height;
    int w = new_width;

    for (y = 0; y < h; y++)
    {
        oy = (int)(y * yFactor);
        for (x = 0; x < w; x++)
        {
            ox = (int)(x * xFactor);
            d_im[y][x] = p_im[oy][ox];
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void IMTFilterSFRP (IMTpixel** p_im, IMTpixel** d_im, unsigned height, unsigned width, int smode)
{
    unsigned y, x, yf, xf, ys, xs, ysf, xsf, ym, xm, yr, xr, d, k, l, bpp = 3, r = 1, rd = (2 * r + 1);
    unsigned widths, heights, imdm, imdt, imdms = (512 * bpp * (rd * rd + 4));
    int im, ims, imd;
    IMTpixel **s_im, **r_im;

    if (smode < 2)
    {
        smode = 2;
    }
    heights = (height + smode - 1) / smode;
    widths = (width + smode - 1) / smode;
    s_im = IMTalloc(heights, widths);
    r_im = IMTalloc(rd, rd);
    IMTFilterSReduce(p_im, s_im, height, width, smode);

    for (y = 0; y < height; y++)
    {
        yr = y * smode;
        for (x = 0; x < width; x++)
        {
            xr = x * smode;
            for (k = 0; k < rd; k++)
            {
                yf = IndexClamp(((int)(y + k) - r), (height - 1));
                for (l = 0; l < rd; l++)
                {
                    xf = IndexClamp(((int)(x + l) - r), (width - 1));
                    r_im[k][l] = p_im[yf][xf];
                }
            }
            ym = 0;
            xm = 0;
            imdm = imdms;
            for (ys = 0; ys < heights; ys++)
            {
                for (xs = 0; xs < widths; xs++)
                {
                    imdt = 0;
                    im = (int)p_im[y][x].s;
                    ims = (int)s_im[ys][xs].s;
                    imd = (im > ims) ? (im - ims) : (ims - im);
                    imdt += imd;
                    imdt += imd;
                    imdt += imd;
                    imdt += imd;
                    for (d = 0; d < bpp; d++)
                    {
                        im = (int)p_im[y][x].c[d];
                        ims = (int)s_im[ys][xs].c[d];
                        imd = (im > ims) ? (im - ims) : (ims - im);
                        imdt += imd;
                        imdt += imd;
                        imdt += imd;
                        imdt += imd;
                    }
                    for (k = 0; k < rd; k++)
                    {
                        ysf = IndexClamp(((int)(ys + k) - r), (heights - 1));
                        for (l = 0; l < rd; l++)
                        {
                            xsf = IndexClamp(((int)(xs + l) - r), (widths - 1));
                            im = (int)r_im[k][l].s;
                            ims = (int)s_im[ysf][xsf].s;
                            imd = (im > ims) ? (im - ims) : (ims - im);
                            imdt += imd;
                            for (d = 0; d < bpp; d++)
                            {
                                im = (int)r_im[k][l].c[d];
                                ims = (int)s_im[ysf][xsf].c[d];
                                imd = (im > ims) ? (im - ims) : (ims - im);
                                imdt += imd;
                            }
                        }
                    }
                    if (imdt < imdm)
                    {
                        ym = ys;
                        xm = xs;
                        imdm = imdt;
                    }
                }
            }
            ym *= smode;
            xm *= smode;
            for (k =  0; k < smode; k++)
            {
                yf = IndexClamp((ym + k), (height - 1));
                for (l =  0; l < smode; l++)
                {
                    xf = IndexClamp((xm + l), (width - 1));
                    d_im[yr + k][xr + l] = p_im[yf][xf];
                }
            }
        }
    }
    IMTfree(r_im, rd);
    IMTfree(s_im, heights);
}

////////////////////////////////////////////////////////////////////////////////

void IMTReduceBW (BYTE** m_im, BYTE** g_im, unsigned height, unsigned width, unsigned heightg, unsigned widthg, unsigned kred, BYTE preval, BYTE result)
{
    unsigned y, x, y0, x0, y1, x1, i, j;
    BYTE val;
    for (y = 0; y < heightg; y++)
    {
        y0 = y * kred;
        y1 = (((y0 + kred) < height) ? (y0 + kred) : height);
        for (x = 0; x < widthg; x++)
        {
            x0 = x * kred;
            x1 = (((x0 + kred) < width) ? (x0 + kred) : width);
            val = preval;
            for (i = y0; i < y1; i++)
            {
                for (j = x0; j < x1; j++)
                {
                    if (m_im[i][j] == 0)
                    {
                        val = result;
                    }
                }
            }
            g_im[y][x] = val;
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
