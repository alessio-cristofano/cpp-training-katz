#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "../helpers.hpp"

using std::abs;
using std::atan;
using std::cin;
using std::cos;
using std::cout;
using std::endl;
using std::pow;
using std::sin;
using std::sqrt;
using std::vector;

// ----------------------------------------------------------------------
void vor2d(float &, float &, float &, float &, const float &, float &, float &, float &);
void decomp(const int &, vector<int> &, vector<vector<float>> &);
void solver(const int &, vector<vector<float>> &, vector<float> &, vector<int> &);

int main()
{
    const int n = 10;
    vector<float> gamma(n);
    vector<float> xc(n);
    vector<float> zc(n);
    vector<float> x(n);
    vector<float> z(n);
    vector<float> enx(n);
    vector<float> enz(n);
    vector<int> ip(n);
    vector<float> dl(n);
    vector<float> dcp(n);
    vector<float> dcp1(n);
    vector<vector<float>> a(n, vector<float>(n));

    float c = 1;
    float epsilon = 0.1 * c;
    float alpha_deg = 10;
    float alpha_rad = deg2rad(alpha_deg);
    float rho = 1;
    float v = 1;
    float u_inf = v * cos(alpha_rad);
    float w_inf = v * sin(alpha_rad);
    float que = 0.5 * rho * v * v;
    //
    float dx = c / n;
    //
    // Grid Generation

    for (int i = 0; i < n; i++)
    {
        xc[i] = dx * (i + 1 - 0.25);
        zc[i] = 4 * epsilon * xc[i] / c * (1 - xc[i] / c);
        x[i] = dx * (i + 1 - 0.75);
        z[i] = 4 * epsilon * x[i] / c * (1 - x[i] / c);

        float detadx = 4. * epsilon / c * (1 - 2 * xc[i] / c);
        float sq = sqrt(1 + detadx * detadx);
        enx[i] = -detadx / sq;
        enz[i] = 1 / sq;
    }
    // //
    // // Influence Coefficients
    // //
    float gamma_strength = 1.0;

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            vor2d(xc[i], zc[i], x[j], z[j], gamma_strength, enx[i], enz[i], a[i][j]);
        }
        gamma[i] = -u_inf * enx[i] - w_inf * enz[i];
    }

    decomp(n, ip, a);
    solver(n, a, gamma, ip);

    //
    //  Compute Aerodynamic Loads
    //
    float bl = 0;
    for (int i = 0; i < n; i++)
    {
        dl[i] = rho * v * gamma[i];
        dcp[i] = dl[i] / dx / que;
        float dd = 32 * epsilon / c * sqrt(x[i] / c * (1 - x[i] / c));
        dcp1[i] = 4 * sqrt((c - x[i]) / x[i]) * alpha_rad + dd;
        bl += dl[i];
    }
    float cl = bl / (que * c);
    float cl1 = 2 * PI * (alpha_rad + 2 * epsilon / c);
    // cl1 and dcp1 are Analytical solutions
    //
    //  Output Printing
    //
    std::ofstream file("FILE.DAT");
    file << " THIN AIRFOIL WITH ELLIPTIC CAMBER " << "\n";
    file << " V=" << v << ", CL=" << cl << ", CL(EXACT)=" << cl1 << ", N=" << n << ", ALPHA=" << alpha_deg << "\n";
    for (int i = 0; i < n; i++)
    {
        file << " I=" << i << ", X=" << x[i] << ", DCP=" << dcp[i] << ", DCP (EXACT)=" << dcp1[i] << "\n";
    }

    return 0;
}

void vor2d(float &x, float &z, float &x1, float &z1, const float &gamma, float &enx_i, float &enz_i, float &a_ij)
{
    float rx = x - x1;
    float rz = z - z1;
    float r = sqrt(rx * rx + rz * rz);
    if (r < 0.001)
    {
        a_ij = 0;
        return;
    }
    float u = 0.5 / PI * gamma / r * rz / r;
    float w = -0.5 / PI * gamma / r * rx / r;
    a_ij = u * enx_i + w * enz_i;
    return;
}

void decomp(const int &n, vector<int> &ip, vector<vector<float>> &a)
{
    ip[n - 1] = 1;
    for (int k = 0; k < n; k++)
    {
        if (k == n - 1)
        {
            if (a[k][k] == 0.0)
                ip[n - 1] = 0;
            break; // Last row, nothing left to eliminate
        }
        int kp1 = k + 1;
        int m = k;
        for (int i = kp1; i < n; i++)
        {
            if (abs(a[i][k]) > abs(a[m][k]))
            {
                m = i;
            }
        }
        ip[k] = m;
        if (m != k)
        {
            ip[n - 1] = -ip[n - 1];
        }
        float t = a[m][k];
        a[m][k] = a[k][k];
        a[k][k] = t;
        if (t == 0.0 && a[k][k] == 0.0)
        {
            ip[n - 1] = 0;
            continue;
        }
        for (int i = kp1; i < n; i++)
        {
            a[i][k] = -a[i][k] / t;
        }
        for (int j = kp1; j < n; j++)
        {
            t = a[m][j];
            a[m][j] = a[k][j];
            a[k][j] = t;
            if (t == 0.0)
                continue;
            for (int i = kp1; i < n; i++)
            {
                a[i][j] = a[i][j] + a[i][k] * t;
            }
        }
        if (a[k][k] == 0.0)
            ip[n - 1] = 0;
    }
}

void solver(const int &n, vector<vector<float>> &a, vector<float> &b, vector<int> &ip)
{

    if (n == 1)
    {
        b[0] /= a[0][0];
        return;
    }

    int nm1 = n - 1;
    for (int k = 0; k < nm1; k++)
    {
        int kp1 = k + 1;
        int m = ip[k];
        float t = b[m];
        b[m] = b[k];
        b[k] = t;
        for (int i = kp1; i < n; i++)
        {
            b[i] += a[i][k] * t;
        }
    }

    for (int k = n - 1; k > 0; k--) // da k = 9 (ultimo) a k = 1 (penultimo)
    {
        int km1 = k - 1;
        b[k] /= a[k][k];
        float t = -b[k];
        for (int i = 0; i < k; i++)
        {
            b[i] += a[i][k] * t;
        }
    }
    b[0] /= a[0][0];

    return;
}
