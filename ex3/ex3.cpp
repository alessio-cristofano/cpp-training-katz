#include "../utilities_lib.hpp"

// ----------------------------------------------------------------------
//  Function declaration
//
int matrix(vector<vector<float>> &, int &, vector<float> &);
// ----------------------------------------------------------------------
int main()
{
    vector<vector<float>> ep;
    vector<vector<float>> ept;
    vector<vector<float>> pt1;
    vector<vector<float>> pt2;
    vector<vector<float>> co;
    vector<float> th;
    // --------------------------------------
    bool matrix_red = true; // perform\ matrix reduction or not

    ifstream file("AFOIL2.DAT");
    if (!file.good())
    {
        cout << "File AFOIL2.DAT was not found" << endl;
        return 1;
    }

    int n = 0;
    string file_row;
    while (getline(file, file_row))
    {
        vector<string> points_xy = split_string(file_row, ",");
        ept.push_back({stof(points_xy[0]), stof(points_xy[1])});
        n++;
    }
    file.close();
    //
    const int m = n - 1; // get number of panels

    vector<vector<float>> a(m, vector<float>(n));
    vector<vector<float>> b(m, vector<float>(m));
    vector<float> g(m);

    //
    // Convert paneling to clockwise
    //
    for (int i = 0; i < n; i++)
    {
        ep.push_back({ept[n - i - 1][0], ept[n - i - 1][1]});
    }
    // Establish coordinates of panel end points
    //
    for (int i = 0; i < m; i++)
    {
        pt1.push_back({ep[i][0], ep[i][1]});
        pt2.push_back({ep[i + 1][0], ep[i + 1][1]});
        float dz = pt2[i][1] - pt1[i][1];
        float dx = pt2[i][0] - pt1[i][0];
        th.push_back(atan2(dz, dx));
        // Collocation points
        float co_x = (pt2[i][0] - pt1[i][0]) / 2 + pt1[i][0];
        float co_z = (pt2[i][1] - pt1[i][1]) / 2 + pt1[i][1];
        co.push_back({co_x, co_z});
    }

    // Establish Influence Coefficients
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < m; j++)
        {
            float xt = co[i][0] - pt1[j][0];
            float zt = co[i][1] - pt1[j][1];
            float x2t = pt2[j][0] - pt1[j][0];
            float z2t = pt2[j][1] - pt1[j][1];
            float x = xt * cos(th[j]) + zt * sin(th[j]);
            float z = -xt * sin(th[j]) + zt * cos(th[j]);
            float x2 = x2t * cos(th[j]) + z2t * sin(th[j]);
            float z2 = 0;

            float r1 = sqrt(x * x + z * z);
            float r2 = sqrt((x - x2) * (x - x2) + z * z);
            float th1 = atan2(z, x);
            float th2 = atan2(z, x - x2);

            float ul = 0;
            float wl = 0.5;
            if (i != j)
            {
                ul = 1 / (2 * PI) * log(r1 / r2);
                wl = 1 / (2 * PI) * (th2 - th1);
            }

            float u = ul * cos(th[j]) - wl * sin(th[j]);
            float w = ul * sin(th[j]) + wl * cos(th[j]);

            a[i][j] = -u * sin(th[i]) + w * cos(th[i]);
            b[i][j] = u * cos(th[i]) + w * sin(th[i]);
        }
        a[i][n - 1] = sin(th[i]);
    }
    ofstream out_file("CPS_cpp.DAT");

    if (matrix_red == true)
    {
        int solution_found = matrix(a, n, g);
        cout << "Matrix Reduction" << endl;
        if (solution_found != 0)
            return 1;
    }

    for (int i = 0; i < m; i++)
    {
        float vel = 0;
        for (int j = 0; j < m; j++)
        {
            vel += b[i][j] * g[j];
        }
        float cp = 1 - pow((vel + cos(th[i])), 2);
        out_file << co[i][0] << " ," << cp << endl;
    }

    return 0;
}

int matrix(vector<vector<float>> &a, int &n, vector<float> &g)
{
    for (int i = 0; i < n - 1; i++)
    {
        while (abs(a[i][i]) < 0.0000001)
        {
            if (i != n - 2)
            {
                for (int j = 0; j < n; j++)
                {
                    float tmp = a[i][j];
                    a[i][j] = a[i + 1][j];
                    a[i + 1][j] = tmp;
                }
            }
            else
            {
                cout << "Unable to find solution. Matrix is singular" << endl;
                return 1;
            }
        }
        // Convert matrix A to upper triangular
        float p = a[i][i];
        for (int j = i; j < n; j++)
            a[i][j] /= p;
        for (int k = i + 1; k < n - 1; k++)
        {
            float p2 = a[k][i];
            for (int l = i; l < n; l++)
                a[k][l] -= p2 * a[i][l];
        }
    }
    // Back substitution
    for (int i = n - 2; i >= 0; i--)
    {
        g[i] = a[i][n - 1];
        for (int j = 0; j < n - 1; j++)
        {
            a[i][i] = 0;
            g[i] = g[i] - a[i][j] * g[j];
        }
    }
    return 0;
}
