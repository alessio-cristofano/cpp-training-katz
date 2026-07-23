#include "../utilities_lib.hpp"
// ----------------------------------------------------------------------
//  Function declaration
//
int matrix(vector<vector<float>>&, vector<float>&);
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
    float alpha = deg2rad(0);
    // 
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
        ept.push_back({ stof(points_xy[0]), stof(points_xy[1]) });
        n++;
    }
    file.close();
    ofstream out_file("CPLV_cpp.DAT");
    //
    const int m = n - 1; // get number of panels

    vector<vector<float>> a(n, vector<float>(n + 1));
    vector<vector<float>> b(n, vector<float>(n));
    vector<float> g(n);
    vector<float> dl(m);

    //
    // Convert paneling to clockwise
    //
    for (int i = 0; i < n; i++)
    {
        ep.push_back({ ept[n - i - 1][0], ept[n - i - 1][1] });
    }
    // Establish coordinates of panel end points & collocation points
    //
    for (int i = 0; i < m; i++)
    {
        pt1.push_back({ ep[i][0], ep[i][1] });
        pt2.push_back({ ep[i + 1][0], ep[i + 1][1] });
        float dz = pt2[i][1] - pt1[i][1];
        float dx = pt2[i][0] - pt1[i][0];
        th.push_back(atan2(dz, dx));
        // Collocation points
        float co_x = (pt2[i][0] - pt1[i][0]) / 2 + pt1[i][0];
        float co_z = (pt2[i][1] - pt1[i][1]) / 2 + pt1[i][1];
        co.push_back({ co_x, co_z });
    }

    float hold_a;
    float hold_b;
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

            if (i == 0) dl.at(j) = x2;

            float r1 = sqrt(x * x + z * z);
            float r2 = sqrt((x - x2) * (x - x2) + z * z);
            float th1 = atan2(z, x);
            float th2 = atan2(z, x - x2);

            float u1l = -(z * log(r2 / r1) + x * (th2 - th1) - x2 * (th2 - th1)) / (2 * x2 * PI);
            float u2l = (z * log(r2 / r1) + x * (th2 - th1)) / (2 * x2 * PI);
            float w1l = -((x2 - z * (th2 - th1)) - x * log(r1 / r2) + x2 * log(r1 / r2)) / (2 * x2 * PI);
            float w2l = ((x2 - z * (th2 - th1)) - x * log(r1 / r2)) / (2 * x2 * PI);

            if (i == j) {
                u1l = -0.5 * (x - x2) / x2;
                u2l = 0.5 * x / x2;
                w1l = -1 / (2 * PI);
                w2l = -w1l;
            }

            float u1 = u1l * cos(th[j]) - w1l * sin(th[j]);
            float u2 = u2l * cos(th[j]) - w2l * sin(th[j]);
            float w1 = u1l * sin(th[j]) + w1l * cos(th[j]);
            float w2 = u2l * sin(th[j]) + w2l * cos(th[j]);

            if (j == 0) {
                a.at(i).at(0) = -u1 * sin(th.at(i)) + w1 * cos(th.at(i));
                hold_a = -u2 * sin(th.at(i)) + w2 * cos(th.at(i));
                b.at(i).at(0) = u1 * cos(th.at(i)) + w1 * sin(th.at(i));
                hold_b = u2 * cos(th.at(i)) + w2 * sin(th.at(i));
            }
            else if (j == m - 1) {
                a.at(i).at(m - 1) = -u1 * sin(th.at(i)) + w1 * cos(th.at(i)) + hold_a;
                a.at(i).at(m) = -u2 * sin(th.at(i)) + w2 * cos(th.at(i));
                b.at(i).at(m - 1) = u1 * cos(th.at(i)) + w1 * sin(th.at(i)) + hold_b;
                b.at(i).at(m) = u2 * cos(th.at(i)) + w2 * sin(th.at(i));
            }
            else {
                a.at(i).at(j) = -u1 * sin(th.at(i)) + w1 * cos(th.at(i)) + hold_a;
                hold_a = -u2 * sin(th.at(i)) + w2 * cos(th.at(i));
                b.at(i).at(j) = u1 * cos(th.at(i)) + w1 * sin(th.at(i)) + hold_b;
                hold_b = u2 * cos(th.at(i)) + w2 * sin(th.at(i));
            }
        }

        a.at(i).at(n) = cos(alpha) * sin(th.at(i)) - sin(alpha) * cos(th.at(i));
    }
    // Kutta Condition
    //
    a.at(n - 1).at(0) = 1.0;
    a.at(n - 1).at(n - 1) = 1.0;

    // Solve with Gauss Elimination Method
    // 
    matrix(a, g);
    float cl = 0;


    for (int i = 0;i < m; i++) {
        float temp = 0;
        for (int j = 0;j < n; j++) {
            temp += b.at(i).at(j) * g.at(j);
        }
        float v = temp + cos(alpha) * cos(th.at(i)) + sin(alpha) * sin(th.at(i));
        cl += v * dl.at(i);
        float cp = 1 - v * v;
        out_file << co.at(i).at(0) << " ," << cp << endl;
    }
    out_file << "LIFT COEFFICIENT = " << cl << endl;


    return 0;
}

int matrix(vector<vector<float>>& a, vector<float>& g)
{
    int n = a.at(0).size();
    for (int i = 0; i < n - 1; i++)
    {
        while (fabs(a[i][i]) < 0.000000001)
        {
            if (i != n - 2)
            {
                for (int j = 0; j < n; j++)
                {
                    float tmp = a[i][j];
                    cout << "swapping zeros" << endl;
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
