#include "../utilities_lib.hpp"
// ----------------------------------------------------------------------
//  Function declaration
//
int matrix(vector<vector<double>>&, vector<double>&);
// ----------------------------------------------------------------------
int main()
{
    vector<vector<double>> ep;
    vector<vector<double>> ept;
    vector<vector<double>> pt1;
    vector<vector<double>> pt2;
    vector<vector<double>> co;
    vector<double> th;
    // --------------------------------------
    double alpha = deg2rad(10);
    bool smooth_distribution = false;
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
    //
    const int m = n - 1; // get number of panels

    vector<vector<double>> a(m, vector<double>(n));
    vector<vector<double>> b(m, vector<double>(m));
    vector<double> g(m);
    vector<double> dl(m);
    vector<double> vel(m);

    //
    // Convert paneling to clockwise
    //
    for (int i = 0; i < n; i++)
    {
        ep.push_back({ ept[n - i - 1][0], ept[n - i - 1][1] });
    }
    // Establish coordinates of panel end points
    //
    for (int i = 0; i < m; i++)
    {
        pt1.push_back({ ep[i][0], ep[i][1] });
        pt2.push_back({ ep[i + 1][0], ep[i + 1][1] });
        double dz = pt2[i][1] - pt1[i][1];
        double dx = pt2[i][0] - pt1[i][0];
        th.push_back(atan2(dz, dx));
        // Collocation points
        double co_x = (pt2[i][0] - pt1[i][0]) / 2 + pt1[i][0];
        double co_z = (pt2[i][1] - pt1[i][1]) / 2 + pt1[i][1];
        co.push_back({ co_x, co_z });
    }

    // Establish Influence Coefficients
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < m; j++)
        {
            double xt = co[i][0] - pt1[j][0];
            double zt = co[i][1] - pt1[j][1];
            double x2t = pt2[j][0] - pt1[j][0];
            double z2t = pt2[j][1] - pt1[j][1];
            double x = xt * cos(th[j]) + zt * sin(th[j]);
            double z = -xt * sin(th[j]) + zt * cos(th[j]);
            double x2 = x2t * cos(th[j]) + z2t * sin(th[j]);
            double z2 = 0;

            //  Save Panel lengths for later use
            if (i == 0) dl.at(j) = x2;

            double r1 = sqrt(x * x + z * z);
            double r2 = sqrt((x - x2) * (x - x2) + z * z);
            double th1 = atan2(z, x);
            double th2 = atan2(z, x - x2);

            double ul = 0.15916 * (th2 - th1);
            double wl = 0.15916 * log(r2 / r1);
            if (i == j) {
                ul = 0.5;
                wl = 0;
            }
            double u = ul * cos(th[j]) - wl * sin(th[j]);
            double w = ul * sin(th[j]) + wl * cos(th[j]);

            a[i][j] = -u * sin(th[i]) + w * cos(th[i]);
            b[i][j] = u * cos(th[i]) + w * sin(th[i]);
        }

        a.at(i).at(n - 1) = cos(alpha) * sin(th.at(i)) - sin(alpha) * cos(th.at(i));
    }
    //  Kutta condition

    // for (int j = i;j < m + 1;j++)
    // {
    // }
    a.at(m / 4 - 1).at(m) = 0;
    a.at(m / 4 - 1).at(0) = 1;
    a.at(m / 4 - 1).at(m - 1) = 1;
    // Solve with Gauss Elimination Method
    // 
    matrix(a, g);

    ofstream out_file("CPV_cpp.DAT");
    double cl = 0;
    for (int i = 0;i < m; i++) {
        double temp = 0;
        for (int j = 0;j < m; j++) {
            temp += b.at(i).at(j) * g.at(j);
        }
        vel.at(i) = temp + cos(alpha) * cos(th.at(i)) + sin(alpha) * sin(th.at(i));
        cl += vel.at(i) * dl.at(i);
    }

    for (int i = 1; i < m; i++) {
        if (smooth_distribution == true) {
            double cp = 1 - pow(0.5 * (vel.at(i) + vel.at(i - 1)), 2);
            out_file << pt2.at(i - 1).at(0) << " ," << cp << endl;
        }
        else {
            double cp = 1 - vel.at(i) * vel.at(i);
            out_file << co.at(i).at(0) << " ," << cp << endl;

        }
    }
    out_file << endl << "LIFT COEFFICIENT = " << cl << endl;
    return 0;
}

int matrix(vector<vector<double>>& a, vector<double>& g)
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
                    double tmp = a[i][j];
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
        double p = a[i][i];
        for (int j = i; j < n; j++)
            a[i][j] /= p;
        for (int k = i + 1; k < n - 1; k++)
        {
            double p2 = a[k][i];
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
