#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "../helpers.hpp"

using std::atan;
using std::cin;
using std::cos;
using std::cout;
using std::pow;
using std::sin;
using std::sqrt;
// ----------------------------------------------------------------------
void read_input_e(float &);
void read_input_k(float &);
void read_input_aoa(float &);
void read_input_panels(int &);

int main()
{
    const int tl = 1;
    float e_coeff; // Thickness Coefficient
    float k_coeff; // Angle Coefficient
    float aoa;     // Angle of attack in degrees
    int panels;    // Number of panels of airfoil
    float x;
    float y;
    float cp;

    read_input_e(e_coeff);
    read_input_k(k_coeff);
    read_input_aoa(aoa);
    read_input_panels(panels);

    while (360 % panels != 0)
    {
        cout << "(NUMBER OF PANELS MUST BE AN EVEN FACTOR OF 360)\n"
             << "ENTER ANOTHER NUMBER";
        cin >> panels;
    }

    int i_theta = 360 / panels;
    float a = 2 * tl * pow(e_coeff + 1, k_coeff - 1) / pow(2, k_coeff);

    std::ofstream airfoil_file("AFOIL2.DAT");
    std::ofstream cp_file("CP.DAT");
    //
    for (int i = 0; i <= 360; i += i_theta)
    {
        if (i == 0 || i == 360)
        {
            x = 1;
            y = 0;
            cp = 1;
            write_to_file(airfoil_file, {x, y});
            if (k_coeff != 2)
            {
                write_to_file(cp_file, {x, cp});
            }
        }
        else
        {
            float th = deg2rad(i);
            float r1 = sqrt(pow(a * (cos(th) - 1), 2) + pow(a * sin(th), 2));
            float r2 = sqrt(pow(a * (cos(th) - e_coeff), 2) + pow(a * sin(th), 2));
            float th1;
            float th2;
            if (th == 0)
            {
                th1 = PI / 2;
            }
            else
            {
                th1 = (atan((a * sin(th)) / (a * (cos(th) - 1)))) + PI;
            }
            if (cos(th) - e_coeff < 0 && sin(th) != 0)
            {
                th2 = atan((a * sin(th)) / (a * (cos(th) - e_coeff))) + PI;
            }
            else if (cos(th) - e_coeff > 0 && sin(th) < 0)
            {
                th2 = (atan((a * sin(th)) / (a * (cos(th) - e_coeff)))) + 2 * PI;
            }
            else
            {
                th2 = (atan((a * sin(th)) / (a * (cos(th) - e_coeff))));
            }
            //
            float com1;
            com1 = (pow(r1, k_coeff) / pow(r2, (k_coeff - 1))) / (pow(cos((k_coeff - 1) * th2), 2) + pow(sin((k_coeff - 1) * th2), 2));
            x = com1 * (cos(k_coeff * th1) * cos((k_coeff - 1) * th2) + sin(k_coeff * th1) * sin((k_coeff - 1) * th2)) + tl;
            y = com1 * (sin(k_coeff * th1) * cos((k_coeff - 1) * th2) - cos(k_coeff * th1) * sin((k_coeff - 1) * th2));
            write_to_file(airfoil_file, {x, y});
            //
            float a1 = cos((k_coeff - 1) * th1) * cos(k_coeff * th2) + sin((k_coeff - 1) * th1) * sin(k_coeff * th2);
            float b1 = sin((k_coeff - 1) * th1) * cos(k_coeff * th2) - cos((k_coeff - 1) * th1) * sin(k_coeff * th2);
            float c1 = pow(cos(k_coeff * th2), 2) + pow(sin(k_coeff * th2), 2);
            float p = a * (1 - k_coeff + k_coeff * e_coeff);
            float d1 = a1 * (a * cos(th) - p) - b1 * a * sin(th);
            float d2 = a1 * a * sin(th) + b1 * (a * cos(th) - p);
            float al = deg2rad(aoa);
            float temp = 2 * c1 * (sin(al) - sin(al - th)) / (pow(d1, 2) + pow(d2, 2));
            float com2 = temp * pow(r2, k_coeff) / pow(r1, (k_coeff - 1));
            float vx = d1 * sin(th) + d2 * cos(th);
            float vy = -(d1 * cos(th) - d2 * sin(th));
            cp = 1 - pow(com2, 2) * (pow(vx, 2) + pow(vy, 2));
            write_to_file(cp_file, {x, cp});
        }
    }

    return 0;
}

void read_input_e(float &e)
{
    cout << "READY TO START VAN DE VOOREN TRANSFORMATION\n";
    cout << "ENTER THICKNESS COEFF. E\n";
    cin >> e;
}
void read_input_k(float &k)
{
    cout << "ENTER T.E. ANGLE COEFF. K\n";
    cin >> k;
}
void read_input_aoa(float &aoa)
{
    cout << "ENTER THE ANGLE OF ATTACK IN DEGREES\n";
    cin >> aoa;
}
void read_input_panels(int &panels)
{
    cout << "ENTER NUMBER OF AIRFOIL PANELS,\n";
    cout << "WITH WHICH TO MODEL THE AIRFOIL,\n";
    cout << "(NOTE THAT M SHOULD BE AN EVEN FACTOR OF 360)\n";
    cin >> panels;
};
