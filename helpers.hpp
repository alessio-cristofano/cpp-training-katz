#pragma once
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>

using std::string;
using std::vector;
// -----------------------------------------------
// CONSTANTS
// -----------------------------------------------
const double PI = 3.14159265358979323846;

// -----------------------------------------------
// FUNCTIONS
// -----------------------------------------------
float deg2rad(float angle_deg)
{
    return angle_deg / 57.2958;
}

float rad2deg(float angle_rad)
{
    return angle_rad * 57.2958;
}

void write_to_file(std::ofstream &filename, const vector<float> &var_array, string separator = ",")
{
    for (int i = 0; i < var_array.size() - 1; i++)
    {
        filename << var_array[i] << separator;
    }
    filename << var_array[var_array.size() - 1] << "\n";
}

vector<string> split_string(string word, string delimiter)
{
    int index = 0;
    vector<string> splitted;
    string rest = word;
    while (index != string::npos)
    {
        index = word.find(delimiter);
        if (index == string::npos)
        {
            rest.erase(std::remove(rest.begin(), rest.end(), '\n'), rest.cend());
            splitted.push_back(rest);
            continue;
        }
        string sep = word.substr(0, index);
        rest = word.substr(index + delimiter.size(), word.size());
        splitted.push_back(sep);
        word = rest;
    }
    return splitted;
}
