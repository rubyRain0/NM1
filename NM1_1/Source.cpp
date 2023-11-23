﻿#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <random>

const double E = 0.0000001;

double generateRandomNumber(double range_min, double range_max) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(range_min, range_max);
    return dis(gen);
}

double roundError(double error) {
    // Находим порядок погрешности
    int power = std::floor(std::log10(std::abs(error)));
    // Округляем погрешность до 3 значащих цифр согласно условию
    double roundedError = std::round(error / std::pow(10, power - 2)) * std::pow(10, power - 2);
    return roundedError;
}

void printArrays(int N, long double** a, long double** b, long double** c, long double* p, long double* q, long double* f)
{
    std::cout << "p: ";
    for (int i = 0; i < N; ++i) {
        std::cout << p[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "q: ";
    for (int i = 0; i < N; ++i) {
        std::cout << q[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "b: ";
    for (int i = 0; i < N; ++i) {
        std::cout << *b[i] << " ";
    }
    std::cout << std::endl;

    /*std::cout << "b_t: ";
    for (int i = 0; i < N; ++i) {
        std::cout << b_t[i] << " ";
    }
    std::cout << std::endl;*/

    std::cout << "c: ";
    for (int i = 0; i < N - 1; ++i) {
        std::cout << *c[i] << " ";
    }
    std::cout << std::endl;

    /*std::cout << "c_t: ";
    for (int i = 0; i < N - 1; ++i) {
        std::cout << c_t[i] << " ";
    }
    std::cout << std::endl;*/

    std::cout << "a: ";
    for (int i = 0; i < N - 1; ++i) {
        std::cout << *a[i] << " ";
    }
    std::cout << std::endl;

    /*std::cout << "a_t: ";
    for (int i = 0; i < N - 1; ++i) {
        std::cout << a_t[i] << " ";
    }
    std::cout << std::endl;*/

    std::cout << "f: ";
    for (int i = 0; i < N; ++i) {
        std::cout << f[i] << " ";
    }
    std::cout << std::endl;
}

void inputMatrixFromFile(int N, long double** a, 
    long double* a_t,
    long double** b, 
    long double* b_t,
    long double** c,
    long double* c_t,
    long double* f,
    long double* ft,
    long double* p, 
    long double* q)
{
    std::ifstream file("data.txt");
    if (!file) {
        std::cout << "Error opening the file." << std::endl;
        return;
    }
    long double el;
    //Not Overlapping Elements:
    for (size_t i = 0; i < N; ++i)
    {
        for (size_t j = 0; j < N; ++j)
        {   
            file >> el;
            if (el != 0 || j==0 || j == N-1 || i+j == N-2)
            {   
                if (j == 0)
                    p[i] = el;
                else if (j == N - 1)
                    q[i] = el;
                else if (i + j == N - 1)
                    b_t[i] = el;
                else if (i + j == N - 2)
                    c_t[i] = el;
                else if (i + j == N)
                {
                    if (i >= 1)
                    {
                        a_t[i - 1] = el;
                    }
                    else
                        a_t[i] = el;
                }
            }
        }
    }
    //Overlapping Elements:
    b[0] = &q[0];
    a[0] = &q[1];
    

    b[N - 1] = &p[N - 1];
    
    c[N-2] = &p[N - 2];
    for (size_t i = 0; i < N; i++)
    {
        file >> f[i];
    }

    //from _t vectors to original
    for (size_t i = 1; i < N-1; i++)
    {
        b[i] = &b_t[i];
    }
    for (size_t i = 0; i < N-2; i++)
    {
        c[i] = &c_t[i];
    }
    for (size_t i = 1; i < N-1; i++)
    {
        a[i] = &a_t[i];
    }
}

void WriteMatrixToFile(int N, long double** a, long double** b, long double** c, long double* p, long double* q, long double* f, std::string filename)
{
    std::ofstream ofile(filename);
    for (size_t i = 0; i < N; ++i)
    {
        for (size_t j = 0; j < N; ++j)
        {
            if (j == 0)
                ofile << p[i] << "\t\t";
            else if (i + j == N - 1)
                ofile << *b[i] << "\t\t";
            else if (i + j == N - 2)
            {
                ofile << *c[i] << "\t\t";
            }
            else if (i + j == N)
                ofile << *a[i-1] << "\t\t";
            else if (j == N - 1)
                ofile << q[i];
            else
                ofile << 0 << "\t\t";

        }
        ofile << "\n";
    }
    ofile << "\n\n";
    for (size_t i = 0; i < N; i++)
    {
        ofile << f[i] << '\n';
    }
}

bool isCorrectToBeSolved(int N, long double* f, long double* p)
{
    if ((f[N - 1] > E || f[N - 1] < -E) && p[N - 1] == 0)
        return false;
    return true;
}

void solution(int N, long double** a, long double** b, long double** c, long double* p, long double* q, long double* f, long double* ft, std::string filename)
{   
    //Step1: Clear bottom diag.
    for (size_t i = 1; i < N-2; ++i)
    {
        p[i] /= *b[i];
        *c[i] /= *b[i];
        q[i] /= *b[i];
        f[i] /= *b[i];
        *b[i] = 1;

        p[i + 1] += p[i] * -*a[i];
        *b[i + 1] += *c[i] * -*a[i];
        q[i + 1] += q[i] * -*a[i];
        f[i + 1] += f[i] * -*a[i];

        *a[i] = 0; //?    
    }

    //N-2 before-last row

    p[N - 2] /= *b[N - 2];
    q[N - 2] /= *b[N - 2];
    f[N - 2] /= *b[N - 2];
    *b[N - 2] = 1;

    p[N - 1] += p[N - 2] * -*a[N-2];
    q[N - 1] += q[N - 2] * -*a[N-2];
    f[N - 1] += f[N - 2] * -*a[N-2];
    *a[N - 2] = 0;
    
    q[N - 1] /= *b[N - 1];
    f[N - 1] /= *b[N - 1];
    p[N - 1] /= *b[N - 1];

    //Step2: (Clear first col)
    for (int i = N-2; i >= 0; i--)
    {   
        if (p[i]!= 0)
        {
            q[i] += q[N - 1] * -p[i];
            f[i] += f[N - 1] * -p[i];
            p[i] += p[N - 1] * -p[i];
         
        }
    }

    //Step3: (Clear above diag.)
    for (size_t i = N-2; i >= 1; --i)
    {
       q[i - 1] += q[i] * -*c[i - 1];
       f[i - 1] += f[i] * -*c[i - 1];
       *c[i - 1] += *b[i] * -*c[i - 1];
    }
    f[0] /= *b[0];
    q[0] = 1;

    //Step4: (Clear last col)
    for (size_t i = 0; i < N-1; ++i)
    {
        if (q[i + 1] != 0)
        {
            f[i + 1] += f[0] * -q[i + 1];
            q[i + 1] += q[0] * -q[i + 1];
        }
    }
    WriteMatrixToFile(N, a, b, c, p, q, f, filename);
    //printArrays(N, a, b, c, p, q, f);
    long double* x = new long double[N];
    if (isCorrectToBeSolved)
    {
        std::cout << "Solution vector: \n";
        for (size_t i = 0; i < N; ++i)
        {
            x[i] = f[i];
            std::cout << x[i] << std::endl;
        }
    }

}

void generateMatrixTestAccuracy(int N, int minValue, int maxValue)
{
    long double** a = new long double*[N - 1]; //under
    long double* a_t = new long double[N - 1]; //under

    long double** b = new long double*[N]; //mid
    long double* b_t = new long double[N - 1]; //mid

    long double** c = new long double*[N - 1]; //upper
    long double* c_t = new long double[N - 1]; //upper

    long double* f = new long double[N];
    long double* ft = new long double[N];
    long double* p = new long double[N];
    long double* q = new long double[N];

    for (size_t i = 0; i < N; ++i)
    {
        *b[i] = generateRandomNumber(minValue, maxValue);
        p[i] = generateRandomNumber(minValue, maxValue);
        q[i] = generateRandomNumber(minValue, maxValue);
        f[i] = generateRandomNumber(minValue, maxValue);

    }
    for (size_t i = 0; i < N-1; ++i)
    {
        *a[i] = generateRandomNumber(minValue, maxValue);
        *c[i] = generateRandomNumber(minValue, maxValue);
    }
    WriteMatrixToFile(N, a, b, c, p, q, f, "randomSystemTest.txt");
    solution(N, a, b, c, p, q, f, ft, "randomSystemResult.txt");
}


int main() {
    std::ifstream file("data.txt");
    if (!file) {
        std::cout << "Error opening the file." << std::endl;
        return 1;
    }
    int N = 10;
    long double** a = new long double*[N-1]; //under
    long double* a_t = new long double[N - 1]; //under

    long double** b = new long double*[N]; //mid
    long double* b_t = new long double[N - 1]; //mid

    long double** c = new long double*[N-1]; //upper
    long double* c_t = new long double[N - 1]; //upper

    long double* f = new long double[N];
    long double* ft = new long double[N];
    long double* p = new long double[N];
    long double* q = new long double[N];
  
    inputMatrixFromFile(N, a, a_t, b, b_t, c, c_t, f, ft, p, q);
    printArrays(N, a, b, c, p, q, f);
    solution(N, a, b, c, p, q, f, ft, "resultSystem.txt");
    //generateMatrixTestAccuracy(100, -10, 10);

    return 0;
}