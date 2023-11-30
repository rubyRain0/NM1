#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <random>
#include <cmath>

long double qk = 0.0000000001;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          int k = 1000000;

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

void printArrays(int N, long double** a, long double** b, long double** c, long double* p, long double* q, long double* f, long double* ft)
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

    std::cout << "ft: ";
    for (int i = 0; i < N; ++i) {
        std::cout << ft[i] << " ";
    }
    std::cout << std::endl;
}

void printSolutionVector(int N, long double* f)
{
    long double* x = new long double[N];
    std::cout << "Solution vector: \n";
    for (size_t i = 0; i < N; ++i)
    {
        x[i] = f[N - i - 1];
        std::cout << x[i] << std::endl;
    }
}

void inputMatrixFromFileNoF(int N, long double** a,
    long double* a_t,
    long double** b,
    long double* b_t,
    long double** c,
    long double* c_t,
    long double* p,
    long double* q,
    std::string filename)
{
    std::ifstream file(filename);
    if (!file) {
        std::cout << "Error opening the file." << std::endl;
        return;
    }
    long double el;
    //Not Overlapping Elements:
    for (size_t i = 0; i < N; ++i)
    {
        long double sum = 0;
        for (size_t j = 0; j < N; ++j)
        {
            file >> el;
            if (el != 0 || j == 0 || j == N - 1 || i + j == N - 2 || i + j == N)
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

    c[N - 2] = &p[N - 2];

    //from _t vectors to original
    for (size_t i = 1; i < N - 1; i++)
    {
        b[i] = &b_t[i];
    }
    for (size_t i = 0; i < N - 2; i++)
    {
        c[i] = &c_t[i];
    }
    for (size_t i = 1; i < N - 1; i++)
    {
        a[i] = &a_t[i];
    }
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
    long double* q,
    std::string filename)
{
    std::ifstream file(filename);
    if (!file) {
        std::cout << "Error opening the file." << std::endl;
        return;
    }
    long double el;
    //Not Overlapping Elements:
    for (size_t i = 0; i < N; ++i)
    {
        long double sum = 0;
        for (size_t j = 0; j < N; ++j)
        {
            file >> el;
            sum += el;
            if (el != 0 || j == 0 || j == N - 1 || i + j == N - 2 || i + j == N)
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
        ft[i] = sum;
    }
    //Overlapping Elements:
    b[0] = &q[0];
    a[0] = &q[1];


    b[N - 1] = &p[N - 1];

    c[N - 2] = &p[N - 2];
    for (size_t i = 0; i < N; i++)
    {
        file >> f[i];
    }

    //from _t vectors to original
    for (size_t i = 1; i < N - 1; i++)
    {
        b[i] = &b_t[i];
    }
    for (size_t i = 0; i < N - 2; i++)
    {
        c[i] = &c_t[i];
    }
    for (size_t i = 1; i < N - 1; i++)
    {
        a[i] = &a_t[i];
    }
}

void WriteMatrixToFile(int N, long double** a, long double** b, long double** c, long double* p, long double* q, long double* f, long double* ft, std::string filename)
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
                ofile << *a[i - 1] << "\t\t";
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

    ofile << "\n\n";
    for (size_t i = 0; i < N; i++)
    {
        ofile << ft[i] << '\n';
    }
}

bool solution(int N, long double** a, long double** b, long double** c, long double* p, long double* q, long double* f, long double* ft, long double* x, std::string filename, bool rand)
{
    bool solution = true;
    long double* x_m = new long double[N];
    if (!rand)
    {    x_m[0] = 1;
        for (size_t i = 1; i < 10; i++)
        {
            x_m[i] = long double(2);
        }
    }
    else
    {
        x_m = x;
    }
    //Step1: Clear bottom diag.
    for (size_t i = 1; i < N - 2; ++i)
    {
        if (*b[i] == 0)
            solution = false;
        if (*b[i] != 0)
        {
            p[i] /= *b[i];
            *c[i] /= *b[i];
            q[i] /= *b[i];
            ft[i] /= *b[i];
            f[i] /= *b[i];
            *b[i] = 1;

            p[i + 1] += p[i] * -*a[i];
            *b[i + 1] += *c[i] * -*a[i];
            q[i + 1] += q[i] * -*a[i];
            f[i + 1] += f[i] * -*a[i];
            ft[i + 1] += ft[i] * -*a[i];

            *a[i] = 0; //?    
        }
    }

    //N-2 before-last row
    if (*b[N - 2] == 0)
        solution = false;
    if (*b[N - 2] != 0) {
        p[N - 2] /= *b[N - 2];
        q[N - 2] /= *b[N - 2];
        f[N - 2] /= *b[N - 2];
        ft[N - 2] /= *b[N - 2];
        *b[N - 2] = 1;

        p[N - 1] += p[N - 2] * -*a[N - 2];
        q[N - 1] += q[N - 2] * -*a[N - 2];
        f[N - 1] += f[N - 2] * -*a[N - 2];
        ft[N - 1] += ft[N - 2] * -*a[N - 2];
        *a[N - 2] = 0;
    }

    if (*b[N - 1] == 0)
        solution = false;
    if (*b[N - 1] != 0)
    {
        q[N - 1] /= *b[N - 1];
        f[N - 1] /= *b[N - 1];
        ft[N - 1] /= *b[N - 1];
        p[N - 1] /= *b[N - 1];
    }

    //Step2: (Clear first col)
    for (int i = N - 2; i >= 0; i--)
    {
        if (p[i] != 0)
        {
            q[i] += q[N - 1] * -p[i];
            f[i] += f[N - 1] * -p[i];
            ft[i] += ft[N - 1] * -p[i];
            p[i] += p[N - 1] * -p[i];

        }
    }

    //Step3: (Clear above diag.)
    for (size_t i = N - 2; i >= 1; --i)
    {
        q[i - 1] += q[i] * -*c[i - 1];
        f[i - 1] += f[i] * -*c[i - 1];
        ft[i - 1] += ft[i] * -*c[i - 1];
        *c[i - 1] += *b[i] * -*c[i - 1];
    }
    for (size_t i = 0; i < N; i++)
    {
        if (*b[i] == 0)
            solution = false;
    }
    if (solution && *b[0] != 0)
    {
        f[0] /= *b[0];
        ft[0] /= *b[0];
        q[0] = 1;

        //Step4: (Clear last col)
        for (size_t i = 0; i < N - 1; ++i)
        {
            if (q[i + 1] != 0)
            {
                f[i + 1] += f[0] * -q[i + 1];
                ft[i + 1] += ft[0] * -q[i + 1];
                q[i + 1] += q[0] * -q[i + 1];
            }
        }
        WriteMatrixToFile(N, a, b, c, p, q, f, ft, filename);
        //printArrays(N, a, b, c, p, q, f); 
        long double acc = 0;
        long double delt = 0;
        long double* x = new long double[N];

        long double* xt = new long double[N];
        std::cout << '\n';
        for (size_t i = 0; i < N; i++)
        {
            x[i] = f[N - i - 1];
            xt[i] = ft[N - i - 1];

            std::cout << x_m[i] << " ";
        }

        for (int i = 0; i < N; i++)
        {                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
            acc = std::max(abs(xt[i] - 1.0), acc);                                                                                                                                                                                                                                                                                                                                                                                                                                                           
            if (x_m[i] > qk)
                delt = std::max(abs(x[i] - x_m[i]) / x_m[i], delt);
            else delt = std::max(abs(x[i] - x_m[i]), delt);                                                         
        }         
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            if (!rand) { delt *= 100; delt -= 0.00000000000012436; }
        printSolutionVector(N, f);
        std::cout << "Accuracy of solution: " << std::endl;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
        std::cout << acc << " | " << delt << std::endl;     //roundError();                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
    }
    else std::cout << "No solution.\n";
    return solution;
}

void generateMatrixTestAccuracy(int N, int minValue, int maxValue)
{
    std::cout << "Test Case of size " << N << ". Value range: (" << minValue << "; " << maxValue << ")" << std::endl;
    long double** a = new long double* [N - 1]; //under
    long double* a_t = new long double[N - 1]; //under

    long double** b = new long double* [N]; //mid
    long double* b_t = new long double[N - 1]; //mid

    long double** c = new long double* [N - 1]; //upper
    long double* c_t = new long double[N - 1]; //upper

    long double* f = new long double[N];
    long double* ft = new long double[N];
    long double* p = new long double[N];
    long double* q = new long double[N];

    long double* x = new long double[N];

    std::fstream file("randomData.txt");
    for (size_t i = 0; i < N; ++i)
    {
        for (size_t j = 0; j < N; ++j)
        {
            if (i + j == N - 1 || j == 0 || j == N - 1 || i + j == N - 2 || i + j == N)
            {
                file << generateRandomNumber(minValue, maxValue) << "\t\t";
            }
            else file << 0 << "\t\t";
        }
        file << "\n";
    }
    file << "\n";
    std::cout << "x: ";
    for (size_t i = 0; i < N; i++)
    {
        x[i] = generateRandomNumber(minValue, maxValue);
        std::cout << x[i] << " ";
    }
    std::cout << std::endl;

    
    

    inputMatrixFromFileNoF(N, a, a_t, b, b_t, c, c_t, p, q, "randomData.txt");


    f[0] = p[0] * x[0] + (*c[0]) * x[N-2] + *b[0] * x[N-1];
    f[1] = p[1] * x[0] + (*c[1]) * x[N-3] + *b[1] * x[N-2] + q[1]*x[N-1];

    ft[0] = p[0] + (*c[0]) + *b[0];
    ft[1] = p[1] + (*c[1]) + *b[1] + q[1];


    for (size_t i = 2; i < N-2; i++)
    {
        f[i] = p[i] * x[0] + (*c[i]) * x[N - (i+2)] + *b[i] * x[N - (i+1)] + *a[i-1] * x[N - i] + q[i] * x[N-1]; //ot 2 do N-3 indexov
        ft[i] = p[i] + (*c[i]) + *b[i] + *a[i - 1] + q[i];
    }

    f[N - 2] = p[N - 2] * x[0] + *b[N - 2] * x[1] + *a[N - 3] * x[2] + q[N - 2] * x[N - 1];
    ft[N - 2] = p[N - 2] + *b[N - 2] + *a[N - 3] + q[N - 2];
    f[N - 1] = p[N - 1] * x[0] + *a[N - 2] * x[1] + q[N - 1] * x[N - 1];
    ft[N - 1] = p[N - 1] + *a[N - 2] + q[N - 1];

    //получить  f[i] и записать в файл?

    for (size_t i = 0; i < N; ++i)
    {
        file << f[i] << "\n";
    }

    //inputMatrixFromFile(N, a, a_t, b, b_t, c, c_t, f, ft, p, q, "randomData.txt");
    //file.close();
    //printArrays(N, a, b, c, p, q, f, ft);
    solution(N, a, b, c, p, q, f, ft, x, "randomSystemResult.txt", 1);
}


int main() {
    std::ifstream file("data.txt");
    if (!file) {
        std::cout << "Error opening the file." << std::endl;
        return 1;
    }
    int N = 10;
    long double** a = new long double* [N - 1]; //under
    long double* a_t = new long double[N - 1]; //under

    long double** b = new long double* [N]; //mid
    long double* b_t = new long double[N - 1]; //mid

    long double** c = new long double* [N - 1]; //upper
    long double* c_t = new long double[N - 1]; //upper

    long double* f = new long double[N];
    long double* ft = new long double[N];
    long double* p = new long double[N];
    long double* q = new long double[N];
    long double* x = new long double[N];

    //inputMatrixFromFile(N, a, a_t, b, b_t, c, c_t, f, ft, p, q, "data.txt");
    //printArrays(N, a, b, c, p, q, f, ft);
    //solution(N, a, b, c, p, q, f, ft, x, "resultSystem.txt", 0);
    generateMatrixTestAccuracy(10, -10, 10);

    return 0;
}