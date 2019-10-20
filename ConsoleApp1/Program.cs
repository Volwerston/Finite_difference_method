using System;
using System.Collections.Generic;
using Accord.Math;

namespace ConsoleApp1
{
    class Program
    {
        private const int A = 0;

        private const int B = 1;

        private const int N = 10;

        private static double h = (B - A) / (double) N;

        private const int l2 = 4;

        private static double CalculateDefiniteIntegral(double[] y1, double[] y2)
        {
            double integralSum = 0;

            for (int i = 0; i < N; ++i)
            {
                integralSum += h * y1[i] * y2[i];
            }

            return integralSum;
        }

        private static double[,] BuildMatrix()
        {
            double[,] matrix = new double[N+1, N+1];

            matrix[0, 0] = 1;
            for (int i = 1; i <= N; ++i)
            {
                matrix[0, N] = 0;
            }

            matrix[N, N - 1] = -1;
            matrix[N, N] = 1;
            for (int i = 0; i < N - 1; ++i)
            {
                matrix[N, i] = 0;
            }

            for (int i = 1; i <= N-1; ++i)
            {
                for (int j = 0; j <= N; ++j)
                {
                    if (j == i - 1)
                    {
                        matrix[i, j] = -(1/Math.Pow(h, 2))*Math.Pow(A + i*h + 1, 2);
                    }
                    else if (j == i)
                    {
                        matrix[i, j] = (2 / Math.Pow(h, 2)) * Math.Pow(A + i * h + 1, 2);
                    }
                    else if (j == i + 1)
                    {
                        matrix[i, j] = -(1 / Math.Pow(h, 2)) * Math.Pow(A + i * h + 1, 2);
                    }
                    else
                    {
                        matrix[i, j] = 0;
                    }
                }
            }

            return matrix;
        }

        static void Main()
        {
            Console.WriteLine("Enter number of iterations: ");
            int iterations = Convert.ToInt32(Console.ReadLine());

            // Setting initial approximation
            Func<double, double> f0 = x => x - Math.Pow(x, 2)/2.0;
            double[] f0Val = new double[N + 1];
            for (int i = 1; i < N; ++i)
            {
                f0Val[i] = f0(A + i * h);
            }

            f0Val[0] = 0;
            f0Val[N] = h;

            double a0 = CalculateDefiniteIntegral(f0Val, f0Val);
            var matrix = BuildMatrix();
            var inversedMatrix = matrix.Inverse();

            double[] f1Val;
            double[] f2Val;

            double a_k_minus_1 = a0;
            double a_k;
            double a_k_plus_1;

            LinkedList<double> lowerBounds = new LinkedList<double>();
            LinkedList<double> upperBounds = new LinkedList<double>();
            LinkedList<double> mus = new LinkedList<double>();

            for (int k = 1; k <= iterations; ++k)
            {
                f1Val = inversedMatrix.Dot(f0Val);

                f2Val = inversedMatrix.Dot(f1Val);

                a_k = CalculateDefiniteIntegral(f0Val, f1Val);
                a_k_plus_1 = CalculateDefiniteIntegral(f1Val, f2Val);

                double mu_k = a_k_minus_1 / a_k;
                double mu_k_plus_1 = a_k / a_k_plus_1;
                mus.AddLast(mu_k_plus_1);

                upperBounds.AddFirst(mu_k_plus_1);

                if (mu_k_plus_1 < l2)
                {
                    double lowerBound = mu_k_plus_1 - ((mu_k - mu_k_plus_1) * mu_k_plus_1) / (l2 - mu_k_plus_1);
                    lowerBounds.AddLast(lowerBound);
                }

                f0Val = f1Val;
                a_k_minus_1 = a_k;
            }

            // Print results: 

            Console.WriteLine("Lower bounds");
            foreach (var elem in lowerBounds)
            {
                Console.WriteLine(elem);
            }

            Console.WriteLine();
            Console.WriteLine("Upper bounds");
            foreach (var elem in upperBounds)
            {
                Console.WriteLine(elem);
            }

            // Mu-s
            //foreach (var elem in mus)
            //{
            //    Console.WriteLine(elem);
            //}

            Console.ReadLine();
        }
    }
}
