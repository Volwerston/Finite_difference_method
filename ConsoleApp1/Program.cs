using System;
using Accord.Math;

namespace ConsoleApp1
{
    class Program
    {
        private const int A = 0;
        private const int B = 1;

        private const int N = 5;
        private static double h = (B - A) / (double)N;

        private const int l2 = 4;

        private static double CalculateDefiniteIntegral(double[] y1, double[] y2)
        {
            double integralSum = 0;

            for (int i = 0; i < N; ++i)
            {
                integralSum += h * y1[i] * y2[i] / Math.Pow(A + i * h + 1, 2);
            }

            return integralSum;
        }

        private static double[,] BuildMatrix()
        {
            double[,] matrix = new double[N + 1, N + 1];

            matrix[0, 0] = h;
            matrix[N, N - 1] = -1;
            matrix[N, N] = 1;

            for (int i = 1; i <= N - 1; ++i)
            {
                matrix[i, i - 1] = -1 / Math.Pow(h, 2) * Math.Pow(A + i * h + 1 , 2);
                matrix[i, i] = 2 / Math.Pow(h, 2) * Math.Pow(A + i * h + 1, 2);
                matrix[i, i + 1] = -1 / Math.Pow(h, 2) * Math.Pow(A + i * h + 1, 2);
            }

            return matrix;
        }

        public static void RunKollatz(int iterations)
        {
            Console.WriteLine($"-----Running Kollatz-Temple ({iterations} iterations)-----");

            Func<double, double> f0 = x => x - Math.Pow(x, 2) / 2.0;
            double[] f0Val = new double[N + 1];
            for (int i = 1; i < N; ++i)
            {
                f0Val[i] = f0(A + i * h);
            }
            f0Val[N] = h;


            var matrix = BuildMatrix();

            for (int k = 1; k <= iterations; ++k)
            {
                double[] f1Val = Matrix.Solve(matrix, f0Val);

                double a_k_minus_1 = CalculateDefiniteIntegral(f0Val, f0Val);
                double a_k = CalculateDefiniteIntegral(f0Val, f1Val);
                double a_k_plus_1 = CalculateDefiniteIntegral(f1Val, f1Val);

                double mu_k = a_k_minus_1 / a_k;
                double mu_k_plus_1 = a_k / a_k_plus_1;

                if (mu_k_plus_1 < l2)
                {
                    double lowerBound = mu_k_plus_1 - ((mu_k - mu_k_plus_1) * mu_k_plus_1) / (l2 - mu_k_plus_1);
                    Console.WriteLine("Lower bound: " + lowerBound);
                }

                Console.WriteLine("Upper bound: " + mu_k_plus_1);
                Console.WriteLine("==========");

                f0Val = f1Val;
            }

            Console.WriteLine("\n\n");
        }

        public static void RunKrylov(int iterations)
        {
            Console.WriteLine($"-----Running Krylov-Bogoliubov ({iterations} iterations)-----");

            Func<double, double> f0 = x => x - Math.Pow(x, 2) / 2.0;
            double[] f0Val = new double[N + 1];
            for (int i = 1; i < N; ++i)
            {
                f0Val[i] = f0(A + i * h);
            }
            f0Val[N] = h;

           
            var matrix = BuildMatrix();

            for (int k = 1; k <= iterations; ++k)
            {
                double[] f1Val = Matrix.Solve(matrix, f0Val);

                double a_k_minus_1 = CalculateDefiniteIntegral(f0Val, f0Val);
                double a_k = CalculateDefiniteIntegral(f0Val, f1Val);
                double a_k_plus_1 = CalculateDefiniteIntegral(f1Val, f1Val);

                double mu_k = a_k_minus_1 / a_k;
                double mu_k_plus_1 = a_k / a_k_plus_1;

                Console.WriteLine("Lower bound: " + (mu_k_plus_1 - Math.Sqrt(Math.Abs(mu_k - mu_k_plus_1) * mu_k_plus_1)));
                Console.WriteLine("Upper bound: " + (mu_k_plus_1 + Math.Sqrt(Math.Abs(mu_k - mu_k_plus_1) * mu_k_plus_1)));
                Console.WriteLine("==========");

                f0Val = f1Val;
            }
        }

        static void Main()
        {
            Console.WriteLine("Enter number of iterations: ");
            int iterations = Convert.ToInt32(Console.ReadLine());
            
            RunKollatz(iterations);
            RunKrylov(iterations);

            Console.ReadLine();
        }
    }
}
