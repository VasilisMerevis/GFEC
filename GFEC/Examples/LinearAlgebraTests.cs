using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Timers;
using System.Diagnostics;
using System.Threading;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;

namespace GFEC
{
    public static class LinearAlgebraTests
    {
        public static void SolveExample()
        {
            Stopwatch watch = new Stopwatch();
            Stopwatch watch2 = new Stopwatch();
            Stopwatch watch3 = new Stopwatch();
            Stopwatch watchS1 = new Stopwatch();
            double[,] matrixA = new double[16000, 16000];
            double[] vectorA = VectorOperations.CreateRandomVector(16000);
            MatrixOperations.FillMatrixWithDoubleValue(matrixA, 1.0);
            
            
 
            watch.Start();
            double[,] matrixC = MatrixOperations.MatrixAddition(matrixA, matrixA);
            watch.Stop();
            long time = watch.ElapsedMilliseconds;

            Thread.Sleep(1000);
            watch2.Start();
            double[,] matrixD = MatrixOperations.MatrixAdditionParallel2(matrixA, matrixA);
            watch2.Stop();
            long time2 = watch2.ElapsedMilliseconds;

            //watchS1.Start();
            //LUFactorization solver = new LUFactorization();
            //solver.Solve(matrixA, vectorA);
            //watchS1.Stop();
           

            Matrix<double> matrixE = Matrix<double>.Build.DenseOfArray(matrixA);
            Matrix<double> matrixF = Matrix<double>.Build.DenseOfArray(matrixA);
            Vector<double> vectorB = Vector<double>.Build.DenseOfArray(vectorA);
            watch3.Start();
            Matrix<double> matrixG = matrixE.Add(matrixF);
            watch3.Stop();
            long time3 = watch3.ElapsedMilliseconds;

            var x = matrixE.Solve(vectorB);
        }

    }
}
