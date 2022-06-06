using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Timers;
using System.Diagnostics;
using System.Threading;

namespace GFEC
{
    public static class LinearAlgebraTests
    {
        public static void SolveExample()
        {
            Stopwatch watch = new Stopwatch();
            Stopwatch watch2 = new Stopwatch();
            double[,] matrixA = new double[16000, 16000];
            double[,] matrixB = new double[16000, 16000];
            watch.Start();
            double[,] matrixC = MatrixOperations.MatrixAddition(matrixA, matrixB);
            watch.Stop();
            long time = watch.ElapsedMilliseconds;

            Thread.Sleep(1000);
            watch2.Start();
            double[,] matrixD = MatrixOperations.MatrixAdditionParallel2(matrixA, matrixB);
            watch2.Stop();
            long time2 = watch2.ElapsedMilliseconds;
        }

    }
}
