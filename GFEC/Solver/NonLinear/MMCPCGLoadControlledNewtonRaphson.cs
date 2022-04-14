using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;

namespace GFEC
{
    public class MMCPCGLoadControlledNewtonRaphson : NonLinearSolution
    {
        private double[] localSolutionVector;

        public MMCPCGLoadControlledNewtonRaphson()
        {

        }
        public MMCPCGLoadControlledNewtonRaphson(double[] exSolution)
        {
            localSolutionVector = exSolution;
        }
        private double[] MMCPCGLoadControlledNR(double[] forceVector)
        {
            ////timers
            //List<double> solutionTimerList = new List<double>();
            //List<double> SystemAssembleTimerList = new List<double>();
            ////
            //
            List<double> IterationsCounters = new List<double>();
            int countIterations = new int();
            //discretization.SeperateContactDoF();
            double[] incrementDf = VectorOperations.VectorScalarProductNew(forceVector, lambda);
            double[] solutionVector = localSolutionVector;
            double[] incrementalExternalForcesVector = new double[forceVector.Length];
            //double[] incrementalExternalForcesVectorModified = new double[forceVector.Length];
            double[] tempSolutionVector = new double[solutionVector.Length];
            double[] deltaU = new double[solutionVector.Length];
            double[] internalForcesTotalVector;
            double[] dU;
            double[] residual;
            double residualNorm;
            //double[] oldInternalForces;
            //double relativeInternalForcesNorm;
            int fillInDegree = 10;
            discretization.MMCPCGUpdateDisplacements(solutionVector);
            var stiffnessMatrices = discretization.MMCPCGCreateTotalStiffnessMatrix();
            double[,] K = stiffnessMatrices.Item2;
            int[,] fillInLevels = MatrixOperations.ICholLevels(K, fillInDegree);
            discretization.InitializeContactTangentialProperties();
            discretization.InitializeContactSurfaceVectors();
            for (int i = 0; i < numberOfLoadSteps; i++)
            {
                //////add timer
                ////var timer = System.Diagnostics.Stopwatch.StartNew();

                //
                //int count = 0;
                //int countIncrements = 1;
                //
                incrementalExternalForcesVector = VectorOperations.VectorVectorAddition(incrementalExternalForcesVector,
                    incrementDf);
                discretization.UpdateContactSurfaceVectors();
                internalForcesTotalVector = discretization.MMCPCGCreateTotalInternalForcesVectors();
                //stiffnessMatrices = discretization.MMCPCGCreateTotalStiffnessMatrix();
                //double[,] totalStifnessMatrix = stiffnessMatrices.Item1;
                //double[,] K = stiffnessMatrices.Item2;
                double[,] B = stiffnessMatrices.Item3;
                double[,] C = stiffnessMatrices.Item4;
                Tuple<double[], double[]> incrementForces = discretization.MMCPCGSeperateReducedExternalForcesVectors(incrementDf);
                double[] incrementDfM = incrementForces.Item1;
                double[] incrementDfC = incrementForces.Item2;

                ////timer.Stop();
                ////var AssembleElapsedMs = (double)timer.ElapsedMilliseconds;
                ////SystemAssembleTimerList.Add(AssembleElapsedMs);
                //////timer
                //////print matrices to files
                ////if(i == 0)
                ////{
                ////    MatrixOperations.PrintMatrixToFile2(K, @"C:\Users\Public\Documents\" + "K0.dat");
                ////    MatrixOperations.PrintMatrixToFile2(B, @"C:\Users\Public\Documents\" + "B0.dat");
                ////    MatrixOperations.PrintMatrixToFile2(C, @"C:\Users\Public\Documents\" + "C0.dat");
                ////}
                ////else if(i == 2)
                ////{
                ////    MatrixOperations.PrintMatrixToFile2(K, @"C:\Users\Public\Documents\" + "K2.dat");
                ////    MatrixOperations.PrintMatrixToFile2(B, @"C:\Users\Public\Documents\" + "B2.dat");
                ////    MatrixOperations.PrintMatrixToFile2(C, @"C:\Users\Public\Documents\" + "C2.dat");
                ////}
                ////else if(i == 20)
                ////{
                ////    MatrixOperations.PrintMatrixToFile2(K, @"C:\Users\Public\Documents\" + "K20.dat");
                ////    MatrixOperations.PrintMatrixToFile2(B, @"C:\Users\Public\Documents\" + "B20.dat");
                ////    MatrixOperations.PrintMatrixToFile2(C, @"C:\Users\Public\Documents\" + "C20.dat");
                ////}
                ////else if (i == 30)
                ////{
                ////    MatrixOperations.PrintMatrixToFile2(K, @"C:\Users\Public\Documents\" + "K30.dat");
                ////    MatrixOperations.PrintMatrixToFile2(B, @"C:\Users\Public\Documents\" + "B30.dat");
                ////    MatrixOperations.PrintMatrixToFile2(C, @"C:\Users\Public\Documents\" + "C30.dat");
                ////}
                //////print matrices to files
                //////add timer
                ////var watch = System.Diagnostics.Stopwatch.StartNew();
                ///

                var solverResults = linearSolver.Solve(K, B, C, incrementDfM, incrementDfC, fillInLevels, fillInDegree, true);
                dU = solverResults.Item1;
                countIterations = solverResults.Item2;
                IterationsCounters.Add((double)(countIterations));
                ////watch.Stop();
                ////var solutionElapsedMs = (double)watch.ElapsedMilliseconds;
                ////solutionTimerList.Add(solutionElapsedMs);
                //////timer
                discretization.MMCPCGUpdateElementsIncrementalDisplacements(dU);
                solutionVector = VectorOperations.VectorVectorAddition(solutionVector, dU);
                residual = VectorOperations.VectorVectorSubtraction(internalForcesTotalVector, incrementalExternalForcesVector);
                residualNorm = VectorOperations.VectorNorm2(residual);
                //oldInternalForces = internalForcesTotalVector;
                int iteration = 0;
                Array.Clear(deltaU, 0, deltaU.Length);
                while (residualNorm > Tolerance && iteration < MaxIterations)
                {
                    //////add timer
                    ////var iterationTimer = System.Diagnostics.Stopwatch.StartNew();
                    Tuple<double[], double[]> residuals = discretization.MMCPCGSeperateReducedExternalForcesVectors(residual);
                    double[] residualM = residuals.Item1;
                    double[] residualC = residuals.Item2;
                    //
                    stiffnessMatrices = discretization.MMCPCGCreateTotalStiffnessMatrix();
                    //totalStifnessMatrix = stiffnessMatrices.Item1;
                    K = stiffnessMatrices.Item2;
                    B = stiffnessMatrices.Item3;
                    C = stiffnessMatrices.Item4;
                    //
                    //
                    //if (residualNorm >= 10.00)
                    //{
                    //    stiffnessMatrices = discretization.MMCPCGCreateTotalStiffnessMatrix();
                    //    K = stiffnessMatrices.Item2;
                    //    B = stiffnessMatrices.Item3;
                    //    C = stiffnessMatrices.Item4; count += 1;
                    //    countIncrements += 1;
                    //}
                    //else if (residualNorm >= 1.00)
                    //{
                    //    if (countIncrements % 2 == 0)
                    //    {
                    //        stiffnessMatrices = discretization.MMCPCGCreateTotalStiffnessMatrix();
                    //        K = stiffnessMatrices.Item2;
                    //        B = stiffnessMatrices.Item3;
                    //        C = stiffnessMatrices.Item4; count += 1;
                    //    }
                    //    countIncrements += 1;
                    //}
                    //else if (1000.0 * Tolerance <= residualNorm && residualNorm < 1.00)
                    //{
                    //    if (countIncrements % 5 == 0)
                    //    {
                    //        stiffnessMatrices = discretization.MMCPCGCreateTotalStiffnessMatrix();
                    //        K = stiffnessMatrices.Item2;
                    //        B = stiffnessMatrices.Item3;
                    //        C = stiffnessMatrices.Item4; count += 1;
                    //    }
                    //    countIncrements += 1;
                    //}
                    //else if (10.0 * Tolerance <= residualNorm && residualNorm < 1000.0 * Tolerance)
                    //{
                    //    if (countIncrements % 10 == 0)
                    //    {
                    //        stiffnessMatrices = discretization.MMCPCGCreateTotalStiffnessMatrix();
                    //        K = stiffnessMatrices.Item2;
                    //        B = stiffnessMatrices.Item3;
                    //        C = stiffnessMatrices.Item4; count += 1;
                    //    }
                    //    countIncrements += 1;
                    //}
                    //
                    //////add timer
                    ////var iterationWatch = System.Diagnostics.Stopwatch.StartNew();
                    //////
                    ///
                    var localSolverResults = linearSolver.Solve(K, B, C, residualM, residualC, fillInLevels, fillInDegree, true);
                    countIterations = localSolverResults.Item2;
                    IterationsCounters.Add((double)(countIterations));
                    deltaU = VectorOperations.VectorVectorSubtraction(deltaU, localSolverResults.Item1);
                    //deltaU = VectorOperations.VectorVectorSubtraction(deltaU, linearSolver.Solve(K, B, C, residualM, residualC, fillInLevels, fillInDegree));

                    //////
                    ////iterationWatch.Stop();
                    ////solutionElapsedMs = (double)iterationWatch.ElapsedMilliseconds;
                    ////solutionTimerList.Add(solutionElapsedMs);
                    //////timer
                    //////
                    tempSolutionVector = VectorOperations.VectorVectorAddition(solutionVector, deltaU);
                    discretization.MMCPCGUpdateElementsIncrementalDisplacements(VectorOperations.VectorVectorAddition(dU, deltaU));
                    discretization.MMCPCGUpdateDisplacements(tempSolutionVector);
                    internalForcesTotalVector = discretization.MMCPCGCreateTotalInternalForcesVectors();
                    residual = VectorOperations.VectorVectorSubtraction(internalForcesTotalVector, incrementalExternalForcesVector);
                    residualNorm = VectorOperations.VectorNorm2(residual);
                    if (residualNorm <= Tolerance)
                    {
                        OnConvergenceResult(new ConvergenceValues() { ResidualNorm = residualNorm, LoadStep = i, Iteration = iteration, Tolerance = Tolerance, ConvergenceResult = true });
                    }
                    else
                    {
                        OnConvergenceResult(new ConvergenceValues() { ResidualNorm = residualNorm, LoadStep = i, Iteration = iteration, Tolerance = Tolerance, ConvergenceResult = false });
                    }
                    iteration = iteration + 1;
                    //(Application.Current.Windows[0] as MainWindow).LogTool.Text = "ok"; 
                    //OnConvergenceResult("Newton-Raphson: Solution not converged at load step" + iteration);

                    //////
                    ////iterationTimer.Stop();
                    ////AssembleElapsedMs = (double)iterationTimer.ElapsedMilliseconds;
                    ////SystemAssembleTimerList.Add(AssembleElapsedMs);
                    //////timer
                }
                double[] internalForcesTotalVectorRearranged = discretization.MMCPCGRearrange(internalForcesTotalVector);
                InternalForces.Add(i + 1, internalForcesTotalVectorRearranged);
                solutionVector = VectorOperations.VectorVectorAddition(solutionVector, deltaU);
                double[] solutionVectorRearranged = discretization.MMCPCGRearrange(solutionVector);
                Solutions.Add(i + 1, solutionVectorRearranged);
                discretization.UpdateContactTangentialProperties();
                if (iteration >= MaxIterations)
                {
                    OnConvergenceResult(new ConvergenceValues() { ResidualNorm = residualNorm, LoadStep = i, Iteration = iteration, Tolerance = Tolerance, ConvergenceResult = false });
                    LoadStepConvergence.Add("Solution not converged.");
                    break;
                }
                LoadStepConvergence.Add("Solution converged.");

                string exportPath = @"C:\Users\Public\Documents\";
                MatrixOperations.PrintMatrixToFile(stiffnessMatrices.Item1, exportPath + "Ktotal" + i+".dat");
                MatrixOperations.PrintMatrixToFile(stiffnessMatrices.Item2, exportPath + "K" + i + ".dat");
                MatrixOperations.PrintMatrixToFile(stiffnessMatrices.Item3, exportPath + "B" + i + ".dat");
                MatrixOperations.PrintMatrixToFile(stiffnessMatrices.Item4, exportPath + "C" + i + ".dat");
                discretization.MMCPCGUpdateDisplacements(solutionVector);
                stiffnessMatrices = discretization.MMCPCGCreateTotalStiffnessMatrix();
                K = stiffnessMatrices.Item2;
            }
            solutionVector = discretization.MMCPCGRearrange(solutionVector);
            //////
            ////VectorOperations.PrintVectorToFile(solutionTimerList.ToArray(), @"C:\Users\Public\Documents\" + "SolutionTime");
            ////VectorOperations.PrintVectorToFile(SystemAssembleTimerList.ToArray(), @"C:\Users\Public\Documents\" + "AssembleTime");
            //////
            VectorOperations.PrintVectorToFile(IterationsCounters.ToArray(), @"C:\Users\Public\Documents\" + "MMPCGIterationsExm2.dat");

            return solutionVector;
        }

        public override double[] Solve(IAssembly assembly, ILinearSolution linearScheme, double[] extForceVector)
        {
            InternalForces = new Dictionary<int, double[]>();
            Solutions = new Dictionary<int, double[]>();
            LoadStepConvergence = new List<string>();
            if (localSolutionVector == null)
            {
                localSolutionVector = new double[extForceVector.Length];
            }
            discretization = assembly;
            linearSolver = new MMCPCGSolver();
            lambda = 1.0 / numberOfLoadSteps;
            //double[] solution = null;

            //Thread tcore1 = new Thread(() => LoadControlledNR(forceVector));
            //tcore1.Start();
            //tcore1.Join();
            double[] solution = MMCPCGLoadControlledNR(extForceVector);
            return solution;
        }

    }
}

