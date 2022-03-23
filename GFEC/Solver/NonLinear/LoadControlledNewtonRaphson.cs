using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;

namespace GFEC
{
    public class LoadControlledNewtonRaphson : NonLinearSolution
    {
        private double[] localSolutionVector;

        public LoadControlledNewtonRaphson()
        {

        }
        public LoadControlledNewtonRaphson(double[] exSolution)
        {
            localSolutionVector = exSolution;
        }
        private double[] LoadControlledNR(double[] forceVector)
        {
            //
            //List<double> IterationsCounters = new List<double>();
            //int countIterations = new int();
            //
            ////timers
            //List<double> solutionTimerList = new List<double>();
            //List<double> SystemAssembleTimerList = new List<double>();
            ////
            double[] incrementDf = VectorOperations.VectorScalarProductNew(forceVector, lambda);
            double[] solutionVector = localSolutionVector;
            double[] incrementalExternalForcesVector = new double[forceVector.Length];
            double[] tempSolutionVector = new double[solutionVector.Length];
            double[] deltaU = new double[solutionVector.Length];
            double[] internalForcesTotalVector;
            double[] dU;
            double[] residual;
            double residualNorm;
            discretization.InitializeContactTangentialProperties();
            discretization.InitializeContactSurfaceVectors();
            discretization.InitializeEASParameters();//Added EAS
            discretization.CalculateEASMatrices();//added EAS
            for (int i = 0; i < numberOfLoadSteps; i++)
            {
                ////add timer
                //var timer = System.Diagnostics.Stopwatch.StartNew();
                incrementalExternalForcesVector = VectorOperations.VectorVectorAddition(incrementalExternalForcesVector, incrementDf);
                discretization.UpdateDisplacements(solutionVector);
                //discretization.UpdateEASParameters(solutionVector);//added EAS
                discretization.StoreFinalStepDisplacementVector(solutionVector);//Added EAS

                discretization.UpdateContactSurfaceVectors();
                internalForcesTotalVector = discretization.CreateTotalInternalForcesVector();
                double[,] stiffnessMatrix = discretization.CreateTotalStiffnessMatrix();
                //-------------------------------------------------------------------------------------------------
                ////timer.Stop();
                ////var AssembleElapsedMs = (double)timer.ElapsedMilliseconds;
                ////SystemAssembleTimerList.Add(AssembleElapsedMs);
                ////timer
                //OnConvergenceResult("Newton-Raphson: Solution not converged at load step" + i); 

                ////add timer
                //var watch = System.Diagnostics.Stopwatch.StartNew();
                //
                //
                //var solverResults = linearSolver.Solve(stiffnessMatrix, incrementDf, true);
                //dU = solverResults.Item1;
                //countIterations = solverResults.Item2;
                //IterationsCounters.Add((double)(countIterations));
                //
                //
                ////watch.Stop();
                ////var solutionElapsedMs = (double)watch.ElapsedMilliseconds;
                ////solutionTimerList.Add(solutionElapsedMs);
                ////timer
                //-------------------------------------------------------------------------------------------------
                dU = linearSolver.Solve(stiffnessMatrix, incrementDf);
                discretization.UpdateElementsIncrementalDisplacements(dU);
                solutionVector = VectorOperations.VectorVectorAddition(solutionVector, dU);
                residual = VectorOperations.VectorVectorSubtraction(internalForcesTotalVector, incrementalExternalForcesVector);
                residualNorm = VectorOperations.VectorNorm2(residual);
                int iteration = 0;
                Array.Clear(deltaU, 0, deltaU.Length);
                while (residualNorm > Tolerance && iteration < MaxIterations)
                {
                    stiffnessMatrix = discretization.CreateTotalStiffnessMatrix();
                    deltaU = VectorOperations.VectorVectorSubtraction(deltaU, linearSolver.Solve(stiffnessMatrix, residual));
                    tempSolutionVector = VectorOperations.VectorVectorAddition(solutionVector, deltaU);
                    discretization.UpdateElementsIncrementalDisplacements(VectorOperations.VectorVectorAddition(dU, deltaU));
                    discretization.UpdateDisplacements(tempSolutionVector);
                    discretization.UpdateEASParameters(tempSolutionVector);//added EAS
                    discretization.CalculateEASMatrices();//added EAS
                    internalForcesTotalVector = discretization.CreateTotalInternalForcesVector();
                    residual = VectorOperations.VectorVectorSubtraction(internalForcesTotalVector, incrementalExternalForcesVector);
                    residualNorm = VectorOperations.VectorNorm2(residual);
                    if (residualNorm <= Tolerance)
                    {
                        OnConvergenceResult("Newton-Raphson: Load Step "+i+ " - Solution converged at iteration " + iteration + " - Residual Norm = " +residualNorm);
                    }
                    else
                    {
                        OnConvergenceResult("Newton-Raphson: Load Step " + i + " - Solution not converged at iteration " + iteration + " - Residual Norm = " + residualNorm);
                    }
                    iteration += 1;
                }
                InternalForces.Add(i + 1, internalForcesTotalVector);
                solutionVector = VectorOperations.VectorVectorAddition(solutionVector, deltaU);
                Solutions.Add(i + 1, solutionVector);
                discretization.UpdateContactTangentialProperties();
                //discretization.NextStepFrictionUpdate();
                if (iteration >= MaxIterations)
                {
                    OnConvergenceResult("Newton-Raphson did not converge at Load Step " + i + ". Exiting solution.");
                    LoadStepConvergence.Add("Solution not converged.");
                    break;
                }

                LoadStepConvergence.Add("Solution converged.");
            }
            ////
            //VectorOperations.PrintVectorToFile(solutionTimerList.ToArray(), @"C:\Users\Public\Documents\" + "SolutionTime");
            //VectorOperations.PrintVectorToFile(SystemAssembleTimerList.ToArray(), @"C:\Users\Public\Documents\" + "AssembleTime");
            ////
            ///
            //VectorOperations.PrintVectorToFile(IterationsCounters.ToArray(), @"C:\Users\Public\Documents\" + "PCGIterationsExm2.dat");
            //
            return solutionVector;
        }
        
        public override double[] Solve(IAssembly assembly, ILinearSolution linearScheme, double[] forceVector)
        {
            InternalForces = new Dictionary<int, double[]>();
            Solutions = new Dictionary<int, double[]>();
            LoadStepConvergence = new List<string>();
            if (localSolutionVector == null)
            {
                localSolutionVector = new double[forceVector.Length];
            }
            discretization = assembly;
            linearSolver = linearScheme;
            lambda = 1.0 / numberOfLoadSteps;
            //double[] solution = null;

            //Thread tcore1 = new Thread(() => LoadControlledNR(forceVector));
            //tcore1.Start();
            //tcore1.Join();
            double[] solution = LoadControlledNR(forceVector);
            return solution;
        }

    }
}

