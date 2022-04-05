using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GFEC
{
    class MMCPCGSolver : LinearSolution
    {
        int maxIterations = 1000;
        double tolerance = 1e-9;
        int fillInDegree = 0;

        //    private double[] MMCPCG(double[,] stiffnessMatrixK, double[,] stiffnessMatrixB, double[,] stiffnessMatrixC,
        //double[] forceVectorM, double[] forceVectorC, int[,] fillInLevels, int fillLevel)
        //    {
        //        //var gauss = new GaussSolver();
        //        double[,] stiffnessMatrixBT = MatrixOperations.Transpose(stiffnessMatrixB);
        //        double[] solutionVector = new double[forceVectorM.Length + forceVectorC.Length];
        //        double[] forceVector = new double[forceVectorM.Length + forceVectorC.Length];
        //        forceVector = VectorOperations.CreateFullVectorFromTwoVectors(forceVectorM, forceVectorC);
        //        double[,] preconditionerLowerMatrix = new double[forceVectorM.Length, forceVectorM.Length];
        //        //var isSymmetric = MatrixOperations.CheckIfSymmetric(stiffnessMatrixK);
        //        //stiffnessMatrixK = MatrixOperations.CorrectSymmetry(stiffnessMatrixK);
        //        //string StiffnessMatrixname = "K" + ".dat";
        //        //MatrixOperations.PrintMatrixToFile2(stiffnessMatrixK, @"C:\Users\Public\Documents\" + StiffnessMatrixname);
        //        preconditionerLowerMatrix = MatrixOperations.IChol(stiffnessMatrixK, fillInLevels, fillLevel);
        //        double[,] preconditionerUpperMatrix = new double[forceVectorM.Length, forceVectorM.Length];
        //        preconditionerUpperMatrix = MatrixOperations.Transpose(preconditionerLowerMatrix);
        //        //
        //        //double[,] cLowerMatrix = new double[forceVectorC.Length, forceVectorC.Length];
        //        //cLowerMatrix = MatrixOperations.CChol(stiffnessMatrixC);
        //        //double[,] cUpperMatrix = new double[forceVectorC.Length, forceVectorC.Length];
        //        //cUpperMatrix = MatrixOperations.Transpose(cLowerMatrix);
        //        //
        //        var cMatrices = MatrixOperations.LUfactorizedMatrices(stiffnessMatrixC);
        //        double[,] cLowerMatrix = cMatrices.Item1;
        //        double[,] cUpperMatrix = cMatrices.Item2;

        //        double[,] totalStiffnessMatrix = new double[forceVectorM.Length + forceVectorC.Length, forceVectorM.Length + forceVectorC.Length];
        //        totalStiffnessMatrix = MatrixOperations.CreateBlockMatrix(stiffnessMatrixK, stiffnessMatrixB, stiffnessMatrixBT, stiffnessMatrixC);
        //        double[] residualM = new double[forceVectorM.Length];
        //        residualM = forceVectorM;
        //        double[] residualC = new double[forceVectorC.Length];
        //        residualC = forceVectorC;
        //        double[] residual = new double[forceVectorM.Length + forceVectorC.Length];
        //        double relativeResidual = 0;
        //        double[] z = new double[solutionVector.Length];
        //        double[] p = new double[solutionVector.Length];
        //        double[] q = new double[solutionVector.Length];
        //        double ihta = new double();

        //        double[] vmIntermediateVector = new double[forceVectorM.Length];
        //        double[] vm = new double[forceVectorM.Length];
        //        double[] w = new double[forceVectorC.Length];
        //        double[] zc = new double[forceVectorC.Length];
        //        double[] zcIntermediateVector = new double[forceVectorC.Length];
        //        double[] zmIntermediateVector = new double[forceVectorM.Length];
        //        double[] zm = new double[forceVectorM.Length];
        //        //
        //        //solve z
        //        vmIntermediateVector = ForwardSubstitution(preconditionerLowerMatrix, residualM);
        //        vm = BackSubstitution(preconditionerUpperMatrix, vmIntermediateVector);

        //        w = VectorOperations.VectorVectorSubtraction(residualC,
        //            VectorOperations.MatrixVectorProduct(stiffnessMatrixBT, vm));

        //        //zc = gauss.Solve(stiffnessMatrixC, w);
        //        zcIntermediateVector = ForwardSubstitution(cLowerMatrix, w);
        //        zc = BackSubstitution(cUpperMatrix, zcIntermediateVector);


        //        vm = VectorOperations.VectorVectorSubtraction(residualM, VectorOperations.MatrixVectorProduct(stiffnessMatrixB, zc));
        //        zmIntermediateVector = ForwardSubstitution(preconditionerLowerMatrix, vm);
        //        zm = BackSubstitution(preconditionerUpperMatrix, zmIntermediateVector);
        //        z = VectorOperations.CreateFullVectorFromTwoVectors(zm, zc);
        //        residual = VectorOperations.CreateFullVectorFromTwoVectors(residualM, residualC);
        //        double initialResidualNorm2 = VectorOperations.VectorNorm2(residual);
        //        //
        //        p = z;
        //        q = VectorOperations.MatrixVectorProduct(totalStiffnessMatrix, p);
        //        ihta = VectorOperations.VectorDotProduct(p, residual) / VectorOperations.VectorDotProduct(p, q);

        //        for (int iter = 0; iter < maxIterations; iter++)
        //        {
        //            solutionVector = VectorOperations.VectorVectorAddition(solutionVector, VectorOperations.VectorScalarProductNew(p, ihta));
        //            double residualDotOld = VectorOperations.VectorDotProduct(z, residual);
        //            //residual = VectorOperations.VectorVectorSubtraction(residual, VectorOperations.VectorScalarProductNew(q, ihta));
        //            residual = VectorOperations.VectorVectorSubtraction(forceVector, VectorOperations.MatrixVectorProduct(totalStiffnessMatrix, solutionVector));
        //            double currentResidualNorm2 = VectorOperations.VectorNorm2(residual);
        //            relativeResidual = currentResidualNorm2 / initialResidualNorm2;
        //            //double m = VectorOperations.VectorDotProduct(residual, residual);
        //            if (relativeResidual < tolerance)
        //            {
        //                break;
        //            }
        //            var residuals = VectorOperations.SeperateVectorToTwoVectors(residual, residualM.Length);
        //            residualM = residuals.Item1;
        //            residualC = residuals.Item2;

        //            vmIntermediateVector = ForwardSubstitution(preconditionerLowerMatrix, residualM);
        //            vm = BackSubstitution(preconditionerUpperMatrix, vmIntermediateVector);
        //            w = VectorOperations.VectorVectorSubtraction(residualC,
        //                VectorOperations.MatrixVectorProduct(stiffnessMatrixBT, vm));
        //            //zc = gauss.Solve(stiffnessMatrixC, w);

        //            zcIntermediateVector = ForwardSubstitution(cLowerMatrix, w);
        //            zc = BackSubstitution(cUpperMatrix, zcIntermediateVector);

        //            vm = VectorOperations.VectorVectorSubtraction(residualM, VectorOperations.MatrixVectorProduct(stiffnessMatrixB, zc));
        //            zmIntermediateVector = ForwardSubstitution(preconditionerLowerMatrix, vm);
        //            zm = BackSubstitution(preconditionerUpperMatrix, zmIntermediateVector);
        //            z = VectorOperations.CreateFullVectorFromTwoVectors(zm, zc);
        //            p = VectorOperations.VectorVectorAddition(z, VectorOperations.VectorScalarProductNew(p, VectorOperations.VectorDotProduct(z, residual) / residualDotOld));
        //            q = VectorOperations.MatrixVectorProduct(totalStiffnessMatrix, p);
        //            ihta = VectorOperations.VectorDotProduct(z, residual) / VectorOperations.VectorDotProduct(p, q);
        //        }
        //        return solutionVector;
        //    }

        //    public override double[] Solve(double[,] stiffnessMatrixK, double[,] stiffnessMatrixB, double[,] stiffnessMatrixC,
        //        double[] forceVectorM, double[] forceVectorC, int[,] fillInLevels, int fillLevel)
        //    {
        //        double[] solution = MMCPCG(stiffnessMatrixK, stiffnessMatrixB, stiffnessMatrixC,
        //            forceVectorM, forceVectorC, fillInLevels, fillLevel);
        //        return solution;
        //    }
        private Tuple<double[], int> MMCPCG(double[,] stiffnessMatrixK, double[,] stiffnessMatrixB, double[,] stiffnessMatrixC,
double[] forceVectorM, double[] forceVectorC, int[,] fillInLevels, int fillLevel)
        {
            int numberOfIterations = new int();
            //var gauss = new GaussSolver();
            double[,] stiffnessMatrixBT = MatrixOperations.Transpose(stiffnessMatrixB);
            double[] solutionVector = new double[forceVectorM.Length + forceVectorC.Length];
            double[] forceVector = new double[forceVectorM.Length + forceVectorC.Length];
            forceVector = VectorOperations.CreateFullVectorFromTwoVectors(forceVectorM, forceVectorC);
            double[,] preconditionerLowerMatrix = new double[forceVectorM.Length, forceVectorM.Length];
            //var isSymmetric = MatrixOperations.CheckIfSymmetric(stiffnessMatrixK);
            //stiffnessMatrixK = MatrixOperations.CorrectSymmetry(stiffnessMatrixK);
            //string StiffnessMatrixname = "K" + ".dat";
            //MatrixOperations.PrintMatrixToFile2(stiffnessMatrixK, @"C:\Users\Public\Documents\" + StiffnessMatrixname);
            preconditionerLowerMatrix = MatrixOperations.IChol(stiffnessMatrixK, fillInLevels, fillLevel);
            double[,] preconditionerUpperMatrix = new double[forceVectorM.Length, forceVectorM.Length];
            preconditionerUpperMatrix = MatrixOperations.Transpose(preconditionerLowerMatrix);
            //
            //double[,] cLowerMatrix = new double[forceVectorC.Length, forceVectorC.Length];
            //cLowerMatrix = MatrixOperations.CChol(stiffnessMatrixC);
            //double[,] cUpperMatrix = new double[forceVectorC.Length, forceVectorC.Length];
            //cUpperMatrix = MatrixOperations.Transpose(cLowerMatrix);
            //
            var cMatrices = MatrixOperations.LUfactorizedMatrices(stiffnessMatrixC);
            double[,] cLowerMatrix = cMatrices.Item1;
            double[,] cUpperMatrix = cMatrices.Item2;

            double[,] totalStiffnessMatrix = new double[forceVectorM.Length + forceVectorC.Length, forceVectorM.Length + forceVectorC.Length];
            totalStiffnessMatrix = MatrixOperations.CreateBlockMatrix(stiffnessMatrixK, stiffnessMatrixB, stiffnessMatrixBT, stiffnessMatrixC);
            double[] residualM = new double[forceVectorM.Length];
            residualM = forceVectorM;
            double[] residualC = new double[forceVectorC.Length];
            residualC = forceVectorC;
            double[] residual = new double[forceVectorM.Length + forceVectorC.Length];
            double relativeResidual = 0;
            double[] z = new double[solutionVector.Length];
            double[] p = new double[solutionVector.Length];
            double[] q = new double[solutionVector.Length];
            double ihta = new double();

            double[] vmIntermediateVector = new double[forceVectorM.Length];
            double[] vm = new double[forceVectorM.Length];
            double[] w = new double[forceVectorC.Length];
            double[] zc = new double[forceVectorC.Length];
            double[] zcIntermediateVector = new double[forceVectorC.Length];
            double[] zmIntermediateVector = new double[forceVectorM.Length];
            double[] zm = new double[forceVectorM.Length];
            //
            //solve z
            vmIntermediateVector = ForwardSubstitution(preconditionerLowerMatrix, residualM);
            vm = BackSubstitution(preconditionerUpperMatrix, vmIntermediateVector);

            w = VectorOperations.VectorVectorSubtraction(residualC,
                VectorOperations.MatrixVectorProduct(stiffnessMatrixBT, vm));

            //zc = gauss.Solve(stiffnessMatrixC, w);
            zcIntermediateVector = ForwardSubstitution(cLowerMatrix, w);
            zc = BackSubstitution(cUpperMatrix, zcIntermediateVector);


            vm = VectorOperations.VectorVectorSubtraction(residualM, VectorOperations.MatrixVectorProduct(stiffnessMatrixB, zc));
            zmIntermediateVector = ForwardSubstitution(preconditionerLowerMatrix, vm);
            zm = BackSubstitution(preconditionerUpperMatrix, zmIntermediateVector);
            z = VectorOperations.CreateFullVectorFromTwoVectors(zm, zc);
            residual = VectorOperations.CreateFullVectorFromTwoVectors(residualM, residualC);
            double initialResidualNorm2 = VectorOperations.VectorNorm2(residual);
            //
            p = z;
            q = VectorOperations.MatrixVectorProduct(totalStiffnessMatrix, p);
            ihta = VectorOperations.VectorDotProduct(p, residual) / VectorOperations.VectorDotProduct(p, q);

            for (int iter = 0; iter < maxIterations; iter++)
            {
                solutionVector = VectorOperations.VectorVectorAddition(solutionVector, VectorOperations.VectorScalarProductNew(p, ihta));
                double residualDotOld = VectorOperations.VectorDotProduct(z, residual);
                //residual = VectorOperations.VectorVectorSubtraction(residual, VectorOperations.VectorScalarProductNew(q, ihta));
                residual = VectorOperations.VectorVectorSubtraction(forceVector, VectorOperations.MatrixVectorProduct(totalStiffnessMatrix, solutionVector));
                double currentResidualNorm2 = VectorOperations.VectorNorm2(residual);
                relativeResidual = currentResidualNorm2 / initialResidualNorm2;
                //double m = VectorOperations.VectorDotProduct(residual, residual);
                if (relativeResidual < tolerance)
                {
                    numberOfIterations = iter + 1;
                    break;
                }
                var residuals = VectorOperations.SeperateVectorToTwoVectors(residual, residualM.Length);
                residualM = residuals.Item1;
                residualC = residuals.Item2;

                vmIntermediateVector = ForwardSubstitution(preconditionerLowerMatrix, residualM);
                vm = BackSubstitution(preconditionerUpperMatrix, vmIntermediateVector);
                w = VectorOperations.VectorVectorSubtraction(residualC,
                    VectorOperations.MatrixVectorProduct(stiffnessMatrixBT, vm));
                //zc = gauss.Solve(stiffnessMatrixC, w);

                zcIntermediateVector = ForwardSubstitution(cLowerMatrix, w);
                zc = BackSubstitution(cUpperMatrix, zcIntermediateVector);

                vm = VectorOperations.VectorVectorSubtraction(residualM, VectorOperations.MatrixVectorProduct(stiffnessMatrixB, zc));
                zmIntermediateVector = ForwardSubstitution(preconditionerLowerMatrix, vm);
                zm = BackSubstitution(preconditionerUpperMatrix, zmIntermediateVector);
                z = VectorOperations.CreateFullVectorFromTwoVectors(zm, zc);
                p = VectorOperations.VectorVectorAddition(z, VectorOperations.VectorScalarProductNew(p, VectorOperations.VectorDotProduct(z, residual) / residualDotOld));
                q = VectorOperations.MatrixVectorProduct(totalStiffnessMatrix, p);
                ihta = VectorOperations.VectorDotProduct(z, residual) / VectorOperations.VectorDotProduct(p, q);
            }
            return new Tuple<double[], int>(solutionVector, numberOfIterations);
        }

        public override Tuple<double[],int> Solve(double[,] stiffnessMatrixK, double[,] stiffnessMatrixB, double[,] stiffnessMatrixC,
            double[] forceVectorM, double[] forceVectorC, int[,] fillInLevels, int fillLevel, bool countIterations)
        {
            var solution = MMCPCG(stiffnessMatrixK, stiffnessMatrixB, stiffnessMatrixC,
                forceVectorM, forceVectorC, fillInLevels, fillLevel);
            return solution;
        }
    }
}