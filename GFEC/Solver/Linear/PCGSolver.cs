﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GFEC
{
    class PCGSolver : LinearSolution
    {
        int maxIterations = 1000;
        double tolerance = 1e-6;

        private double[] PCG(double[,] stiffnessMatrix, double[] forceVector)
        {
            double[] solutionVector = new double[forceVector.Length];
            double[,] preconditioner = new double[stiffnessMatrix.GetLength(0), stiffnessMatrix.GetLength(1)];
            for (int i = 0; i < preconditioner.GetLength(0); i++)
            {
                preconditioner[i, i] = 1 / stiffnessMatrix[i, i];
            }
            double[] residual = VectorOperations.VectorVectorSubtraction(
                forceVector,
                VectorOperations.MatrixVectorProduct(stiffnessMatrix, solutionVector)
                );
            double[] preconVector = VectorOperations.MatrixVectorProduct(preconditioner, residual);
            for (int iter = 0; iter < maxIterations; iter++)
            {
                double[] u = VectorOperations.MatrixVectorProduct(stiffnessMatrix, preconVector);
                double residualDotOld = VectorOperations.VectorDotProduct(residual, preconVector);
                double alpha = residualDotOld / VectorOperations.VectorDotProduct(preconVector, u);
                solutionVector = VectorOperations.VectorVectorAddition(solutionVector, VectorOperations.VectorScalarProductNew(preconVector, alpha));
                residual = VectorOperations.VectorVectorSubtraction(residual, VectorOperations.VectorScalarProductNew(u, alpha));
                if (VectorOperations.VectorDotProduct(residual, residual) < tolerance)
                {
                    break;
                }
                double residualDotNew = VectorOperations.VectorDotProduct(residual, VectorOperations.MatrixVectorProduct(preconditioner, residual));
                double beta = residualDotNew / residualDotOld;
                preconVector = VectorOperations.VectorVectorAddition(
                    VectorOperations.MatrixVectorProduct(preconditioner, residual),
                    VectorOperations.VectorScalarProductNew(preconVector, beta)
                    );
            }
            return solutionVector;
        }

        public override double[] Solve(double[,] stiffnessMatrix, double[] forceVector)
        {
            double[] solution = PCG(stiffnessMatrix, forceVector);
            return solution;
        }
    }
}
