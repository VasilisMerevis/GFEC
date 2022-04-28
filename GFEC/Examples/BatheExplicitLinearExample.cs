using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GFEC
{
    
    public static class BatheExplicitLinearExample
    {
        public static void SolveExample()
        {
            double[,] M = new double[,]
            {
                { 2.0 , 0.0 },
                { 0.0 , 1.0 }
            };

            double[,] K = new double[,]
            {
                { 6.0 , -2.0 },
                { -2.0 , 4.0 }
            };

            double[] F = new double[] { 0.0, 10.0 };

            InitialConditions initialValues = new InitialConditions();
            initialValues.InitialAccelerationVector = new double[] { 0.0, 10.0 };
            initialValues.InitialDisplacementVector = new double[] { 0.0, 0.0 };
            initialValues.InitialVelocityVector = new double[] { 0.0, 0.0 };
            initialValues.InitialTime = 0.0;

            BatheExplicit solver = new BatheExplicit(new LUFactorization(), initialValues, 10.0, 10, K, M, F);
            //solver.InitialValues = initialValues;
            //solver.CustomMassMatrix = M;
            //solver.CustomStiffnessMatrix = K;
            //solver.CustomDampingMatrix = new double[2, 2];
            //solver.ExternalForcesVector = F;
            //solver.ActivateNonLinearSolution = false;
            //solver.LinearSolver = new LUFactorization();
            solver.SolveBatheExplicit();
            //solver.PrintExplicitSolution();
        }

        

        public static Results RunStaticExample()
        {
            return new Results();
        }
    }
}
