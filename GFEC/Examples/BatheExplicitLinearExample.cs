using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GFEC
{
    
    public static class BatheExplicitLinearExample
    {
        private static Results finalresutls;
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

            BatheExplicit solver = new BatheExplicit(new LUFactorization(), initialValues, 2.8, 10, K, M, F);
            
            solver.SolveBatheExplicit();

            ExplicitSolver solver2 = new ExplicitSolver(2.8, 10);
            solver2.InitialValues = initialValues;
            solver2.CustomMassMatrix = M;
            solver2.CustomStiffnessMatrix = K;
            solver2.CustomDampingMatrix = new double[2, 2];
            solver2.ExternalForcesVector = F;
            solver2.ActivateNonLinearSolution = false;
            solver2.LinearSolver = new LUFactorization();
            solver2.SolveExplicit();
            //solver.PrintExplicitSolution();
            finalresutls = new Results() { DynamicSolution = solver.displacement, TimeSteps = solver.TimeAtEachStep, SelectedDOF = 1, SelectedInterval = 1, SolutionType = "Dynamic" };
        }

        

        public static Results RunStaticExample()
        {
            return finalresutls;
        }
    }
}