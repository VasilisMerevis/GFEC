﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GFEC
{
    //Example Based on Bathe's Example 9.1
    public static class ExplicitLinearExample
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

            ExplicitSolver solver = new ExplicitSolver(2.8, 9);
            solver.InitialValues = initialValues;
            solver.CustomMassMatrix = M;
            solver.CustomStiffnessMatrix = K;
            solver.CustomDampingMatrix = new double[2, 2];
            solver.ExternalForcesVector = F;
            solver.ActivateNonLinearSolution = false;
            solver.LinearSolver = new LUFactorization();
            solver.SolveExplicit();
            finalresutls = new Results() { DynamicSolution = solver.displacement, TimeSteps = solver.TimeAtEachStep, SelectedDOF = 1, SelectedInterval = 1, SolutionType = "Dynamic" };
            //solver.PrintExplicitSolution();
        }

        public static void SolveNewmarkExample()
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

            ExplicitSolver solver = new ExplicitSolver(2.8, 10);
            solver.InitialValues = initialValues;
            solver.CustomMassMatrix = M;
            solver.CustomStiffnessMatrix = K;
            solver.CustomDampingMatrix = new double[2, 2];
            solver.ExternalForcesVector = F;
            solver.ActivateNonLinearSolution = false;
            solver.LinearSolver = new LUFactorization();
            solver.SolveNewmark();
            
            //solver.PrintExplicitSolution();
        }

        public static Results RunStaticExample()
        {
            return finalresutls;
        }

        public static Results RunStaticExample()
        {
            return new Results();
        }
    }
}
