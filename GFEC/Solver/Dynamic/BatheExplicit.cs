using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GFEC
{
    public class BatheExplicit
    {
        private double totalTime, timeStep;
        private int timeStepsNumber;
        private Dictionary<int, double[]> displacement = new Dictionary<int, double[]>();
        private Dictionary<int, double[]> velocity = new Dictionary<int, double[]>();
        private Dictionary<int, double[]> acceleration = new Dictionary<int, double[]>();
        private Dictionary<int, double[]> exForces = new Dictionary<int, double[]>();
        int totalDOFs;
        double[,] massMatrix, dampingMatrix;
        double[,] stiffnessMatrix;
        double[] externalForcesVector;
        double a0, a1, a2, a3, a4, a5, a6, a7;
        double q0, q1, q2, p;
        double[] initialDisplacementVector, initialVelocityVector, initialAccelerationVector;
        double initialTime;
        private ILinearSolution linearSolver;

        public BatheExplicit(ILinearSolution linearSolver, InitialConditions initialValues, double totalTime, int timeStepsNumber, double[,] stiffnessMatrix, double[,] massMatrix, double[] externalForcesVector)
        {
            totalDOFs = stiffnessMatrix.GetLength(0);
            this.totalTime = totalTime;
            this.timeStepsNumber = timeStepsNumber;
            timeStep = totalTime / timeStepsNumber;
            this.massMatrix = massMatrix;
            this.stiffnessMatrix = stiffnessMatrix;
            dampingMatrix = new double[totalDOFs, totalDOFs];
            this.externalForcesVector = externalForcesVector;
            initialDisplacementVector = initialValues.InitialDisplacementVector;
            initialVelocityVector = initialValues.InitialVelocityVector;
            initialAccelerationVector = initialValues.InitialAccelerationVector;
            initialTime = initialValues.InitialTime;
            this.linearSolver = linearSolver;
            p = 0.54;
        }

        private void Calculate_qValues()
        {
            q1 = (1.0 - 2.0 * p) / (2.0 * p * (1.0 - p));
            q2 = 0.5 - p * q1;
            q0 = -q1 - q2 + 0.5;
        }

        private void Calculate_aValues()
        {
            double deltat = timeStep;
            a0 = p * deltat;
            a1 = 0.5 * Math.Pow(p * deltat, 2);
            a2 = a0 / 2.0;
            a3 = (1.0 - p) * deltat;
            a4 = 0.5 * Math.Pow((1.0 - p) * deltat, 2);
            a5 = q0 * a3;
            a6 = (0.5 + q1) * a3;
            a7 = q2 * a3;
        }

        private double[] U_middle(double[] u_previous, double[] du_previous, double[] ddu_previous, double a0, double a1)
        {
            double[] a0du = VectorOperations.VectorScalarProductNew(du_previous, a0);
            double[] a1ddu = VectorOperations.VectorScalarProductNew(ddu_previous, a1);
            double[] u_middle = VectorOperations.VectorVectorAddition(u_previous,
                VectorOperations.VectorVectorAddition(a0du, a1ddu));
            return u_middle;
        }

        private double[] R_hat_middle(double[] r_previous, double[] r_current, double p)
        {
            double[] part1 = VectorOperations.VectorScalarProductNew(r_previous, 1.0 - p);
            double[] part2 = VectorOperations.VectorScalarProductNew(r_current, p);
            double[] r_hat_middle = VectorOperations.VectorVectorAddition(part1, part2);
            return r_hat_middle;
        }

        private double[] R_roundhat_middle(double[] r_hat_middle, double[] u_middle, double[] du_previous, double[] ddu_previous,
            double[,] kMatrix, double[,] cMatrix, double a0)
        {
            double[] stifPart = VectorOperations.MatrixVectorProduct(kMatrix, u_middle);
            double[] duplusa0ddu = VectorOperations.VectorVectorAddition(du_previous,
                VectorOperations.VectorScalarProductNew(ddu_previous, a0));
            double[] dampPart = VectorOperations.MatrixVectorProduct(cMatrix, duplusa0ddu);
            double[] r_roundhat1 = VectorOperations.VectorVectorSubtraction(r_hat_middle, stifPart);
            double[] r_roundhat2 = VectorOperations.VectorVectorSubtraction(r_roundhat1, dampPart);
            return r_roundhat2;
        }

        private double[] DDU_middle(double[,] mMatrix, double[] r_roundhatmiddle)
        {
            ILinearSolution solver = new CholeskyFactorization();
            double[] ddu_middle = solver.Solve(mMatrix, r_roundhatmiddle);
            return ddu_middle;
        }

        private double[] DU_middle(double[] du_previous, double[] ddu_previous, double[] ddu_middle, double a2)
        {
            double[] parenthesis = VectorOperations.VectorVectorAddition(ddu_previous, ddu_middle);
            double[] part2 = VectorOperations.VectorScalarProductNew(parenthesis, a2);
            double[] du_middle = VectorOperations.VectorVectorAddition(du_previous, part2);
            return du_middle;
        }

        private double[] U_current(double[] u_middle, double[] du_middle, double[] ddu_middle, double a3, double a4)
        {
            double[] part2 = VectorOperations.VectorScalarProductNew(du_middle, a3);
            double[] part3 = VectorOperations.VectorScalarProductNew(ddu_middle, a4);
            double[] u_current = VectorOperations.VectorVectorAddition(u_middle,
                VectorOperations.VectorVectorAddition(part2, part3));
            return u_current;
        }

        private double[] R_roundhat_current(double[] r_current, double[,] kMatrix, double[] u_current, double[,] cMatrix,
            double[] du_middle, double[] ddu_middle, double a3)
        {
            double[] part2 = VectorOperations.MatrixVectorProduct(kMatrix, u_current);
            double[] parenthesis = VectorOperations.VectorVectorAddition(du_middle,
                VectorOperations.VectorScalarProductNew(ddu_middle, a3));
            double[] part3 = VectorOperations.MatrixVectorProduct(cMatrix, parenthesis);
            double[] r_roundhat_current1 = VectorOperations.VectorVectorSubtraction(r_current, part2);
            double[] r_roundhat_current2 = VectorOperations.VectorVectorSubtraction(r_roundhat_current1, part3);
            return r_roundhat_current2;
        }

        private double[] DDU_current(double[,] mMatrix, double[] r_roundhat_current)
        {
            ILinearSolution solver = new CholeskyFactorization();
            double[] ddu_current = solver.Solve(mMatrix, r_roundhat_current);
            return ddu_current;
        }

        private double[] DU_current(double[] du_middle, double[] ddu_previous, double[] ddu_middle,
            double[] ddu_current, double a5, double a6, double a7)
        {
            double[] part2 = VectorOperations.VectorScalarProductNew(ddu_previous, a5);
            double[] part3 = VectorOperations.VectorScalarProductNew(ddu_middle, a6);
            double[] part4 = VectorOperations.VectorScalarProductNew(ddu_current, a7);
            double[] du_current = VectorOperations.VectorVectorAddition(du_middle,
                VectorOperations.VectorVectorAddition(part2,
                VectorOperations.VectorVectorAddition(part3, part4)));
            return du_current;
        }

        private void CreateRForAllSteps()
        {
            exForces.Add(0, new double[externalForcesVector.Length]);
            for (int i = 1; i < timeStepsNumber; i++)
            {
                double[] forceStep = VectorOperations.VectorScalarProductNew(externalForcesVector, 1.0 / timeStepsNumber);
                exForces.Add(i, VectorOperations.VectorScalarProductNew(forceStep, i));
            }
        }
        public void SolveBatheExplicit()
        {

            displacement.Add(0, initialDisplacementVector);
            velocity.Add(0, initialVelocityVector);
            acceleration.Add(0, initialAccelerationVector);
            Calculate_qValues();
            Calculate_aValues();
            for (int i = 1; i < timeStepsNumber; i++)
            {
                double[] u_middle = U_middle(displacement[i - 1], velocity[i - 1], acceleration[i - 1], a0, a1);
                double[] r_hat_middle = R_hat_middle(exForces[i - 1], exForces[i], p);
                double[] r_roundhat_middle = R_roundhat_middle(r_hat_middle, u_middle, velocity[i - 1], acceleration[i - 1], stiffnessMatrix, dampingMatrix, a0);
                double[] ddu_middle = DDU_middle(massMatrix, r_roundhat_middle);
                double[] du_middle = DU_middle(velocity[i - 1], acceleration[i - 1], ddu_middle, a2);

                double[] u_current = U_current(u_middle, du_middle, ddu_middle, a3, a4);
                double[] r_roundhat_current = R_roundhat_current(exForces[i], stiffnessMatrix, u_current, dampingMatrix, du_middle, ddu_middle, a3);
                double[] ddu_current = DDU_current(massMatrix, r_roundhat_current);
                double[] du_current = DU_current(du_middle, acceleration[i - 1], ddu_middle, ddu_current, a5, a6, a7);

                displacement.Add(i, u_current);
                velocity.Add(i, du_current);
                acceleration.Add(i, ddu_current);
            }
        }
    }
}
