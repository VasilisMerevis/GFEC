﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GFEC
{
    class ContactNtS2Df : IElement
    {
        public Dictionary<int, INode> Nodes { get; }
        public IElementProperties Properties { get; set; }
        public Dictionary<int, bool[]> ElementFreedomSignature { get; } = new Dictionary<int, bool[]>();
        public List<int> ElementFreedomList { get; set; }
        public double[] DisplacementVector { get; set; }
        public double[] AccelerationVector { get; set; }
        private double PenaltyFactor { get; set; }
        private double TangentPenaltyFactor { get; set; }
        private double FrictionCoef { get; set; }
        private double mhid;
        private double Ksi1Initial { get; set; }
        private double Ksi1Current { get; set; }
        private int counter;
        public void InitializeTangentialProperties()
        {
            throw new Exception("Needs to be removed. Has beeb used only for testing purposes");
        }
        public void UpdateTangentialProperties()
        {
            throw new Exception("Needs to be removed. Has beeb used only for testing purposes");
        }
        public void InitializeContactSurfaceGeometry()
        {
            throw new Exception("Needs to be removed. Has beeb used only for testing purposes");
        }
        public void UpdateContactSurfaceGeometry()
        {
            throw new Exception("Needs to be removed. Has beeb used only for testing purposes");
        }
        public void UpdateIncrementalDisplacements(double[] deltaU)
        {
            throw new Exception("Needs to be removed. Has beeb used only for testing purposes");
        }
        public List<double[]> GetStressFromElements(List<double[]> parametricCoordinatesVector)
        {
            List<double[]> StessVectorsList = new List<double[]>();
            StessVectorsList.Add(new double[] { 0.0, 0.0 });
            //double[,] E = CalculateStressStrainMatrix(Properties.YoungMod, Properties.PoissonRatio);
            //int count = parametricCoordinatesVector.Count;           
            //for (int i = 0; i < count; i++)
            //{
            //    double[] nodalParamCoord = parametricCoordinatesVector[i];
            //    Dictionary<string, double[]> localdN = CalculateShapeFunctionsLocalDerivatives(nodalParamCoord);
            //    double[,] J = CalculateJacobian(localdN);
            //    double[,] invJ = CalculateInverseJacobian(J).Item1;
            //    Dictionary<int, double[]> globaldN = CalculateShapeFunctionsGlobalDerivatives(localdN, invJ);
            //    double[,] B = CalculateBMatrix(globaldN);
            //    double[] strainVector = CalculateStrainsVector(B);
            //    double[] stressVector = CalculateStressVector(E, strainVector);
            //    StessVectorsList.Add(stressVector);
            //}
            return StessVectorsList;
        }
        public List<double[]> GetphysicalCoordinatesFromElements(List<double[]> parametricCoordinatesVector)
        {
            List<double[]> PositionVectorsList = new List<double[]>();
            //double[] xUpdated = UpdateNodalCoordinates(DisplacementVector);
            PositionVectorsList.Add(new double[] { 0.0, 0.0 });
            //int count = parametricCoordinatesVector.Count;
            //for (int i = 0; i < count; i++)
            //{
            //    double[] parametricCoordinatesVec = parametricCoordinatesVector[i];
            //    double[] positionVector = VectorOperations.MatrixVectorProduct(CalculateShapeFunctionMatrix(parametricCoordinatesVec[0], parametricCoordinatesVec[1], parametricCoordinatesVec[2]), xUpdated);
            //    PositionVectorsList.Add(positionVector);
            //}
            return PositionVectorsList;
        }
        public ContactNtS2Df(IElementProperties properties, Dictionary<int, INode> nodes)
        {
            Properties = properties;
            this.Nodes = nodes;
            ElementFreedomSignature[1] = new bool[] { true, true, false, false, false, false };
            ElementFreedomSignature[2] = new bool[] { true, true, false, false, false, false };
            ElementFreedomSignature[3] = new bool[] { true, true, false, false, false, false };
            DisplacementVector = new double[6];
            PenaltyFactor = properties.YoungMod * 10.0;
            TangentPenaltyFactor = properties.YoungMod * 10.0;
            counter = 1;
            mhid = 0.2;
            FrictionCoef = 0.2;
        }

        public void CalculateElementEASMatrices()
        {
            throw new Exception("This method is to be used only for EAS method elements");
        }
        public void InitializeElementEASParameters()
        {
            throw new Exception("Needs to be removed. Has beeb used only for testing purposes");
        }
        public void UpdateElementEASParameters(double[] solutionVector)
        {
            throw new Exception("Needs to be removed. Has beeb used only for testing purposes");
        }
        public void StoreElementFinalStepDisplacementVector(double[] solutionVector)
        {
            throw new Exception("Needs to be removed. Has beeb used only for testing purposes");
        }

        //private double[] CalculateNormalUnitVector(double detm)
        //{
        //    double X1 = Nodes[1].XCoordinate;
        //    double Y1 = Nodes[1].YCoordinate;
        //    double X2 = Nodes[2].XCoordinate;
        //    double Y2 = Nodes[2].YCoordinate;
        //    double[] vector = new double[] { Y2 - Y1, X2 - X1 };
        //    double scalarCoef = -1 / (2 * Math.Sqrt(detm));
        //    double[] normalUnitVec = VectorOperations.VectorScalarProductNew(vector, scalarCoef);
        //    return normalUnitVec;
        //}

        public Dictionary<int, INode> NodesAtFinalState()
        {
            throw new Exception("Method not implemenented");
        }
        public List<double[]> GetStressVector()
        {
            throw new Exception("Needs to be removed. Has beeb used only for testing purposes");
        }
        public List<double[]> GetStrainVector()
        {
            throw new Exception("Needs to be removed. Has beeb used only for testing purposes");
        }

        public List<double[]> GetGaussPointsInPhysicalSpace()
        {
            List<double[]> l = new List<double[]>();
            l.Add(new double[] { 0.0, 0.0 });
            return l;
        }
        public List<double[]> GetStressFromElementsNodes()
        {
            List<double[]> l = new List<double[]>();
            l.Add(new double[] { 0.0, 0.0, 0.0 });
            return l;
        }
        public List<double[]> GetStrainFromElementsNodes()
        {
            List<double[]> l = new List<double[]>();
            l.Add(new double[] { 0.0, 0.0, 0.0 });
            return l;
        }
        public double ClosestPointProjection()
        {
            double Xm1 = Nodes[1].XCoordinate + DisplacementVector[0];
            double Ym1 = Nodes[1].YCoordinate + DisplacementVector[1];
            double Xm2 = Nodes[2].XCoordinate + DisplacementVector[2];
            double Ym2 = Nodes[2].YCoordinate + DisplacementVector[3];
            double Xs = Nodes[3].XCoordinate + DisplacementVector[4];
            double Ys = Nodes[3].YCoordinate + DisplacementVector[5];

            double ksi1 = (2.0 * (Xs * (Xm2 - Xm1) + Ys * (Ym2 - Ym1)) - Math.Pow(Xm2, 2) - Math.Pow(Ym2, 2) + Math.Pow(Xm1, 2) + Math.Pow(Ym1, 2)) / (Math.Pow(Xm2 - Xm1, 2) + Math.Pow(Ym2 - Ym1, 2));
            return ksi1;
        }

        private Tuple<double[,], double[,]> CalculatePositionMatrix(double ksi1)
        {
            double N1 = 1.0 / 2.0 * (1.0 - ksi1);
            double N2 = 1.0 / 2.0 * (1.0 + ksi1);
            double dN1 = -1.0 / 2.0;
            double dN2 = 1.0 / 2.0;
            double[,] aMatrix = new double[,]
                {
                    { -N1 ,0.0 ,-N2 ,0.0 ,1.0 ,0.0 },
                    {0.0, -N1 , 0.0 ,-N2, 0.0, 1.0 }
                };

            double[,] daMatrix = new double[,]
                {
                    { -dN1 ,0.0 ,-dN2 ,0.0 ,0.0 ,0.0 },
                    {0.0, -dN1 , 0.0 ,-dN2, 0.0, 0.0 }
                };
            return new Tuple<double[,], double[,]>(aMatrix, daMatrix);
        }

        private Tuple<double[], double, double[], double[], double> SurfaceGeometry(double[,] daMatrix)
        {
            double Xm1 = Nodes[1].XCoordinate + DisplacementVector[0];
            double Ym1 = Nodes[1].YCoordinate + DisplacementVector[1];
            double Xm2 = Nodes[2].XCoordinate + DisplacementVector[2];
            double Ym2 = Nodes[2].YCoordinate + DisplacementVector[3];
            double Xs = Nodes[3].XCoordinate + DisplacementVector[4];
            double Ys = Nodes[3].YCoordinate + DisplacementVector[5];

            double[] xupd = new double[] { -Xm1, -Ym1, -Xm2, -Ym2, -Xs, -Ys };
            double[] surfaceVector = VectorOperations.MatrixVectorProduct(daMatrix, xupd);
            double detm = VectorOperations.VectorDotProduct(surfaceVector, surfaceVector);
            double m11 = 1.0 / detm;
            double[] vector = new double[] { Ym2 - Ym1, Xm1 - Xm2 };
            double scalarCoef = -1.0 / (2.0 * Math.Sqrt(detm));
            double[] normalUnitVec = VectorOperations.VectorScalarProductNew(vector, scalarCoef);
            double[] tVector = VectorOperations.VectorScalarProductNew(surfaceVector, 1.0 / Math.Sqrt(detm));
            return new Tuple<double[], double, double[], double[], double>(surfaceVector, m11, normalUnitVec, tVector, detm);
        }

        private double CalculateNormalGap(double[,] aMatrix, double[] n)
        {
            double[,] AT = MatrixOperations.Transpose(aMatrix);
            double[] AT_n = VectorOperations.MatrixVectorProduct(AT, n);
            double[] xupd = new double[] {
                Nodes[1].XCoordinate + DisplacementVector[0],
                Nodes[1].YCoordinate + DisplacementVector[1],
                Nodes[2].XCoordinate + DisplacementVector[2],
                Nodes[2].YCoordinate + DisplacementVector[3],
                Nodes[3].XCoordinate + DisplacementVector[4],
                Nodes[3].YCoordinate + DisplacementVector[5]
            };
            double normalGap = VectorOperations.VectorDotProduct(xupd, AT_n);
            return normalGap;
        }

        private double CalculateTangentialVelocity(double ksi1New, double ksi1Initial)
        {
            double deltaKsi = ksi1New - ksi1Initial;
            return deltaKsi;
        }

        private double CalculateTangentialTraction(double deltaKsi, double m11)
        {
            double tr1 = -TangentPenaltyFactor * m11 * deltaKsi;
            return tr1;
        }

        private double[,] CalculateMainStiffnessPart(double ksi1, double[] n)
        {
            double[,] mainStiffnessMatrix = new double[6, 6];
            double N1 = 1.0 / 2.0 * (1.0 - ksi1);
            double N2 = 1.0 / 2.0 * (1.0 + ksi1);
            Tuple<double[,], double[,]> positionMatrices = CalculatePositionMatrix(ksi1);
            double[,] A = positionMatrices.Item1;
            double[,] nxn = VectorOperations.VectorVectorTensorProduct(n, n);
            double[,] nxn_A = MatrixOperations.MatrixProduct(nxn, A);
            double[,] AT_nxn_A = MatrixOperations.MatrixProduct(MatrixOperations.Transpose(A), nxn_A);
            mainStiffnessMatrix = MatrixOperations.ScalarMatrixProductNew(PenaltyFactor, AT_nxn_A);
            //mainStiffnessMatrix[0, 0] = N1 * N1 * n[0] * n[0] * PenaltyFactor;
            //mainStiffnessMatrix[0, 1] = N1 * N1 * n[0] * n[1] * PenaltyFactor;
            //mainStiffnessMatrix[0, 2] = N1 * N2 * n[0] * n[0] * PenaltyFactor;
            //mainStiffnessMatrix[0, 3] = N1 * N1 * n[0] * n[1] * PenaltyFactor;
            //mainStiffnessMatrix[0, 4] = -N1 * n[0] * n[0] * PenaltyFactor;
            //mainStiffnessMatrix[0, 5] = -N1 * n[0] * n[1] * PenaltyFactor;

            //mainStiffnessMatrix[1, 0] = mainStiffnessMatrix[0, 1];
            //mainStiffnessMatrix[1, 1] = N1 * N1 * n[1] * n[1] * PenaltyFactor;
            //mainStiffnessMatrix[1, 2] = N1 * N2 * n[0] * n[1] * PenaltyFactor;
            //mainStiffnessMatrix[1, 3] = N1 * N2 * n[1] * n[1] * PenaltyFactor;
            //mainStiffnessMatrix[1, 4] = -N1 * n[0] * n[1] * PenaltyFactor;
            //mainStiffnessMatrix[1, 5] = -N1 * n[1] * n[1] * PenaltyFactor;

            //mainStiffnessMatrix[2, 0] = mainStiffnessMatrix[0, 2];
            //mainStiffnessMatrix[2, 1] = mainStiffnessMatrix[1, 2];
            //mainStiffnessMatrix[2, 2] = N2 * N2 * n[0] * n[0] * PenaltyFactor;
            //mainStiffnessMatrix[2, 3] = N2 * N2 * n[0] * n[1] * PenaltyFactor;
            //mainStiffnessMatrix[2, 4] = -N2 * n[0] * n[0] * PenaltyFactor;
            //mainStiffnessMatrix[2, 5] = -N2 * n[0] *n[1] * PenaltyFactor;

            //mainStiffnessMatrix[3, 0] = mainStiffnessMatrix[0, 3];
            //mainStiffnessMatrix[3, 1] = mainStiffnessMatrix[1, 3];
            //mainStiffnessMatrix[3, 2] = mainStiffnessMatrix[2, 3];
            //mainStiffnessMatrix[3, 3] = N2 * N2 * n[1] * n[1] * PenaltyFactor;
            //mainStiffnessMatrix[3, 4] = -N2 * n[0] * n[1] * PenaltyFactor;
            //mainStiffnessMatrix[3, 5] = -N2 * n[1] * n[1] * PenaltyFactor;

            //mainStiffnessMatrix[4, 0] = mainStiffnessMatrix[0, 4];
            //mainStiffnessMatrix[4, 1] = mainStiffnessMatrix[1, 4];
            //mainStiffnessMatrix[4, 2] = mainStiffnessMatrix[2, 4];
            //mainStiffnessMatrix[4, 3] = mainStiffnessMatrix[3, 4];
            //mainStiffnessMatrix[4, 4] = n[0] * n[0] * PenaltyFactor;
            //mainStiffnessMatrix[4, 5] = n[0] * n[1] * PenaltyFactor;

            //mainStiffnessMatrix[5, 0] = mainStiffnessMatrix[0, 5];
            //mainStiffnessMatrix[5, 1] = mainStiffnessMatrix[1, 5];
            //mainStiffnessMatrix[5, 2] = mainStiffnessMatrix[2, 5];
            //mainStiffnessMatrix[5, 3] = mainStiffnessMatrix[3, 5];
            //mainStiffnessMatrix[5, 4] = mainStiffnessMatrix[4, 5];
            //mainStiffnessMatrix[5, 5] = n[1] * n[1] * PenaltyFactor;
            return mainStiffnessMatrix;
        }

        private double[,] CalculateRotationalStiffnessPart(double[,] A, double[,] dA, double[] n, double ksi3, double m11, double[] dRho)
        {
            double coef = PenaltyFactor * ksi3 * m11;
            double[,] rotationalPart;
            double[,] n_x_dRho = VectorOperations.VectorVectorTensorProduct(n, dRho);
            double[,] dRho_x_n = VectorOperations.VectorVectorTensorProduct(dRho, n);
            double[,] firstTerm = MatrixOperations.MatrixProduct(
                                                                    MatrixOperations.Transpose(dA),
                                                                    MatrixOperations.MatrixProduct(n_x_dRho, A)
                                                                    );
            double[,] secondTerm = MatrixOperations.MatrixProduct(
                                                                    MatrixOperations.Transpose(A),
                                                                    MatrixOperations.MatrixProduct(dRho_x_n, dA)
                                                                    );
            rotationalPart = MatrixOperations.ScalarMatrixProductNew(
                                                                        coef,
                                                                        MatrixOperations.MatrixAddition(firstTerm, secondTerm)
                                                                        );
            //double[,] rotationalPart = new double[6, 6];
            //double N1 = 1 / 2 * (1 - ksi1);
            //double N2 = 1 / 2 * (1 + ksi1);
            //double coef = PenaltyFactor * ksi3 * m11;
            //rotationalPart[0, 0] = -coef * N1 * drho[0] * n[0];
            //rotationalPart[0, 1] = -coef * (N1 * drho[0] * n[1] / 2) - coef * (N1 * drho[1] * n[0] / 2);
            //rotationalPart[0, 2] = coef * (N1 * drho[0] * n[0] / 2) - coef * (N2 * drho[0] * n[0] / 2);
            //rotationalPart[0, 3] = coef * (N1 * drho[0] * n[1] / 2) - coef * (N2 * drho[1] * n[0] / 2);
            //rotationalPart[0, 4] = coef * (drho[0] * n[0] / 2);
            //rotationalPart[0, 5] = coef * (drho[1] * n[0] / 2) - coef * (N1 * drho[0] * n[1]);

            //rotationalPart[1, 0] = rotationalPart[0, 1];
            //rotationalPart[1, 1] = -coef * N1 * drho[1] * n[1];
            //rotationalPart[1, 2] = coef * (N1 * drho[1] * n[0] / 2) - coef * (N2 * drho[0] * n[1] / 2);
            //rotationalPart[1, 3] = coef * (N1 * drho[1] * n[1] / 2) - coef * (N2 * drho[1] * n[1] / 2);
            //rotationalPart[1, 4] = coef * drho[0] * n[1] / 2;
            //rotationalPart[1, 5] = coef * (drho[1] * n[1] / 2) - coef * (N1 * drho[1] * n[1]);
            return rotationalPart;
        }

        public double[,] CreateGlobalStiffnessMatrix()
        {
            //double ksi1 = ClosestPointProjection();
            if (counter == 1) { Ksi1Initial = Ksi1Current; }
            counter = counter + 1;
            if (Math.Abs(Ksi1Current) <= 1.05)
            {
                Tuple<double[,], double[,]> positionMatrices = CalculatePositionMatrix(Ksi1Current);
                double[,] aMatrix = positionMatrices.Item1;
                double[,] daMatrix = positionMatrices.Item2;

                Tuple<double[], double, double[], double[], double> surfaceCharacteristics = SurfaceGeometry(daMatrix);
                double m11 = surfaceCharacteristics.Item2;
                double[] dRho = surfaceCharacteristics.Item1;
                double[] n = surfaceCharacteristics.Item3;
                double[] tVector = surfaceCharacteristics.Item4;

                double detM = surfaceCharacteristics.Item5;

                double ksi3 = CalculateNormalGap(aMatrix, n);
                if (ksi3 <= 0)
                {
                    double[,] sN = CalculateMainStiffnessPart(Ksi1Current, n);
                    double deltaKsi = CalculateTangentialVelocity(Ksi1Current, Ksi1Initial);
                    double Tr1 = CalculateTangentialTraction(deltaKsi, detM);
                    double phi = Math.Sqrt(Tr1 * Tr1 * m11) - FrictionCoef * PenaltyFactor * Math.Abs(ksi3);
                    if (phi <= 0.0)
                    {
                        double T1 = Tr1;
                        double[,] sT1 = MatrixOperations.ScalarMatrixProductNew(TangentPenaltyFactor,
                            MatrixOperations.MatrixProduct(
                                MatrixOperations.Transpose(aMatrix),
                                MatrixOperations.MatrixProduct(
                                VectorOperations.VectorVectorTensorProduct(tVector, tVector), aMatrix)));
                        double[,] sT2 = MatrixOperations.ScalarMatrixProductNew(T1 * m11,
                            MatrixOperations.MatrixProduct(
                                MatrixOperations.Transpose(daMatrix),
                                MatrixOperations.MatrixProduct(
                                VectorOperations.VectorVectorTensorProduct(tVector, tVector), aMatrix)));
                        double[,] sT = MatrixOperations.MatrixAddition(
                            MatrixOperations.ScalarMatrixProductNew(-1.0, sT1),
                            MatrixOperations.MatrixAddition(sT2,
                            MatrixOperations.Transpose(sT2)));
                        double[,] globalStiffnessMatrix = MatrixOperations.MatrixAddition(sN, sT);
                        return globalStiffnessMatrix;
                    }
                    else
                    {
                        double T1 = (Tr1 / Math.Abs(Tr1)) * mhid * PenaltyFactor * Math.Abs(ksi3) * Math.Sqrt(detM);
                        double[,] sT1 = MatrixOperations.ScalarMatrixProductNew(mhid * PenaltyFactor * (Tr1 / Math.Abs(Tr1)),
                            MatrixOperations.MatrixProduct(
                                MatrixOperations.Transpose(aMatrix),
                                MatrixOperations.MatrixProduct(
                                VectorOperations.VectorVectorTensorProduct(tVector, n), aMatrix)));
                        double[,] sT2 = MatrixOperations.ScalarMatrixProductNew(mhid * PenaltyFactor * Math.Abs(ksi3) * (Tr1 / Math.Abs(Tr1)) * Math.Sqrt(m11),
                            MatrixOperations.MatrixProduct(
                                MatrixOperations.Transpose(daMatrix),
                                MatrixOperations.MatrixProduct(
                                VectorOperations.VectorVectorTensorProduct(tVector, tVector), aMatrix)));
                        double[,] sT = MatrixOperations.MatrixAddition(
                            MatrixOperations.ScalarMatrixProductNew(-1.0, sT1),
                            MatrixOperations.MatrixAddition(sT2,
                            MatrixOperations.Transpose(sT2)));
                        double[,] globalStiffnessMatrix = MatrixOperations.MatrixAddition(sN, sT);
                        return globalStiffnessMatrix;
                    }
                    //double[,] rotationalPart = CalculateRotationalStiffnessPart(aMatrix, daMatrix, n, ksi3, m11, dRho);
                    //double[,] globalStiffnessMatrix = MatrixOperations.MatrixAddition(mainPart, rotationalPart);


                }
                else
                {
                    double[,] globalStifnessMatrix = new double[6, 6];
                    return globalStifnessMatrix;
                }
            }
            else
            {
                double[,] globalStifnessMatrix = new double[6, 6];
                return globalStifnessMatrix;
            }
        }

        public double[] CreateInternalGlobalForcesVector()
        {
            //double ksi1 = ClosestPointProjection();
            Ksi1Current = ClosestPointProjection();
            if (counter == 1) { Ksi1Initial = Ksi1Current; }
            counter = counter + 1;
            if (Math.Abs(Ksi1Current) <= 1.05)
            {
                Tuple<double[,], double[,]> positionMatrices = CalculatePositionMatrix(Ksi1Current);
                double[,] aMatrix = positionMatrices.Item1;
                double[,] daMatrix = positionMatrices.Item2;

                Tuple<double[], double, double[], double[], double> surfaceCharacteristics = SurfaceGeometry(daMatrix);
                double m11 = surfaceCharacteristics.Item2;
                double[] n = surfaceCharacteristics.Item3;
                double[] tVector = surfaceCharacteristics.Item4;

                double detM = surfaceCharacteristics.Item5;

                double ksi3 = CalculateNormalGap(aMatrix, n);
                if (ksi3 <= 0)
                {
                    double[,] AT = MatrixOperations.Transpose(aMatrix);
                    double[] AT_n = VectorOperations.MatrixVectorProduct(AT, n);
                    double[] AT_t = VectorOperations.MatrixVectorProduct(AT, tVector);

                    double deltaKsi = CalculateTangentialVelocity(Ksi1Current, Ksi1Initial);
                    //Console.WriteLine(Ksi1Initial);
                    double Tr1 = CalculateTangentialTraction(deltaKsi, detM);
                    double phi = Math.Sqrt(Tr1 * Tr1 * m11) - FrictionCoef * PenaltyFactor * Math.Abs(ksi3);
                    if (phi <= 0.0)
                    {
                        double T1 = Tr1;
                        double[] internalGlobalForcesVector = VectorOperations.VectorVectorAddition(
                                                VectorOperations.VectorScalarProductNew(AT_n, PenaltyFactor * ksi3),
                                                VectorOperations.VectorScalarProductNew(AT_t, T1 * Math.Sqrt(m11)));
                        return internalGlobalForcesVector;
                    }
                    else
                    {
                        double T1 = (Tr1 / Math.Abs(Tr1)) * mhid * PenaltyFactor * Math.Abs(ksi3) * Math.Sqrt(detM);
                        double[] internalGlobalForcesVector = VectorOperations.VectorVectorAddition(
                                                VectorOperations.VectorScalarProductNew(AT_n, PenaltyFactor * ksi3),
                                                VectorOperations.VectorScalarProductNew(AT_t, T1 * Math.Sqrt(m11)));
                        return internalGlobalForcesVector;
                    }
                }
                else
                {
                    double[] internalGlobalForcesVector = new double[6];
                    return internalGlobalForcesVector;
                }
            }
            else
            {
                double[] internalGlobalForcesVector = new double[6];
                return internalGlobalForcesVector;
            }
        }

        public double[,] CreateMassMatrix()
        {
            return new double[6, 6];
        }

        public double[,] CreateDampingMatrix()
        {
            return new double[6, 6];
        }
    }
}


