using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GFEC
{
    class ContactNtS3Df : IElement
    {
        public Dictionary<int, INode> Nodes { get; }
        public IElementProperties Properties { get; set; }
        public Dictionary<int, bool[]> ElementFreedomSignature { get; } = new Dictionary<int, bool[]>();
        public List<int> ElementFreedomList { get; set; }
        public double[] DisplacementVector { get ; set ; }
        private double[] DisplacementVectorPrevious { get; set; }
        public double[] AccelerationVector { get; set; }
        private double PenaltyFactor { get; set; }
        private double PenaltyFactorT { get; set; }
        double[] lastKsiVector;
        private bool intForceSwitch, stiffMatrixSwitch;

        private double[] ksiVectorPrevious { get; set;}
        private double[] ksiVectorCurrent { get; set; }
        private int counter;
        private double FrictionCoef { get; set; }
        private double[,] metricTensorPrevious;
        private double[] tTractionPrevious;
        private List<double[]> dRhoPrevious;

        public ContactNtS3Df(IElementProperties properties, Dictionary<int, INode> nodes)
        {
            Properties = properties;
            this.Nodes = nodes;
            ElementFreedomSignature[1] = new bool[] { true, true, true, false, false, false };
            ElementFreedomSignature[2] = new bool[] { true, true, true, false, false, false };
            ElementFreedomSignature[3] = new bool[] { true, true, true, false, false, false };
            ElementFreedomSignature[4] = new bool[] { true, true, true, false, false, false };
            ElementFreedomSignature[5] = new bool[] { true, true, true, false, false, false };
            DisplacementVector = new double[15];
            DisplacementVectorPrevious = new double[15];
            PenaltyFactor = properties.YoungMod * 2.0 ;
            PenaltyFactorT = properties.YoungMod * 2.0 ;
            lastKsiVector = new double[2];
            counter = 1;
            intForceSwitch = false;
            stiffMatrixSwitch = false;
        }

        public Dictionary<int, INode> NodesAtFinalState()
        {
            Dictionary<int, INode> finalNodes = new Dictionary<int, INode>();
            finalNodes[1] = new Node(Nodes[1].XCoordinate + DisplacementVector[0], Nodes[1].YCoordinate + DisplacementVector[1], Nodes[1].ZCoordinate + DisplacementVector[2]);
            finalNodes[2] = new Node(Nodes[2].XCoordinate + DisplacementVector[3], Nodes[2].YCoordinate + DisplacementVector[4], Nodes[2].ZCoordinate + DisplacementVector[5]);
            finalNodes[3] = new Node(Nodes[3].XCoordinate + DisplacementVector[6], Nodes[3].YCoordinate + DisplacementVector[7], Nodes[3].ZCoordinate + DisplacementVector[8]);
            finalNodes[4] = new Node(Nodes[4].XCoordinate + DisplacementVector[9], Nodes[4].YCoordinate + DisplacementVector[10], Nodes[4].ZCoordinate + DisplacementVector[11]);
            finalNodes[5] = new Node(Nodes[5].XCoordinate + DisplacementVector[12], Nodes[5].YCoordinate + DisplacementVector[13], Nodes[5].ZCoordinate + DisplacementVector[14]);
            return finalNodes;
        }

        private double[] xUpdatedVector()
        {
            Dictionary<int, INode> finalNodes = NodesAtFinalState();
            double[] xVectorUpdated = new double[15];
            for (int i = 0; i < 5; i++)
            {
                xVectorUpdated[3 * i] = finalNodes[i + 1].XCoordinate;
                xVectorUpdated[3 * i + 1] = finalNodes[i + 1].YCoordinate;
                xVectorUpdated[3 * i + 2] = finalNodes[i + 1].ZCoordinate;
            }
            return xVectorUpdated;
        }

        private double[] xUpdatedVector(double[] displacementVector)
        {
            Dictionary<int, INode> finalNodes = new Dictionary<int, INode>();
            finalNodes[1] = new Node(Nodes[1].XCoordinate + displacementVector[0], Nodes[1].YCoordinate + displacementVector[1], Nodes[1].ZCoordinate + displacementVector[2]);
            finalNodes[2] = new Node(Nodes[2].XCoordinate + displacementVector[3], Nodes[2].YCoordinate + displacementVector[4], Nodes[2].ZCoordinate + displacementVector[5]);
            finalNodes[3] = new Node(Nodes[3].XCoordinate + displacementVector[6], Nodes[3].YCoordinate + displacementVector[7], Nodes[3].ZCoordinate + displacementVector[8]);
            finalNodes[4] = new Node(Nodes[4].XCoordinate + displacementVector[9], Nodes[4].YCoordinate + displacementVector[10], Nodes[4].ZCoordinate + displacementVector[11]);
            finalNodes[5] = new Node(Nodes[5].XCoordinate + displacementVector[12], Nodes[5].YCoordinate + displacementVector[13], Nodes[5].ZCoordinate + displacementVector[14]);

            double[] xVectorUpdated = new double[15];
            for (int i = 0; i < 5; i++)
            {
                xVectorUpdated[3 * i] = finalNodes[i + 1].XCoordinate;
                xVectorUpdated[3 * i + 1] = finalNodes[i + 1].YCoordinate;
                xVectorUpdated[3 * i + 2] = finalNodes[i + 1].ZCoordinate;
            }
            return xVectorUpdated;
        }

        private Tuple<double[,], double[,], double[,]> CalculatePositionMatrix(double ksi1, double ksi2)
        {
            double N1 = 1.0 / 4.0 * (1.0 - ksi1) * (1.0 - ksi2);
            double N2 = 1.0 / 4.0 * (1.0 + ksi1) * (1.0 - ksi2);
            double N3 = 1.0 / 4.0 * (1.0 + ksi1) * (1.0 + ksi2);
            double N4 = 1.0 / 4.0 * (1.0 - ksi1) * (1.0 + ksi2);

            double dN11 = -1.0 / 4.0 * (1.0 - ksi2);
            double dN21 = 1.0 / 4.0 * (1.0 - ksi2);
            double dN31 = 1.0 / 4.0 * (1.0 + ksi2);
            double dN41 = -1.0 / 4.0 * (1.0 + ksi2);

            double dN12 = -1.0 / 4.0 * (1.0 - ksi1);
            double dN22 = -1.0 / 4.0 * (1.0 + ksi1);
            double dN32 = 1.0 / 4.0 * (1.0 + ksi1);
            double dN42 = 1.0 / 4.0 * (1.0 - ksi1);

            double[,] aMatrix = new double[,]
                {
                    { -N1 ,0.0, 0.0 ,-N2 ,0.0 ,0.0 ,-N3, 0.0, 0.0, -N4, 0.0, 0.0, 1.0, 0.0, 0.0 },
                    { 0.0, -N1 ,0.0, 0.0 ,-N2 ,0.0 ,0.0 ,-N3, 0.0, 0.0, -N4, 0.0, 0.0, 1.0, 0.0 },
                    { 0.0, 0.0, -N1 ,0.0, 0.0 ,-N2 ,0.0 ,0.0 ,-N3, 0.0, 0.0, -N4, 0.0, 0.0, 1.0 }
                };

            double[,] da1Matrix = new double[,]
                {
                    { -dN11 ,0.0, 0.0 ,-dN21 ,0.0 ,0.0 ,-dN31, 0.0, 0.0, -dN41, 0.0, 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, -dN11 ,0.0, 0.0 ,-dN21 ,0.0 ,0.0 ,-dN31, 0.0, 0.0, -dN41, 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, 0.0, -dN11 ,0.0, 0.0 ,-dN21 ,0.0 ,0.0 ,-dN31, 0.0, 0.0, -dN41, 0.0, 0.0, 0.0 }
                };

            double[,] da2Matrix = new double[,]
                {
                    { -dN12 ,0.0, 0.0 ,-dN22 ,0.0 ,0.0 ,-dN32, 0.0, 0.0, -dN42, 0.0, 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, -dN12 ,0.0, 0.0 ,-dN22 ,0.0 ,0.0 ,-dN32, 0.0, 0.0, -dN42, 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, 0.0, -dN12 ,0.0, 0.0 ,-dN22 ,0.0 ,0.0 ,-dN32, 0.0, 0.0, -dN42, 0.0, 0.0, 0.0 }
                };
            return new Tuple<double[,], double[,], double[,]>(aMatrix, da1Matrix, da2Matrix);
        }

        private List<double[]> SurfaceVectors(double ksi1, double ksi2)
        {
            Tuple<double[,], double[,], double[,]> positionsMatrices = CalculatePositionMatrix(ksi1, ksi2);
            double[,] da1 = positionsMatrices.Item2;
            double[,] da2 = positionsMatrices.Item3;
            double[] xUpdated = xUpdatedVector();

            List<double[]> dRho = new List<double[]>();
            dRho.Add(VectorOperations.VectorScalarProductNew(
                VectorOperations.MatrixVectorProduct(da1, xUpdated), -1.0));
            dRho.Add(VectorOperations.VectorScalarProductNew(
                VectorOperations.MatrixVectorProduct(da2, xUpdated), -1.0));
            return dRho;
        }

        private double[,] MetricTensor(List<double[]> dRho)
        {
            double[,] m = new double[2, 2];
            m[0, 0] = VectorOperations.VectorDotProduct(dRho[0], dRho[0]);
            m[0, 1] = VectorOperations.VectorDotProduct(dRho[0], dRho[1]);
            m[1, 0] = VectorOperations.VectorDotProduct(dRho[1], dRho[0]);
            m[1, 1] = VectorOperations.VectorDotProduct(dRho[1], dRho[1]);
            return m;
        }

        private double MetricTensorDet(double[,] m)
        {
            double detm = m[0, 0] * m[1, 1] - m[1, 0] * m[0, 1];

            return detm;
        }

        private double[,] InverseMetricTensor(double[,] m)
        {
            double detm = MetricTensorDet(m);

            double[,] mInv = MatrixOperations.ScalarMatrixProductNew(1.0 / detm,
                new double[,] {
                    { m[1, 1], -m[0, 1] },
                    {-m[1,0], m[0,0] }
                });
            return mInv;
        }
        private double[] NormalVector(double[,] m, List<double[]> dRho)
        {
            double detm = MetricTensorDet(m);
            double[] drhoBydrho = VectorOperations.VectorCrossProduct(dRho[0], dRho[1]);
            double[] n = VectorOperations.VectorScalarProductNew(drhoBydrho, 1.0 / (Math.Sqrt(detm)));
            return n;
        }

        private double CalculatePenetration(double ksi1, double ksi2)
        {
            double[] xUpdated = xUpdatedVector();
            List<double[]> dRho = SurfaceVectors(ksi1, ksi2);
            double[,] m = MetricTensor(dRho);
            double[] n = NormalVector(m, dRho);
            Tuple<double[,], double[,], double[,]> aMatrices = CalculatePositionMatrix(ksi1, ksi2);
            double[,] a = aMatrices.Item1;
            double[,] aT = MatrixOperations.Transpose(a);
            double ksi3 = VectorOperations.VectorDotProduct(
                xUpdated, VectorOperations.MatrixVectorProduct(
                    aT, n));
            return ksi3;
        }

        private double[] Calculate_f(List<double[]> dRho, double[,] aMatrix, double[] xUpdated)
        {
            double[] f = new double[2];
            f[0] = VectorOperations.VectorDotProduct(
                dRho[0], VectorOperations.MatrixVectorProduct(
                    aMatrix, xUpdated));
            f[1] = VectorOperations.VectorDotProduct(
                dRho[1], VectorOperations.MatrixVectorProduct(
                    aMatrix, xUpdated));
            return f;
        }

        private double Calculate_e(double[,] aMatrix, double[] xUpdated)
        {

            double[] dRho12 = new double[] {
                0.25*(xUpdated[0] - xUpdated[3] + xUpdated[6] - xUpdated[9]),
                0.25*(xUpdated[1] - xUpdated[4] + xUpdated[7] - xUpdated[10]),
                0.25*(xUpdated[2] - xUpdated[5] + xUpdated[8] - xUpdated[11])
            };

            double e = VectorOperations.VectorDotProduct(
                dRho12, VectorOperations.MatrixVectorProduct(
                    aMatrix, xUpdated));
            return e;
        }

        private double[] CalculateDeltaKsi(double detm, double[,] mTensor, double[] fVector, double e)
        {
            double scalar = 1.0 / (detm - Math.Pow(e, 2) + 2.0 * e * mTensor[0, 1]);
            double[,] matrix = new double[,]
            {
                {mTensor[1,1], e-mTensor[0,1] },
                {e-mTensor[1,0], mTensor[0,0] }
            };

            double[] deltaKsi = VectorOperations.VectorScalarProductNew(
                VectorOperations.MatrixVectorProduct(matrix, fVector),
                scalar);

            return deltaKsi;
        }

        private double[] Project(double[] ksiVector)
        {
            double[] xUpdated = xUpdatedVector();
            Tuple<double[,], double[,], double[,]> aMatrices = CalculatePositionMatrix(ksiVector[0], ksiVector[1]);
            List<double[]> dRho = SurfaceVectors(ksiVector[0], ksiVector[1]);
            double[] f = Calculate_f(dRho, aMatrices.Item1, xUpdated);
            double e = Calculate_e(aMatrices.Item2, xUpdated);
            double[,] m = MetricTensor(dRho);
            double detm = MetricTensorDet(m);
            double[] deltaKsi = CalculateDeltaKsi(detm, m, f, e);
            ksiVector = VectorOperations.VectorVectorAddition(ksiVector, deltaKsi);
            return ksiVector;
        }

        private double CalculatePenetration(double[] normalVector, double[,] aMatrix, double[] xUpdated)
        {
            double ksi3 = VectorOperations.VectorDotProduct(
                xUpdated, VectorOperations.MatrixVectorProduct(
                    MatrixOperations.Transpose(aMatrix), normalVector));
            return ksi3;
        }

        private double[] CalculateTangentialTraction(double phi, double N, double[] trialT)
        {
            double[] tractionT = new double[2];
            if (phi <= 0)
            {
                tractionT = trialT;
            }
            else
            {
                tractionT[0] = FrictionCoef * N * trialT[0] / VectorOperations.VectorNorm2(trialT);
                tractionT[1] = FrictionCoef * N * trialT[1] / VectorOperations.VectorNorm2(trialT);
            }
            return tractionT;
        }

        private double CalculatePhi(double[] trialT, double[,] metricTensor, double N)
        {
            double phi = 0;
            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    phi = phi + Math.Sqrt(trialT[i] * trialT[j] * metricTensor[i, j]);
                }
            }
            phi = phi + FrictionCoef * N;
            return phi;
        }

        private double[] CalculateRhoC(double ksi1, double ksi2, double[] xUpdated)
        {
            double N1 = 1.0 / 4.0 * (1.0 - ksi1) * (1.0 - ksi2);
            double N2 = 1.0 / 4.0 * (1.0 + ksi1) * (1.0 - ksi2);
            double N3 = 1.0 / 4.0 * (1.0 + ksi1) * (1.0 + ksi2);
            double N4 = 1.0 / 4.0 * (1.0 - ksi1) * (1.0 + ksi2);

            double[] x1 = new double[] { xUpdated[0], xUpdated[1], xUpdated[2] };
            double[] x2 = new double[] { xUpdated[3], xUpdated[4], xUpdated[5] };
            double[] x3 = new double[] { xUpdated[6], xUpdated[7], xUpdated[8] };
            double[] x4 = new double[] { xUpdated[9], xUpdated[10], xUpdated[11] };

            double[] N1x1 = VectorOperations.VectorScalarProductNew(x1, N1);
            double[] N2x2 = VectorOperations.VectorScalarProductNew(x2, N2);
            double[] N3x3 = VectorOperations.VectorScalarProductNew(x3, N3);
            double[] N4x4 = VectorOperations.VectorScalarProductNew(x4, N4);

            double[] rhoC = VectorOperations.VectorVectorAddition(
                N1x1, VectorOperations.VectorVectorAddition(
                    N2x2, VectorOperations.VectorVectorAddition(
                        N3x3, N4x4)));
            return rhoC;
        }

        private double[] CalculateUC(double ksi1, double ksi2, double[] displacementVector)
        {
            double N1 = 1.0 / 4.0 * (1.0 - ksi1) * (1.0 - ksi2);
            double N2 = 1.0 / 4.0 * (1.0 + ksi1) * (1.0 - ksi2);
            double N3 = 1.0 / 4.0 * (1.0 + ksi1) * (1.0 + ksi2);
            double N4 = 1.0 / 4.0 * (1.0 - ksi1) * (1.0 + ksi2);

            double[] u1 = new double[] { displacementVector[0], displacementVector[1], displacementVector[2] };
            double[] u2 = new double[] { displacementVector[3], displacementVector[4], displacementVector[5] };
            double[] u3 = new double[] { displacementVector[6], displacementVector[7], displacementVector[8] };
            double[] u4 = new double[] { displacementVector[9], displacementVector[10], displacementVector[11] };

            double[] N1u1 = VectorOperations.VectorScalarProductNew(u1, N1);
            double[] N2u2 = VectorOperations.VectorScalarProductNew(u2, N2);
            double[] N3u3 = VectorOperations.VectorScalarProductNew(u3, N3);
            double[] N4u4 = VectorOperations.VectorScalarProductNew(u4, N4);

            double[] UC = VectorOperations.VectorVectorAddition(
                N1u1, VectorOperations.VectorVectorAddition(
                    N2u2, VectorOperations.VectorVectorAddition(
                        N3u3, N4u4)));
            return UC;
        }
        private double[] CalculateTrialTraction(double[] ksiVectorCurrent, double[] ksiVectorPrevious, double[] tPrevious, double[,] aPrevious, List<double[]> dRhoPrevious, List<double[]> dRhoCurrent)
        {
            double[] xUpdatedPrevious = xUpdatedVector(DisplacementVectorPrevious);
            double[] UCPrevious = CalculateUC(ksiVectorPrevious[0], ksiVectorPrevious[1], DisplacementVectorPrevious);
            double[] rhoCPrevious = CalculateRhoC(ksiVectorPrevious[0], ksiVectorPrevious[1], xUpdatedPrevious);
            double[] rhoCCurrent = CalculateRhoC(ksiVectorCurrent[0], ksiVectorCurrent[1], xUpdatedVector());
            double[] deltaRho = VectorOperations.VectorVectorSubtraction(rhoCCurrent, VectorOperations.VectorVectorAddition(rhoCPrevious, UCPrevious));
            double[] tTrial = new double[2];
            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    for (int k = 0; k < 2; k++)
                    {
                        tTrial[i] = tPrevious[k] * aPrevious[k, j] * (VectorOperations.VectorDotProduct(dRhoPrevious[j], dRhoCurrent[i]))
                            - PenaltyFactorT * (VectorOperations.VectorDotProduct(deltaRho, dRhoCurrent[i]));
                    }
                    
                }
            }
            return tTrial;
        }

        

        private double[,] Calculate_Kn(double penaltyFactor, double[] normalVector, double[,] A, List<double[,]> dA, List<double[]> dRho, double ksi3, double[,] alpha)
        {
            double normalTraction = penaltyFactor * ksi3;
            double[,] AT = MatrixOperations.Transpose(A);
            List<double[,]> dAT = new List<double[,]>() { MatrixOperations.Transpose(dA[0]), MatrixOperations.Transpose(dA[1]) };
            double[,] Kn;

            double[,] nxn = VectorOperations.VectorVectorTensorProduct(normalVector, normalVector);
            double[,] nxna = MatrixOperations.MatrixProduct(nxn, A);
            double[,] aTnxna = MatrixOperations.MatrixProduct(AT, nxna);
            double[,] Kn1 = MatrixOperations.ScalarMatrixProductNew(penaltyFactor, aTnxna);

            Kn = Kn1;

            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    double[,] Kn2a = 
                        MatrixOperations.ScalarMatrixProductNew
                        (
                            normalTraction * alpha[i, j], MatrixOperations.MatrixProduct
                            (
                                dAT[j], MatrixOperations.MatrixProduct
                                (
                                    VectorOperations.VectorVectorTensorProduct(normalVector,dRho[i]), A
                                    )
                                )
                            );
                    Kn = MatrixOperations.MatrixAddition(Kn, Kn2a);

                    double[,] Kn2b =
                        MatrixOperations.ScalarMatrixProductNew
                        (
                            normalTraction * alpha[i, j], MatrixOperations.MatrixProduct
                            (
                                AT, MatrixOperations.MatrixProduct
                                (
                                    VectorOperations.VectorVectorTensorProduct(dRho[j], normalVector), dA[i]
                                    )
                                )
                            );
                    Kn = MatrixOperations.MatrixAddition(Kn, Kn2b);

                    //double[,] Kn3 =
                    //   MatrixOperations.ScalarMatrixProductNew
                    //   (
                    //       normalTraction * h[i, j], MatrixOperations.MatrixProduct
                    //       (
                    //           AT, MatrixOperations.MatrixProduct
                    //           (
                    //               VectorOperations.VectorVectorTensorProduct(dRho[i], dRho[j]), A
                    //               )
                    //           )
                    //       );
                    //Kn = MatrixOperations.MatrixSubtraction(Kn, Kn3);
                }
            }
           
            return Kn;
        }

        public double[,] Calculate_Kst(double et, double[,] A, List<double[,]> dA, double[,] alpha, List<double[]> dRho, double[] T, double[] normalVector)
        {
            double[,] Kst = new double[15, 15];
            double[,] AT = MatrixOperations.Transpose(A);
            List<double[,]> dAT = new List<double[,]>() { MatrixOperations.Transpose(dA[0]), MatrixOperations.Transpose(dA[1]) };
            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    double[,] Kst1 =
                               MatrixOperations.ScalarMatrixProductNew
                               (
                                   -et * alpha[i, j], MatrixOperations.MatrixProduct
                                   (
                               AT, MatrixOperations.MatrixProduct
                               (
                                   VectorOperations.VectorVectorTensorProduct(dRho[i], dRho[j]), A
                                   )
                               )
                           );

                    //double[,] Kst3 =
                    //            MatrixOperations.ScalarMatrixProductNew
                    //            (
                    //                T[i] * h[i, j], MatrixOperations.MatrixProduct
                    //                (
                    //                    AT, MatrixOperations.MatrixProduct
                    //                    (
                    //                        MatrixOperations.MatrixAddition
                    //                        (
                    //                            VectorOperations.VectorVectorTensorProduct(dRho[j], normalVector), VectorOperations.VectorVectorTensorProduct(normalVector, dRho[j])
                    //                            ), A
                    //                        )
                    //                    )
                    //                );

                    Kst = Kst1;//MatrixOperations.MatrixAddition(Kst, MatrixOperations.MatrixAddition(Kst1, Kst3));

                    for (int k = 0; k < 2; k++)
                    {
                        for (int l = 0; l < 2; l++)
                        {
                           

                            double[,] Kst2a =
                                MatrixOperations.ScalarMatrixProductNew
                                (
                                    T[i] * alpha[i, l] * alpha[j, k], MatrixOperations.MatrixProduct
                                    (
                                        AT, MatrixOperations.MatrixProduct
                                        (
                                            VectorOperations.VectorVectorTensorProduct(dRho[k], dRho[l]), dA[j]
                                            )
                                        )
                                    );

                            double[,] Kst2b =
                                MatrixOperations.ScalarMatrixProductNew
                                (
                                    T[i] * alpha[i, k] * alpha[j, l], MatrixOperations.MatrixProduct
                                    (
                                        dAT[j], MatrixOperations.MatrixProduct
                                        (
                                            VectorOperations.VectorVectorTensorProduct(dRho[k], dRho[l]), A
                                            )
                                        )
                                    );

                            Kst = MatrixOperations.MatrixAddition(Kst, MatrixOperations.MatrixAddition(Kst2a, Kst2b));
                        }
                    }
                }
            }
            return Kst;
        }

        public double[,] Calculate_Ksl(double mhi, double ksi3, double et, double en, double[,] A, List<double[,]> dA, double[,] alpha, List<double[]> dRho, double[] T, double[] normalVector)
        {
            double[,] Ksl = new double[15, 15];
            double[,] AT = MatrixOperations.Transpose(A);
            double normT = VectorOperations.VectorNorm2(T);
            double normalTraction = en * ksi3;
            List<double[,]> dAT = new List<double[,]>() { MatrixOperations.Transpose(dA[0]), MatrixOperations.Transpose(dA[1]) };

            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    double[,] Ksl1 =
                               MatrixOperations.ScalarMatrixProductNew
                               (
                                   -en * mhi * (T[i]/normT) * alpha[i, j], MatrixOperations.MatrixProduct
                                   (
                               AT, MatrixOperations.MatrixProduct
                               (
                                   VectorOperations.VectorVectorTensorProduct(dRho[j], normalVector), A
                                   )
                               )
                           );

                    double[,] Ksl2 =
                               MatrixOperations.ScalarMatrixProductNew
                               (
                                   -et * mhi * (Math.Abs(normalTraction) / normT) * alpha[i, j], MatrixOperations.MatrixProduct
                                   (
                               AT, MatrixOperations.MatrixProduct
                               (
                                   VectorOperations.VectorVectorTensorProduct(dRho[i], dRho[j]), A
                                   )
                               )
                           );

                    Ksl = MatrixOperations.MatrixAddition(Ksl1, Ksl2);

                    //double[,] Ksl5a =
                    //           MatrixOperations.ScalarMatrixProductNew
                    //           (
                    //               mhi * (Math.Abs(normalTraction) * T[i] / normT) * h[i, j], MatrixOperations.MatrixProduct
                    //               (
                    //           AT, MatrixOperations.MatrixProduct
                    //           (
                    //               VectorOperations.VectorVectorTensorProduct(dRho[j], normalVector), A
                    //               )
                    //           )
                    //       );

                    //double[,] Ksl5b =
                    //           MatrixOperations.ScalarMatrixProductNew
                    //           (
                    //               mhi * (Math.Abs(normalTraction) * T[i] / normT) * h[i, j], MatrixOperations.MatrixProduct
                    //               (
                    //           AT, MatrixOperations.MatrixProduct
                    //           (
                    //               VectorOperations.VectorVectorTensorProduct(normalVector, dRho[j]), A
                    //               )
                    //           )
                    //       );

                    for (int k = 0; k < 2; k++)
                    {
                        for (int l = 0; l < 2; l++)
                        {
                            double[,] Ksl3 =
                               MatrixOperations.ScalarMatrixProductNew
                               (
                                   et * mhi * (Math.Abs(normalTraction) * T[i] * T[j] / Math.Pow(normT,3)) * alpha[i, k] * alpha[j,l], MatrixOperations.MatrixProduct
                                   (
                               AT, MatrixOperations.MatrixProduct
                               (
                                   VectorOperations.VectorVectorTensorProduct(dRho[k], dRho[l]), A
                                   )
                               )
                           );

                            double[,] Ksl4a =
                               MatrixOperations.ScalarMatrixProductNew
                               (
                                   mhi * (Math.Abs(normalTraction) * T[i] / normT) * alpha[i, l] * alpha[j, k], MatrixOperations.MatrixProduct
                                   (
                               AT, MatrixOperations.MatrixProduct
                               (
                                   VectorOperations.VectorVectorTensorProduct(dRho[k], dRho[l]), dA[j]
                                   )
                               )
                           );

                            double[,] Ksl4b =
                               MatrixOperations.ScalarMatrixProductNew
                               (
                                   mhi * (Math.Abs(normalTraction) * T[i] / normT) * alpha[i, k] * alpha[j, l], MatrixOperations.MatrixProduct
                                   (
                               dAT[j], MatrixOperations.MatrixProduct
                               (
                                   VectorOperations.VectorVectorTensorProduct(dRho[k], dRho[l]), A
                                   )
                               )
                           );

                            Ksl = MatrixOperations.MatrixAddition(Ksl, Ksl3);
                            Ksl = MatrixOperations.MatrixAddition(Ksl, Ksl4a);
                            Ksl = MatrixOperations.MatrixAddition(Ksl, Ksl4b);
                        }
                    }
                }
            }
            return Ksl;
        }

        private void SetAsPreviousValues(bool intFSwitch, bool stifMSwitch, double[] currentDispVector, double[] ksiVectorCurrent, List<double[]> dRhoCurrent, double[,] metricTensorCurrent, double[] tTractionCurrent)
        {
            if (intForceSwitch == false && stiffMatrixSwitch == false)
            {
                intForceSwitch = intFSwitch;
                stiffMatrixSwitch = stifMSwitch;
                return;
            }
            else
            {
                DisplacementVectorPrevious = currentDispVector;
                ksiVectorPrevious = ksiVectorCurrent;
                dRhoPrevious = dRhoCurrent;
                metricTensorPrevious = metricTensorCurrent;
                tTractionPrevious = tTractionCurrent;
            }
            
        }

        private void CalculateInitialValues(double[] displacementVector)
        {

        }

        public double[,] CreateGlobalStiffnessMatrix()
        {
            ksiVectorCurrent = Project(new double[2]);
            if (counter == 1)
            {
                ksiVectorPrevious = ksiVectorCurrent;
                metricTensorPrevious = new double[2, 2];
                tTractionPrevious = new double[2];
                dRhoPrevious = new List<double[]>() { new double[3], new double[3] };
            }
            counter = counter + 1;
            //double[] ksiVector = Project(new double[2]);


            if (Math.Abs(ksiVectorCurrent[0]) <= 1.05 && ksiVectorCurrent[1] <= 1.05)
            {
                Tuple<double[,], double[,], double[,]> aMatrices = CalculatePositionMatrix(ksiVectorCurrent[0], ksiVectorCurrent[1]);
                List<double[]> dRho = SurfaceVectors(ksiVectorCurrent[0], ksiVectorCurrent[1]);
                double[,] m = MetricTensor(dRho);
                double[] n = NormalVector(m, dRho);
                double[] xUpdated = xUpdatedVector();
                double ksi3 = CalculatePenetration(n, aMatrices.Item1, xUpdated);

                if (ksi3 <= 0)
                {
                    double normalTraction = PenaltyFactor * ksi3;
                    double[] tTrial = CalculateTrialTraction(ksiVectorCurrent, ksiVectorPrevious, tTractionPrevious, metricTensorPrevious, dRhoPrevious, dRho);
                    double phi = CalculatePhi(tTrial, m, normalTraction);
                    double[] tangentialTraction = CalculateTangentialTraction(phi, normalTraction, tTrial);
                    double[,] AT = MatrixOperations.Transpose(aMatrices.Item1);
                    double[] AT_n = VectorOperations.MatrixVectorProduct(AT, n);
                    List<double[,]> dA = new List<double[,]>() { aMatrices.Item2, aMatrices.Item3 };
                    double[,] K;

                    if (phi<=0)
                    {
                        double[,] Kn = Calculate_Kn(PenaltyFactor, n, aMatrices.Item1, dA, dRho, ksi3, m);
                        double[,] Kst = Calculate_Kst(PenaltyFactorT, aMatrices.Item1, dA, m, dRho, tangentialTraction, n);
                        K = MatrixOperations.MatrixAddition(Kn, Kst);
                        SetAsPreviousValues(false, true, DisplacementVector, ksiVectorCurrent, dRho, m, tangentialTraction);
                        return K;
                    }
                    else
                    {
                        double[,] Kn = Calculate_Kn(PenaltyFactor, n, aMatrices.Item1, dA, dRho, ksi3, m);
                        double[,] Ksl = Calculate_Kst(PenaltyFactorT, aMatrices.Item1, dA, m, dRho, tangentialTraction, n);
                        K = MatrixOperations.MatrixAddition(Kn, Ksl);
                        SetAsPreviousValues(false, true, DisplacementVector, ksiVectorCurrent, dRho, m, tangentialTraction);
                        return K;
                    }
                }
                else
                {
                    return new double[15, 15];
                }
            }
            else
            {
                return new double[15, 15];
            }
        }

        public double[] CreateInternalGlobalForcesVector()
        {
            ksiVectorCurrent = Project(new double[2]);
            if (counter == 1)
            {
                ksiVectorPrevious = ksiVectorCurrent;
                metricTensorPrevious = new double[2, 2];
                tTractionPrevious = new double[2];
                dRhoPrevious = new List<double[]>() { new double[3], new double[3] };
            }
            counter = counter + 1;
            if (Math.Abs(ksiVectorCurrent[0]) <= 1.05 && ksiVectorCurrent[1] <= 1.05)
            {
                Tuple<double[,], double[,], double[,]> aMatrices = CalculatePositionMatrix(ksiVectorCurrent[0], ksiVectorCurrent[1]);
                List<double[]> dRho = SurfaceVectors(ksiVectorCurrent[0], ksiVectorCurrent[1]);
                double[,] m = MetricTensor(dRho);
                double[] n = NormalVector(m, dRho);
                double[] xUpdated = xUpdatedVector();
                double ksi3 = CalculatePenetration(n, aMatrices.Item1, xUpdated);

                if (ksi3 <= 0)
                {
                    double normalTraction = PenaltyFactor * ksi3;
                    double[] tTrial = CalculateTrialTraction(ksiVectorCurrent, ksiVectorPrevious, tTractionPrevious, metricTensorPrevious, dRhoPrevious, dRho);
                    double phi = CalculatePhi(tTrial, m, normalTraction);
                    double[] tangentialTraction = CalculateTangentialTraction(phi, normalTraction, tTrial);
                    double[,] AT = MatrixOperations.Transpose(aMatrices.Item1);
                    double[] AT_n = VectorOperations.MatrixVectorProduct(AT, n);
                    double[] normalPart = VectorOperations.VectorScalarProductNew(AT_n, PenaltyFactor * ksi3);

                    double[] tangentialPart = new double[15];
                    for (int i = 0; i < 2; i++)
                    {
                        for (int j = 0; j < 2; j++)
                        {
                            tangentialPart = VectorOperations.VectorVectorAddition(tangentialPart, VectorOperations.VectorScalarProductNew(
                                VectorOperations.MatrixVectorProduct(AT, dRho[i]),
                                tangentialTraction[j] * m[i, j]));
                        }
                    }

                    double[] internalGlobalForcesVector = VectorOperations.VectorVectorAddition(normalPart, tangentialPart);
                    SetAsPreviousValues(true, false, DisplacementVector, ksiVectorCurrent, dRho, m, tangentialTraction);
                    return internalGlobalForcesVector;
                }
                else
                {
                    return new double[15];
                }
            }
            else
            {
                return new double[15];
            }
        }

        public double ClosestPointProjection()
        {
            throw new Exception("Alternative method <Project> has been used for 3D contact");
        }

        public double[,] CreateMassMatrix()
        {
            return new double[15, 15];
        }

        public double[,] CreateDampingMatrix()
        {
            return new double[15, 15];
        }
    }
}
