using System;
using System.Collections.Generic;

namespace GFEC
{
    class ANSSolidShell8LEAS7 : IElement
    {
        public Dictionary<int, INode> Nodes { get; }
        public IElementProperties Properties { get; set; }
        public Dictionary<int, bool[]> ElementFreedomSignature { get; } = new Dictionary<int, bool[]>();
        public List<int> ElementFreedomList { get; set; }
        public double[] DisplacementVector { get; set; }
        public double[] AccelerationVector { get; set; }
        public double poisson { get; set; }
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
        public ANSSolidShell8LEAS7(IElementProperties properties, Dictionary<int, INode> nodes)
        {
            Properties = properties;
            this.Nodes = nodes;
            ElementFreedomSignature[1] = new bool[] { true, true, true, false, false, false };
            ElementFreedomSignature[2] = new bool[] { true, true, true, false, false, false };
            ElementFreedomSignature[3] = new bool[] { true, true, true, false, false, false };
            ElementFreedomSignature[4] = new bool[] { true, true, true, false, false, false };
            ElementFreedomSignature[5] = new bool[] { true, true, true, false, false, false };
            ElementFreedomSignature[6] = new bool[] { true, true, true, false, false, false };
            ElementFreedomSignature[7] = new bool[] { true, true, true, false, false, false };
            ElementFreedomSignature[8] = new bool[] { true, true, true, false, false, false };
            DisplacementVector = new double[24];
            Properties.EASVector = new double[7];
            Properties.EASFEnhancedVector = new double[7];
            Properties.EASLMatrix = new double[7, 24];
            Properties.EASDMatrix = new double[7, 7];
            Properties.DisplacementVectorPreviousIncrement = new double[24];
        }
        public void CalculateElementEASMatrices()
        {
            throw new Exception("This method is to be used only for EAS method elements");
        }
        public void InitializeElementEASParameters()
        {
            Properties.EASVector[0] = 0.0;
            Properties.EASVector[1] = 0.0;
            Properties.EASVector[2] = 0.0;
            Properties.EASVector[3] = 0.0;
            Properties.EASVector[4] = 0.0;
            Properties.EASVector[5] = 0.0;
            Properties.EASVector[6] = 0.0;
        }
        public void UpdateElementEASParameters(double[] totalU)
        {
            //double[] aNext = EasInternalParametersIterations(100, 0.00001);
            //Properties.EASVector = aNext;
            //Properties.dU = deltaU;
            double[] deltaU = VectorOperations.VectorVectorSubtraction(totalU, Properties.DisplacementVectorPreviousIncrement);
            double[] aPrevious = Properties.EASVector;
            double[] Pe = Properties.EASFEnhancedVector;
            double[,] Le = Properties.EASLMatrix;
            double[,] DeInv = Properties.EASDMatrix;
            //double[] deltaAplha = VectorOperations.VectorScalarProductNew(
            //    VectorOperations.MatrixVectorProduct(DeInv, VectorOperations.VectorVectorAddition(
            //        VectorOperations.MatrixVectorProduct(Le, deltaU), Pe)), -1.0);
            //Properties.EASVector = VectorOperations.VectorVectorAddition(aPrevious, deltaAplha);
            Properties.EASVector = VectorOperations.VectorVectorAddition(aPrevious,
                VectorOperations.VectorScalarProductNew(
                VectorOperations.MatrixVectorProduct(DeInv, VectorOperations.VectorVectorAddition(
                VectorOperations.MatrixVectorProduct(Le, deltaU), Pe)), -1.0));
            Properties.DisplacementVectorPreviousIncrement = totalU;
        }
        //private double[] EasInternalParametersIterations(int maxIterations, double tolerance)
        //{
        //    double[] a = Properties.EASVector;
        //    double[,] E = CalculateStressStrainMatrix(Properties.YoungMod, Properties.PoissonRatio);
        //    double[,] DeMatrix = Properties.EASDMatrix;
        //    //double[] PeVector = new double[1];
        //    double[,] LeMatrix = Properties.EASLMatrix;
        //    Dictionary<int, double[]> xNodal = NodalCartesianCoordinatesInitial();
        //    //Dictionary<int, double[]> xNodal = NodalCartesianCoordinatesCurrent(DisplacementVector);
        //    Dictionary<int, double[]> ksiNodal = NodalNaturalCoordinates();
        //    Dictionary<int, double[]> hVectors = PositionVectors();
        //    double[] centerNaturalCoordinates = new double[3];
        //    var J0Matrices = JacobianMatrixDecomposed(xNodal, ksiNodal, hVectors, centerNaturalCoordinates);
        //    var BMatricesNaturalSystem = CalculateNaturalSystemBANSMatricesDecomposed(xNodal, ksiNodal, hVectors);
        //    var JInverse0 = CalculateInverseJacobian(J0Matrices["complete"]);
        //    double detJ0 = JInverse0.Item2;
        //    double[,] invJacobian0 = JInverse0.Item1;
        //    var transformationMatrices = StrainTransformationMatrixDecomposed(invJacobian0, J0Matrices["zero"],
        //        J0Matrices["ksi"], J0Matrices["ihta"], J0Matrices["zita"]);
        //    var BMatricesCartesianSystem = CalculateCartesianCoordinateSystemBANSMatricesDecomposed(BMatricesNaturalSystem, transformationMatrices);
        //    double[,] transformationMat0 = MatrixOperations.BlockMatrixInversion6X6(TransformationMatrixTransposed(J0Matrices["complete"]));
        //    int iteration = 0;
        //    while (iteration < maxIterations)
        //    {
        //        //DeMatrix[0, 0] = 0.0;
        //        //PeVector[0] = 0.0;
        //        //for (int k = 0; k < 24; k++)
        //        //{
        //        //    LeMatrix[0, k] = 0.0;
        //        //}
        //        //double[,] DeMatrix = new double[1, 1];
        //        double[] PeVector = new double[1];
        //        //double[,] LeMatrix = new double[1, 24];
        //        for (int i = 0; i < 2; i++)
        //        {
        //            for (int j = 0; j < 2; j++)
        //            {
        //                for (int k = 0; k < 2; k++)
        //                {
        //                    double[] gP = GaussPoints(i, j, k).Item1;
        //                    double[] gW = GaussPoints(i, j, k).Item2;
        //                    double[,] J =
        //                    MatrixOperations.MatrixAddition(
        //                    MatrixOperations.MatrixAddition(
        //                    MatrixOperations.MatrixAddition(
        //                    MatrixOperations.MatrixAddition(
        //                    MatrixOperations.MatrixAddition(
        //                    MatrixOperations.MatrixAddition(J0Matrices["zero"],
        //                    MatrixOperations.ScalarMatrixProductNew(gP[2], J0Matrices["zita"])),
        //                    MatrixOperations.ScalarMatrixProductNew(gP[0], J0Matrices["ksi"])),
        //                    MatrixOperations.ScalarMatrixProductNew(gP[1], J0Matrices["ihta"])),
        //                    MatrixOperations.ScalarMatrixProductNew(gP[0] * gP[1], J0Matrices["ksiIhta"])),
        //                    MatrixOperations.ScalarMatrixProductNew(gP[0] * gP[2], J0Matrices["zitaKsi"]))
        //                    , MatrixOperations.ScalarMatrixProductNew(gP[1] * gP[2], J0Matrices["ihtaZita"]));
        //                    double detJ = CalculateInverseJacobian(J).Item2;
        //                    double[,] B =
        //                    MatrixOperations.MatrixAddition(
        //                    MatrixOperations.MatrixAddition(
        //                    MatrixOperations.MatrixAddition(
        //                    MatrixOperations.MatrixAddition(
        //                    MatrixOperations.MatrixAddition(
        //                    MatrixOperations.MatrixAddition(BMatricesCartesianSystem["zero"],
        //                    MatrixOperations.ScalarMatrixProductNew(gP[2], BMatricesCartesianSystem["zita"])),
        //                    MatrixOperations.ScalarMatrixProductNew(gP[0], BMatricesCartesianSystem["ksi"])),
        //                    MatrixOperations.ScalarMatrixProductNew(gP[1], BMatricesCartesianSystem["ihta"])),
        //                    MatrixOperations.ScalarMatrixProductNew(gP[0] * gP[1], BMatricesCartesianSystem["ksiIhta"])),
        //                    MatrixOperations.ScalarMatrixProductNew(gP[0] * gP[2], BMatricesCartesianSystem["zitaKsi"])),
        //                    MatrixOperations.ScalarMatrixProductNew(gP[1] * gP[2], BMatricesCartesianSystem["ihtaZita"]));
        //                    //double[,] Gamma = CalculateEnhancedStrainMatrixGamma(transformationMatrices["zero"], CreateEnhancedStrainsInterpolationMatrix(gP));
        //                    double[,] Gamma = CalculateEnhancedStrainMatrixGamma(transformationMat0, CreateEnhancedStrainsInterpolationMatrix(gP),
        //                        detJ0, detJ);
        //                    double[] alphaVector = a;
        //                    double[] EnhStrainVector = VectorOperations.MatrixVectorProduct(Gamma, alphaVector);
        //                    double[] modifiedStrainVector = VectorOperations.VectorVectorAddition(CalculateStrainsVector(B), EnhStrainVector);
        //                    double[] modifiedStressVector = CalculateStressVector(E, modifiedStrainVector);
        //                    //DeMatrix = MatrixOperations.MatrixAddition(DeMatrix, MatrixOperations.ScalarMatrixProductNew(detJ * gW[0] * gW[1] * gW[2],
        //                    //    MatrixOperations.MatrixProduct(
        //                    //    MatrixOperations.Transpose(Gamma), MatrixOperations.MatrixProduct(E, Gamma))));
        //                    PeVector = VectorOperations.VectorVectorAddition(PeVector,
        //                   VectorOperations.VectorScalarProductNew(
        //                   VectorOperations.MatrixVectorProduct(MatrixOperations.Transpose(Gamma), modifiedStressVector),
        //                   detJ * gW[0] * gW[1] * gW[2]));
        //                    //LeMatrix = MatrixOperations.MatrixAddition(LeMatrix, MatrixOperations.ScalarMatrixProductNew(detJ * gW[0] * gW[1] * gW[2],
        //                    //MatrixOperations.MatrixProduct(
        //                    //MatrixOperations.Transpose(Gamma), MatrixOperations.MatrixProduct(E, B))));
        //                }
        //            }
        //        }
        //        //a = VectorOperations.VectorScalarProductNew(PeVector, -1.0 / DeMatrix[0, 0]);
        //        //double norm = Math.Abs(VectorOperations.VectorDotProduct(a, PeVector));
        //        double[] deltaAlpha = VectorOperations.VectorScalarProductNew(PeVector, -DeMatrix[0, 0]);
        //        a = VectorOperations.VectorVectorAddition(a, deltaAlpha);
        //        double norm = Math.Abs(VectorOperations.VectorDotProduct(deltaAlpha, PeVector));
        //        if (norm < tolerance) break;
        //        iteration += 1;
        //    }
        //    if (iteration >= maxIterations)
        //    {
        //        throw new Exception("EAS internal parameters did not converge");
        //    }
        //    else
        //    {
        //        //double[,] DeMatrixInv = new double[1, 1];
        //        //DeMatrixInv[0, 0] = 1.0 / DeMatrix[0, 0];
        //        //Properties.EASLMatrix = LeMatrix;
        //        //Properties.EASDMatrix = DeMatrixInv;
        //        //Properties.EASFEnhancedVector = PeVector;
        //        return a;
        //    }
        //}
        public void StoreElementFinalStepDisplacementVector(double[] solutionVector)
        {
            //Properties.DisplacementVectorPreviousStep = solutionVector;
            Properties.DisplacementVectorPreviousIncrement = solutionVector;
            //double[] Fh = Properties.HourglassInternalForceVector;
            //Properties.HourglassInternalForceVectorPreviousConvergedSolution = Fh;
        }
        public double ClosestPointProjection()
        {
            throw new Exception("Needs to be removed. Has beeb used only for testing purposes");
        }

        public Dictionary<int, INode> NodesAtFinalState()
        {
            throw new Exception("Method not implemenented");
        }
        public List<double[]> GetStressVector()
        {
            throw new Exception("Method to be added later in element ANSSolidShell8LEAS1RI");
        }
        public List<double[]> GetStrainVector()
        {
            throw new Exception("Method to be added later in element ANSSolidShell8LEAS1RI");
        }
        public List<double[]> GetGaussPointsInPhysicalSpace()
        {
            List<double[]> GpointsPhysicalCoordinates = new List<double[]>();
            double[] xUpdated = UpdateNodalCoordinates(DisplacementVector);
            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    for (int k = 0; k < 2; k++)
                    {
                        double[] gP = GaussPoints(i, j, k).Item1;
                        double[] gaussPoint = VectorOperations.MatrixVectorProduct(CalculateShapeFunctionMatrix(gP[0], gP[1], gP[2]), xUpdated);
                        GpointsPhysicalCoordinates.Add(gaussPoint);
                    }
                }
            }
            return GpointsPhysicalCoordinates;
        }
        public List<double[]> GetStressFromElementsNodes()
        {
            throw new Exception("Method to be added later in element ANSSolidShell8LEAS1RI");

        }
        public List<double[]> GetStrainFromElementsNodes()
        {
            throw new Exception("Method to be added later in element ANSSolidShell8LEAS1RI");

        }
        public List<double[]> GetStressFromElements(List<double[]> parametricCoordinatesVector)
        {
            throw new Exception("Method to be added later in element ANSSolidShell8LEAS1RI");
        }
        public List<double[]> GetphysicalCoordinatesFromElements(List<double[]> parametricCoordinatesVector)
        {
            List<double[]> PositionVectorsList = new List<double[]>();
            double[] xUpdated = UpdateNodalCoordinates(DisplacementVector);
            //PositionVectorsList.Add(new double[] { 0.0, 0.0 });
            int count = parametricCoordinatesVector.Count;
            for (int i = 0; i < count; i++)
            {
                double[] parametricCoordinatesVec = parametricCoordinatesVector[i];
                double[] positionVector = VectorOperations.MatrixVectorProduct(CalculateShapeFunctionMatrix(parametricCoordinatesVec[0], parametricCoordinatesVec[1], parametricCoordinatesVec[2]), xUpdated);
                PositionVectorsList.Add(positionVector);
            }
            return PositionVectorsList;
        }
        private double[] UpdateNodalCoordinates(double[] displacementVector)
        {
            double[] updatedCoor = new double[24];
            for (int i = 1; i <= 8; i++)
            {
                updatedCoor[3 * i - 3] = Nodes[i].XCoordinate + displacementVector[3 * i - 3];
                updatedCoor[3 * i - 2] = Nodes[i].YCoordinate + displacementVector[3 * i - 2];
                updatedCoor[3 * i - 1] = Nodes[i].ZCoordinate + displacementVector[3 * i - 1];
            }
            return updatedCoor;
        }
        private double[] InitialNodalCoordinates()
        {
            double[] initialCoordinates = new double[24];
            for (int i = 1; i <= 8; i++)
            {
                initialCoordinates[3 * i - 3] = Nodes[i].XCoordinate;
                initialCoordinates[3 * i - 2] = Nodes[i].YCoordinate;
                initialCoordinates[3 * i - 1] = Nodes[i].ZCoordinate;
            }
            return initialCoordinates;
        }
        Dictionary<int, double[]> NodalCartesianCoordinatesCurrent(double[] displacementVector)
        {
            Dictionary<int, double[]> updatedCoor = new Dictionary<int, double[]>();
            for (int i = 0; i < 8; i++)
            {
                double[] nodalPositionVector = new double[]
                {
                    Nodes[i + 1].XCoordinate + displacementVector[3 * i],
                    Nodes[i + 1].YCoordinate + displacementVector[3 * i + 1],
                    Nodes[i + 1].ZCoordinate + displacementVector[3 * i + 2]
                };
                updatedCoor.Add(i, nodalPositionVector);
            }
            return updatedCoor;
        }
        Dictionary<int, double[]> NodalCartesianCoordinatesInitial()
        {
            Dictionary<int, double[]> updatedCoor = new Dictionary<int, double[]>();
            for (int i = 0; i < 8; i++)
            {
                double[] nodalPositionVector = new double[]
                {
                    Nodes[i + 1].XCoordinate,
                    Nodes[i + 1].YCoordinate,
                    Nodes[i + 1].ZCoordinate
                };
                updatedCoor.Add(i, nodalPositionVector);
            }
            return updatedCoor;
        }
        Dictionary<int, double[]> NodalNaturalCoordinates()
        {
            Dictionary<int, double[]> naturalCoor = new Dictionary<int, double[]>();
            naturalCoor.Add(0, new double[] { 1.0, 1.0, 1.0 });
            naturalCoor.Add(1, new double[] { -1.0, 1.0, 1.0 });
            naturalCoor.Add(2, new double[] { -1.0, -1.0, 1.0 });
            naturalCoor.Add(3, new double[] { 1.0, -1.0, 1.0 });
            naturalCoor.Add(4, new double[] { 1.0, 1.0, -1.0 });
            naturalCoor.Add(5, new double[] { -1.0, 1.0, -1.0 });
            naturalCoor.Add(6, new double[] { -1.0, -1.0, -1.0 });
            naturalCoor.Add(7, new double[] { 1.0, -1.0, -1.0 });
            return naturalCoor;
        }
        Dictionary<int, double[]> PositionVectors()
        {
            Dictionary<int, double[]> positionVectors = new Dictionary<int, double[]>();
            positionVectors.Add(0, new double[] { 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 1.0, 1.0 });
            positionVectors.Add(1, new double[] { 1.0, -1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0 });
            positionVectors.Add(2, new double[] { 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0 });
            positionVectors.Add(3, new double[] { 1.0, -1.0, 1.0, -1.0, -1.0, 1.0, -1.0, 1.0 });
            return positionVectors;
        }
        private Dictionary<int, double> CalculateShapeFunctions(double ksi, double ihta, double mhi)
        {
            Dictionary<int, double> shapeFunctions = new Dictionary<int, double>();
            double N1 = 1.0 / 8.0 * (1 + ksi) * (1 + ihta) * (1 + mhi); shapeFunctions.Add(1, N1);
            double N2 = 1.0 / 8.0 * (1 - ksi) * (1 + ihta) * (1 + mhi); shapeFunctions.Add(2, N2);
            double N3 = 1.0 / 8.0 * (1 - ksi) * (1 - ihta) * (1 + mhi); shapeFunctions.Add(3, N3);
            double N4 = 1.0 / 8.0 * (1 + ksi) * (1 - ihta) * (1 + mhi); shapeFunctions.Add(4, N4);
            double N5 = 1.0 / 8.0 * (1 + ksi) * (1 + ihta) * (1 - mhi); shapeFunctions.Add(5, N5);
            double N6 = 1.0 / 8.0 * (1 - ksi) * (1 + ihta) * (1 - mhi); shapeFunctions.Add(6, N6);
            double N7 = 1.0 / 8.0 * (1 - ksi) * (1 - ihta) * (1 - mhi); shapeFunctions.Add(7, N7);
            double N8 = 1.0 / 8.0 * (1 + ksi) * (1 - ihta) * (1 - mhi); shapeFunctions.Add(8, N8);
            return shapeFunctions;
        }
        private double[,] CalculateShapeFunctionMatrix(double ksi, double ihta, double zita)
        {
            double N1 = 1.0 / 8.0 * (1 + ksi) * (1 + ihta) * (1 + zita);
            double N2 = 1.0 / 8.0 * (1 - ksi) * (1 + ihta) * (1 + zita);
            double N3 = 1.0 / 8.0 * (1 - ksi) * (1 - ihta) * (1 + zita);
            double N4 = 1.0 / 8.0 * (1 + ksi) * (1 - ihta) * (1 + zita);
            double N5 = 1.0 / 8.0 * (1 + ksi) * (1 + ihta) * (1 - zita);
            double N6 = 1.0 / 8.0 * (1 - ksi) * (1 + ihta) * (1 - zita);
            double N7 = 1.0 / 8.0 * (1 - ksi) * (1 - ihta) * (1 - zita);
            double N8 = 1.0 / 8.0 * (1 + ksi) * (1 - ihta) * (1 - zita);
            double[,] shapeFunctionsMat = new double[,] {
                {N1, 0.0, 0.0, N2, 0.0, 0.0, N3, 0.0, 0.0, N4, 0.0, 0.0, N5, 0.0, 0.0, N6, 0.0, 0.0, N7, 0.0, 0.0, N8, 0.0, 0.0 },
                {0.0, N1, 0.0, 0.0, N2, 0.0, 0.0, N3, 0.0, 0.0, N4, 0.0, 0.0, N5, 0.0, 0.0, N6, 0.0, 0.0, N7, 0.0, 0.0, N8, 0.0 },
                {0.0, 0.0, N1, 0.0, 0.0, N2, 0.0, 0.0, N3, 0.0, 0.0, N4, 0.0, 0.0, N5, 0.0, 0.0, N6, 0.0, 0.0, N7, 0.0, 0.0, N8 },
            };
            return shapeFunctionsMat;
        }

        private Dictionary<string, double[]> CalculateShapeFunctionsLocalDerivatives(double[] naturalCoordinates)
        {
            double ksi = naturalCoordinates[0];
            double ihta = naturalCoordinates[1];
            double mhi = naturalCoordinates[2];

            double[] dN_ksi = new double[]
            {
                (1.0/8.0*(1+ihta)*(1+mhi)),
                (-1.0/8.0*(1+ihta)*(1+mhi)),
                (-1.0/8.0*(1-ihta)*(1+mhi)),
                (1.0/8.0*(1-ihta)*(1+mhi)),
                (1.0/8.0*(1+ihta)*(1-mhi)),
                (-1.0/8.0*(1+ihta)*(1-mhi)),
                (-1.0/8.0*(1-ihta)*(1-mhi)),
                (1.0/8.0*(1-ihta)*(1-mhi))
            };

            double[] dN_ihta = new double[]
            {
                (1.0/8.0*(1+ksi)*(1+mhi)),
                (1.0/8.0*(1-ksi)*(1+mhi)),
                (-1.0/8.0*(1-ksi)*(1+mhi)),
                (-1.0/8.0*(1+ksi)*(1+mhi)),
                (1.0/8.0*(1+ksi)*(1-mhi)),
                (1.0/8.0*(1-ksi)*(1-mhi)),
                (-1.0/8.0*(1-ksi)*(1-mhi)),
                (-1.0/8.0*(1+ksi)*(1-mhi))
            };

            double[] dN_mhi = new double[]
            {
                (1.0/8.0*(1+ksi)*(1+ihta)),
                (1.0/8.0*(1-ksi)*(1+ihta)),
                (1.0/8.0*(1-ksi)*(1-ihta)),
                (1.0/8.0*(1+ksi)*(1-ihta)),
                (-1.0/8.0*(1+ksi)*(1+ihta)),
                (-1.0/8.0*(1-ksi)*(1+ihta)),
                (-1.0/8.0*(1-ksi)*(1-ihta)),
                (-1.0/8.0*(1+ksi)*(1-ihta))
            };

            Dictionary<string, double[]> dN = new Dictionary<string, double[]>();
            dN.Add("ksi", dN_ksi);
            dN.Add("ihta", dN_ihta);
            dN.Add("mhi", dN_mhi);
            return dN;
        }
        private double[,] CalculateFullJacobianMatrix(double[] xUpdated, Dictionary<string, double[]> dN)
        {
            double[,] jacobianMatrix = new double[3, 3];

            int k = 0;
            for (int i = 0; i < 8; i++)
            {
                jacobianMatrix[0, 0] = jacobianMatrix[0, 0] + xUpdated[k] * dN["ksi"][i];
                k = k + 3;
            }
            k = 1;
            for (int i = 0; i < 8; i++)
            {
                jacobianMatrix[0, 1] = jacobianMatrix[0, 1] + xUpdated[k] * dN["ksi"][i];
                k = k + 3;
            }
            k = 2;
            for (int i = 0; i < 8; i++)
            {
                jacobianMatrix[0, 2] = jacobianMatrix[0, 2] + xUpdated[k] * dN["ksi"][i];
                k = k + 3;
            }

            k = 0;
            for (int i = 0; i < 8; i++)
            {
                jacobianMatrix[1, 0] = jacobianMatrix[1, 0] + xUpdated[k] * dN["ihta"][i];
                k = k + 3;
            }
            k = 1;
            for (int i = 0; i < 8; i++)
            {
                jacobianMatrix[1, 1] = jacobianMatrix[1, 1] + xUpdated[k] * dN["ihta"][i];
                k = k + 3;
            }
            k = 2;
            for (int i = 0; i < 8; i++)
            {
                jacobianMatrix[1, 2] = jacobianMatrix[1, 2] + xUpdated[k] * dN["ihta"][i];
                k = k + 3;
            }

            k = 0;
            for (int i = 0; i < 8; i++)
            {
                jacobianMatrix[2, 0] = jacobianMatrix[2, 0] + xUpdated[k] * dN["mhi"][i];
                k = k + 3;
            }
            k = 1;
            for (int i = 0; i < 8; i++)
            {
                jacobianMatrix[2, 1] = jacobianMatrix[2, 1] + xUpdated[k] * dN["mhi"][i];
                k = k + 3;
            }
            k = 2;
            for (int i = 0; i < 8; i++)
            {
                jacobianMatrix[2, 2] = jacobianMatrix[2, 2] + xUpdated[k] * dN["mhi"][i];
                k = k + 3;
            }

            return jacobianMatrix;
        }
        private Dictionary<string, double[,]> JacobianMatrixDecomposed(Dictionary<int, double[]> xNodal,
                Dictionary<int, double[]> ksiNodal,
                Dictionary<int, double[]> hVectors,
                double[] naturalCoordinates)
        {
            double[,] jacobianMatrixZero = new double[3, 3];
            double[,] jacobianMatrixKsi = new double[3, 3];
            double[,] jacobianMatrixIhta = new double[3, 3];
            double[,] jacobianMatrixZita = new double[3, 3];
            double[,] jacobianMatrixKsiIhta = new double[3, 3];
            double[,] jacobianMatrixIhtaZita = new double[3, 3];
            double[,] jacobianMatrixZitaKsi = new double[3, 3];

            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    double sum0 = 0.0;
                    if (i == 0)
                    {
                        double sumIhta = 0.0;
                        double sumZita = 0.0;
                        double sumIhtaZita = 0.0;
                        for (int k = 0; k < 8; k++)
                        {
                            sum0 += ksiNodal[k][i] * xNodal[k][j];
                            sumIhta += hVectors[2][k] * xNodal[k][j];
                            sumZita += hVectors[1][k] * xNodal[k][j];
                            sumIhtaZita += hVectors[3][k] * xNodal[k][j];
                        }
                        jacobianMatrixIhta[i, j] = sumIhta / 8.0;
                        jacobianMatrixZita[i, j] = sumZita / 8.0;
                        jacobianMatrixIhtaZita[i, j] = sumIhtaZita / 8.0;
                    }
                    else if (i == 1)
                    {
                        double sumKsi = 0.0;
                        double sumZita = 0.0;
                        double sumZitaKsi = 0.0;
                        for (int k = 0; k < 8; k++)
                        {
                            sum0 += ksiNodal[k][i] * xNodal[k][j];
                            sumKsi += hVectors[2][k] * xNodal[k][j];
                            sumZita += hVectors[0][k] * xNodal[k][j];
                            sumZitaKsi += hVectors[3][k] * xNodal[k][j];
                        }
                        jacobianMatrixKsi[i, j] = sumKsi / 8.0;
                        jacobianMatrixZita[i, j] = sumZita / 8.0;
                        jacobianMatrixZitaKsi[i, j] = sumZitaKsi / 8.0;
                    }
                    else
                    {
                        double sumKsi = 0.0;
                        double sumIhta = 0.0;
                        double sumKsiIta = 0.0;
                        for (int k = 0; k < 8; k++)
                        {
                            sum0 += ksiNodal[k][i] * xNodal[k][j];
                            sumKsi += hVectors[1][k] * xNodal[k][j];
                            sumIhta += hVectors[0][k] * xNodal[k][j];
                            sumKsiIta += hVectors[3][k] * xNodal[k][j];
                        }
                        jacobianMatrixKsi[i, j] = sumKsi / 8.0;
                        jacobianMatrixIhta[i, j] = sumIhta / 8.0;
                        jacobianMatrixKsiIhta[i, j] = sumKsiIta / 8.0;
                    }
                    jacobianMatrixZero[i, j] = sum0 / 8.0;
                }
            }
            double[,] jacobianMatrix = MatrixOperations.MatrixAddition(
                MatrixOperations.MatrixAddition(
                MatrixOperations.MatrixAddition(
                MatrixOperations.MatrixAddition(
                MatrixOperations.MatrixAddition(
                MatrixOperations.MatrixAddition(jacobianMatrixZero,
                MatrixOperations.ScalarMatrixProductNew(naturalCoordinates[0], jacobianMatrixKsi)),
                MatrixOperations.ScalarMatrixProductNew(naturalCoordinates[1], jacobianMatrixIhta)),
                MatrixOperations.ScalarMatrixProductNew(naturalCoordinates[2], jacobianMatrixZita)),
                MatrixOperations.ScalarMatrixProductNew(naturalCoordinates[0] * naturalCoordinates[1],
                jacobianMatrixKsiIhta)),
                MatrixOperations.ScalarMatrixProductNew(naturalCoordinates[1] * naturalCoordinates[2],
                jacobianMatrixIhtaZita)),
                MatrixOperations.ScalarMatrixProductNew(naturalCoordinates[0] * naturalCoordinates[2],
                jacobianMatrixZitaKsi));

            Dictionary<string, double[,]> jacobianMatricesDecomposed = new Dictionary<string, double[,]>();

            jacobianMatricesDecomposed.Add("complete", jacobianMatrix);
            jacobianMatricesDecomposed.Add("zero", jacobianMatrixZero);
            jacobianMatricesDecomposed.Add("ksi", jacobianMatrixKsi);
            jacobianMatricesDecomposed.Add("ihta", jacobianMatrixIhta);
            jacobianMatricesDecomposed.Add("zita", jacobianMatrixZita);
            jacobianMatricesDecomposed.Add("ksiIhta", jacobianMatrixKsiIhta);
            jacobianMatricesDecomposed.Add("ihtaZita", jacobianMatrixIhtaZita);
            jacobianMatricesDecomposed.Add("zitaKsi", jacobianMatrixZitaKsi);

            return jacobianMatricesDecomposed;
        }
        private double CalculateJacobianDet(double[,] jacobianMatrix)
        {
            double[,] jacobianInverseMatrix = new double[3, 3];

            jacobianInverseMatrix[0, 0] = jacobianMatrix[1, 1] * jacobianMatrix[2, 2] - jacobianMatrix[1, 2] * jacobianMatrix[2, 1];
            jacobianInverseMatrix[1, 1] = jacobianMatrix[2, 2] * jacobianMatrix[0, 0] - jacobianMatrix[2, 0] * jacobianMatrix[2, 2];
            jacobianInverseMatrix[2, 2] = jacobianMatrix[0, 0] * jacobianMatrix[1, 1] - jacobianMatrix[0, 1] * jacobianMatrix[1, 0];

            jacobianInverseMatrix[0, 1] = jacobianMatrix[1, 2] * jacobianMatrix[2, 0] - jacobianMatrix[1, 0] * jacobianMatrix[2, 2];
            jacobianInverseMatrix[1, 2] = jacobianMatrix[2, 0] * jacobianMatrix[0, 1] - jacobianMatrix[2, 1] * jacobianMatrix[0, 0];
            jacobianInverseMatrix[2, 0] = jacobianMatrix[0, 1] * jacobianMatrix[1, 2] - jacobianMatrix[0, 2] * jacobianMatrix[1, 1];

            jacobianInverseMatrix[1, 0] = jacobianMatrix[2, 1] * jacobianMatrix[0, 2] - jacobianMatrix[0, 1] * jacobianMatrix[2, 2];
            jacobianInverseMatrix[2, 1] = jacobianMatrix[0, 2] * jacobianMatrix[1, 0] - jacobianMatrix[1, 2] * jacobianMatrix[0, 0];
            jacobianInverseMatrix[0, 2] = jacobianMatrix[1, 0] * jacobianMatrix[1, 1] - jacobianMatrix[2, 0] * jacobianMatrix[1, 1];

            double detj = jacobianMatrix[0, 0] * jacobianInverseMatrix[0, 0] + jacobianMatrix[0, 1] * jacobianInverseMatrix[1, 0] + jacobianMatrix[0, 2] * jacobianInverseMatrix[2, 0];
            return detj;
        }
        private double[,] CalculateTransformationMatrix(double[,] jacobianInverseMatrix)
        {
            double[,] T = new double[6, 6];
            T[0, 0] = jacobianInverseMatrix[0, 0] * jacobianInverseMatrix[0, 0];
            T[0, 1] = jacobianInverseMatrix[0, 1] * jacobianInverseMatrix[0, 1];
            T[0, 2] = jacobianInverseMatrix[0, 2] * jacobianInverseMatrix[0, 2];
            T[0, 3] = jacobianInverseMatrix[0, 0] * jacobianInverseMatrix[0, 1];
            T[0, 4] = jacobianInverseMatrix[0, 1] * jacobianInverseMatrix[0, 2];
            T[0, 5] = jacobianInverseMatrix[0, 2] * jacobianInverseMatrix[0, 0];

            T[1, 0] = jacobianInverseMatrix[1, 0] * jacobianInverseMatrix[1, 0];
            T[1, 1] = jacobianInverseMatrix[1, 1] * jacobianInverseMatrix[1, 1];
            T[1, 2] = jacobianInverseMatrix[1, 2] * jacobianInverseMatrix[1, 2];
            T[1, 3] = jacobianInverseMatrix[1, 0] * jacobianInverseMatrix[1, 1];
            T[1, 4] = jacobianInverseMatrix[1, 1] * jacobianInverseMatrix[1, 2];
            T[1, 5] = jacobianInverseMatrix[1, 2] * jacobianInverseMatrix[1, 0];

            T[2, 0] = jacobianInverseMatrix[2, 0] * jacobianInverseMatrix[2, 0];
            T[2, 1] = jacobianInverseMatrix[2, 1] * jacobianInverseMatrix[2, 1];
            T[2, 2] = jacobianInverseMatrix[2, 2] * jacobianInverseMatrix[2, 2];
            T[2, 3] = jacobianInverseMatrix[2, 0] * jacobianInverseMatrix[2, 1];
            T[2, 4] = jacobianInverseMatrix[2, 1] * jacobianInverseMatrix[2, 2];
            T[2, 5] = jacobianInverseMatrix[2, 2] * jacobianInverseMatrix[2, 0];

            T[3, 0] = 2 * jacobianInverseMatrix[0, 0] * jacobianInverseMatrix[1, 0];
            T[3, 1] = 2 * jacobianInverseMatrix[0, 1] * jacobianInverseMatrix[1, 1];
            T[3, 2] = 2 * jacobianInverseMatrix[0, 2] * jacobianInverseMatrix[1, 2];

            T[3, 3] = jacobianInverseMatrix[0, 0] * jacobianInverseMatrix[1, 1] +
                                jacobianInverseMatrix[0, 1] * jacobianInverseMatrix[1, 0];

            T[3, 4] = jacobianInverseMatrix[0, 1] * jacobianInverseMatrix[1, 2] +
                                jacobianInverseMatrix[0, 2] * jacobianInverseMatrix[1, 1];

            T[3, 5] = jacobianInverseMatrix[0, 2] * jacobianInverseMatrix[1, 0] +
                                jacobianInverseMatrix[0, 0] * jacobianInverseMatrix[2, 1];

            T[4, 0] = 2 * jacobianInverseMatrix[1, 0] * jacobianInverseMatrix[2, 0];
            T[4, 1] = 2 * jacobianInverseMatrix[1, 1] * jacobianInverseMatrix[2, 1];
            T[4, 2] = 2 * jacobianInverseMatrix[1, 2] * jacobianInverseMatrix[2, 2];

            T[4, 3] = jacobianInverseMatrix[1, 0] * jacobianInverseMatrix[2, 1] +
                                jacobianInverseMatrix[1, 1] * jacobianInverseMatrix[2, 0];

            T[4, 4] = jacobianInverseMatrix[1, 1] * jacobianInverseMatrix[2, 2] +
                                jacobianInverseMatrix[1, 2] * jacobianInverseMatrix[2, 1];

            T[4, 5] = jacobianInverseMatrix[1, 2] * jacobianInverseMatrix[2, 0] +
                                jacobianInverseMatrix[1, 0] * jacobianInverseMatrix[2, 2];

            T[5, 0] = 2 * jacobianInverseMatrix[2, 0] * jacobianInverseMatrix[0, 0];
            T[5, 1] = 2 * jacobianInverseMatrix[2, 1] * jacobianInverseMatrix[0, 1];
            T[5, 2] = 2 * jacobianInverseMatrix[2, 2] * jacobianInverseMatrix[0, 2];

            T[5, 3] = jacobianInverseMatrix[2, 0] * jacobianInverseMatrix[0, 1] +
                                jacobianInverseMatrix[2, 1] * jacobianInverseMatrix[0, 0];

            T[5, 4] = jacobianInverseMatrix[2, 1] * jacobianInverseMatrix[0, 2] +
                                jacobianInverseMatrix[2, 2] * jacobianInverseMatrix[0, 1];

            T[5, 5] = jacobianInverseMatrix[2, 2] * jacobianInverseMatrix[0, 0] +
                                jacobianInverseMatrix[2, 0] * jacobianInverseMatrix[0, 2];
            return T;
        }
        private Dictionary<string, double[,]> StrainTransformationMatrixDecomposed(double[,] JacobianInverseMatrix,
                double[,] JacobianMatrix0,
                double[,] JacobianMatrixKsi,
                double[,] JacobianMatrixIhta,
                double[,] JacobianMatrixZita)
        {
            double[,] JacobianMatrix0Inverse = CalculateInverseJacobian(JacobianMatrix0).Item1;
            double[,] JacobianMatrixKsiInverse = MatrixOperations.ScalarMatrixProductNew(-1.0, MatrixOperations.MatrixProduct(
                MatrixOperations.MatrixProduct(JacobianMatrix0Inverse, JacobianMatrixKsi), JacobianMatrix0Inverse));
            double[,] JacobianMatrixIhtaInverse = MatrixOperations.ScalarMatrixProductNew(-1.0, MatrixOperations.MatrixProduct(
                MatrixOperations.MatrixProduct(JacobianMatrix0Inverse, JacobianMatrixIhta), JacobianMatrix0Inverse));
            double[,] JacobianMatrixZitaInverse = MatrixOperations.ScalarMatrixProductNew(-1.0, MatrixOperations.MatrixProduct(
                MatrixOperations.MatrixProduct(JacobianMatrix0Inverse, JacobianMatrixZita), JacobianMatrix0Inverse));

            double[,] Tzero = CalculateTransformationMatrix(JacobianInverseMatrix);
            double[,] TKsi = CalculateTransformationMatrix(JacobianMatrixKsiInverse);
            double[,] TIhta = CalculateTransformationMatrix(JacobianMatrixIhtaInverse);
            double[,] TZita = CalculateTransformationMatrix(JacobianMatrixZitaInverse);

            Dictionary<string, double[,]> StrainTransformationMatrices = new Dictionary<string, double[,]>();
            StrainTransformationMatrices.Add("zero", Tzero);
            StrainTransformationMatrices.Add("ksi", TKsi);
            StrainTransformationMatrices.Add("ihta", TIhta);
            StrainTransformationMatrices.Add("zita", TZita);
            return StrainTransformationMatrices;
        }
        private Tuple<double[,], double> CalculateInverseJacobian(double[,] jacobianMatrix)
        {
            double[,] jacobianInverseMatrix = new double[3, 3];

            jacobianInverseMatrix[0, 0] = jacobianMatrix[1, 1] * jacobianMatrix[2, 2] - jacobianMatrix[1, 2] * jacobianMatrix[2, 1];
            jacobianInverseMatrix[1, 1] = jacobianMatrix[2, 2] * jacobianMatrix[0, 0] - jacobianMatrix[2, 0] * jacobianMatrix[2, 2];
            jacobianInverseMatrix[2, 2] = jacobianMatrix[0, 0] * jacobianMatrix[1, 1] - jacobianMatrix[0, 1] * jacobianMatrix[1, 0];

            jacobianInverseMatrix[0, 1] = jacobianMatrix[1, 2] * jacobianMatrix[2, 0] - jacobianMatrix[1, 0] * jacobianMatrix[2, 2];
            jacobianInverseMatrix[1, 2] = jacobianMatrix[2, 0] * jacobianMatrix[0, 1] - jacobianMatrix[2, 1] * jacobianMatrix[0, 0];
            jacobianInverseMatrix[2, 0] = jacobianMatrix[0, 1] * jacobianMatrix[1, 2] - jacobianMatrix[0, 2] * jacobianMatrix[1, 1];

            jacobianInverseMatrix[1, 0] = jacobianMatrix[2, 1] * jacobianMatrix[0, 2] - jacobianMatrix[0, 1] * jacobianMatrix[2, 2];
            jacobianInverseMatrix[2, 1] = jacobianMatrix[0, 2] * jacobianMatrix[1, 0] - jacobianMatrix[1, 2] * jacobianMatrix[0, 0];
            jacobianInverseMatrix[0, 2] = jacobianMatrix[1, 0] * jacobianMatrix[1, 1] - jacobianMatrix[2, 0] * jacobianMatrix[1, 1];

            double detj = jacobianMatrix[0, 0] * jacobianInverseMatrix[0, 0] + jacobianMatrix[0, 1] * jacobianInverseMatrix[1, 0] + jacobianMatrix[0, 2] * jacobianInverseMatrix[2, 0];

            jacobianInverseMatrix[0, 0] = jacobianInverseMatrix[0, 0] / detj;
            jacobianInverseMatrix[1, 1] = jacobianInverseMatrix[1, 1] / detj;
            jacobianInverseMatrix[2, 2] = jacobianInverseMatrix[2, 2] / detj;

            jacobianInverseMatrix[0, 1] = jacobianInverseMatrix[0, 1] / detj;
            jacobianInverseMatrix[1, 2] = jacobianInverseMatrix[1, 2] / detj;
            jacobianInverseMatrix[2, 0] = jacobianInverseMatrix[2, 0] / detj;

            jacobianInverseMatrix[1, 0] = jacobianInverseMatrix[1, 0] / detj;
            jacobianInverseMatrix[2, 1] = jacobianInverseMatrix[2, 1] / detj;
            jacobianInverseMatrix[0, 2] = jacobianInverseMatrix[0, 2] / detj;

            return new Tuple<double[,], double>(jacobianInverseMatrix, detj);
        }
        private double[] CalculateStrainsVector(double[,] Bmatrix)
        {
            double[] strains = VectorOperations.MatrixVectorProduct(Bmatrix, DisplacementVector);
            return strains;
        }
        private double[] CalculateIncrementStrainsVector(double[,] Bmatrix)
        {
            double[] deltaU = VectorOperations.VectorVectorSubtraction(DisplacementVector, Properties.DisplacementVectorPreviousStep);
            double[] strainsIncrement = VectorOperations.MatrixVectorProduct(Bmatrix, deltaU);
            return strainsIncrement;
        }
        private Dictionary<string, double[,]> CalculateNaturalSystemBANSMatricesDecomposed(Dictionary<int, double[]> xNodal,
                Dictionary<int, double[]> ksiNodal,
                Dictionary<int, double[]> hVectors)
        {
            double[,] BZero = new double[6, 24];
            double[,] BKsi = new double[6, 24];
            double[,] BIhta = new double[6, 24];
            double[,] BZita = new double[6, 24];
            double[,] BKsiIhta = new double[6, 24];
            double[,] BIhtaZita = new double[6, 24];
            double[,] BZitaKsi = new double[6, 24];

            for (int node = 0; node < 8; node++)
            {
                for (int i = 0; i < 6; i++)
                {
                    for (int j = 0; j < 3; j++)
                    {
                        double sum0 = 0.0;
                        if (i == 0)
                        {
                            double sumIhta = 0.0;
                            double sumZita = 0.0;
                            double sumIhtaZita = 0.0;
                            for (int k = 0; k < 8; k++)
                            {
                                sum0 += ksiNodal[k][i] * xNodal[k][j] * ksiNodal[node][i];
                                //-------------------------------------------------------------
                                sumIhta += ksiNodal[k][0] * xNodal[k][j] * hVectors[2][node];
                                sumIhta += hVectors[2][k] * xNodal[k][j] * ksiNodal[node][0];
                                //-------------------------------------------------------------
                                sumZita += ksiNodal[k][0] * xNodal[k][j] * hVectors[1][node];
                                sumZita += hVectors[1][k] * xNodal[k][j] * ksiNodal[node][0];
                                //-------------------------------------------------------------
                                sumIhtaZita += ksiNodal[k][0] * xNodal[k][j] * hVectors[3][node];
                                sumIhtaZita += hVectors[3][k] * xNodal[k][j] * ksiNodal[node][0];
                                sumIhtaZita += hVectors[1][k] * xNodal[k][j] * hVectors[2][node];
                                sumIhtaZita += hVectors[2][k] * xNodal[k][j] * hVectors[1][node];
                            }
                            BIhta[i, node * 3 + j] = sumIhta / 64.0;
                            BZita[i, node * 3 + j] = sumZita / 64.0;
                            BIhtaZita[i, node * 3 + j] = sumIhtaZita / 64.0;
                        }
                        else if (i == 1)
                        {
                            double sumKsi = 0.0;
                            double sumZita = 0.0;
                            double sumZitaKsi = 0.0;
                            for (int k = 0; k < 8; k++)
                            {
                                sum0 += ksiNodal[k][i] * xNodal[k][j] * ksiNodal[node][i];
                                //-------------------------------------------------------------
                                sumKsi += ksiNodal[k][1] * xNodal[k][j] * hVectors[2][node];
                                sumKsi += hVectors[2][k] * xNodal[k][j] * ksiNodal[node][1];
                                //-------------------------------------------------------------
                                sumZita += ksiNodal[k][1] * xNodal[k][j] * hVectors[0][node];
                                sumZita += hVectors[0][k] * xNodal[k][j] * ksiNodal[node][1];
                                //-------------------------------------------------------------
                                sumZitaKsi += ksiNodal[k][1] * xNodal[k][j] * hVectors[3][node];
                                sumZitaKsi += hVectors[3][k] * xNodal[k][j] * ksiNodal[node][1];
                                sumZitaKsi += hVectors[2][k] * xNodal[k][j] * hVectors[0][node];
                                sumZitaKsi += hVectors[0][k] * xNodal[k][j] * hVectors[2][node];
                            }
                            BKsi[i, node * 3 + j] = sumKsi / 64.0;
                            BZita[i, node * 3 + j] = sumZita / 64.0;
                            BZitaKsi[i, node * 3 + j] = sumZitaKsi / 64.0;
                        }
                        else if (i == 2)
                        {
                            double sumKsi = 0.0;
                            double sumIhta = 0.0;
                            double sumKsiIhta = 0.0;
                            for (int k = 0; k < 8; k++)
                            {
                                sum0 += ksiNodal[k][i] * xNodal[k][j] * ksiNodal[node][i];
                                sum0 += hVectors[0][k] * xNodal[k][j] * hVectors[0][node];
                                sum0 += hVectors[1][k] * xNodal[k][j] * hVectors[1][node];
                                sum0 += hVectors[3][k] * xNodal[k][j] * hVectors[3][node];
                                //-------------------------------------------------------------
                                sumKsi += ksiNodal[k][2] * xNodal[k][j] * hVectors[1][node];
                                sumKsi += hVectors[1][k] * xNodal[k][j] * ksiNodal[node][2];
                                sumKsi += hVectors[0][k] * xNodal[k][j] * hVectors[3][node];
                                sumKsi += hVectors[3][k] * xNodal[k][j] * hVectors[0][node];
                                //-------------------------------------------------------------
                                sumIhta += ksiNodal[k][2] * xNodal[k][j] * hVectors[0][node];
                                sumIhta += hVectors[0][k] * xNodal[k][j] * ksiNodal[node][2];
                                sumIhta += hVectors[1][k] * xNodal[k][j] * hVectors[3][node];
                                sumIhta += hVectors[3][k] * xNodal[k][j] * hVectors[1][node];
                                //-------------------------------------------------------------
                                sumKsiIhta += ksiNodal[k][2] * xNodal[k][j] * hVectors[3][node];
                                sumKsiIhta += hVectors[3][k] * xNodal[k][j] * ksiNodal[node][2];
                                sumKsiIhta += hVectors[0][k] * xNodal[k][j] * hVectors[1][node];
                                sumKsiIhta += hVectors[1][k] * xNodal[k][j] * hVectors[0][node];
                            }
                            BKsi[i, node * 3 + j] = sumKsi / 64.0;
                            BIhta[i, node * 3 + j] = sumIhta / 64.0;
                            BKsiIhta[i, node * 3 + j] = sumKsiIhta / 64.0;
                        }
                        else if (i == 3)
                        {
                            double sumKsi = 0.0;
                            double sumIhta = 0.0;
                            double sumZita = 0.0;
                            double sumKsiIhta = 0.0;
                            double sumIhtaZita = 0.0;
                            double sumZitaKsi = 0.0;
                            for (int k = 0; k < 8; k++)
                            {
                                sum0 += ksiNodal[k][0] * xNodal[k][j] * ksiNodal[node][1];
                                sum0 += ksiNodal[k][1] * xNodal[k][j] * ksiNodal[node][0];
                                //-------------------------------------------------------------
                                sumKsi += ksiNodal[k][0] * xNodal[k][j] * hVectors[2][node];
                                sumKsi += hVectors[2][k] * xNodal[k][j] * ksiNodal[node][0];
                                //-------------------------------------------------------------
                                sumIhta += ksiNodal[k][1] * xNodal[k][j] * hVectors[2][node];
                                sumIhta += hVectors[2][k] * xNodal[k][j] * ksiNodal[node][1];
                                //-------------------------------------------------------------
                                sumZita += ksiNodal[k][0] * xNodal[k][j] * hVectors[0][node];
                                sumZita += hVectors[0][k] * xNodal[k][j] * ksiNodal[node][0];
                                sumZita += ksiNodal[k][1] * xNodal[k][j] * hVectors[1][node];
                                sumZita += hVectors[1][k] * xNodal[k][j] * ksiNodal[node][1];
                                //-------------------------------------------------------------
                                sumKsiIhta += 2.0 * hVectors[2][k] * xNodal[k][j] * hVectors[2][node];
                                //-------------------------------------------------------------
                                sumIhtaZita += ksiNodal[k][1] * xNodal[k][j] * hVectors[3][node];
                                sumIhtaZita += hVectors[3][k] * xNodal[k][j] * ksiNodal[node][1];
                                sumIhtaZita += hVectors[2][k] * xNodal[k][j] * hVectors[0][node];
                                sumIhtaZita += hVectors[0][k] * xNodal[k][j] * hVectors[2][node];
                                //-------------------------------------------------------------
                                sumZitaKsi += ksiNodal[k][0] * xNodal[k][j] * hVectors[3][node];
                                sumZitaKsi += hVectors[3][k] * xNodal[k][j] * ksiNodal[node][0];
                                sumZitaKsi += hVectors[1][k] * xNodal[k][j] * hVectors[2][node];
                                sumZitaKsi += hVectors[2][k] * xNodal[k][j] * hVectors[1][node];
                            }
                            BKsi[i, node * 3 + j] = sumKsi / 64.0;
                            BIhta[i, node * 3 + j] = sumIhta / 64.0;
                            BZita[i, node * 3 + j] = sumZita / 64.0;
                            BKsiIhta[i, node * 3 + j] = sumKsiIhta / 64.0;
                            BIhtaZita[i, node * 3 + j] = sumIhtaZita / 64.0;
                            BZitaKsi[i, node * 3 + j] = sumZitaKsi / 64.0;
                        }
                        else if (i == 4)
                        {
                            double sumKsi = 0.0;
                            double sumZita = 0.0;
                            double sumZitaKsi = 0.0;
                            for (int k = 0; k < 8; k++)
                            {
                                sum0 += ksiNodal[k][1] * xNodal[k][j] * ksiNodal[node][2];
                                sum0 += ksiNodal[k][2] * xNodal[k][j] * ksiNodal[node][1];
                                sum0 += hVectors[1][k] * xNodal[k][j] * hVectors[2][node];
                                sum0 += hVectors[2][k] * xNodal[k][j] * hVectors[1][node];
                                //-------------------------------------------------------------
                                sumKsi += ksiNodal[k][1] * xNodal[k][j] * hVectors[1][node];
                                sumKsi += hVectors[1][k] * xNodal[k][j] * ksiNodal[node][1];
                                sumKsi += ksiNodal[k][2] * xNodal[k][j] * hVectors[2][node];
                                sumKsi += hVectors[2][k] * xNodal[k][j] * ksiNodal[node][2];
                                //-------------------------------------------------------------
                                sumZita += ksiNodal[k][2] * xNodal[k][j] * hVectors[0][node];
                                sumZita += hVectors[0][k] * xNodal[k][j] * ksiNodal[node][2];
                                sumZita += hVectors[1][k] * xNodal[k][j] * hVectors[3][node];
                                sumZita += hVectors[3][k] * xNodal[k][j] * hVectors[1][node];
                                //-------------------------------------------------------------
                                sumZitaKsi += ksiNodal[k][2] * xNodal[k][j] * hVectors[3][node];
                                sumZitaKsi += hVectors[3][k] * xNodal[k][j] * ksiNodal[node][2];
                                sumZitaKsi += hVectors[0][k] * xNodal[k][j] * hVectors[1][node];
                                sumZitaKsi += hVectors[1][k] * xNodal[k][j] * hVectors[0][node];
                            }
                            BKsi[i, node * 3 + j] = sumKsi / 64.0;
                            BZita[i, node * 3 + j] = sumZita / 64.0;
                            BZitaKsi[i, node * 3 + j] = sumZitaKsi / 64.0;
                        }
                        else if (i == 5)
                        {
                            double sumIhta = 0.0;
                            double sumZita = 0.0;
                            double sumIhtaZita = 0.0;
                            for (int k = 0; k < 8; k++)
                            {
                                sum0 += ksiNodal[k][2] * xNodal[k][j] * ksiNodal[node][0];
                                sum0 += ksiNodal[k][0] * xNodal[k][j] * ksiNodal[node][2];
                                sum0 += hVectors[2][k] * xNodal[k][j] * hVectors[0][node];
                                sum0 += hVectors[0][k] * xNodal[k][j] * hVectors[2][node];
                                //-------------------------------------------------------------
                                sumIhta += ksiNodal[k][2] * xNodal[k][j] * hVectors[2][node];
                                sumIhta += hVectors[2][k] * xNodal[k][j] * ksiNodal[node][2];
                                sumIhta += ksiNodal[k][0] * xNodal[k][j] * hVectors[0][node];
                                sumIhta += hVectors[0][k] * xNodal[k][j] * ksiNodal[node][0];
                                //-------------------------------------------------------------
                                sumZita += ksiNodal[k][2] * xNodal[k][j] * hVectors[1][node];
                                sumZita += hVectors[1][k] * xNodal[k][j] * ksiNodal[node][2];
                                sumZita += hVectors[0][k] * xNodal[k][j] * hVectors[3][node];
                                sumZita += hVectors[3][k] * xNodal[k][j] * hVectors[0][node];
                                //-------------------------------------------------------------
                                sumIhtaZita += ksiNodal[k][2] * xNodal[k][j] * hVectors[3][node];
                                sumIhtaZita += hVectors[3][k] * xNodal[k][j] * ksiNodal[node][2];
                                sumIhtaZita += hVectors[0][k] * xNodal[k][j] * hVectors[1][node];
                                sumIhtaZita += hVectors[1][k] * xNodal[k][j] * hVectors[0][node];
                            }
                            BIhta[i, node * 3 + j] = sumIhta / 64.0;
                            BZita[i, node * 3 + j] = sumZita / 64.0;
                            BIhtaZita[i, node * 3 + j] = sumIhtaZita / 64.0;
                        }
                        BZero[i, node * 3 + j] = sum0 / 64.0;
                    }
                }
            }
            Dictionary<string, double[,]> BMatrices = new Dictionary<string, double[,]>();
            BMatrices.Add("zero", BZero);
            BMatrices.Add("ksi", BKsi);
            BMatrices.Add("ihta", BIhta);
            BMatrices.Add("zita", BZita);
            BMatrices.Add("ksiIhta", BKsiIhta);
            BMatrices.Add("ihtaZita", BIhtaZita);
            BMatrices.Add("zitaKsi", BZitaKsi);
            return BMatrices;
        }
        private Dictionary<string, double[,]> CalculateCartesianCoordinateSystemBANSMatricesDecomposed(Dictionary<string, double[,]> BANSMatrices,
                Dictionary<string, double[,]> TransformationMatrices)
        {
            double[,] BZero = MatrixOperations.MatrixProduct(TransformationMatrices["zero"], BANSMatrices["zero"]);
            double[,] BKsi = MatrixOperations.MatrixAddition(
                MatrixOperations.MatrixProduct(TransformationMatrices["zero"], BANSMatrices["ksi"]),
                MatrixOperations.MatrixProduct(TransformationMatrices["ksi"], BANSMatrices["zero"])
                );
            double[,] BIhta = MatrixOperations.MatrixAddition(
                MatrixOperations.MatrixProduct(TransformationMatrices["zero"], BANSMatrices["ihta"]),
                MatrixOperations.MatrixProduct(TransformationMatrices["ihta"], BANSMatrices["zero"])
                );
            double[,] BZita = MatrixOperations.MatrixAddition(
                MatrixOperations.MatrixProduct(TransformationMatrices["zero"], BANSMatrices["zita"]),
                MatrixOperations.MatrixProduct(TransformationMatrices["zita"], BANSMatrices["zero"])
                );
            double[,] BKsiIhta = MatrixOperations.MatrixAddition(
                MatrixOperations.MatrixAddition(
                MatrixOperations.MatrixProduct(TransformationMatrices["zero"], BANSMatrices["ksiIhta"]),
                MatrixOperations.MatrixProduct(TransformationMatrices["ksi"], BANSMatrices["ihta"])
                ),
                MatrixOperations.MatrixProduct(TransformationMatrices["ihta"], BANSMatrices["ksi"])
                );
            double[,] BIhtaZita = MatrixOperations.MatrixAddition(
                MatrixOperations.MatrixAddition(
                MatrixOperations.MatrixProduct(TransformationMatrices["zero"], BANSMatrices["ihtaZita"]),
                MatrixOperations.MatrixProduct(TransformationMatrices["ihta"], BANSMatrices["zita"])
                ),
                MatrixOperations.MatrixProduct(TransformationMatrices["zita"], BANSMatrices["ihta"])
                );
            double[,] BZitaKsi = MatrixOperations.MatrixAddition(
                MatrixOperations.MatrixAddition(
                MatrixOperations.MatrixProduct(TransformationMatrices["zero"], BANSMatrices["zitaKsi"]),
                MatrixOperations.MatrixProduct(TransformationMatrices["zita"], BANSMatrices["ksi"])
                ),
                MatrixOperations.MatrixProduct(TransformationMatrices["ksi"], BANSMatrices["zita"])
                );
            Dictionary<string, double[,]> matricesB = new Dictionary<string, double[,]>();
            matricesB.Add("zero", BZero);
            matricesB.Add("ksi", BKsi);
            matricesB.Add("ihta", BIhta);
            matricesB.Add("zita", BZita);
            matricesB.Add("ksiIhta", BKsiIhta);
            matricesB.Add("ihtaZita", BIhtaZita);
            matricesB.Add("zitaKsi", BZitaKsi);
            return matricesB;
        }
        private Dictionary<string, double[,]> HourglassComponentModification(Dictionary<string, double[,]> Bmatrices)
        {
            Dictionary<string, double[,]> matricesB = new Dictionary<string, double[,]>();
            double[,] modificationMatrix = new double[6, 6];
            modificationMatrix[0, 0] = 1.0;
            modificationMatrix[0, 1] = 1.0;
            modificationMatrix[0, 2] = 1.0;
            modificationMatrix[1, 0] = 1.0;
            modificationMatrix[1, 1] = 1.0;
            modificationMatrix[1, 2] = 1.0;
            modificationMatrix[2, 0] = 1.0;
            modificationMatrix[2, 1] = 1.0;
            modificationMatrix[2, 2] = 1.0;
            double[,] BKsi = MatrixOperations.MatrixSubtraction(Bmatrices["ksi"],
                             MatrixOperations.ScalarMatrixProductNew(1.0 / 3.0, MatrixOperations.MatrixProduct(modificationMatrix, Bmatrices["ksi"])));
            double[,] BIhta = MatrixOperations.MatrixSubtraction(Bmatrices["ihta"],
                 MatrixOperations.ScalarMatrixProductNew(1.0 / 3.0, MatrixOperations.MatrixProduct(modificationMatrix, Bmatrices["ihta"])));
            double[,] BIhtaZita = MatrixOperations.MatrixSubtraction(Bmatrices["ihtaZita"],
                 MatrixOperations.ScalarMatrixProductNew(1.0 / 3.0, MatrixOperations.MatrixProduct(modificationMatrix, Bmatrices["ihtaZita"])));
            double[,] BZitaKsi = MatrixOperations.MatrixSubtraction(Bmatrices["zitaKsi"],
                 MatrixOperations.ScalarMatrixProductNew(1.0 / 3.0, MatrixOperations.MatrixProduct(modificationMatrix, Bmatrices["zitaKsi"])));
            BKsi = MatrixOperations.ΕliminateMatrixLine(BKsi, 4);
            BKsi = MatrixOperations.ΕliminateMatrixLine(BKsi, 5);
            BIhta = MatrixOperations.ΕliminateMatrixLine(BIhta, 4);
            BIhta = MatrixOperations.ΕliminateMatrixLine(BIhta, 5);
            matricesB.Add("ksi", BKsi);
            matricesB.Add("ihta", BIhta);
            matricesB.Add("ihtaZita", BIhtaZita);
            matricesB.Add("zitaKsi", BZitaKsi);
            return matricesB;
        }
        private double[,] CalculateStressStrainMatrix(double E, double v)
        {
            double[,] Ematrix = new double[6, 6];
            double Ehat = E / ((1.0 - 2.0 * v) * (1.0 + v));
            double G = (1.0 / 2.0) * (E / (1.0 + v));

            Ematrix[0, 0] = Ehat * (1.0 - v);
            Ematrix[0, 1] = Ehat * v;
            Ematrix[0, 2] = Ehat * v;
            Ematrix[1, 0] = Ehat * v;
            Ematrix[1, 1] = Ehat * (1.0 - v);
            Ematrix[1, 2] = Ehat * v;
            Ematrix[2, 0] = Ehat * v;
            Ematrix[2, 1] = Ehat * v;
            Ematrix[2, 2] = Ehat * (1.0 - v);
            Ematrix[3, 3] = G;
            Ematrix[4, 4] = G;
            Ematrix[5, 5] = G;
            return Ematrix;
        }
        private double[,] CalculateHourglassConstitutiveMatrix(double E, double v)
        {
            double[,] ConstitutiveMatrix = new double[6, 6];
            double G = (1.0 / 2.0) * (E / (1.0 + v));

            ConstitutiveMatrix[0, 0] = 4.0 / 3.0;
            ConstitutiveMatrix[0, 1] = -2.0 / 3.0;
            ConstitutiveMatrix[0, 2] = -2.0 / 3.0;
            ConstitutiveMatrix[1, 0] = -2.0 / 3.0;
            ConstitutiveMatrix[1, 1] = 4.0 / 3.0;
            ConstitutiveMatrix[1, 2] = -2.0 / 3.0;
            ConstitutiveMatrix[2, 0] = -2.0 / 3.0;
            ConstitutiveMatrix[2, 1] = -2.0 / 3.0;
            ConstitutiveMatrix[2, 2] = 4.0 / 3.0;
            ConstitutiveMatrix[3, 3] = 1.0;
            ConstitutiveMatrix[4, 4] = 1.0;
            ConstitutiveMatrix[5, 5] = 1.0;
            ConstitutiveMatrix = MatrixOperations.ScalarMatrixProductNew(G, ConstitutiveMatrix);
            return ConstitutiveMatrix;
        }
        private Tuple<double[], double[]> GaussPoints(int i, int j, int k)
        {
            double[] gaussPoints = new double[] { -1.0 / Math.Sqrt(3), 1.0 / Math.Sqrt(3) };
            double[] gaussWeights = new double[] { 1.0, 1.0 };
            //double[] gaussPointsThickness = new double[] { -Math.Sqrt(3.0 / 5.0), 0.0, Math.Sqrt(3.0 / 5.0) };
            //double[] gaussWeightsThickness = new double[] { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };
            double[] vectorWithPoints = new double[] { gaussPoints[i], gaussPoints[j], gaussPoints[k] };
            double[] vectorWithWeights = new double[] { gaussWeights[i], gaussWeights[j], gaussWeights[k] };
            return new Tuple<double[], double[]>(vectorWithPoints, vectorWithWeights);
        }
        private Tuple<double[], double[]> LobattoPoints(int i, int j, int k)
        {
            double[] gaussPoints = new double[] { -1.0, 1.0 };
            double[] gaussWeights = new double[] { 1.0, 1.0 };

            double[] vectorWithPoints = new double[] { gaussPoints[i], gaussPoints[j], gaussPoints[k] };
            double[] vectorWithWeights = new double[] { gaussWeights[i], gaussWeights[j], gaussWeights[k] };
            return new Tuple<double[], double[]>(vectorWithPoints, vectorWithWeights);
        }
        private double[,] CreateEnhancedStrainsInterpolationMatrix(double[] ksi)
        {
            double[,] M = new double[6, 7];
            M[0, 0] = ksi[0];
            M[1, 1] = ksi[1];
            M[2, 2] = ksi[2];
            M[2, 5] = ksi[0] * ksi[2];
            M[2, 6] = ksi[1] * ksi[2];
            M[3, 3] = ksi[0];
            M[3, 4] = ksi[1];
            return M;
        }
        private double[,] CalculateEnhancedStrainMatrixGamma(double[,] transformationMat0, double[,] M,
            double detJ0, double detJ)
        {
            double scalar = detJ0 / detJ;
            double[,] gamma = MatrixOperations.ScalarMatrixProductNew(scalar, MatrixOperations.MatrixProduct(transformationMat0, M));
            return gamma;
        }
        private double[,] TransformationMatrixTransposed(double[,] jacobianMatrix)
        {
            double[,] TransposedT = new double[6, 6];
            TransposedT[0, 0] = jacobianMatrix[0, 0] * jacobianMatrix[0, 0];
            TransposedT[0, 1] = jacobianMatrix[0, 1] * jacobianMatrix[0, 1];
            TransposedT[0, 2] = jacobianMatrix[0, 2] * jacobianMatrix[0, 2];
            TransposedT[0, 3] = jacobianMatrix[0, 0] * jacobianMatrix[0, 1];
            TransposedT[0, 4] = jacobianMatrix[0, 0] * jacobianMatrix[0, 2];
            TransposedT[0, 5] = jacobianMatrix[0, 1] * jacobianMatrix[0, 2];

            TransposedT[1, 0] = jacobianMatrix[1, 0] * jacobianMatrix[1, 0];
            TransposedT[1, 1] = jacobianMatrix[1, 1] * jacobianMatrix[1, 1];
            TransposedT[1, 2] = jacobianMatrix[1, 2] * jacobianMatrix[1, 2];
            TransposedT[1, 3] = jacobianMatrix[1, 0] * jacobianMatrix[1, 1];
            TransposedT[1, 4] = jacobianMatrix[1, 0] * jacobianMatrix[1, 2];
            TransposedT[1, 5] = jacobianMatrix[1, 1] * jacobianMatrix[1, 2];

            TransposedT[2, 0] = jacobianMatrix[2, 0] * jacobianMatrix[2, 0];
            TransposedT[2, 1] = jacobianMatrix[2, 1] * jacobianMatrix[2, 1];
            TransposedT[2, 2] = jacobianMatrix[2, 2] * jacobianMatrix[2, 2];
            TransposedT[2, 3] = jacobianMatrix[2, 0] * jacobianMatrix[2, 1];
            TransposedT[2, 4] = jacobianMatrix[2, 0] * jacobianMatrix[2, 2];
            TransposedT[2, 5] = jacobianMatrix[2, 1] * jacobianMatrix[2, 2];

            TransposedT[3, 0] = 2 * jacobianMatrix[0, 0] * jacobianMatrix[1, 0];
            TransposedT[3, 1] = 2 * jacobianMatrix[0, 1] * jacobianMatrix[1, 1];
            TransposedT[3, 2] = 2 * jacobianMatrix[0, 2] * jacobianMatrix[1, 2];

            TransposedT[3, 3] = jacobianMatrix[0, 0] * jacobianMatrix[1, 1] +
                                jacobianMatrix[0, 1] * jacobianMatrix[1, 0];

            TransposedT[3, 4] = jacobianMatrix[0, 0] * jacobianMatrix[1, 2] +
                                jacobianMatrix[0, 2] * jacobianMatrix[1, 0];

            TransposedT[3, 5] = jacobianMatrix[0, 1] * jacobianMatrix[1, 2] +
                                jacobianMatrix[0, 2] * jacobianMatrix[1, 1];

            TransposedT[4, 0] = 2 * jacobianMatrix[0, 0] * jacobianMatrix[2, 0];
            TransposedT[4, 1] = 2 * jacobianMatrix[0, 1] * jacobianMatrix[2, 1];
            TransposedT[4, 2] = 2 * jacobianMatrix[0, 2] * jacobianMatrix[2, 2];

            TransposedT[4, 3] = jacobianMatrix[0, 0] * jacobianMatrix[2, 1] +
                                jacobianMatrix[0, 1] * jacobianMatrix[2, 0];

            TransposedT[4, 4] = jacobianMatrix[0, 0] * jacobianMatrix[2, 2] +
                                jacobianMatrix[0, 2] * jacobianMatrix[2, 0];

            TransposedT[4, 5] = jacobianMatrix[0, 1] * jacobianMatrix[2, 2] +
                                jacobianMatrix[0, 2] * jacobianMatrix[2, 1];

            TransposedT[5, 0] = 2 * jacobianMatrix[1, 0] * jacobianMatrix[2, 0];
            TransposedT[5, 1] = 2 * jacobianMatrix[1, 1] * jacobianMatrix[2, 1];
            TransposedT[5, 2] = 2 * jacobianMatrix[1, 2] * jacobianMatrix[2, 2];

            TransposedT[5, 3] = jacobianMatrix[1, 0] * jacobianMatrix[2, 1] +
                                jacobianMatrix[1, 1] * jacobianMatrix[2, 0];

            TransposedT[5, 4] = jacobianMatrix[1, 0] * jacobianMatrix[2, 2] +
                                jacobianMatrix[1, 2] * jacobianMatrix[2, 0];

            TransposedT[5, 5] = jacobianMatrix[1, 1] * jacobianMatrix[2, 2] +
                                jacobianMatrix[1, 2] * jacobianMatrix[2, 1];
            return TransposedT;
        }
        public double[,] CreateGlobalStiffnessMatrix()
        {
            double[,] ansK = new double[24, 24];
            double[,] LeMatrix = new double[7, 24];
            double[,] DeMatrix = new double[7, 7];
            double[,] E = CalculateStressStrainMatrix(Properties.YoungMod, Properties.PoissonRatio);
            Dictionary<int, double[]> xNodalInitial = NodalCartesianCoordinatesInitial();
            Dictionary<int, double[]> ksiNodal = NodalNaturalCoordinates();
            Dictionary<int, double[]> hVectors = PositionVectors();
            double[] centerNaturalCoordinates = new double[3];
            var J0Matrices = JacobianMatrixDecomposed(xNodalInitial, ksiNodal, hVectors, centerNaturalCoordinates);
            var BMatricesNaturalSystem = CalculateNaturalSystemBANSMatricesDecomposed(xNodalInitial, ksiNodal, hVectors);
            var JInverse0 = CalculateInverseJacobian(J0Matrices["complete"]);
            double detJ0 = JInverse0.Item2;
            double[,] invJacobian0 = JInverse0.Item1;
            var transformationMatrices = StrainTransformationMatrixDecomposed(invJacobian0, J0Matrices["zero"],
                J0Matrices["ksi"], J0Matrices["ihta"], J0Matrices["zita"]);
            var BMatricesCartesianSystem = CalculateCartesianCoordinateSystemBANSMatricesDecomposed(BMatricesNaturalSystem, transformationMatrices);
            var BMatricesModified = HourglassComponentModification(BMatricesCartesianSystem);
            double[,] transformationMat0 = MatrixOperations.BlockMatrixInversion6X6(TransformationMatrixTransposed(J0Matrices["complete"]));
            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    for (int k = 0; k < 2; k++)
                    {
                        double[] gP = GaussPoints(i, j, k).Item1;
                        double[] gW = GaussPoints(i, j, k).Item2;
                        double[,] J = MatrixOperations.MatrixAddition(
                        MatrixOperations.MatrixAddition(
                        MatrixOperations.MatrixAddition(
                        MatrixOperations.MatrixAddition(
                        MatrixOperations.MatrixAddition(
                        MatrixOperations.MatrixAddition(J0Matrices["zero"],
                        MatrixOperations.ScalarMatrixProductNew(gP[2], J0Matrices["zita"])),
                        MatrixOperations.ScalarMatrixProductNew(gP[0], J0Matrices["ksi"])),
                        MatrixOperations.ScalarMatrixProductNew(gP[1], J0Matrices["ihta"])),
                        MatrixOperations.ScalarMatrixProductNew(gP[0] * gP[1], J0Matrices["ksiIhta"])),
                        MatrixOperations.ScalarMatrixProductNew(gP[0] * gP[2], J0Matrices["zitaKsi"])),
                        MatrixOperations.ScalarMatrixProductNew(gP[1] * gP[2], J0Matrices["ihtaZita"]));
                        double detJ = CalculateInverseJacobian(J).Item2;
                        double[,] B = MatrixOperations.MatrixAddition(
                        MatrixOperations.MatrixAddition(
                        MatrixOperations.MatrixAddition(
                        MatrixOperations.MatrixAddition(
                        MatrixOperations.MatrixAddition(MatrixOperations.MatrixAddition(BMatricesCartesianSystem["zero"],
                        MatrixOperations.ScalarMatrixProductNew(gP[2], BMatricesCartesianSystem["zita"])),
                        MatrixOperations.ScalarMatrixProductNew(gP[0], BMatricesCartesianSystem["ksi"])),
                        MatrixOperations.ScalarMatrixProductNew(gP[1], BMatricesCartesianSystem["ihta"])),
                        MatrixOperations.ScalarMatrixProductNew(gP[0] * gP[1], BMatricesCartesianSystem["ksiIhta"])),
                        MatrixOperations.ScalarMatrixProductNew(gP[0] * gP[2], BMatricesCartesianSystem["zitaKsi"])),
                        MatrixOperations.ScalarMatrixProductNew(gP[1] * gP[2], BMatricesCartesianSystem["ihtaZita"]));
                        ansK = MatrixOperations.MatrixAddition(ansK, MatrixOperations.ScalarMatrixProductNew(detJ * gW[0] * gW[1] * gW[2],
                            MatrixOperations.MatrixProduct(MatrixOperations.Transpose(B), MatrixOperations.MatrixProduct(E, B))));
                        double[,] Gamma = CalculateEnhancedStrainMatrixGamma(transformationMat0, CreateEnhancedStrainsInterpolationMatrix(gP),
                        detJ0, detJ);
                        LeMatrix = MatrixOperations.MatrixAddition(LeMatrix, MatrixOperations.ScalarMatrixProductNew(detJ * gW[0] * gW[1] * gW[2],
                            MatrixOperations.MatrixProduct(
                            MatrixOperations.Transpose(Gamma), MatrixOperations.MatrixProduct(E, B))));
                        DeMatrix = MatrixOperations.MatrixAddition(DeMatrix, MatrixOperations.ScalarMatrixProductNew(detJ * gW[0] * gW[1] * gW[2],
                            MatrixOperations.MatrixProduct(
                            MatrixOperations.Transpose(Gamma), MatrixOperations.MatrixProduct(E, Gamma))));
                    }
                }
            }
            double[,] DeMatrixInv = MatrixOperations.CalculateInverseMatrixGaussJordanMethod(DeMatrix);
            //DeMatrixInv[0, 0] = 1.0 / DeMatrix[0, 0];
            double[,] tangentMatrix =
                MatrixOperations.MatrixSubtraction(
                ansK,
                MatrixOperations.MatrixProduct(
                MatrixOperations.MatrixProduct(MatrixOperations.Transpose(LeMatrix), DeMatrixInv), LeMatrix));
            Properties.EASLMatrix = LeMatrix;
            Properties.EASDMatrix = DeMatrixInv;
            return tangentMatrix;
        }
        public double[,] CreateMassMatrix()
        {
            double[,] M = new double[24, 24];
            double[] nodalX = UpdateNodalCoordinates(DisplacementVector);
            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    for (int k = 0; k < 2; k++)
                    {
                        double[] lP = LobattoPoints(i, j, k).Item1;
                        double[] lW = LobattoPoints(i, j, k).Item2;
                        Dictionary<string, double[]> localdN = CalculateShapeFunctionsLocalDerivatives(lP);
                        double[,] J = CalculateFullJacobianMatrix(nodalX, localdN);
                        //double[,] invJ = CalculateInverseJacobian(J).Item1;
                        double detJ = CalculateInverseJacobian(J).Item2;
                        double[,] N = CalculateShapeFunctionMatrix(lP[0], lP[1], lP[2]);
                        M = MatrixOperations.MatrixAddition(M, MatrixOperations.ScalarMatrixProductNew(detJ * lW[0] * lW[1] * lW[2] * Properties.Density,
                            MatrixOperations.MatrixProduct(MatrixOperations.Transpose(N), N)));
                    }
                }
            }
            return M;
        }
        public double[,] CreateDampingMatrix()
        {
            return new double[24, 24];
        }
        private double[] CalculateStressVector(double[,] E, double[] strain)
        {
            double[] stressVector = VectorOperations.MatrixVectorProduct(E, strain);
            return stressVector;
        }
        public double[] CreateInternalGlobalForcesVector()
        {
            double[] VectorRe = new double[24];
            double[] PeVector = new double[7];
            double[,] matrixLe = Properties.EASLMatrix;
            double[,] matrixDeInv = Properties.EASDMatrix;
            double[,] E = CalculateStressStrainMatrix(Properties.YoungMod, Properties.PoissonRatio);
            Dictionary<int, double[]> xNodalInitial = NodalCartesianCoordinatesInitial();
            Dictionary<int, double[]> ksiNodal = NodalNaturalCoordinates();
            Dictionary<int, double[]> hVectors = PositionVectors();
            double[] centerNaturalCoordinates = new double[3];
            var J0Matrices = JacobianMatrixDecomposed(xNodalInitial, ksiNodal, hVectors, centerNaturalCoordinates);
            var BMatricesNaturalSystem = CalculateNaturalSystemBANSMatricesDecomposed(xNodalInitial, ksiNodal, hVectors);
            var JInverse0 = CalculateInverseJacobian(J0Matrices["complete"]);
            double detJ0 = JInverse0.Item2;
            double[,] invJacobian0 = JInverse0.Item1;
            var transformationMatrices = StrainTransformationMatrixDecomposed(invJacobian0, J0Matrices["zero"],
                J0Matrices["ksi"], J0Matrices["ihta"], J0Matrices["zita"]);
            var BMatricesCartesianSystem = CalculateCartesianCoordinateSystemBANSMatricesDecomposed(BMatricesNaturalSystem, transformationMatrices);
            double[,] transformationMat0 = MatrixOperations.BlockMatrixInversion6X6(TransformationMatrixTransposed(J0Matrices["complete"]));
            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    for (int k = 0; k < 2; k++)
                    {
                        double[] gP = GaussPoints(i, j, k).Item1;
                        double[] gW = GaussPoints(i, j, k).Item2;
                        double[,] J = MatrixOperations.MatrixAddition(
                            MatrixOperations.MatrixAddition(
                            MatrixOperations.MatrixAddition(
                            MatrixOperations.MatrixAddition(
                            MatrixOperations.MatrixAddition(
                            MatrixOperations.MatrixAddition(J0Matrices["zero"],
                            MatrixOperations.ScalarMatrixProductNew(gP[2], J0Matrices["zita"])),
                            MatrixOperations.ScalarMatrixProductNew(gP[0], J0Matrices["ksi"])),
                            MatrixOperations.ScalarMatrixProductNew(gP[1], J0Matrices["ihta"])),
                            MatrixOperations.ScalarMatrixProductNew(gP[0] * gP[1], J0Matrices["ksiIhta"])),
                            MatrixOperations.ScalarMatrixProductNew(gP[0] * gP[2], J0Matrices["zitaKsi"])),
                            MatrixOperations.ScalarMatrixProductNew(gP[1] * gP[2], J0Matrices["ihtaZita"]));
                        double detJ = CalculateInverseJacobian(J).Item2;
                        double[,] B = MatrixOperations.MatrixAddition(
                            MatrixOperations.MatrixAddition(
                            MatrixOperations.MatrixAddition(
                            MatrixOperations.MatrixAddition(
                            MatrixOperations.MatrixAddition(MatrixOperations.MatrixAddition(BMatricesCartesianSystem["zero"],
                            MatrixOperations.ScalarMatrixProductNew(gP[2], BMatricesCartesianSystem["zita"])),
                            MatrixOperations.ScalarMatrixProductNew(gP[0], BMatricesCartesianSystem["ksi"])),
                            MatrixOperations.ScalarMatrixProductNew(gP[1], BMatricesCartesianSystem["ihta"])),
                            MatrixOperations.ScalarMatrixProductNew(gP[0] * gP[1], BMatricesCartesianSystem["ksiIhta"])),
                            MatrixOperations.ScalarMatrixProductNew(gP[0] * gP[2], BMatricesCartesianSystem["zitaKsi"])),
                            MatrixOperations.ScalarMatrixProductNew(gP[1] * gP[2], BMatricesCartesianSystem["ihtaZita"]));
                        double[,] Gamma = CalculateEnhancedStrainMatrixGamma(transformationMat0, CreateEnhancedStrainsInterpolationMatrix(gP),
                        detJ0, detJ);
                        double[] alphaVector = Properties.EASVector;
                        double[] EnhStrainVector = VectorOperations.MatrixVectorProduct(Gamma, alphaVector);
                        double[] modifiedStrainVector = VectorOperations.VectorVectorAddition(CalculateStrainsVector(B), EnhStrainVector);
                        double[] stressVectorModified = CalculateStressVector(E, modifiedStrainVector);
                        VectorRe = VectorOperations.VectorVectorAddition(VectorRe, VectorOperations.VectorScalarProductNew(
                            VectorOperations.MatrixVectorProduct(MatrixOperations.Transpose(B), stressVectorModified), detJ * gW[0] * gW[1] * gW[2]));
                        PeVector = VectorOperations.VectorVectorAddition(PeVector,
                           VectorOperations.VectorScalarProductNew(
                           VectorOperations.MatrixVectorProduct(MatrixOperations.Transpose(Gamma), stressVectorModified),
                           detJ * gW[0] * gW[1] * gW[2]));
                    }
                }
            }
            double[] fInt = VectorOperations.VectorVectorSubtraction(VectorRe,
            VectorOperations.MatrixVectorProduct(MatrixOperations.MatrixProduct(MatrixOperations.Transpose(matrixLe), matrixDeInv), PeVector));
            Properties.EASFEnhancedVector = PeVector;
            return fInt;
        }
    }
}