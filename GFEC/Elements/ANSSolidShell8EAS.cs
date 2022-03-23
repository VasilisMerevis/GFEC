using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GFEC
{
    class ANSSolidShell8EAS : IElement
    {
        public Dictionary<int, INode> Nodes { get; }
        public IElementProperties Properties { get; set; }
        public Dictionary<int, bool[]> ElementFreedomSignature { get; } = new Dictionary<int, bool[]>();
        public List<int> ElementFreedomList { get; set; }
        public double[] DisplacementVector { get; set; }
        public double[] DisplacementVectorPreviousStep { get; set; }
        public double[] EASVector { get; set; }
        public double[] EASFEnhancedVector { get; set; }
        public double[,] EASLMatrix { get; set; }
        public double[,] EASDMatrix { get; set; }
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
        public ANSSolidShell8EAS(IElementProperties properties, Dictionary<int, INode> nodes)
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
            Properties.DisplacementVectorPreviousStep = new double[24];
            Properties.EASVector = new double[4];
            Properties.EASFEnhancedVector = new double[4];
            Properties.EASLMatrix = new double[4, 24];
            Properties.EASDMatrix = new double[4, 4];
        }
        public void CalculateElementEASMatrices()
        {
            double[,] LeMatrix = new double[4, 24];
            double[,] DeMatrix = new double[4, 4];
            double[] PeVector = new double[4];
            double[,] E = CalculateStressStrainMatrix(Properties.YoungMod, Properties.PoissonRatio);
            //double[] nodalX = UpdateNodalCoordinates(DisplacementVector);
            double[] nodalXInitial = InitialNodalCoordinates();
            double[] center = new double[3];
            double[,] J0 = CalculateJacobian(nodalXInitial, CalculateShapeFunctionsLocalDerivatives(center));
            double detJ0 = CalculateJacobianDet(J0);
            double[,] transformationMat0 = MatrixOperations.BlockMatrixInversion6X6(TransformationMatrixTransposed(J0));
            Dictionary<int, Dictionary<string, double[]>> aStrainPoints = AssumedStrainSamplingPoints(nodalXInitial);
            Dictionary<int, Dictionary<string, double[]>> aStrainPointsDU = AssumedStrainSamplingPointsDU();

            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    for (int k = 0; k < 2; k++)
                    {
                        double[] gP = GaussPoints(i, j, k).Item1;
                        double[] gW = GaussPoints(i, j, k).Item2;
                        Dictionary<string, double[]> localdN = CalculateShapeFunctionsLocalDerivatives(gP);
                        double[,] J = CalculateJacobian(nodalXInitial, localdN);
                        double[,] transformationMat = MatrixOperations.BlockMatrixInversion6X6(TransformationMatrixTransposed(J));
                        double detJ = CalculateInverseJacobian(J).Item2;
                        Dictionary<int, double[]> gVectors = CalculategVectors(nodalXInitial, localdN);
                        Dictionary<int, double[]> gVectorsdU = CalculatedgdUVectors(localdN);

                        var BMatrices = CalculateBANSMatrix(gP, localdN, transformationMat, gVectors, gVectorsdU,
                            aStrainPoints, aStrainPointsDU);
                        double[,] B = BMatrices.Item1;
                        double[,] Gamma = CalculateEnhancedStrainMatrixGamma(transformationMat0, CreateEnhancedStrainsInterpolationMatrix(gP),
                            detJ0, detJ);
                        LeMatrix = MatrixOperations.MatrixAddition(LeMatrix, MatrixOperations.ScalarMatrixProductNew(detJ * gW[0] * gW[1] * gW[2],
                            MatrixOperations.MatrixProduct(
                            MatrixOperations.Transpose(Gamma), MatrixOperations.MatrixProduct(E, B))));
                        DeMatrix = MatrixOperations.MatrixAddition(DeMatrix, MatrixOperations.ScalarMatrixProductNew(detJ * gW[0] * gW[1] * gW[2],
                            MatrixOperations.MatrixProduct(
                            MatrixOperations.Transpose(Gamma), MatrixOperations.MatrixProduct(E, Gamma))));

                        double[] alphaVector = Properties.EASVector;
                        double[] EnhStrainVector = VectorOperations.MatrixVectorProduct(Gamma, alphaVector);
                        double[] modifiedStrainVector = VectorOperations.VectorVectorAddition(CalculateStrainsVector(B), EnhStrainVector);
                        double[] modifiedStressVector = CalculateStressVector(E, modifiedStrainVector);
                        PeVector = VectorOperations.VectorVectorAddition(PeVector,
                                   VectorOperations.VectorScalarProductNew(
                                   VectorOperations.MatrixVectorProduct(MatrixOperations.Transpose(Gamma), modifiedStressVector),
                                   detJ * gW[0] * gW[1] * gW[2]));
                    }
                }
            }
            double[,] DeMatrixInv = MatrixOperations.BlockMatrixInversion4X4(DeMatrix);
            Properties.EASLMatrix = LeMatrix;
            Properties.EASDMatrix = DeMatrixInv;
            Properties.EASFEnhancedVector = PeVector;
        }
        public void InitializeElementEASParameters()
        {
            Properties.EASVector[0] = 0.0;
            Properties.EASVector[1] = 0.0;
            Properties.EASVector[2] = 0.0;
            Properties.EASVector[3] = 0.0;
        }
        public void UpdateElementEASParameters(double[] solutionVector)
        {
            double[] aPrevious = Properties.EASVector;
            double[] vectorPe = Properties.EASFEnhancedVector;
            double[] uPrevious = Properties.DisplacementVectorPreviousStep;
            double[,] matrixLe = Properties.EASLMatrix;
            double[,] matrixDeInv = Properties.EASDMatrix;
            double[] dU = VectorOperations.VectorVectorSubtraction(solutionVector, uPrevious);
            double[] aNext = VectorOperations.VectorVectorSubtraction(aPrevious,
                            VectorOperations.MatrixVectorProduct(matrixDeInv,
                            VectorOperations.VectorVectorAddition(
                            VectorOperations.MatrixVectorProduct(matrixLe, dU), vectorPe)));
            Properties.EASVector = aNext;
        }
        public void StoreElementFinalStepDisplacementVector(double[] solutionVector)
        {
            Properties.DisplacementVectorPreviousStep = solutionVector;
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
            //List<double[]> GpointsStress = new List<double[]>();
            //double[,] E = CalculateStressStrainMatrix(Properties.YoungMod, Properties.PoissonRatio);
            //for (int i = 0; i < 2; i++)
            //{
            //    for (int j = 0; j < 2; j++)
            //    {
            //        for (int k = 0; k < 2; k++)
            //        {
            //            double[] gP = GaussPoints(i, j, k).Item1;
            //            double[] gW = GaussPoints(i, j, k).Item2;
            //            Dictionary<string, double[]> localdN = CalculateShapeFunctionsLocalDerivatives(gP);
            //            double[,] J = CalculateJacobian(localdN);
            //            double[,] invJ = CalculateInverseJacobian(J).Item1;
            //            double detJ = CalculateInverseJacobian(J).Item2;
            //            Dictionary<int, double[]> globaldN = CalculateShapeFunctionsGlobalDerivatives(localdN, invJ);
            //            double[,] B = CalculateBANSMatrix(globaldN);
            //            double[] strainVector = CalculateStrainsVector(B);
            //            double[] stressVector = CalculateStressVector(E, strainVector);
            //            GpointsStress.Add(stressVector);
            //        }
            //    }
            //}
            //return GpointsStress;
            throw new Exception("Method not implemenented");
        }
        public List<double[]> GetStrainVector()
        {
            //List<double[]> GpointsDeformation = new List<double[]>();
            //for (int i = 0; i < 2; i++)
            //{
            //    for (int j = 0; j < 2; j++)
            //    {
            //        for (int k = 0; k < 2; k++)
            //        {
            //            double[] gP = GaussPoints(i, j, k).Item1;
            //            double[] gW = GaussPoints(i, j, k).Item2;
            //            Dictionary<string, double[]> localdN = CalculateShapeFunctionsLocalDerivatives(gP);
            //            double[,] J = CalculateJacobian(localdN);
            //            double[,] invJ = CalculateInverseJacobian(J).Item1;
            //            double detJ = CalculateInverseJacobian(J).Item2;
            //            Dictionary<int, double[]> globaldN = CalculateShapeFunctionsGlobalDerivatives(localdN, invJ);
            //            double[,] B = CalculateBANSMatrix(globaldN);
            //            double[] strainVector = CalculateStrainsVector(B);
            //            GpointsDeformation.Add(strainVector);
            //        }
            //    }
            //}
            //return GpointsDeformation;
            throw new Exception("Method not implemenented");
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
        public List<double[]> GetStressFromElements(List<double[]> parametricCoordinatesVector)
        {
            //List<double[]> StessVectorsList = new List<double[]>();
            //double[,] E = CalculateStressStrainMatrix(Properties.YoungMod, Properties.PoissonRatio);
            ////double[] xUpdated = UpdateNodalCoordinates(DisplacementVector);
            //for (int i = 0; i < parametricCoordinatesVector.Count; i++)
            //{
            //    Dictionary<string, double[]> localdN = CalculateShapeFunctionsLocalDerivatives(parametricCoordinatesVector[i]);
            //    double[,] J = CalculateJacobian(localdN);
            //    double[,] invJ = CalculateInverseJacobian(J).Item1;
            //    Dictionary<int, double[]> globaldN = CalculateShapeFunctionsGlobalDerivatives(localdN, invJ);
            //    double[,] B = CalculateBANSMatrix(globaldN);
            //    double[] strainVector = CalculateStrainsVector(B);
            //    double[] stressVector = CalculateStressVector(E, strainVector);
            //    StessVectorsList.Add(stressVector);
            //}
            //return StessVectorsList;
            throw new Exception("Method not implemenented");
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
            double[] InitialCoor = new double[24];
            for (int i = 1; i <= 8; i++)
            {
                InitialCoor[3 * i - 3] = Nodes[i].XCoordinate;
                InitialCoor[3 * i - 2] = Nodes[i].YCoordinate;
                InitialCoor[3 * i - 1] = Nodes[i].ZCoordinate;
            }
            return InitialCoor;
        }
        private Dictionary<int, double> CalculateShapeFunctions(double ksi, double ihta, double mhi)
        {
            Dictionary<int, double> shapeFunctions = new Dictionary<int, double>();
            double N1 = 1.0 / 8.0 * (1 - ksi) * (1 - ihta) * (1 - mhi); shapeFunctions.Add(1, N1);
            double N2 = 1.0 / 8.0 * (1 + ksi) * (1 - ihta) * (1 - mhi); shapeFunctions.Add(2, N2);
            double N3 = 1.0 / 8.0 * (1 + ksi) * (1 + ihta) * (1 - mhi); shapeFunctions.Add(3, N3);
            double N4 = 1.0 / 8.0 * (1 - ksi) * (1 + ihta) * (1 - mhi); shapeFunctions.Add(4, N4);
            double N5 = 1.0 / 8.0 * (1 - ksi) * (1 - ihta) * (1 + mhi); shapeFunctions.Add(5, N5);
            double N6 = 1.0 / 8.0 * (1 + ksi) * (1 - ihta) * (1 + mhi); shapeFunctions.Add(6, N6);
            double N7 = 1.0 / 8.0 * (1 + ksi) * (1 + ihta) * (1 + mhi); shapeFunctions.Add(7, N7);
            double N8 = 1.0 / 8.0 * (1 - ksi) * (1 + ihta) * (1 + mhi); shapeFunctions.Add(8, N8);

            return shapeFunctions;
        }
        private double[,] CalculateShapeFunctionMatrix(double ksi, double ihta, double zita)
        {
            double N1 = 1.0 / 8.0 * (1 - ksi) * (1 - ihta) * (1 - zita);
            double N2 = 1.0 / 8.0 * (1 + ksi) * (1 - ihta) * (1 - zita);
            double N3 = 1.0 / 8.0 * (1 + ksi) * (1 + ihta) * (1 - zita);
            double N4 = 1.0 / 8.0 * (1 - ksi) * (1 + ihta) * (1 - zita);
            double N5 = 1.0 / 8.0 * (1 - ksi) * (1 - ihta) * (1 + zita);
            double N6 = 1.0 / 8.0 * (1 + ksi) * (1 - ihta) * (1 + zita);
            double N7 = 1.0 / 8.0 * (1 + ksi) * (1 + ihta) * (1 + zita);
            double N8 = 1.0 / 8.0 * (1 - ksi) * (1 + ihta) * (1 + zita);
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
                (-1.0/8.0*(1-ihta)*(1-mhi)),
                (1.0/8.0*(1-ihta)*(1-mhi)),
                (1.0/8.0*(1+ihta)*(1-mhi)),
                (-1.0/8.0*(1+ihta)*(1-mhi)),
                (-1.0/8.0*(1-ihta)*(1+mhi)),
                (1.0/8.0*(1-ihta)*(1+mhi)),
                (1.0/8.0*(1+ihta)*(1+mhi)),
                (-1.0/8.0*(1+ihta)*(1+mhi))
            };

            double[] dN_ihta = new double[]
            {
                (-1.0/8.0*(1-ksi)*(1-mhi)),
                (-1.0/8.0*(1+ksi)*(1-mhi)),
                (1.0/8.0*(1+ksi)*(1-mhi)),
                (1.0/8.0*(1-ksi)*(1-mhi)),
                (-1.0/8.0*(1-ksi)*(1+mhi)),
                (-1.0/8.0*(1+ksi)*(1+mhi)),
                (1.0/8.0*(1+ksi)*(1+mhi)),
                (1.0/8.0*(1-ksi)*(1+mhi))
            };

            double[] dN_mhi = new double[]
            {
                (-1.0/8.0*(1-ksi)*(1-ihta)),
                (-1.0/8.0*(1+ksi)*(1-ihta)),
                (-1.0/8.0*(1+ksi)*(1+ihta)),
                (-1.0/8.0*(1-ksi)*(1+ihta)),
                (1.0/8.0*(1-ksi)*(1-ihta)),
                (1.0/8.0*(1+ksi)*(1-ihta)),
                (1.0/8.0*(1+ksi)*(1+ihta)),
                (1.0/8.0*(1-ksi)*(1+ihta))
            };

            Dictionary<string, double[]> dN = new Dictionary<string, double[]>();
            dN.Add("ksi", dN_ksi);
            dN.Add("ihta", dN_ihta);
            dN.Add("mhi", dN_mhi);
            return dN;
        }

        private double[,] CalculateJacobian(double[] xNodal, Dictionary<string, double[]> dN)
        {
            double[,] jacobianMatrix = new double[3, 3];
            //DisplacementVector = new double[24];
            //double[] xInitial = UpdateNodalCoordinates(DisplacementVector);

            int k = 0;
            for (int i = 0; i < 8; i++)
            {
                jacobianMatrix[0, 0] = jacobianMatrix[0, 0] + xNodal[k] * dN["ksi"][i];
                k = k + 3;
            }
            k = 1;
            for (int i = 0; i < 8; i++)
            {
                jacobianMatrix[0, 1] = jacobianMatrix[0, 1] + xNodal[k] * dN["ksi"][i];
                k = k + 3;
            }
            k = 2;
            for (int i = 0; i < 8; i++)
            {
                jacobianMatrix[0, 2] = jacobianMatrix[0, 2] + xNodal[k] * dN["ksi"][i];
                k = k + 3;
            }

            k = 0;
            for (int i = 0; i < 8; i++)
            {
                jacobianMatrix[1, 0] = jacobianMatrix[1, 0] + xNodal[k] * dN["ihta"][i];
                k = k + 3;
            }
            k = 1;
            for (int i = 0; i < 8; i++)
            {
                jacobianMatrix[1, 1] = jacobianMatrix[1, 1] + xNodal[k] * dN["ihta"][i];
                k = k + 3;
            }
            k = 2;
            for (int i = 0; i < 8; i++)
            {
                jacobianMatrix[1, 2] = jacobianMatrix[1, 2] + xNodal[k] * dN["ihta"][i];
                k = k + 3;
            }

            k = 0;
            for (int i = 0; i < 8; i++)
            {
                jacobianMatrix[2, 0] = jacobianMatrix[2, 0] + xNodal[k] * dN["mhi"][i];
                k = k + 3;
            }
            k = 1;
            for (int i = 0; i < 8; i++)
            {
                jacobianMatrix[2, 1] = jacobianMatrix[2, 1] + xNodal[k] * dN["mhi"][i];
                k = k + 3;
            }
            k = 2;
            for (int i = 0; i < 8; i++)
            {
                jacobianMatrix[2, 2] = jacobianMatrix[2, 2] + xNodal[k] * dN["mhi"][i];
                k = k + 3;
            }

            return jacobianMatrix;
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
        private Dictionary<int, double[]> CalculateShapeFunctionsGlobalDerivatives(Dictionary<string, double[]> dN, double[,] Jinv)
        {
            Dictionary<int, double[]> dNg = new Dictionary<int, double[]>();

            for (int i = 0; i < 8; i++)
            {
                double[] dNlocal = new double[] { dN["ksi"][i], dN["ihta"][i], dN["mhi"][i] };
                double[] dNglobal = VectorOperations.MatrixVectorProduct(Jinv, dNlocal);
                dNg.Add(i, dNglobal);
            }
            return dNg;
        }

        private double[] CalculateStrainsVector(double[,] Bmatrix)
        {
            double[] strains = VectorOperations.MatrixVectorProduct(Bmatrix, DisplacementVector);
            return strains;
        }
        private Dictionary<int, double[]> CalculategVectors(double[] nodalXInitial, Dictionary<string, double[]> localdN)
        {
            Dictionary<int, double[]> gVectors = new Dictionary<int, double[]>();
            double g11 = new double();
            double g12 = new double();
            double g13 = new double();
            double g21 = new double();
            double g22 = new double();
            double g23 = new double();
            double g31 = new double();
            double g32 = new double();
            double g33 = new double();
            for (int i = 0; i < 8; i++)
            {
                g11 += localdN["ksi"][i] * (nodalXInitial[3 * i] + DisplacementVector[3 * i]);
                g21 += localdN["ihta"][i] * (nodalXInitial[3 * i] + DisplacementVector[3 * i]);
                g31 += localdN["mhi"][i] * (nodalXInitial[3 * i] + DisplacementVector[3 * i]);
                g12 += localdN["ksi"][i] * (nodalXInitial[3 * i + 1] + DisplacementVector[3 * i + 1]);
                g22 += localdN["ihta"][i] * (nodalXInitial[3 * i + 1] + DisplacementVector[3 * i + 1]);
                g32 += localdN["mhi"][i] * (nodalXInitial[3 * i + 1] + DisplacementVector[3 * i + 1]);
                g13 += localdN["ksi"][i] * (nodalXInitial[3 * i + 2] + DisplacementVector[3 * i + 2]);
                g23 += localdN["ihta"][i] * (nodalXInitial[3 * i + 2] + DisplacementVector[3 * i + 2]);
                g33 += localdN["mhi"][i] * (nodalXInitial[3 * i + 2] + DisplacementVector[3 * i + 2]);
            }
            double[] g1 = new double[] { g11, g12, g13 };
            double[] g2 = new double[] { g21, g22, g23 };
            double[] g3 = new double[] { g31, g32, g33 };
            gVectors.Add(1, g1);
            gVectors.Add(2, g2);
            gVectors.Add(3, g3);
            return gVectors;
        }
        private Dictionary<int, double[]> CalculatedgdUVectors(Dictionary<string, double[]> localdN)
        {
            Dictionary<int, double[]> dgVectors = new Dictionary<int, double[]>();
            double dg11 = new double();
            double dg12 = new double();
            double dg13 = new double();
            double dg21 = new double();
            double dg22 = new double();
            double dg23 = new double();
            double dg31 = new double();
            double dg32 = new double();
            double dg33 = new double();
            for (int i = 0; i < 8; i++)
            {
                dg11 += localdN["ksi"][i];
                dg21 += localdN["ihta"][i];
                dg31 += localdN["mhi"][i];
                dg12 += localdN["ksi"][i];
                dg22 += localdN["ihta"][i];
                dg32 += localdN["mhi"][i];
                dg13 += localdN["ksi"][i];
                dg23 += localdN["ihta"][i];
                dg33 += localdN["mhi"][i];
            }
            double[] dg1 = new double[] { dg11, dg12, dg13 };
            double[] dg2 = new double[] { dg21, dg22, dg23 };
            double[] dg3 = new double[] { dg31, dg32, dg33 };
            dgVectors.Add(1, dg1);
            dgVectors.Add(2, dg2);
            dgVectors.Add(3, dg3);
            return dgVectors;
        }
        private Tuple<double[,], double[,]> CalculateBANSMatrix(double[] GaussPointNaturalCoordinates, Dictionary<string, double[]> localdNGaussPoint,
                                              double[,] transformationMatrix,Dictionary<int, double[]> gVectors,
                                              Dictionary<int, double[]> gVectorsDU,
                                              Dictionary<int, Dictionary<string, double[]>> StrainSamplingPoints,
                                              Dictionary<int, Dictionary<string, double[]>> StrainSamplingPointsDU)
        {
            double[,] BANSmatrix = new double[6, 24];
            double[,] BANSmatrixDU = new double[6, 24];
            for (int i = 0; i < 8; i++)
            {
                double[,] nodalB = new double[6, 3];
                double[,] nodalBdU = new double[6, 3];

                nodalB[0, 0] = localdNGaussPoint["ksi"][i] * gVectors[1][0];
                nodalB[1, 0] = localdNGaussPoint["ihta"][i] * gVectors[2][0];
                nodalB[2, 0] = localdNGaussPoint["mhi"][i] * gVectors[3][0];
                nodalB[0, 1] = localdNGaussPoint["ksi"][i] * gVectors[1][1];
                nodalB[1, 1] = localdNGaussPoint["ihta"][i] * gVectors[2][1];
                nodalB[2, 1] = localdNGaussPoint["mhi"][i] * gVectors[3][1];
                nodalB[0, 2] = localdNGaussPoint["ksi"][i] * gVectors[1][2];
                nodalB[1, 2] = localdNGaussPoint["ihta"][i] * gVectors[2][2];
                nodalB[2, 2] = localdNGaussPoint["mhi"][i] * gVectors[3][2];

                nodalB[3, 0] = localdNGaussPoint["ksi"][i] * gVectors[2][0] +
                               localdNGaussPoint["ihta"][i] * gVectors[1][0];
                nodalB[3, 1] = localdNGaussPoint["ksi"][i] * gVectors[2][1] +
                               localdNGaussPoint["ihta"][i] * gVectors[1][1];
                nodalB[3, 2] = localdNGaussPoint["ksi"][i] * gVectors[2][2] +
                               localdNGaussPoint["ihta"][i] * gVectors[1][2];

                nodalB[4, 0] = (1.0 / 2.0) * (1.0 - GaussPointNaturalCoordinates[1]) *
                               StrainSamplingPoints[i]["A"][0] +
                               (1.0 / 2.0) * (1.0 + GaussPointNaturalCoordinates[1]) *
                               StrainSamplingPoints[i]["C"][0];

                nodalB[4, 1] = (1.0 / 2.0) * (1.0 - GaussPointNaturalCoordinates[1]) *
                               StrainSamplingPoints[i]["A"][1] +
                               (1.0 / 2.0) * (1.0 + GaussPointNaturalCoordinates[1]) *
                               StrainSamplingPoints[i]["C"][1];

                nodalB[4, 2] = (1.0 / 2.0) * (1.0 - GaussPointNaturalCoordinates[1]) *
                               StrainSamplingPoints[i]["A"][2] +
                               (1.0 / 2.0) * (1.0 + GaussPointNaturalCoordinates[1]) *
                               StrainSamplingPoints[i]["C"][2];

                nodalB[5, 0] = (1.0 / 2.0) * (1.0 - GaussPointNaturalCoordinates[0]) *
                               StrainSamplingPoints[i]["D"][0] +
                               (1.0 / 2.0) * (1.0 + GaussPointNaturalCoordinates[0]) *
                               StrainSamplingPoints[i]["B"][0];

                nodalB[5, 1] = (1.0 / 2.0) * (1.0 - GaussPointNaturalCoordinates[0]) *
                               StrainSamplingPoints[i]["D"][1] +
                               (1.0 / 2.0) * (1.0 + GaussPointNaturalCoordinates[0]) *
                               StrainSamplingPoints[i]["B"][1];

                nodalB[5, 2] = (1.0 / 2.0) * (1.0 - GaussPointNaturalCoordinates[0]) *
                               StrainSamplingPoints[i]["D"][2] +
                               (1.0 / 2.0) * (1.0 + GaussPointNaturalCoordinates[0]) *
                               StrainSamplingPoints[i]["B"][2];

                double[,] nodalBTransformed = MatrixOperations.MatrixProduct(transformationMatrix, nodalB);
                BANSmatrix[0, i * 3] = nodalBTransformed[0, 0];
                BANSmatrix[1, i * 3] = nodalBTransformed[1, 0];
                BANSmatrix[2, i * 3] = nodalBTransformed[2, 0];
                BANSmatrix[3, i * 3] = nodalBTransformed[3, 0];
                BANSmatrix[4, i * 3] = nodalBTransformed[4, 0];
                BANSmatrix[5, i * 3] = nodalBTransformed[5, 0];
                BANSmatrix[0, i * 3 + 1] = nodalBTransformed[0, 1];
                BANSmatrix[1, i * 3 + 1] = nodalBTransformed[1, 1];
                BANSmatrix[2, i * 3 + 1] = nodalBTransformed[2, 1];
                BANSmatrix[3, i * 3 + 1] = nodalBTransformed[3, 1];
                BANSmatrix[4, i * 3 + 1] = nodalBTransformed[4, 1];
                BANSmatrix[5, i * 3 + 1] = nodalBTransformed[5, 1];
                BANSmatrix[0, i * 3 + 2] = nodalBTransformed[0, 2];
                BANSmatrix[1, i * 3 + 2] = nodalBTransformed[1, 2];
                BANSmatrix[2, i * 3 + 2] = nodalBTransformed[2, 2];
                BANSmatrix[3, i * 3 + 2] = nodalBTransformed[3, 2];
                BANSmatrix[4, i * 3 + 2] = nodalBTransformed[4, 2];
                BANSmatrix[5, i * 3 + 2] = nodalBTransformed[5, 2];
                //
                //-------------------------------------------------------------------
                //
                nodalBdU[0, 0] = localdNGaussPoint["ksi"][i] * gVectorsDU[1][0];
                nodalBdU[1, 0] = localdNGaussPoint["ihta"][i] * gVectorsDU[2][0];
                nodalBdU[2, 0] = localdNGaussPoint["mhi"][i] * gVectorsDU[3][0];
                nodalBdU[0, 1] = localdNGaussPoint["ksi"][i] * gVectorsDU[1][1];
                nodalBdU[1, 1] = localdNGaussPoint["ihta"][i] * gVectorsDU[2][1];
                nodalBdU[2, 1] = localdNGaussPoint["mhi"][i] * gVectorsDU[3][1];
                nodalBdU[0, 2] = localdNGaussPoint["ksi"][i] * gVectorsDU[1][2];
                nodalBdU[1, 2] = localdNGaussPoint["ihta"][i] * gVectorsDU[2][2];
                nodalBdU[2, 2] = localdNGaussPoint["mhi"][i] * gVectorsDU[3][2];

                nodalBdU[3, 0] = localdNGaussPoint["ksi"][i] * gVectorsDU[2][0] +
                               localdNGaussPoint["ihta"][i] * gVectorsDU[1][0];
                nodalBdU[3, 1] = localdNGaussPoint["ksi"][i] * gVectorsDU[2][1] +
                               localdNGaussPoint["ihta"][i] * gVectorsDU[1][1];
                nodalBdU[3, 2] = localdNGaussPoint["ksi"][i] * gVectorsDU[2][2] +
                               localdNGaussPoint["ihta"][i] * gVectorsDU[1][2];

                nodalBdU[4, 0] = (1.0 / 2.0) * (1.0 - GaussPointNaturalCoordinates[1]) *
                               StrainSamplingPointsDU[i]["A"][0] +
                               (1.0 / 2.0) * (1.0 + GaussPointNaturalCoordinates[1]) *
                               StrainSamplingPointsDU[i]["C"][0];

                nodalBdU[4, 1] = (1.0 / 2.0) * (1.0 - GaussPointNaturalCoordinates[1]) *
                               StrainSamplingPointsDU[i]["A"][1] +
                               (1.0 / 2.0) * (1.0 + GaussPointNaturalCoordinates[1]) *
                               StrainSamplingPointsDU[i]["C"][1];

                nodalBdU[4, 2] = (1.0 / 2.0) * (1.0 - GaussPointNaturalCoordinates[1]) *
                               StrainSamplingPointsDU[i]["A"][2] +
                               (1.0 / 2.0) * (1.0 + GaussPointNaturalCoordinates[1]) *
                               StrainSamplingPointsDU[i]["C"][2];

                nodalBdU[5, 0] = (1.0 / 2.0) * (1.0 - GaussPointNaturalCoordinates[0]) *
                               StrainSamplingPointsDU[i]["D"][0] +
                               (1.0 / 2.0) * (1.0 + GaussPointNaturalCoordinates[0]) *
                               StrainSamplingPointsDU[i]["B"][0];

                nodalBdU[5, 1] = (1.0 / 2.0) * (1.0 - GaussPointNaturalCoordinates[0]) *
                               StrainSamplingPointsDU[i]["D"][1] +
                               (1.0 / 2.0) * (1.0 + GaussPointNaturalCoordinates[0]) *
                               StrainSamplingPointsDU[i]["B"][1];

                nodalBdU[5, 2] = (1.0 / 2.0) * (1.0 - GaussPointNaturalCoordinates[0]) *
                               StrainSamplingPointsDU[i]["D"][2] +
                               (1.0 / 2.0) * (1.0 + GaussPointNaturalCoordinates[0]) *
                               StrainSamplingPointsDU[i]["B"][2];

                double[,] nodalBDUTransformed = MatrixOperations.MatrixProduct(transformationMatrix, nodalBdU);
                BANSmatrixDU[0, i * 3] = nodalBDUTransformed[0, 0];
                BANSmatrixDU[1, i * 3] = nodalBDUTransformed[1, 0];
                BANSmatrixDU[2, i * 3] = nodalBDUTransformed[2, 0];
                BANSmatrixDU[3, i * 3] = nodalBDUTransformed[3, 0];
                BANSmatrixDU[4, i * 3] = nodalBDUTransformed[4, 0];
                BANSmatrixDU[5, i * 3] = nodalBDUTransformed[5, 0];
                BANSmatrixDU[0, i * 3 + 1] = nodalBDUTransformed[0, 1];
                BANSmatrixDU[1, i * 3 + 1] = nodalBDUTransformed[1, 1];
                BANSmatrixDU[2, i * 3 + 1] = nodalBDUTransformed[2, 1];
                BANSmatrixDU[3, i * 3 + 1] = nodalBDUTransformed[3, 1];
                BANSmatrixDU[4, i * 3 + 1] = nodalBDUTransformed[4, 1];
                BANSmatrixDU[5, i * 3 + 1] = nodalBDUTransformed[5, 1];
                BANSmatrixDU[0, i * 3 + 2] = nodalBDUTransformed[0, 2];
                BANSmatrixDU[1, i * 3 + 2] = nodalBDUTransformed[1, 2];
                BANSmatrixDU[2, i * 3 + 2] = nodalBDUTransformed[2, 2];
                BANSmatrixDU[3, i * 3 + 2] = nodalBDUTransformed[3, 2];
                BANSmatrixDU[4, i * 3 + 2] = nodalBDUTransformed[4, 2];
                BANSmatrixDU[5, i * 3 + 2] = nodalBDUTransformed[5, 2];
            }
            return new Tuple<double[,], double[,]>(BANSmatrix, BANSmatrixDU);
        }
        private Dictionary<int, Dictionary<string, double[]>> AssumedStrainSamplingPoints(double[] nodalXInitial)
        {
            Dictionary<int, Dictionary<string, double[]>> assumedStrainPoints = new Dictionary<int, 
                Dictionary<string, double[]>>();

            double[] A = new double[] { 0.0, -1.0, 0.0 };
            double[] B = new double[] { 1.0, 0.0, 0.0 };
            double[] C = new double[] { 0.0, 1.0, 0.0 };
            double[] D = new double[] { -1.0, 0.0, 0.0 };

            Dictionary<string, double[]> localdNA = CalculateShapeFunctionsLocalDerivatives(A);
            Dictionary<string, double[]> localdNB = CalculateShapeFunctionsLocalDerivatives(B);
            Dictionary<string, double[]> localdNC = CalculateShapeFunctionsLocalDerivatives(C);
            Dictionary<string, double[]> localdND = CalculateShapeFunctionsLocalDerivatives(D);
            var gA = CalculategVectors(nodalXInitial, localdNA);
            var gB = CalculategVectors(nodalXInitial, localdNB);
            var gC = CalculategVectors(nodalXInitial, localdNC);
            var gD = CalculategVectors(nodalXInitial, localdND);
            for(int i = 0; i < 8; i++)
            {
                double[] e13A = VectorOperations.VectorVectorAddition(
                    VectorOperations.VectorScalarProductNew(gA[3], localdNA["ksi"][i]),
                    VectorOperations.VectorScalarProductNew(gA[1], localdNA["mhi"][i]));

                double[] e13C = VectorOperations.VectorVectorAddition(
                    VectorOperations.VectorScalarProductNew(gC[3], localdNC["ksi"][i]),
                    VectorOperations.VectorScalarProductNew(gC[1], localdNC["mhi"][i]));

                double[] e23B = VectorOperations.VectorVectorAddition(
                    VectorOperations.VectorScalarProductNew(gB[3], localdNB["ihta"][i]),
                    VectorOperations.VectorScalarProductNew(gB[2], localdNB["mhi"][i]));

                double[] e23D = VectorOperations.VectorVectorAddition(
                    VectorOperations.VectorScalarProductNew(gD[3], localdND["ihta"][i]),
                    VectorOperations.VectorScalarProductNew(gD[2], localdND["mhi"][i]));

                Dictionary<string, double[]> strainPoints = new Dictionary<string, double[]>();
                strainPoints.Add("A", e13A);
                strainPoints.Add("B", e23B);
                strainPoints.Add("C", e13C);
                strainPoints.Add("D", e23D);
                assumedStrainPoints.Add(i, strainPoints);
            }
            return assumedStrainPoints;
        }
        private Dictionary<int, Dictionary<string, double[]>> AssumedStrainSamplingPointsDU()
        {
            Dictionary<int, Dictionary<string, double[]>> assumedStrainPoints = new Dictionary<int,
                Dictionary<string, double[]>>();

            double[] A = new double[] { 0.0, -1.0, 0.0 };
            double[] B = new double[] { 1.0, 0.0, 0.0 };
            double[] C = new double[] { 0.0, 1.0, 0.0 };
            double[] D = new double[] { -1.0, 0.0, 0.0 };

            Dictionary<string, double[]> localdNA = CalculateShapeFunctionsLocalDerivatives(A);
            Dictionary<string, double[]> localdNB = CalculateShapeFunctionsLocalDerivatives(B);
            Dictionary<string, double[]> localdNC = CalculateShapeFunctionsLocalDerivatives(C);
            Dictionary<string, double[]> localdND = CalculateShapeFunctionsLocalDerivatives(D);
            var dgA = CalculatedgdUVectors(localdNA);
            var dgB = CalculatedgdUVectors(localdNB);
            var dgC = CalculatedgdUVectors(localdNC);
            var dgD = CalculatedgdUVectors(localdND);
            for (int i = 0; i < 8; i++)
            {
                double[] e13A = VectorOperations.VectorVectorAddition(
                    VectorOperations.VectorScalarProductNew(dgA[3], localdNA["ksi"][i]),
                    VectorOperations.VectorScalarProductNew(dgA[1], localdNA["mhi"][i]));

                double[] e13C = VectorOperations.VectorVectorAddition(
                    VectorOperations.VectorScalarProductNew(dgC[3], localdNC["ksi"][i]),
                    VectorOperations.VectorScalarProductNew(dgC[1], localdNC["mhi"][i]));

                double[] e23B = VectorOperations.VectorVectorAddition(
                    VectorOperations.VectorScalarProductNew(dgB[3], localdNB["ihta"][i]),
                    VectorOperations.VectorScalarProductNew(dgB[2], localdNB["mhi"][i]));

                double[] e23D = VectorOperations.VectorVectorAddition(
                    VectorOperations.VectorScalarProductNew(dgD[3], localdND["ihta"][i]),
                    VectorOperations.VectorScalarProductNew(dgD[2], localdND["mhi"][i]));

                Dictionary<string, double[]> strainPoints = new Dictionary<string, double[]>();
                strainPoints.Add("A", e13A);
                strainPoints.Add("B", e23B);
                strainPoints.Add("C", e13C);
                strainPoints.Add("D", e23D);
                assumedStrainPoints.Add(i, strainPoints);
            }
            return assumedStrainPoints;
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
        private double[] CalculateStressVector(double[,] E, double[] strain)
        {
            double[] stressVector = VectorOperations.MatrixVectorProduct(E, strain);
            return stressVector;
        }
        private Tuple<double[], double[]> GaussPoints(int i, int j, int k)
        {
            double[] gaussPoints = new double[] { -1.0 / Math.Sqrt(3), 1.0 / Math.Sqrt(3) };
            double[] gaussWeights = new double[] { 1.0, 1.0 };

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
        //
        private double[,] CreateEnhancedStrainsInterpolationMatrix(double[] ksi)
        {
            double[,] M = new double[6, 4];
            M[2, 0] = ksi[2];
            M[2, 1] = ksi[0] * ksi[2];
            M[2, 2] = ksi[1] * ksi[2];
            M[2, 3] = ksi[0] * ksi[1] * ksi[2];
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
            double[,] TransposedT = new double[6,6];
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
        private double[,] CalculateKGeoMatrix(Dictionary<string, double[]> localdN, double[] modifiedStressVector)
        {
            double[,] kGeoMatrix = new double[24, 24];
            double[,] dNMatrix = CalculatedNMatrix(localdN);
            double[,] stressMatrix = CalculateStressMatrix(modifiedStressVector);
            for(int i = 0; i < 8; i++)
            {
                for(int j = 0; j < 8; j++)
                {
                    double s = 0.0;
                    for(int k = 0; k < 3; k++)
                    {
                        for(int l = 0; l < 3; l++)
                        {
                            s += dNMatrix[i, k] * stressMatrix[k, l] * dNMatrix[j, l];
                        }
                    }
                    kGeoMatrix[3 * i, 3 *  j] = s;
                    kGeoMatrix[3 * i + 1, 3 * j + 1] = s;
                    kGeoMatrix[3 * i + 2, 3 * j + 2] = s;

                }
            }
            return kGeoMatrix;
        }
        private double[,] CalculatedNMatrix(Dictionary<string, double[]> localdN)
        {
            double[,] dNMatrix = new double[8, 3];
            for(int i = 0; i < 8; i++)
            {
                dNMatrix[i, 0] = localdN["ksi"][i];
                dNMatrix[i, 1] = localdN["ihta"][i];
                dNMatrix[i, 2] = localdN["mhi"][i];
            }
            return dNMatrix;
        }
        private double[,] CalculateStressMatrix(double[] StressVector)
        {
            double[,] s = new double[3, 3];

            s[0, 0] = StressVector[0];
            s[0, 1] = StressVector[3];
            s[0, 2] = StressVector[4];

            s[1, 0] = StressVector[3];
            s[1, 1] = StressVector[1];
            s[1, 2] = StressVector[5];

            s[2, 0] = StressVector[4];
            s[2, 1] = StressVector[5];
            s[2, 2] = StressVector[2];

            return s;
        }
        public double[,] CreateGlobalStiffnessMatrix()
        {
            double[,] K1 = new double[24, 24];
            double[,] K2 = new double[24, 24];
            double[,] LeMatrix = Properties.EASLMatrix;
            double[,] DeMatrixInv = Properties.EASDMatrix;
            double[,] E = CalculateStressStrainMatrix(Properties.YoungMod, Properties.PoissonRatio);
            double[] nodalXInitial = InitialNodalCoordinates();
            //double[] nodalX = UpdateNodalCoordinates(DisplacementVector);
            double[] center = new double[3];
            double[,] J0 = CalculateJacobian(nodalXInitial, CalculateShapeFunctionsLocalDerivatives(center));
            double detJ0 = CalculateJacobianDet(J0);
            double[,] transformationMat0 = MatrixOperations.BlockMatrixInversion6X6(TransformationMatrixTransposed(J0));
            Dictionary<int, Dictionary<string, double[]>> aStrainPoints = AssumedStrainSamplingPoints(nodalXInitial);
            Dictionary<int, Dictionary<string, double[]>> aStrainPointsDU = AssumedStrainSamplingPointsDU();

            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    for (int k = 0; k < 2; k++)
                    {
                        double[] gP = GaussPoints(i, j, k).Item1;
                        double[] gW = GaussPoints(i, j, k).Item2;
                        Dictionary<string, double[]> localdN = CalculateShapeFunctionsLocalDerivatives(gP);
                        double[,] J = CalculateJacobian(nodalXInitial, localdN);
                        double[,] transformationMat =MatrixOperations.BlockMatrixInversion6X6(TransformationMatrixTransposed(J));
                        //double[,] invJ = CalculateInverseJacobian(J).Item1;
                        double detJ = CalculateInverseJacobian(J).Item2;
                        Dictionary<int, double[]> gVectors = CalculategVectors(nodalXInitial, localdN);
                        Dictionary<int, double[]> gVectorsdU = CalculatedgdUVectors(localdN);

                        var BMatrices = CalculateBANSMatrix(gP, localdN,transformationMat, gVectors, gVectorsdU,
                            aStrainPoints, aStrainPointsDU);
                        double[,] B = BMatrices.Item1;
                        //double[,] dBdU = BMatrices.Item2;
                        double[,] Gamma = CalculateEnhancedStrainMatrixGamma(transformationMat0, CreateEnhancedStrainsInterpolationMatrix(gP),
                            detJ0, detJ);
                        double[] alphaVector = Properties.EASVector;
                        double[] EnhStrainVector = VectorOperations.MatrixVectorProduct(Gamma, alphaVector);
                        double[] modifiedStrainVector = VectorOperations.VectorVectorAddition(CalculateStrainsVector(B), EnhStrainVector);
                        double[] modifiedStressVector = CalculateStressVector(E, modifiedStrainVector);
                        double[,] kGeoMatrix = CalculateKGeoMatrix(localdN, modifiedStressVector);
                        K1 = MatrixOperations.MatrixAddition(K1, MatrixOperations.ScalarMatrixProductNew(detJ * gW[0] * gW[1] * gW[2],
                            MatrixOperations.MatrixProduct(MatrixOperations.Transpose(B), MatrixOperations.MatrixProduct(E, B))));
                        K2 = MatrixOperations.MatrixAddition(K2, MatrixOperations.ScalarMatrixProductNew(detJ * gW[0] * gW[1] * gW[2], kGeoMatrix));
                        //LeMatrix = MatrixOperations.MatrixAddition(LeMatrix, MatrixOperations.ScalarMatrixProductNew(detJ * gW[0] * gW[1] * gW[2],
                        //    MatrixOperations.MatrixProduct(
                        //    MatrixOperations.Transpose(Gamma), MatrixOperations.MatrixProduct(E, B))));
                        //DeMatrix = MatrixOperations.MatrixAddition(DeMatrix, MatrixOperations.ScalarMatrixProductNew(detJ * gW[0] * gW[1] * gW[2],
                        //    MatrixOperations.MatrixProduct(
                        //    MatrixOperations.Transpose(Gamma), MatrixOperations.MatrixProduct(E, Gamma))));
                    }
                }
            }
            //double[,]  = MatrixOperations.BlockMatrixInversion4X4(DeMatrix);
            double[,] K = MatrixOperations.MatrixSubtraction(
                MatrixOperations.MatrixAddition(K1, K2),
                MatrixOperations.MatrixProduct(
                MatrixOperations.MatrixProduct(MatrixOperations.Transpose(LeMatrix), DeMatrixInv), LeMatrix));
            return K;
        }

        public double[,] CreateMassMatrix()
        {
            double[,] M = new double[24, 24];
            double[] nodalXUpdated = UpdateNodalCoordinates(DisplacementVector);
            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    for (int k = 0; k < 2; k++)
                    {
                        double[] lP = LobattoPoints(i, j, k).Item1;
                        double[] lW = LobattoPoints(i, j, k).Item2;
                        Dictionary<string, double[]> localdN = CalculateShapeFunctionsLocalDerivatives(lP);
                        double[,] J = CalculateJacobian(nodalXUpdated, localdN);
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

        public double[] CreateInternalGlobalForcesVector()
        {
            double[] VectorRe = new double[24];
            double[] vectorPe = Properties.EASFEnhancedVector;
            double[,] matrixLe = Properties.EASLMatrix;
            double[,] matrixDeInv = Properties.EASDMatrix;

            double[] nodalXInitial = InitialNodalCoordinates();
            //double[] nodalX = UpdateNodalCoordinates(DisplacementVector);
            double[,] E = CalculateStressStrainMatrix(Properties.YoungMod, Properties.PoissonRatio);
            Dictionary<int, Dictionary<string, double[]>> aStrainPoints = AssumedStrainSamplingPoints(nodalXInitial);
            Dictionary<int, Dictionary<string, double[]>> aStrainPointsDU = AssumedStrainSamplingPointsDU();
            double[] center = new double[3];
            double[,] J0 = CalculateJacobian(nodalXInitial, CalculateShapeFunctionsLocalDerivatives(center));
            double detJ0 = CalculateJacobianDet(J0);
            double[,] transformationMat0 = MatrixOperations.BlockMatrixInversion6X6(TransformationMatrixTransposed(J0));
            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    for (int k = 0; k < 2; k++)
                    {
                        double[] gP = GaussPoints(i, j, k).Item1;
                        double[] gW = GaussPoints(i, j, k).Item2;
                        Dictionary<string, double[]> localdN = CalculateShapeFunctionsLocalDerivatives(gP);
                        Dictionary<int, double[]> gVectors = CalculategVectors(nodalXInitial, localdN);
                        Dictionary<int, double[]> gVectorsdU = CalculatedgdUVectors(localdN);
                        double[,] J = CalculateJacobian(nodalXInitial, localdN);
                        double[,] transformationMat = MatrixOperations.BlockMatrixInversion6X6(TransformationMatrixTransposed(J));
                        double[,] invJ = CalculateInverseJacobian(J).Item1;
                        double detJ = CalculateInverseJacobian(J).Item2;
                        Dictionary<int, double[]> globaldN = CalculateShapeFunctionsGlobalDerivatives(localdN, invJ);
                        var BMatrices = CalculateBANSMatrix(gP, localdN, transformationMat, gVectors, gVectorsdU,
                            aStrainPoints, aStrainPointsDU);
                        double[,] B = BMatrices.Item1;
                        double[,] Gamma = CalculateEnhancedStrainMatrixGamma(transformationMat0, CreateEnhancedStrainsInterpolationMatrix(gP),
                                                                             detJ0, detJ);
                        double[] alphaVector = Properties.EASVector;
                        double[] EnhStrainVector = VectorOperations.MatrixVectorProduct(Gamma, alphaVector);
                        double[] modifiedStrainVector = VectorOperations.VectorVectorAddition(CalculateStrainsVector(B), EnhStrainVector);
                        double[] stressVectorModified = CalculateStressVector(E, modifiedStrainVector);
                        VectorRe = VectorOperations.VectorVectorAddition(VectorRe, VectorOperations.VectorScalarProductNew(
                            VectorOperations.MatrixVectorProduct(MatrixOperations.Transpose(B), stressVectorModified), detJ * gW[0] * gW[1] * gW[2]));
                    }

                }
            }
            double[] fInt = VectorOperations.VectorVectorSubtraction(VectorRe,
                VectorOperations.MatrixVectorProduct(MatrixOperations.MatrixProduct(MatrixOperations.Transpose(matrixLe), matrixDeInv), vectorPe));
            return fInt;
        }
    }
}