using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GFEC
{
    class Shell2DQuadratic4 : IElement
    {
        public Dictionary<int, INode> Nodes { get; }
        public IElementProperties Properties { get; set; }
        public Dictionary<int, bool[]> ElementFreedomSignature { get; } = new Dictionary<int, bool[]>();
        public List<int> ElementFreedomList { get; set; }
        public double[] DisplacementVector { get; set; }
        public double[] AccelerationVector { get; set; }
        public double poisson { get; set; }
        public List<double> DU = new List<double>();
        //private double thickness = 1.0; //To be included in Element Properties
        //private double density = 1.0; //To be included in Element Properties
        public void InitializeTangentialProperties()
        {
            for(int i = 1; i <= 4; i++)
            {
                DU.Add(0);
                DU.Add(0);
                DU.Add(0);
                DU.Add(0);
                DU.Add(0);
                DU.Add(0);
            }
            Properties.DU = DU;
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
            for (int i = 0; i < 4; i++)
            {
                Properties.DU[i * 6] = deltaU[i * 6];
                Properties.DU[i * 6 + 1] = deltaU[i * 6 + 1];
                Properties.DU[i * 6 + 2] = deltaU[i * 6 + 2];
                Properties.DU[i * 6 + 3] = deltaU[i * 6 + 3];
                Properties.DU[i * 6 + 4] = deltaU[i * 6 + 4];
                Properties.DU[i * 6 + 5] = deltaU[i * 6 + 5];
            }
        }
        public Shell2DQuadratic4(IElementProperties properties, Dictionary<int, INode> nodes)
        {
            Properties = properties;
            this.Nodes = nodes;
            ElementFreedomSignature[1] = new bool[] { true, true, true, true, true, true };
            ElementFreedomSignature[2] = new bool[] { true, true, true, true, true, true };
            ElementFreedomSignature[3] = new bool[] { true, true, true, true, true, true };
            ElementFreedomSignature[4] = new bool[] { true, true, true, true, true, true };
            DisplacementVector = new double[24];
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
        public double ClosestPointProjection()
        {
            throw new Exception("Needs to be removed. Has beeb used only for testing purposes");
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
            throw new Exception("Needs to be removed. Has beeb used only for testing purposes");

        }
        public List<double[]> GetStressFromElementsNodes()
        {
            throw new Exception("Needs to be removed. Has beeb used only for testing purposes");

        }
        public List<double[]> GetStrainFromElementsNodes()
        {
            throw new Exception("Needs to be removed. Has beeb used only for testing purposes");

        }
        public List<double[]> GetStressFromElements(List<double[]> parametricCoordinatesVector)
        {
            throw new Exception("Needs to be removed. Has beeb used only for testing purposes");

        }
        public List<double[]> GetphysicalCoordinatesFromElements(List<double[]> parametricCoordinatesVector)
        {
            throw new Exception("Needs to be removed. Has beeb used only for testing purposes");

        }
        private double[] UpdateNodalCoordinates(double[] displacementVector)
        {
            double[] updatedCoor = new double[24];
            for (int i = 1; i <= 4; i++)
            {
                updatedCoor[6 * i - 6] = Nodes[i].XCoordinate + displacementVector[6 * i - 6];
                updatedCoor[6 * i - 5] = Nodes[i].YCoordinate + displacementVector[6 * i - 5];
                updatedCoor[6 * i - 4] = Nodes[i].ZCoordinate + displacementVector[6 * i - 4];
                updatedCoor[6 * i - 3] = Nodes[i].RX + displacementVector[6 * i - 3];
                updatedCoor[6 * i - 2] = Nodes[i].RY + displacementVector[6 * i - 2];
                updatedCoor[6 * i - 1] = Nodes[i].RZ + displacementVector[6 * i - 1];
            }
            return updatedCoor;
        }

        public Dictionary<int, INode> NodesAtFinalState()
        {
            Dictionary<int, INode> finalNodes = new Dictionary<int, INode>();
            finalNodes[1] = new Node(Nodes[1].XCoordinate + DisplacementVector[0], Nodes[1].YCoordinate + DisplacementVector[1],
                Nodes[1].ZCoordinate + DisplacementVector[2], Nodes[1].RX + DisplacementVector[3], Nodes[1].RY + DisplacementVector[4],
                Nodes[1].RZ + DisplacementVector[5]);

            finalNodes[2] = new Node(Nodes[2].XCoordinate + DisplacementVector[6], Nodes[2].YCoordinate + DisplacementVector[7],
                Nodes[2].ZCoordinate + DisplacementVector[8], Nodes[2].RX + DisplacementVector[9], Nodes[2].RY + DisplacementVector[10],
                Nodes[2].RZ + DisplacementVector[11]);

            finalNodes[3] = new Node(Nodes[3].XCoordinate + DisplacementVector[12], Nodes[3].YCoordinate + DisplacementVector[13],
                Nodes[3].ZCoordinate + DisplacementVector[14], Nodes[3].RX + DisplacementVector[15], Nodes[3].RY + DisplacementVector[16],
                Nodes[3].RZ + DisplacementVector[17]);

            finalNodes[4] = new Node(Nodes[4].XCoordinate + DisplacementVector[18], Nodes[4].YCoordinate + DisplacementVector[19],
                Nodes[4].ZCoordinate + DisplacementVector[20], Nodes[4].RX + DisplacementVector[21], Nodes[4].RY + DisplacementVector[22],
                Nodes[4].RZ + DisplacementVector[23]);
            return finalNodes;
        }

        private Dictionary<int, double> CalculateShapeFunctions(double ksi, double ihta)
        {
            Dictionary<int, double> shapeFunctions = new Dictionary<int, double>();
            double N1 = 1.0 / 4.0 * (1 - ksi) * (1 - ihta); shapeFunctions.Add(1, N1);
            double N2 = 1.0 / 4.0 * (1 + ksi) * (1 - ihta); shapeFunctions.Add(2, N2);
            double N3 = 1.0 / 4.0 * (1 + ksi) * (1 + ihta); shapeFunctions.Add(3, N3);
            double N4 = 1.0 / 4.0 * (1 - ksi) * (1 + ihta); shapeFunctions.Add(4, N4);

            return shapeFunctions;
        }

        //private double[,] CalculateShapeFunctionMatrix(double ksi, double ihta)
        //{
        //    Dictionary<int, double> shapeFunctions = CalculateShapeFunctions(ksi, ihta);
        //    double[,] N = new double[,]
        //    {
        //        {shapeFunctions[1], 0, shapeFunctions[2], 0, shapeFunctions[3], 0, shapeFunctions[4], 0 },
        //        {0, shapeFunctions[1], 0, shapeFunctions[2], 0, shapeFunctions[3], 0, shapeFunctions[4] }
        //    };
        //    return N;
        //}

        //private Dictionary<string, double[]> CalculateShapeFunctionsLocalDerivatives(double[] naturalCoordinates)
        //{
        //    double ksi = naturalCoordinates[0];
        //    double ihta = naturalCoordinates[1];

        //    double[] dN_ksi = new double[]
        //    {
        //        (-1.0/4.0*(1-ihta)),
        //        (1.0/4.0*(1-ihta)),
        //        (1.0/4.0*(1+ihta)),
        //        (-1.0/4.0*(1+ihta)),
        //    };

        //    double[] dN_ihta = new double[]
        //    {
        //        (-1.0/4.0*(1-ksi)),
        //        (-1.0/4.0*(1+ksi)),
        //        (1.0/4.0*(1+ksi)),
        //        (1.0/4.0*(1-ksi)),
        //    };

        //    Dictionary<string, double[]> dN = new Dictionary<string, double[]>();
        //    dN.Add("ksi", dN_ksi);
        //    dN.Add("ihta", dN_ihta);
        //    return dN;
        //}

        //private double[,] CalculateJacobian(Dictionary<string, double[]> dN)
        //{
        //    double[,] jacobianMatrix = new double[2, 2];

        //    double[] xUpdated = UpdateNodalCoordinates(DisplacementVector);

        //    int k = 0;
        //    for (int i = 0; i < 4; i++)
        //    {
        //        jacobianMatrix[0, 0] = jacobianMatrix[0, 0] + xUpdated[k] * dN["ksi"][i];
        //        k = k + 2;
        //    }
        //    k = 1;
        //    for (int i = 0; i < 4; i++)
        //    {
        //        jacobianMatrix[0, 1] = jacobianMatrix[0, 1] + xUpdated[k] * dN["ksi"][i];
        //        k = k + 2;
        //    }

        //    k = 0;
        //    for (int i = 0; i < 4; i++)
        //    {
        //        jacobianMatrix[1, 0] = jacobianMatrix[1, 0] + xUpdated[k] * dN["ihta"][i];
        //        k = k + 2;
        //    }
        //    k = 1;
        //    for (int i = 0; i < 4; i++)
        //    {
        //        jacobianMatrix[1, 1] = jacobianMatrix[1, 1] + xUpdated[k] * dN["ihta"][i];
        //        k = k + 2;
        //    }

        //    return jacobianMatrix;
        //}

        //private Tuple<double[,], double> CalculateInverseJacobian(double[,] jacobianMatrix)
        //{
        //    double[,] jacobianInverseMatrix = new double[2, 2];

        //    double detj = jacobianMatrix[0, 0] * jacobianMatrix[1, 1] - jacobianMatrix[0, 1] * jacobianMatrix[1, 0];

        //    jacobianInverseMatrix[0, 0] = jacobianMatrix[1, 1] / detj;
        //    jacobianInverseMatrix[0, 1] = -jacobianMatrix[0, 1] / detj;
        //    jacobianInverseMatrix[1, 0] = -jacobianMatrix[1, 0] / detj;
        //    jacobianInverseMatrix[1, 1] = jacobianMatrix[0, 0] / detj;

        //    return new Tuple<double[,], double>(jacobianInverseMatrix, detj);
        //}

        //private Dictionary<int, double[]> CalculateShapeFunctionsGlobalDerivatives(Dictionary<string, double[]> dN, double[,] Jinv)
        //{
        //    Dictionary<int, double[]> dNg = new Dictionary<int, double[]>();

        //    for (int i = 0; i < 4; i++)
        //    {
        //        double[] dNlocal = new double[] { dN["ksi"][i], dN["ihta"][i] };
        //        double[] dNglobal = VectorOperations.MatrixVectorProduct(Jinv, dNlocal);
        //        dNg.Add(i, dNglobal);
        //    }
        //    return dNg;
        //}

        //private double[] CalculateStrainsVector(double[,] Bmatrix)
        //{
        //    double[] strains = VectorOperations.MatrixVectorProduct(Bmatrix, DisplacementVector);
        //    return strains;
        //}

        //private double[,] CalculateBMatrix(Dictionary<int, double[]> dNglobal)
        //{
        //    double[,] Bmatrix = new double[3, 8];

        //    for (int i = 0; i < 4; i++)
        //    {
        //        Bmatrix[0, i * 2] = dNglobal[i][0];
        //        Bmatrix[1, i * 2 + 1] = dNglobal[i][1];
        //        Bmatrix[2, i * 2] = dNglobal[i][1];
        //        Bmatrix[2, i * 2 + 1] = dNglobal[i][0];
        //    }
        //    return Bmatrix;
        //}

        private double[,] CalculateStressStrainMatrixPlaneStressPart(double E, double v)
        {
            double[,] Ematrix = new double[3, 3];
            //v = 0.30;
            //v = Properties.PoissonRatio;
            double Ehat = E / ((1.0 - Math.Pow(v, 2)));

            Ematrix[0, 0] = Ehat;
            Ematrix[0, 1] = Ehat * v;
            Ematrix[1, 0] = Ehat * v;
            Ematrix[1, 1] = Ehat;
            Ematrix[2, 2] = Ehat * (1.0 / 2.0) * (1.0 - v);

            return Ematrix;
        }

        private double[] CalculateStressVector(double[,] E, double[] strain)
        {
            double[] stressVector = VectorOperations.MatrixVectorProduct(E, strain);
            return stressVector;
        }

        //private Tuple<double[], double[]> GaussPoints(int i, int j)
        //{
        //    double[] gaussPoints = new double[] { -1.0 / Math.Sqrt(3), 1.0 / Math.Sqrt(3) };
        //    double[] gaussWeights = new double[] { 1.0, 1.0 };

        //    double[] vectorWithPoints = new double[] { gaussPoints[i], gaussPoints[j] };
        //    double[] vectorWithWeights = new double[] { gaussWeights[i], gaussWeights[j] };
        //    return new Tuple<double[], double[]>(vectorWithPoints, vectorWithWeights);
        //}

        private double[,] CreateLocalStiffnessMatrixPlaneStressPart(double[] xUpdated)
        {
            double[,] K = new double[8, 8];
            //double[] xUpdated = UpdateNodalCoordinates(DisplacementVector);
            double[] positionVector1 = new double[] { xUpdated[0], xUpdated[1], xUpdated[2] };
            double[] positionVector2 = new double[] { xUpdated[6], xUpdated[7], xUpdated[8] };
            //double[] positionVector3 = new double[] { xUpdated[12], xUpdated[13], xUpdated[14] };
            double[] positionVector4 = new double[] { xUpdated[18], xUpdated[19], xUpdated[20] };

            double alpha = VectorOperations.VectorNorm2(VectorOperations.VectorVectorSubtraction(positionVector2, positionVector1)) / 2.0;
            double beta = VectorOperations.VectorNorm2(VectorOperations.VectorVectorSubtraction(positionVector4, positionVector1)) / 2.0;

            double r = alpha / beta;
            double rInv = beta / alpha;
            double p = (1.0 - Properties.PoissonRatio) / 2.0;
            double mi = 3.0 * (1.0 + Properties.PoissonRatio) / 2.0;
            double lamda = 3.0 * (1.0 - 3.0 * Properties.PoissonRatio) / 2.0;

            K[0, 0] = 4.0 * rInv + 4.0 * p * r;
            K[1, 0] = mi;
            K[2, 0] = -4.0 * rInv + 2.0 * p * r;
            K[3, 0] = -lamda;
            K[4, 0] = -2.0 * rInv - 2.0 * p * r;
            K[5, 0] = -mi;
            K[6, 0] = 2.0 * rInv - 4.0 * p * r;
            K[7, 0] = lamda;

            K[0, 1] = K[1, 0];
            K[0, 2] = K[2, 0];
            K[0, 3] = K[3, 0];
            K[0, 4] = K[4, 0];
            K[0, 5] = K[5, 0];
            K[0, 6] = K[6, 0];
            K[0, 7] = K[7, 0];
            //
            K[1, 1] = 4.0 * r + 4.0 * p * rInv;
            K[2, 1] = lamda;
            K[3, 1] = 2.0 * r - 4.0 * p * rInv;
            K[4, 1] = -mi;
            K[5, 1] = -2.0 * r - 2.0 * p * rInv;
            K[6, 1] = -lamda;
            K[7, 1] = -4.0 * r + 2.0 * p * rInv;

            K[1, 2] = K[2, 1];
            K[1, 3] = K[3, 1];
            K[1, 4] = K[4, 1];
            K[1, 5] = K[5, 1];
            K[1, 6] = K[6, 1];
            K[1, 7] = K[7, 1];
            //
            K[2, 2] = 4.0 * rInv + 4.0 * p * r;
            K[3, 2] = -mi;
            K[4, 2] = 2.0 * rInv - 4.0 * p * r;
            K[5, 2] = -lamda;
            K[6, 2] = -2.0 * rInv - 2.0 * p * r;
            K[7, 2] = mi;

            K[2, 3] = K[3, 2];
            K[2, 4] = K[4, 2];
            K[2, 5] = K[5, 2];
            K[2, 6] = K[6, 2];
            K[2, 7] = K[7, 2];
            //
            K[3, 3] = 4.0 * r + 4.0 * p * rInv;
            K[4, 3] = lamda;
            K[5, 3] = -4.0 * r + 2.0 * p * rInv;
            K[6, 3] = mi;
            K[7, 3] = -2.0 * r - 2.0 * p * rInv;

            K[3, 4] = K[4, 3];
            K[3, 5] = K[5, 3];
            K[3, 6] = K[6, 3];
            K[3, 7] = K[7, 3];
            //
            K[4, 4] = 4.0 * rInv + 4.0 * p * r;
            K[5, 4] = mi;
            K[6, 4] = -4.0 * rInv + 2.0 * p * r;
            K[7, 4] = -lamda;

            K[4, 5] = K[5, 4];
            K[4, 6] = K[6, 4];
            K[4, 7] = K[7, 4];
            //
            K[5, 5] = 4.0 * r + 4.0 * p * rInv;
            K[6, 5] = lamda;
            K[7, 5] = 2.0 * r - 4.0 * p * rInv;

            K[5, 6] = K[6, 5];
            K[5, 7] = K[7, 5];
            //
            K[6, 6] = 4.0 * rInv + 4.0 * p * r;
            K[7, 6] = -mi;

            K[6, 7] = K[7, 6];
            K[7, 7] = 4.0 * r + 4.0 * p * rInv;
            //
            double scalar = Properties.YoungMod * Properties.Thickness / (12.0 * (1.0 - Math.Pow(Properties.PoissonRatio, 2.0)));
            K = MatrixOperations.ScalarMatrixProduct(scalar, K);
            return K;
        }

        private double[,] CreateLocalMassMatrixPlaneStressPart(double[] xUpdated)
        {
            //double[] xUpdated = UpdateNodalCoordinates(DisplacementVector);
            double[] positionVector1 = new double[] { xUpdated[0], xUpdated[1], xUpdated[2] };
            double[] positionVector2 = new double[] { xUpdated[6], xUpdated[7], xUpdated[8] };
            double[] positionVector3 = new double[] { xUpdated[12], xUpdated[13], xUpdated[14] };
            double[] positionVector4 = new double[] { xUpdated[18], xUpdated[19], xUpdated[20] };

            double alpha = VectorOperations.VectorNorm2(VectorOperations.VectorVectorSubtraction(positionVector2, positionVector1));
            double beta = VectorOperations.VectorNorm2(VectorOperations.VectorVectorSubtraction(positionVector4, positionVector1));
            double DLMM = Properties.Density * alpha * beta * Properties.Thickness / 4.0;
            double[,] M = MatrixOperations.CreateDiagonalMatrix(8, DLMM);
            return M;
            //for(int i = 0; i < 4; i++)
            //{
            //    M[i, i] = 4.0;
            //}
            //M[0, 1] = 2.0;
            //M[0, 2] = 1.0;
            //M[0, 3] = 2.0;
            //M[1, 0] = 2.0;
            //M[2, 0] = 1.0;
            //M[3, 0] = 2.0;

            //M[1, 2] = 2.0;
            //M[1, 3] = 1.0;
            //M[2, 1] = 2.0;
            //M[3, 1] = 1.0;

            //M[2, 3] = 2.0;
            //M[3, 2] = 2.0;
            //for(int i = 4; i < 8; i++)
            //{
            //    for(int j = 4; j < 8; j++)
            //    {
            //        M[i, j] = M[i - 4, j - 4];
            //    }
            //}
            //double scalar = Properties.Density * Properties.Thickness * alpha * beta / 36.0;
            //M = MatrixOperations.ScalarMatrixProduct(scalar, M);
            //return M;
        }

        //private double[,] CreateLocalStiffnessMatrixBendingPart(double[] xUpdated)
        //{
        //    double[] positionVector1 = new double[] { xUpdated[0], xUpdated[1], xUpdated[2] };
        //    double[] positionVector2 = new double[] { xUpdated[6], xUpdated[7], xUpdated[8] };
        //    double[] positionVector3 = new double[] { xUpdated[12], xUpdated[13], xUpdated[14] };
        //    double[] positionVector4 = new double[] { xUpdated[18], xUpdated[19], xUpdated[20] };

        //    double alpha = VectorOperations.VectorNorm2(VectorOperations.VectorVectorSubtraction(positionVector2, positionVector1)) / 2.0;
        //    double beta = VectorOperations.VectorNorm2(VectorOperations.VectorVectorSubtraction(positionVector4, positionVector1)) / 2.0;

        //    double[,] K = new double[12, 12];
        //    double[,] K1 = new double[12, 12];
        //    double[,] K2 = new double[12, 12];
        //    double[,] K3 = new double[12, 12];
        //    double[,] K4 = new double[12, 12];
        //    //K1
        //    K1[0, 0] = 6.0;
        //    K1[2, 0] = -6.0 * alpha;
        //    K1[3, 0] = -6.0;
        //    K1[5, 0] = -6.0 * alpha;
        //    K1[6, 0] = -3.0;
        //    K1[8, 0] = -3.0 * alpha;
        //    K1[9, 0] = 3.0;
        //    K1[11, 0] = -3.0 * alpha;
        //    //
        //    K1[0, 2] = K1[2, 0];
        //    K1[0, 3] = K1[3, 0];
        //    K1[0, 5] = K1[5, 0];
        //    K1[0, 6] = K1[6, 0];
        //    K1[0, 8] = K1[8, 0];
        //    K1[0, 9] = K1[9, 0];
        //    K1[0, 11] = K1[11, 0];
        //    //
        //    K1[2, 2] = 8.0 * alpha * alpha;
        //    K1[2, 3] = 6.0 * alpha;
        //    K1[2, 5] = 4.0 * alpha * alpha;
        //    K1[2, 6] = 3.0 * alpha;
        //    K1[2, 8] = 2.0 * alpha * alpha;
        //    K1[2, 9] = -3.0 * alpha;
        //    K1[2, 11] = 4.0 * alpha * alpha;
        //    //
        //    K1[3, 2] = K1[2, 3];
        //    K1[5, 2] = K1[2, 5];
        //    K1[6, 2] = K1[2, 6];
        //    K1[8, 2] = K1[2, 8];
        //    K1[9, 2] = K1[2, 9];
        //    K1[11, 2] = K1[2, 11];
        //    //
        //    K1[3, 3] = 6.0;
        //    K1[5, 3] = 6.0 * alpha;
        //    K1[6, 3] = 3.0;
        //    K1[8, 3] = 3.0 * alpha;
        //    K1[9, 3] = -3.0;
        //    K1[11, 3] = 3.0 * alpha;
        //    //
        //    K1[3, 5] = K1[5, 3];
        //    K1[3, 6] = K1[6, 3];
        //    K1[3, 8] = K1[8, 3];
        //    K1[3, 9] = K1[9, 3];
        //    K1[3, 11] = K1[11, 3];
        //    //
        //    K1[5, 5] = 8.0 * alpha * alpha;
        //    K1[6, 5] = 3.0 * alpha;
        //    K1[8, 5] = 4.0 * alpha * alpha;
        //    K1[9, 5] = -3.0 * alpha;
        //    K1[11, 5] = 2.0 * alpha * alpha;
        //    //
        //    K1[5, 6] = K1[6, 5];
        //    K1[5, 8] = K1[8, 5];
        //    K1[5, 9] = K1[9, 5];
        //    K1[5, 11] = K1[11, 5];
        //    //
        //    K1[6, 6] = 6.0;
        //    K1[8, 6] = 6.0 * alpha;
        //    K1[9, 6] = -6.0;
        //    K1[11, 6] = 6.0 * alpha;
        //    //
        //    K1[6, 8] = K1[8, 6];
        //    K1[6, 9] = K1[9, 6];
        //    K1[6, 11] = K1[11, 6];
        //    //
        //    K1[8, 8] = 8.0 * alpha * alpha;
        //    K1[9, 8] = -6.0 * alpha;
        //    K1[11, 8] = 4.0 * alpha * alpha;
        //    //
        //    K1[8, 9] = K1[9, 8];
        //    K1[8, 11] = K1[11, 8];
        //    //
        //    K1[9, 9] = 6.0;
        //    K1[11, 9] = -6.0 * alpha;
        //    //
        //    K1[9, 11] = K1[11, 9];
        //    K1[11, 11] = 8.0 * alpha * alpha;

        //    //K2
        //    K2[0, 0] = 6.0;
        //    K2[1, 0] = 6.0 * beta;
        //    K2[3, 0] = 3.0;
        //    K2[4, 0] = 3.0 * beta;
        //    K2[6, 0] = -3.0;
        //    K2[7, 0] = 3.0 * beta;
        //    K2[9, 0] = -6.0;
        //    K2[10, 0] = 6.0 * beta;
        //    //
        //    K2[0, 1] = K2[1, 0];
        //    K2[0, 3] = K2[3, 0];
        //    K2[0, 4] = K2[4, 0];
        //    K2[0, 6] = K2[6, 0];
        //    K2[0, 7] = K2[7, 0];
        //    K2[0, 9] = K2[9, 0];
        //    K2[0, 10] = K2[10, 0];
        //    //
        //    K2[1, 1] = 8.0 * beta * beta;
        //    K2[3, 1] = 3.0 * beta;
        //    K2[4, 1] = 4.0 * beta * beta;
        //    K2[6, 1] = -3.0 * beta;
        //    K2[7, 1] = 2.0 * beta * beta;
        //    K2[9, 1] = -6.0 * beta;
        //    K2[10, 1] = 4.0 * beta * beta;
        //    //
        //    K2[1, 3] = K2[3, 1];
        //    K2[1, 4] = K2[4, 1];
        //    K2[1, 6] = K2[6, 1];
        //    K2[1, 7] = K2[7, 1];
        //    K2[1, 9] = K2[9, 1];
        //    K2[1, 10] = K2[10, 1];
        //    //
        //    K2[3, 3] = 6.0;
        //    K2[4, 3] = 6.0 * beta;
        //    K2[6, 3] = -6.0;
        //    K2[7, 3] = 6.0 * beta;
        //    K2[9, 3] = -3.0;
        //    K2[10, 3] = 3.0 * beta;
        //    //
        //    K2[3, 4] = K2[4, 3];
        //    K2[3, 6] = K2[6, 3];
        //    K2[3, 7] = K2[7, 3];
        //    K2[3, 9] = K2[9, 3];
        //    K2[3, 10] = K2[10, 3];
        //    //
        //    K2[4, 4] = 8.0 * beta * beta;
        //    K2[6, 4] = -6.0 * beta;
        //    K2[7, 4] = 4.0 * beta * beta;
        //    K2[9, 4] = -3.0 * beta;
        //    K2[10, 4] = 2.0 * beta * beta;
        //    //
        //    K2[4, 6] = K2[6, 4];
        //    K2[4, 7] = K2[7, 4];
        //    K2[4, 9] = K2[9, 4];
        //    K2[4, 10] = K2[10, 4];
        //    //
        //    K2[6, 6] = 6.0;
        //    K2[7, 6] = -6.0 * beta;
        //    K2[9, 6] = 3.0;
        //    K2[10, 6] = -3.0 * beta;
        //    //
        //    K2[6, 7] = K2[7, 6];
        //    K2[6, 9] = K2[9, 6];
        //    K2[6, 10] = K2[10, 6];
        //    //
        //    K2[7, 7] = 8.0 * beta * beta;
        //    K2[9, 7] = -3.0 * beta;
        //    K2[10, 7] = 4.0 * beta * beta;
        //    //
        //    K2[7, 9] = K2[9, 7];
        //    K2[7, 10] = K2[10, 7];
        //    //
        //    K2[9, 9] = 6.0;
        //    K2[10, 9] = -6.0 * beta;
        //    //
        //    K2[9, 10] = K2[10, 9];
        //    K2[10, 10] = 8.0 * beta * beta;

        //    //K3
        //    K3[0, 0] = 1.0;
        //    K3[1, 0] = beta;
        //    K3[2, 0] = -alpha;
        //    K3[3, 0] = -1.0;
        //    K3[4, 0] = -beta;
        //    K3[6, 0] = 1.0;
        //    K3[9, 0] = -1.0;
        //    K3[11, 0] = alpha;
        //    //
        //    K3[0, 1] = K3[1, 0];
        //    K3[0, 2] = K3[2, 0];
        //    K3[0, 3] = K3[3, 0];
        //    K3[0, 4] = K3[4, 0];
        //    K3[0, 6] = K3[6, 0];
        //    K3[0, 9] = K3[9, 0];
        //    K3[0, 11] = K3[11, 0];
        //    //
        //    K3[2, 1] = -2.0 * alpha * beta;
        //    K3[3, 1] = -beta;
        //    //
        //    K3[1, 2] = K3[2, 1];
        //    K3[1, 3] = K3[3, 1];
        //    //
        //    K3[3, 3] = 1.0;
        //    K3[4, 3] = beta;
        //    K3[5, 3] = alpha;
        //    K3[6, 3] = -1.0;
        //    K3[8, 3] = -alpha;
        //    K3[9, 3] = 1.0;
        //    //
        //    K3[3, 4] = K3[4, 3];
        //    K3[3, 5] = K3[5, 3];
        //    K3[3, 6] = K3[6, 3];
        //    K3[3, 8] = K3[8, 3];
        //    K3[3, 9] = K3[9, 3];
        //    //
        //    K3[4, 5] = 2.0 * alpha * beta;
        //    K3[5, 4] = K3[4, 5];
        //    //
        //    K3[5, 6] = -alpha;
        //    K3[6, 5] = K3[5, 6];
        //    //
        //    K3[6, 6] = 1.0;
        //    K3[7, 6] = -beta;
        //    K3[8, 6] = alpha;
        //    K3[9, 6] = -1.0;
        //    K3[10, 6] = beta;

        //    K3[6, 7] = K3[7, 6];
        //    K3[6, 8] = K3[8, 6];
        //    K3[6, 9] = K3[9, 6];
        //    K3[6, 10] = K3[10, 6];
        //    //
        //    K3[8, 7] = -2.0 * alpha * beta;
        //    K3[9, 7] = beta;
        //    K3[7, 8] = K3[8, 7];
        //    K3[7, 9] = K3[9, 7];
        //    //
        //    K3[9, 9] = 1.0;
        //    K3[10, 9] = -beta;
        //    K3[11, 9] = -alpha;
        //    K3[9, 10] = K3[10, 9];
        //    K3[9, 11] = K3[11, 9];
        //    //
        //    K3[11, 10] = 2.0 * alpha * beta;
        //    K3[10, 11] = K3[11, 10];

        //    //K4
        //    K4[0, 0] = 21.0;
        //    K4[1, 0] = 3.0 * beta;
        //    K4[2, 0] = -3.0 * alpha;
        //    K4[3, 0] = -21.0;
        //    K4[4, 0] = -3.0 * beta;
        //    K4[5, 0] = -3.0 * alpha;
        //    K4[6, 0] = 21.0;
        //    K4[7, 0] = -3.0 * beta;
        //    K4[8, 0] = 3.0 * alpha;
        //    K4[9, 0] = -21.0;
        //    K4[10, 0] = 3.0 * beta;
        //    K4[11, 0] = 3.0 * alpha;
        //    //
        //    K4[0, 1] = K4[1, 0];
        //    K4[0, 2] = K4[2, 0];
        //    K4[0, 3] = K4[3, 0];
        //    K4[0, 4] = K4[4, 0];
        //    K4[0, 5] = K4[5, 0];
        //    K4[0, 6] = K4[6, 0];
        //    K4[0, 7] = K4[7, 0];
        //    K4[0, 8] = K4[8, 0];
        //    K4[0, 9] = K4[9, 0];
        //    K4[0, 10] = K4[10, 0];
        //    K4[0, 11] = K4[11, 0];
        //    //
        //    K4[1, 1] = 8.0 * beta * beta;
        //    K4[3, 1] = -3.0 * beta;
        //    K4[4, 1] = -8.0 * beta * beta;
        //    K4[6, 1] = 3.0 * beta;
        //    K4[7, 1] = 2.0 * beta * beta;
        //    K4[9, 1] = -3.0 * beta;
        //    K4[10, 1] = -2.0 * beta * beta;
        //    //
        //    K4[1, 3] = K4[3, 1];
        //    K4[1, 4] = K4[4, 1];
        //    K4[1, 6] = K4[6, 1];
        //    K4[1, 7] = K4[7, 1];
        //    K4[1, 9] = K4[9, 1];
        //    K4[1, 10] = K4[10, 1];
        //    //
        //    K4[2, 2] = 8.0 * alpha * alpha;
        //    K4[3, 2] = 3.0 * alpha;
        //    K4[5, 2] = -2.0 * alpha * alpha;
        //    K4[6, 2] = -3.0 * alpha;
        //    K4[8, 2] = 2.0 * alpha * alpha;
        //    K4[9, 2] = 3.0 * alpha;
        //    K4[11, 2] = -8.0 * alpha * alpha;
        //    //
        //    K4[2, 3] = K4[3, 2];
        //    K4[2, 5] = K4[5, 2];
        //    K4[2, 6] = K4[6, 2];
        //    K4[2, 8] = K4[8, 2];
        //    K4[2, 9] = K4[9, 2];
        //    K4[2, 11] = K4[11, 2];
        //    //
        //    K4[3, 3] = 21.0;
        //    K4[4, 3] = 3.0 * beta;
        //    K4[5, 3] = 3.0 * alpha;
        //    K4[6, 3] = -21.0;
        //    K4[7, 3] = 3.0 * beta;
        //    K4[8, 3] = -3.0 * alpha;
        //    K4[9, 3] = 21.0;
        //    K4[10, 3] = -3.0 * beta;
        //    K4[11, 3] = -3.0 * alpha;
        //    //
        //    K4[3, 4] = K4[4, 3];
        //    K4[3, 5] = K4[5, 3];
        //    K4[3, 6] = K4[6, 3];
        //    K4[3, 7] = K4[7, 3];
        //    K4[3, 8] = K4[8, 3];
        //    K4[3, 9] = K4[9, 3];
        //    K4[3, 10] = K4[10, 3];
        //    K4[3, 11] = K4[11, 3];
        //    //
        //    K4[4, 4] = 8.0 * beta * beta;
        //    K4[6, 4] = -3.0 * beta;
        //    K4[7, 4] = -2.0 * beta * beta;
        //    K4[9, 4] = 3.0 * beta;
        //    K4[10, 4] = 2.0 * beta * beta;
        //    //
        //    K4[4, 6] = K4[6, 4];
        //    K4[4, 7] = K4[7, 4];
        //    K4[4, 9] = K4[9, 4];
        //    K4[4, 10] = K4[10, 4];
        //    //
        //    K4[5, 5] = 8.0 * alpha * alpha;
        //    K4[6, 5] = -3.0 * alpha;
        //    K4[8, 5] = -8.0 * alpha * alpha;
        //    K4[9, 5] = 3.0 * alpha;
        //    K4[11, 5] = 2.0 * alpha * alpha;
        //    //
        //    K4[5, 6] = K4[6, 5];
        //    K4[5, 8] = K4[8, 5];
        //    K4[5, 9] = K4[9, 5];
        //    K4[5, 11] = K4[11, 5];
        //    //
        //    K4[6, 6] = 21.0;
        //    K4[7, 6] = -3.0 * beta;
        //    K4[8, 6] = 3.0 * alpha;
        //    K4[9, 6] = -21.0;
        //    K4[10, 6] = 3.0 * beta;
        //    K4[11, 6] = 3.0 * alpha;
        //    //
        //    K4[6, 7] = K4[7, 6];
        //    K4[6, 8] = K4[8, 6];
        //    K4[6, 9] = K4[9, 6];
        //    K4[6, 10] = K4[10, 6];
        //    K4[6, 11] = K4[11, 6];
        //    //
        //    K4[7, 7] = 8.0 * beta * beta;
        //    K4[9, 7] = 3.0 * beta;
        //    K4[10, 7] = -8.0 * beta * beta;
        //    //
        //    K4[7, 9] = K4[9, 7];
        //    K4[7, 10] = K4[10, 7];
        //    //
        //    K4[8, 8] = 8.0 * alpha * alpha;
        //    K4[9, 8] = -3.0 * alpha;
        //    K4[11, 8] = -2.0 * alpha * alpha;
        //    //
        //    K4[8, 9] = K4[9, 8];
        //    K4[8, 11] = K4[11, 8];
        //    //
        //    K4[9, 9] = 21.0;
        //    K4[10, 9] = -3.0 * beta;
        //    K4[11, 9] = -3.0 * alpha;
        //    //
        //    K4[9, 10] = K4[10, 9];
        //    K4[9, 11] = K4[11, 9];
        //    //
        //    K4[10, 10] = 8.0 * beta * beta;
        //    K4[11, 11] = 8.0 * alpha * alpha;
        //    //
        //    double scalar = -Properties.YoungMod * Math.Pow(Properties.Thickness, 3.0) / (12.0 * (1.0 - Math.Pow(Properties.PoissonRatio, 2.0)) * 60.0 * alpha * beta);
        //    K = MatrixOperations.MatrixAddition(K1, MatrixOperations.MatrixAddition(K2,
        //        MatrixOperations.MatrixAddition(K3, K4)));
        //    K = MatrixOperations.ScalarMatrixProduct(scalar, K);
        //    return K;
        //}
        //private double[,] CreateLocalStiffnessMatrixBendingPart(double[] xUpdated)
        //{
        //    double[] positionVector1 = new double[] { xUpdated[0], xUpdated[1], xUpdated[2] };
        //    double[] positionVector2 = new double[] { xUpdated[6], xUpdated[7], xUpdated[8] };
        //    double[] positionVector3 = new double[] { xUpdated[12], xUpdated[13], xUpdated[14] };
        //    double[] positionVector4 = new double[] { xUpdated[18], xUpdated[19], xUpdated[20] };

        //    double alpha = VectorOperations.VectorNorm2(VectorOperations.VectorVectorSubtraction(positionVector2, positionVector1));
        //    double beta = VectorOperations.VectorNorm2(VectorOperations.VectorVectorSubtraction(positionVector4, positionVector1));

        //    double[,] K = new double[12, 12];
        //    double[,] K1 = new double[12, 12];
        //    double[,] K2 = new double[12, 12];
        //    double[,] K3 = new double[12, 12];
        //    double[,] K4 = new double[12, 12];
        //    //K1
        //    K1[0, 0] = 60.0;
        //    K1[2, 0] = 30.0;
        //    K1[3, 0] = 30.0;
        //    K1[5, 0] = 15.0;
        //    K1[6, 0] = -60.0;
        //    K1[8, 0] = 30.0;
        //    K1[9, 0] = -30.0;
        //    K1[11, 0] = 15.0;
        //    //
        //    K1[0, 2] = K1[2, 0];
        //    K1[0, 3] = K1[3, 0];
        //    K1[0, 5] = K1[5, 0];
        //    K1[0, 6] = K1[6, 0];
        //    K1[0, 8] = K1[8, 0];
        //    K1[0, 9] = K1[9, 0];
        //    K1[0, 11] = K1[11, 0];
        //    //
        //    K1[2, 2] = 20.0;
        //    K1[2, 3] = 15.0;
        //    K1[2, 5] = 10.0;
        //    K1[2, 6] = -30.0;
        //    K1[2, 8] = 10.0;
        //    K1[2, 9] = -15.0;
        //    K1[2, 11] = 5.0;
        //    //
        //    K1[3, 2] = K1[2, 3];
        //    K1[5, 2] = K1[2, 5];
        //    K1[6, 2] = K1[2, 6];
        //    K1[8, 2] = K1[2, 8];
        //    K1[9, 2] = K1[2, 9];
        //    K1[11, 2] = K1[2, 11];
        //    //
        //    K1[3, 3] = 60.0;
        //    K1[5, 3] = 30.0;
        //    K1[6, 3] = -30.0;
        //    K1[8, 3] = 15.0;
        //    K1[9, 3] = -60.0;
        //    K1[11, 3] = 30.0;
        //    //
        //    K1[3, 5] = K1[5, 3];
        //    K1[3, 6] = K1[6, 3];
        //    K1[3, 8] = K1[8, 3];
        //    K1[3, 9] = K1[9, 3];
        //    K1[3, 11] = K1[11, 3];
        //    //
        //    K1[5, 5] = 20.0;
        //    K1[6, 5] = -15.0;
        //    K1[8, 5] = 5.0;
        //    K1[9, 5] = -30.0;
        //    K1[11, 5] = 10.0;
        //    //
        //    K1[5, 6] = K1[6, 5];
        //    K1[5, 8] = K1[8, 5];
        //    K1[5, 9] = K1[9, 5];
        //    K1[5, 11] = K1[11, 5];
        //    //
        //    K1[6, 6] = 60.0;
        //    K1[8, 6] = -30.0;
        //    K1[9, 6] = 30.0;
        //    K1[11, 6] = -15.0;
        //    //
        //    K1[6, 8] = K1[8, 6];
        //    K1[6, 9] = K1[9, 6];
        //    K1[6, 11] = K1[11, 6];
        //    //
        //    K1[8, 8] = 20.0;
        //    K1[9, 8] = -15.0;
        //    K1[11, 8] = 10.0;
        //    //
        //    K1[8, 9] = K1[9, 8];
        //    K1[8, 11] = K1[11, 8];
        //    //
        //    K1[9, 9] = 60.0;
        //    K1[11, 9] = -30.0;
        //    //
        //    K1[9, 11] = K1[11, 9];
        //    K1[11, 11] = 20.0;
        //    double scal1 = Math.Pow(beta / alpha, 2.0);
        //    K1 = MatrixOperations.ScalarMatrixProduct(scal1, K1);
        //    //K2
        //    K2[0, 0] = 60.0;
        //    K2[1, 0] = -30.0;
        //    K2[3, 0] = -60.0;
        //    K2[4, 0] = -30.0;
        //    K2[6, 0] = 30.0;
        //    K2[7, 0] = -15.0;
        //    K2[9, 0] = -30.0;
        //    K2[10, 0] = -15.0;
        //    //
        //    K2[0, 1] = K2[1, 0];
        //    K2[0, 3] = K2[3, 0];
        //    K2[0, 4] = K2[4, 0];
        //    K2[0, 6] = K2[6, 0];
        //    K2[0, 7] = K2[7, 0];
        //    K2[0, 9] = K2[9, 0];
        //    K2[0, 10] = K2[10, 0];
        //    //
        //    K2[1, 1] = 20.0;
        //    K2[3, 1] = 30.0;
        //    K2[4, 1] = 10.0;
        //    K2[6, 1] = -15.0;
        //    K2[7, 1] = 10.0;
        //    K2[9, 1] = 15.0;
        //    K2[10, 1] = 5.0;
        //    //
        //    K2[1, 3] = K2[3, 1];
        //    K2[1, 4] = K2[4, 1];
        //    K2[1, 6] = K2[6, 1];
        //    K2[1, 7] = K2[7, 1];
        //    K2[1, 9] = K2[9, 1];
        //    K2[1, 10] = K2[10, 1];
        //    //
        //    K2[3, 3] = 60.0;
        //    K2[4, 3] = 30.0;
        //    K2[6, 3] = -30.0;
        //    K2[7, 3] = 15.0;
        //    K2[9, 3] = 30.0;
        //    K2[10, 3] = 15.0;
        //    //
        //    K2[3, 4] = K2[4, 3];
        //    K2[3, 6] = K2[6, 3];
        //    K2[3, 7] = K2[7, 3];
        //    K2[3, 9] = K2[9, 3];
        //    K2[3, 10] = K2[10, 3];
        //    //
        //    K2[4, 4] = 20.0;
        //    K2[6, 4] = -15.0;
        //    K2[7, 4] = 5.0;
        //    K2[9, 4] = 15.0;
        //    K2[10, 4] = 10.0;
        //    //
        //    K2[4, 6] = K2[6, 4];
        //    K2[4, 7] = K2[7, 4];
        //    K2[4, 9] = K2[9, 4];
        //    K2[4, 10] = K2[10, 4];
        //    //
        //    K2[6, 6] = 60.0;
        //    K2[7, 6] = -30.0;
        //    K2[9, 6] = -60.0;
        //    K2[10, 6] = -30.0;
        //    //
        //    K2[6, 7] = K2[7, 6];
        //    K2[6, 9] = K2[9, 6];
        //    K2[6, 10] = K2[10, 6];
        //    //
        //    K2[7, 7] = 20.0;
        //    K2[9, 7] = 30.0;
        //    K2[10, 7] = 10.0;
        //    //
        //    K2[7, 9] = K2[9, 7];
        //    K2[7, 10] = K2[10, 7];
        //    //
        //    K2[9, 9] = 60.0;
        //    K2[10, 9] = 30.0;
        //    //
        //    K2[9, 10] = K2[10, 9];
        //    K2[10, 10] = 20.0;
        //    double scal2 = Math.Pow(alpha / beta, 2.0);
        //    K2 = MatrixOperations.ScalarMatrixProduct(scal2, K2);
        //    //K3
        //    K3[0, 0] = 30.0;
        //    K3[1, 0] = -15.0;
        //    K3[2, 0] = 15.0;
        //    K3[3, 0] = -30.0;
        //    K3[5, 0] = -15.0;
        //    K3[6, 0] = -30.0;
        //    K3[7, 0] = 15.0;
        //    K3[9, 0] = 30.0;
        //    //
        //    K3[0, 1] = K3[1, 0];
        //    K3[0, 2] = K3[2, 0];
        //    K3[0, 3] = K3[3, 0];
        //    K3[0, 5] = K3[5, 0];
        //    K3[0, 6] = K3[6, 0];
        //    K3[0, 7] = K3[7, 0];
        //    K3[0, 9] = K3[9, 0];
        //    //
        //    K3[2, 1] = -15.0;
        //    K3[6, 1] = 15.0;
        //    //
        //    K3[1, 2] = K3[2, 1];
        //    K3[1, 6] = K3[6, 1];
        //    //
        //    K3[3, 2] = -15.0;
        //    K3[2, 3] = K3[3, 2];
        //    //
        //    K3[3, 3] = 30.0;
        //    K3[4, 3] = 15.0;
        //    K3[5, 3] = 15.0;
        //    K3[6, 3] = 30.0;
        //    K3[9, 3] = -30.0;
        //    K3[10, 3] = -15.0;
        //    //
        //    K3[3, 4] = K3[4, 3];
        //    K3[3, 5] = K3[5, 3];
        //    K3[3, 6] = K3[6, 3];
        //    K3[3, 9] = K3[9, 3];
        //    K3[3, 10] = K3[10, 3];
        //    //
        //    K3[4, 5] = 15.0;
        //    K3[5, 4] = K3[4, 5];
        //    K3[9, 4] = -15.0;
        //    K3[4, 9] = K3[9, 4];
        //    //
        //    K3[6, 6] = 30.0;
        //    K3[7, 6] = -15.0;
        //    K3[8, 6] = -15.0;
        //    K3[9, 6] = -30.0;
        //    K3[11, 6] = 15.0;

        //    K3[6, 7] = K3[7, 6];
        //    K3[6, 8] = K3[8, 6];
        //    K3[6, 9] = K3[9, 6];
        //    K3[6, 11] = K3[11, 6];
        //    //
        //    K3[8, 7] = 15.0;
        //    K3[7, 8] = K3[8, 7];
        //    //
        //    K3[9, 9] = 30.0;
        //    K3[10, 9] = 15.0;
        //    K3[11, 9] = -15.0;
        //    K3[9, 10] = K3[10, 9];
        //    K3[9, 11] = K3[11, 9];
        //    //
        //    K3[11, 10] = -15.0;
        //    K3[10, 11] = K3[11, 10];
        //    K3 = MatrixOperations.ScalarMatrixProduct(Properties.PoissonRatio, K3);
        //    //K4
        //    K4[0, 0] = 84.0;
        //    K4[1, 0] = -6.0;
        //    K4[2, 0] = 6.0;
        //    K4[3, 0] = -84.0;
        //    K4[4, 0] = -6.0;
        //    K4[5, 0] = -6.0;
        //    K4[6, 0] = -84.0;
        //    K4[7, 0] = 6.0;
        //    K4[8, 0] = 6.0;
        //    K4[9, 0] = 84.0;
        //    K4[10, 0] = 6.0;
        //    K4[11, 0] = -6.0;
        //    //
        //    K4[0, 1] = K4[1, 0];
        //    K4[0, 2] = K4[2, 0];
        //    K4[0, 3] = K4[3, 0];
        //    K4[0, 4] = K4[4, 0];
        //    K4[0, 5] = K4[5, 0];
        //    K4[0, 6] = K4[6, 0];
        //    K4[0, 7] = K4[7, 0];
        //    K4[0, 8] = K4[8, 0];
        //    K4[0, 9] = K4[9, 0];
        //    K4[0, 10] = K4[10, 0];
        //    K4[0, 11] = K4[11, 0];
        //    //
        //    K4[1, 1] = 8.0;
        //    K4[3, 1] = 6.0;
        //    K4[4, 1] = -2.0;
        //    K4[6, 1] = 6.0;
        //    K4[7, 1] = -8.0;
        //    K4[9, 1] = -6.0;
        //    K4[10, 1] = 2.0;
        //    //
        //    K4[1, 3] = K4[3, 1];
        //    K4[1, 4] = K4[4, 1];
        //    K4[1, 6] = K4[6, 1];
        //    K4[1, 7] = K4[7, 1];
        //    K4[1, 9] = K4[9, 1];
        //    K4[1, 10] = K4[10, 1];
        //    //
        //    K4[2, 2] = 8.0;
        //    K4[3, 2] = -6.0;
        //    K4[5, 2] = -8.0;
        //    K4[6, 2] = -6.0;
        //    K4[8, 2] = -2.0;
        //    K4[9, 2] = 6.0;
        //    K4[11, 2] = 2.0;
        //    //
        //    K4[2, 3] = K4[3, 2];
        //    K4[2, 5] = K4[5, 2];
        //    K4[2, 6] = K4[6, 2];
        //    K4[2, 8] = K4[8, 2];
        //    K4[2, 9] = K4[9, 2];
        //    K4[2, 11] = K4[11, 2];
        //    //
        //    K4[3, 3] = 84.0;
        //    K4[4, 3] = 6.0;
        //    K4[5, 3] = 6.0;
        //    K4[6, 3] = 84.0;
        //    K4[7, 3] = -6.0;
        //    K4[8, 3] = -6.0;
        //    K4[9, 3] = -84.0;
        //    K4[10, 3] = -6.0;
        //    K4[11, 3] = 6.0;
        //    //
        //    K4[3, 4] = K4[4, 3];
        //    K4[3, 5] = K4[5, 3];
        //    K4[3, 6] = K4[6, 3];
        //    K4[3, 7] = K4[7, 3];
        //    K4[3, 8] = K4[8, 3];
        //    K4[3, 9] = K4[9, 3];
        //    K4[3, 10] = K4[10, 3];
        //    K4[3, 11] = K4[11, 3];
        //    //
        //    K4[4, 4] = 8.0;
        //    K4[6, 4] = 6.0;
        //    K4[7, 4] = 2.0;
        //    K4[9, 4] = -6.0;
        //    K4[10, 4] = -8.0;
        //    //
        //    K4[4, 6] = K4[6, 4];
        //    K4[4, 7] = K4[7, 4];
        //    K4[4, 9] = K4[9, 4];
        //    K4[4, 10] = K4[10, 4];
        //    //
        //    K4[5, 5] = 8.0;
        //    K4[6, 5] = 6.0;
        //    K4[8, 5] = 2.0;
        //    K4[9, 5] = -6.0;
        //    K4[11, 5] = -2.0;
        //    //
        //    K4[5, 6] = K4[6, 5];
        //    K4[5, 8] = K4[8, 5];
        //    K4[5, 9] = K4[9, 5];
        //    K4[5, 11] = K4[11, 5];
        //    //
        //    K4[6, 6] = 84.0;
        //    K4[7, 6] = -6.0;
        //    K4[8, 6] = -6.0;
        //    K4[9, 6] = -84.0;
        //    K4[10, 6] = -6.0;
        //    K4[11, 6] = 6.0;
        //    //
        //    K4[6, 7] = K4[7, 6];
        //    K4[6, 8] = K4[8, 6];
        //    K4[6, 9] = K4[9, 6];
        //    K4[6, 10] = K4[10, 6];
        //    K4[6, 11] = K4[11, 6];
        //    //
        //    K4[7, 7] = 8.0;
        //    K4[9, 7] = 6.0;
        //    K4[10, 7] = -2.0;
        //    //
        //    K4[7, 9] = K4[9, 7];
        //    K4[7, 10] = K4[10, 7];
        //    //
        //    K4[8, 8] = 8.0;
        //    K4[9, 8] = 6.0;
        //    K4[11, 8] = -8.0;
        //    //
        //    K4[8, 9] = K4[9, 8];
        //    K4[8, 11] = K4[11, 8];
        //    //
        //    K4[9, 9] = 84.0;
        //    K4[10, 9] = 6.0;
        //    K4[11, 9] = -6.0;
        //    //
        //    K4[9, 10] = K4[10, 9];
        //    K4[9, 11] = K4[11, 9];
        //    //
        //    K4[10, 10] = 8.0;
        //    K4[11, 11] = 8.0;
        //    K4 = MatrixOperations.ScalarMatrixProduct((1.0 - Properties.PoissonRatio) / 2.0, K4);
        //    //
        //    double scalar = Properties.YoungMod * Math.Pow(Properties.Thickness, 3.0) / (12.0 * (1.0 - Math.Pow(Properties.PoissonRatio, 2.0)) * 60.0 * alpha * beta);
        //    K = MatrixOperations.MatrixAddition(K1, MatrixOperations.MatrixAddition(K2,
        //        MatrixOperations.MatrixAddition(K3, K4)));
        //    K = MatrixOperations.ScalarMatrixProduct(scalar, K);
        //    double[,] L = new double[12, 12];
        //    L[0, 0] = 1.0;
        //    L[1, 1] = 2.0 * beta;
        //    L[2, 2] = 2.0 * alpha;
        //    L[3, 3] = 1.0;
        //    L[4, 4] = 2.0 * beta;
        //    L[5, 5] = 2.0 * alpha;
        //    L[6, 6] = 1.0;
        //    L[7, 7] = 2.0 * beta;
        //    L[8, 8] = 2.0 * alpha;
        //    L[9, 9] = 1.0;
        //    L[10, 10] = 2.0 * beta;
        //    L[11, 11] = 2.0 * alpha;
        //    K = MatrixOperations.MatrixProduct(L, K);
        //    K = MatrixOperations.MatrixProduct(K, L);
        //    return K;
        //}
        private double[,] CreateLocalStiffnessMatrixBendingPart(double[] xUpdated)
        {
            double[] positionVector1 = new double[] { xUpdated[0], xUpdated[1], xUpdated[2] };
            double[] positionVector2 = new double[] { xUpdated[6], xUpdated[7], xUpdated[8] };
            double[] positionVector3 = new double[] { xUpdated[12], xUpdated[13], xUpdated[14] };
            double[] positionVector4 = new double[] { xUpdated[18], xUpdated[19], xUpdated[20] };

            double alpha = VectorOperations.VectorNorm2(VectorOperations.VectorVectorSubtraction(positionVector2, positionVector1)) / 2.0;
            double beta = VectorOperations.VectorNorm2(VectorOperations.VectorVectorSubtraction(positionVector4, positionVector1)) / 2.0;

            double[,] K = new double[12, 12];
            double[,] matrix = new double[12, 12];
            double p = Math.Pow(alpha / beta, 2.0);
            double pInv = Math.Pow(beta / alpha, 2.0);

            //
            matrix[0, 0] = 60.0 * p + 60.0 * pInv + 42.0 - 12.0 * Properties.PoissonRatio;
            matrix[1, 0] = beta * (60.0 * p + 6.0 + 24.0 * Properties.PoissonRatio);
            matrix[2, 0] = alpha * (60.0 * pInv + 6.0 + 24.0 * Properties.PoissonRatio);
            matrix[3, 0] = 30.0 * p - 60.0 * pInv - 42.0 + 12.0 * Properties.PoissonRatio;
            matrix[4, 0] = beta * (30.0 * p - 6.0 - 24.0 * Properties.PoissonRatio);
            matrix[5, 0] = -alpha * (-60.0 * pInv - 6.0 + 6.0 * Properties.PoissonRatio);
            matrix[6, 0] = -30.0 * p - 30.0 * pInv + 42.0 - 12.0 * Properties.PoissonRatio;
            matrix[7, 0] = -beta * (-30.0 * p + 6.0 - 6.0 * Properties.PoissonRatio);
            matrix[8, 0] = -alpha * (-30.0 * pInv + 6.0 - 6.0 * Properties.PoissonRatio);
            matrix[9, 0] = -60.0 * p + 30.0 * pInv - 42.0 + 12.0 * Properties.PoissonRatio;
            matrix[10, 0] = -beta * (-60.0 * p - 6.0 + 6.0 * Properties.PoissonRatio);
            matrix[11, 0] = alpha * (30.0 * pInv - 6.0 - 24.0 * Properties.PoissonRatio);
            //
            matrix[0, 1] = matrix[1, 0];
            matrix[0, 2] = matrix[2, 0];
            matrix[0, 3] = matrix[3, 0];
            matrix[0, 4] = matrix[4, 0];
            matrix[0, 5] = matrix[5, 0];
            matrix[0, 6] = matrix[6, 0];
            matrix[0, 7] = matrix[7, 0];
            matrix[0, 8] = matrix[8, 0];
            matrix[0, 9] = matrix[9, 0];
            matrix[0, 10] = matrix[10, 0];
            matrix[0, 11] = matrix[11, 0];
            //
            matrix[1, 1] = beta * beta * (80.0 * p + 16.0 - 16.0 * Properties.PoissonRatio);
            matrix[2, 1] = 60.0 * Properties.PoissonRatio * alpha * beta;
            matrix[3, 1] = beta * (30.0 * p - 6.0 - 24.0 * Properties.PoissonRatio);
            matrix[4, 1] = beta * beta * (40.0 * p - 16.0 + 16.0 * Properties.PoissonRatio);
            matrix[5, 1] = 0.0;
            matrix[6, 1] = beta * (-30.0 * p + 6.0 - 6.0 * Properties.PoissonRatio);
            matrix[7, 1] = beta * beta * (20.0 * p + 4.0 - 4.0 * Properties.PoissonRatio);
            matrix[8, 1] = 0.0;
            matrix[9, 1] = beta * (-60.0 * p - 6.0 + 6.0 * Properties.PoissonRatio);
            matrix[10, 1] = beta * beta * (40.0 * p - 4.0 + 4.0 * Properties.PoissonRatio);
            matrix[11, 1] = 0.0;
            //
            matrix[1, 2] = matrix[2, 1];
            matrix[1, 3] = matrix[3, 1];
            matrix[1, 4] = matrix[4, 1];
            matrix[1, 5] = matrix[5, 1];
            matrix[1, 6] = matrix[6, 1];
            matrix[1, 7] = matrix[7, 1];
            matrix[1, 8] = matrix[8, 1];
            matrix[1, 9] = matrix[9, 1];
            matrix[1, 10] = matrix[10, 1];
            matrix[1, 11] = matrix[11, 1];
            //
            matrix[2, 2] = alpha * alpha * (80.0 * pInv + 16.0 - 16.0 * Properties.PoissonRatio);
            matrix[3, 2] = alpha * (-60.0 * pInv - 6.0 + 6.0 * Properties.PoissonRatio);
            matrix[4, 2] = 0.0;
            matrix[5, 2] = alpha * alpha * (40.0 * pInv - 4.0 + 4.0 * Properties.PoissonRatio);
            matrix[6, 2] = alpha * (-30.0 * pInv + 6.0 - 6.0 * Properties.PoissonRatio);
            matrix[7, 2] = 0.0;
            matrix[8, 2] = alpha * alpha * (20.0 * pInv + 4.0 - 4.0 * Properties.PoissonRatio);
            matrix[9, 2] = alpha * (30.0 * pInv - 6.0 - 24.0 * Properties.PoissonRatio);
            matrix[10, 2] = 0.0;
            matrix[11, 2] = alpha * alpha * (40.0 * pInv - 16.0 + 16.0 * Properties.PoissonRatio);
            //
            matrix[2, 3] = matrix[3, 2];
            matrix[2, 4] = matrix[4, 2];
            matrix[2, 5] = matrix[5, 2];
            matrix[2, 6] = matrix[6, 2];
            matrix[2, 7] = matrix[7, 2];
            matrix[2, 8] = matrix[8, 2];
            matrix[2, 9] = matrix[9, 2];
            matrix[2, 10] = matrix[10, 2];
            matrix[2, 11] = matrix[11, 2];
            //
            matrix[3, 3] = matrix[0, 0];
            matrix[4, 3] = matrix[1, 0];
            matrix[5, 3] = -matrix[2, 0];
            matrix[6, 3] = matrix[9, 0];
            matrix[7, 3] = matrix[10, 0];
            matrix[8, 3] = -matrix[11, 0];
            matrix[9, 3] = matrix[6, 0];
            matrix[10, 3] = matrix[7, 0];
            matrix[11, 3] = -matrix[8, 0];
            //
            matrix[3, 4] = matrix[4, 3];
            matrix[3, 5] = matrix[5, 3];
            matrix[3, 6] = matrix[6, 3];
            matrix[3, 7] = matrix[7, 3];
            matrix[3, 8] = matrix[8, 3];
            matrix[3, 9] = matrix[9, 3];
            matrix[3, 10] = matrix[10, 3];
            matrix[3, 11] = matrix[11, 3];
            //
            matrix[4, 4] = matrix[1, 1];
            matrix[5, 4] = -matrix[2, 1];
            matrix[6, 4] = matrix[9, 1];
            matrix[7, 4] = matrix[10, 1];
            matrix[8, 4] = -matrix[11, 1];
            matrix[9, 4] = matrix[6, 1];
            matrix[10, 4] = matrix[7, 1];
            matrix[11, 4] = -matrix[8, 1];
            //
            matrix[4, 5] = matrix[5, 4];
            matrix[4, 6] = matrix[6, 4];
            matrix[4, 7] = matrix[7, 4];
            matrix[4, 8] = matrix[8, 4];
            matrix[4, 9] = matrix[9, 4];
            matrix[4, 10] = matrix[10, 4];
            matrix[4, 11] = matrix[11, 4];
            //
            matrix[5, 5] = matrix[2, 2];
            matrix[6, 5] = -matrix[9, 2];
            matrix[7, 5] = -matrix[10, 2];
            matrix[8, 5] = matrix[11, 2];
            matrix[9, 5] = -matrix[6, 2];
            matrix[10, 5] = -matrix[7, 2];
            matrix[11, 5] = matrix[8, 2];
            //
            matrix[5, 6] = matrix[6, 5];
            matrix[5, 7] = matrix[7, 5];
            matrix[5, 8] = matrix[8, 5];
            matrix[5, 9] = matrix[9, 5];
            matrix[5, 10] = matrix[10, 5];
            matrix[5, 11] = matrix[11, 5];
            //
            matrix[6, 6] = matrix[0, 0];
            matrix[7, 6] = -matrix[1, 0];
            matrix[8, 6] = -matrix[2, 0];
            matrix[9, 6] = matrix[3, 0];
            matrix[10, 6] = -matrix[4, 0];
            matrix[11, 6] = -matrix[5, 0];
            //
            matrix[6, 7] = matrix[7, 6];
            matrix[6, 8] = matrix[8, 6];
            matrix[6, 9] = matrix[9, 6];
            matrix[6, 10] = matrix[10, 6];
            matrix[6, 11] = matrix[11, 6];
            //
            matrix[7, 7] = matrix[1, 1];
            matrix[8, 7] = matrix[2, 1];
            matrix[9, 7] = -matrix[3, 1];
            matrix[10, 7] = matrix[4, 1];
            matrix[11, 7] = matrix[5, 1];
            //
            matrix[7, 8] = matrix[8, 7];
            matrix[7, 9] = matrix[9, 7];
            matrix[7, 10] = matrix[10, 7];
            matrix[7, 11] = matrix[11, 7];
            //
            matrix[8, 8] = matrix[2, 2];
            matrix[9, 8] = -matrix[3, 2];
            matrix[10, 8] = matrix[4, 2];
            matrix[11, 8] = matrix[5, 2];
            //
            matrix[8, 9] = matrix[9, 8];
            matrix[8, 10] = matrix[10, 8];
            matrix[8, 11] = matrix[11, 8];
            //
            matrix[9, 9] = matrix[0, 0];
            matrix[10, 9] = -matrix[1, 0];
            matrix[11, 9] = matrix[2, 0];
            //
            matrix[9, 10] = matrix[10, 9];
            matrix[9, 11] = matrix[11, 9];
            //
            matrix[10, 10] = matrix[1, 1];
            matrix[11, 10] = -matrix[2, 1];
            matrix[10, 11] = matrix[11, 10];
            //
            matrix[11, 11] = matrix[2, 2];
            //
            double scalar = (1.0 / (60.0 * alpha * beta)) * Properties.YoungMod * Math.Pow(Properties.Thickness, 3.0) /
                (12.0 * (1.0 - Math.Pow(Properties.PoissonRatio, 2.0)));
            K = MatrixOperations.ScalarMatrixProduct(scalar, matrix);
            return K;
        }
        private double[,] CreateLocalMassMatrixBendingPart(double[] xUpdated)
        {
            //double[] xUpdated = UpdateNodalCoordinates(DisplacementVector);
            double[] positionVector1 = new double[] { xUpdated[0], xUpdated[1], xUpdated[2] };
            double[] positionVector2 = new double[] { xUpdated[6], xUpdated[7], xUpdated[8] };
            double[] positionVector3 = new double[] { xUpdated[12], xUpdated[13], xUpdated[14] };
            double[] positionVector4 = new double[] { xUpdated[18], xUpdated[19], xUpdated[20] };

            double alpha = VectorOperations.VectorNorm2(VectorOperations.VectorVectorSubtraction(positionVector2, positionVector1)) / 2.0;
            double beta = VectorOperations.VectorNorm2(VectorOperations.VectorVectorSubtraction(positionVector4, positionVector1)) / 2.0;
            double[,] M = new double[12, 12];
            //
            double totalMass = 4.0 * Properties.Density * Properties.Thickness * alpha * beta;
            double alphaX = 1.0 / (64.0 * Properties.Thickness * Properties.Thickness) -
                1.0 / (6.0 * beta * beta);
            double alphaY = 1.0 / (6.0 * alpha * alpha) +
                1.0 / (64.0 * Properties.Thickness * Properties.Thickness);
            M[0, 0] = 1.0 / 4.0;
            M[1, 1] = alphaX;
            M[2, 2] = alphaY;
            M[3, 3] = 1.0 / 4.0;
            M[4, 4] = alphaX;
            M[5, 5] = alphaY;
            M[6, 6] = 1.0 / 4.0;
            M[7, 7] = alphaX;
            M[8, 8] = alphaY;
            M[9, 9] = 1.0 / 4.0;
            M[10, 10] = alphaX;
            M[11, 11] = alphaY;
            M = MatrixOperations.ScalarMatrixProduct(totalMass, M);
            return M;
            //
            //M[0, 0] = 1727.0;
            //M[1, 0] = 461.0 * beta;
            //M[2, 0] = -461.0 * alpha;
            //M[3, 0] = 613.0;
            //M[4, 0] = 199.0 * beta;
            //M[5, 0] = 274.0 * alpha;
            //M[6, 0] = 197.0;
            //M[7, 0] = -116.0 * beta;
            //M[8, 0] = 116.0 * alpha;
            //M[9, 0] = 613.0;
            //M[10, 0] = -274.0 * beta;
            //M[11, 0] = -199.0 * alpha;

            //M[0, 1] = M[1, 0];
            //M[0, 2] = M[2, 0];
            //M[0, 3] = M[3, 0];
            //M[0, 4] = M[4, 0];
            //M[0, 5] = M[5, 0];
            //M[0, 6] = M[6, 0];
            //M[0, 7] = M[7, 0];
            //M[0, 8] = M[8, 0];
            //M[0, 9] = M[9, 0];
            //M[0, 10] = M[10, 0];
            //M[0, 11] = M[11, 0];
            ////
            //M[1, 1] = 160.0 * beta * beta;
            //M[2, 1] = -126.0 * alpha * beta;
            //M[3, 1] = 199.0 * beta;
            //M[4, 1] = 80.0 * beta * beta;
            //M[5, 1] = 84.0 * alpha * beta;
            //M[6, 1] = 116.0 * beta;
            //M[7, 1] = -60.0 * beta * beta;
            //M[8, 1] = 56.0 * alpha * beta;
            //M[9, 1] = 274.0 * beta;
            //M[10, 1] = -120.0 * beta * beta;
            //M[11, 1] = -84.0 * alpha * beta;

            //M[1, 2] = M[2, 1];
            //M[1, 3] = M[3, 1];
            //M[1, 4] = M[4, 1];
            //M[1, 5] = M[5, 1];
            //M[1, 6] = M[6, 1];
            //M[1, 7] = M[7, 1];
            //M[1, 8] = M[8, 1];
            //M[1, 9] = M[9, 1];
            //M[1, 10] = M[10, 1];
            //M[1, 11] = M[11, 1];
            ////
            //M[2, 2] = 160.0 * alpha * alpha;
            //M[3, 2] = -272.0 * alpha;
            //M[4, 2] = -84.0 * alpha * beta;
            //M[5, 2] = -120.0 * alpha * alpha;
            //M[6, 2] = -116.0 * alpha;
            //M[7, 2] = 56.0 * alpha * beta;
            //M[8, 2] = -60.0 * alpha * alpha;
            //M[9, 2] = -199.0 * alpha;
            //M[10, 2] = 84.0 * alpha * beta;
            //M[11, 2] = 80.0 * alpha * alpha;

            //M[2, 3] = M[3, 2];
            //M[2, 4] = M[4, 2];
            //M[2, 5] = M[5, 2];
            //M[2, 6] = M[6, 2];
            //M[2, 7] = M[7, 2];
            //M[2, 8] = M[8, 2];
            //M[2, 9] = M[9, 2];
            //M[2, 10] = M[10, 2];
            //M[2, 11] = M[11, 2];
            ////
            //M[3, 3] = 1727.0;
            //M[4, 3] = 461.0 * beta;
            //M[5, 3] = 461.0 * alpha;
            //M[6, 3] = 613.0;
            //M[7, 3] = -274.0 * beta;
            //M[8, 3] = 199.0 * alpha;
            //M[9, 3] = 197.0;
            //M[10, 3] = -116.0 * beta;
            //M[11, 3] = -116.0 * alpha;

            //M[3, 4] = M[4, 3];
            //M[3, 5] = M[5, 3];
            //M[3, 6] = M[6, 3];
            //M[3, 7] = M[7, 3];
            //M[3, 8] = M[8, 3];
            //M[3, 9] = M[9, 3];
            //M[3, 10] = M[10, 3];
            //M[3, 11] = M[11, 3];
            ////
            //M[4, 4] = 160.0 * beta * beta;
            //M[5, 4] = 126.0 * alpha * beta;
            //M[6, 4] = 274.0 * beta;
            //M[7, 4] = -120.0 * beta * beta;
            //M[8, 4] = 84.0 * alpha * beta;
            //M[9, 4] = 116.0 * beta;
            //M[10, 4] = -60.0 * beta * beta;
            //M[11, 4] = -56.0 * alpha * beta;

            //M[4, 5] = M[5, 4];
            //M[4, 6] = M[6, 4];
            //M[4, 7] = M[7, 4];
            //M[4, 8] = M[8, 4];
            //M[4, 9] = M[9, 4];
            //M[4, 10] = M[10, 4];
            //M[4, 11] = M[11, 4];

            ////
            //M[5, 5] = 160.0 * alpha * alpha;
            //M[6, 5] = 199.0 * alpha;
            //M[7, 5] = -84.0 * alpha * beta;
            //M[8, 5] = 80.0 * alpha * alpha;
            //M[9, 5] = 116.0 * alpha;
            //M[10, 5] = -56.0 * alpha * beta;
            //M[11, 5] = -60.0 * alpha * alpha;

            //M[5, 6] = M[6, 5];
            //M[5, 7] = M[7, 5];
            //M[5, 8] = M[8, 5];
            //M[5, 9] = M[9, 5];
            //M[5, 10] = M[10, 5];
            //M[5, 11] = M[11, 5];
            ////
            //M[6, 6] = 1727.0;
            //M[7, 6] = -461.0 * beta;
            //M[8, 6] = 461.0 * alpha;
            //M[9, 6] = 613.0;
            //M[10, 6] = -199.0 * beta;
            //M[11, 6] = -274.0 * alpha;

            //M[6, 7] = M[7, 6];
            //M[6, 8] = M[8, 6];
            //M[6, 9] = M[9, 6];
            //M[6, 10] = M[10, 6];
            //M[6, 11] = M[11, 6];
            ////
            //M[7, 7] = 160.0 * beta * beta;
            //M[8, 7] = -126.0 * alpha * beta;
            //M[9, 7] = -199.0 * beta;
            //M[10, 7] = 80.0 * beta * beta;
            //M[11, 7] = 84.0 * alpha * beta;

            //M[7, 8] = M[8, 7];
            //M[7, 9] = M[9, 7];
            //M[7, 10] = M[10, 7];
            //M[7, 11] = M[11, 7];
            ////
            //M[8, 8] = 160.0 * alpha * alpha;
            //M[9, 8] = 274.0 * alpha;
            //M[10, 8] = -84.0 * alpha * beta;
            //M[11, 8] = -120.0 * alpha * alpha;

            //M[8, 9] = M[9, 8];
            //M[8, 10] = M[10, 8];
            //M[8, 11] = M[11, 8];
            ////
            //M[9, 9] = 1727.0;
            //M[10, 9] = -461.0 * beta;
            //M[11, 9] = -461.0 * alpha;

            //M[9, 10] = M[10, 9];
            //M[9, 11] = M[11, 9];
            ////
            //M[10, 10] = 160.0 * beta * beta;
            //M[11, 10] = 126.0 * alpha * beta;

            //M[10, 11] = M[11, 10];
            //M[11, 11] = 160.0 * alpha * alpha;

            //double scalar = Properties.Density * Properties.Thickness * alpha * beta / 3150.0;
            //M = MatrixOperations.ScalarMatrixProduct(scalar, M);
            //return M;
        }
        private double[,] transformationMatrix(double[] xUpdated)
        {
            double[,] trMat = new double[6, 6];
            double[] positionVector1 = new double[] { xUpdated[0], xUpdated[1], xUpdated[2] };
            double[] positionVector2 = new double[] { xUpdated[6], xUpdated[7], xUpdated[8] };
            double[] positionVector3 = new double[] { xUpdated[12], xUpdated[13], xUpdated[14] };
            double[] positionVector4 = new double[] { xUpdated[18], xUpdated[19], xUpdated[20] };

            double[] n1 = VectorOperations.VectorScalarProduct(VectorOperations.VectorVectorSubtraction(positionVector2, positionVector1),
                1.0 / VectorOperations.VectorNorm2(VectorOperations.VectorVectorSubtraction(positionVector2, positionVector1)));
            double[] n2 = VectorOperations.VectorScalarProduct(VectorOperations.VectorVectorSubtraction(positionVector4, positionVector1),
                            1.0 / VectorOperations.VectorNorm2(VectorOperations.VectorVectorSubtraction(positionVector4, positionVector1)));
            double[] n3 = VectorOperations.VectorScalarProduct(VectorOperations.VectorCrossProduct(n1, n2),
                1.0 / VectorOperations.VectorNorm2(VectorOperations.VectorCrossProduct(n1, n2)));
            double phix_X = n1[0];
            double phix_Y = n1[1];
            double phix_Z = n1[2];
            double phiy_X = n2[0];
            double phiy_Y = n2[1];
            double phiy_Z = n2[2];
            double phiz_X = n3[0];
            double phiz_Y = n3[1];
            double phiz_Z = n3[2];
            trMat[0, 0] = phix_X; trMat[0, 1] = phix_Y; trMat[0, 2] = phix_Z;
            trMat[1, 0] = phiy_X; trMat[1, 1] = phiy_Y; trMat[1, 2] = phiy_Z;
            trMat[2, 0] = phiz_X; trMat[2, 1] = phiz_Y; trMat[2, 2] = phiz_Z;
            trMat[3, 3] = phix_X; trMat[3, 4] = phix_Y; trMat[3, 5] = phix_Z;
            trMat[4, 3] = phiy_X; trMat[4, 4] = phiy_Y; trMat[4, 5] = phiy_Z;
            trMat[5, 3] = phiz_X; trMat[5, 4] = phiz_Y; trMat[5, 5] = phiz_Z;

            return trMat;
        }
        public double[,] CreateGlobalStiffnessMatrix()
        {
            double[,] globalK = new double[24, 24];
            double[,] localK = new double[24, 24];
            double[,] transformationMat = new double[24, 24];
            double[] xUpdated = UpdateNodalCoordinates(DisplacementVector);
            double[,] tMat = transformationMatrix(xUpdated);
            double[,] planeStressPart = CreateLocalStiffnessMatrixPlaneStressPart(xUpdated);
            double[,] bendingPart = CreateLocalStiffnessMatrixBendingPart(xUpdated);
            List<int> dofList1 = new List<int>();
            dofList1.Add(0);
            dofList1.Add(1);
            dofList1.Add(6);
            dofList1.Add(7);
            dofList1.Add(12);
            dofList1.Add(13);
            dofList1.Add(18);
            dofList1.Add(19);
            for(int i = 0; i < planeStressPart.GetLength(0); i++)
            {
                for(int j = 0; j < planeStressPart.GetLength(0); j++)
                {
                    int ind1 = dofList1[i];
                    int ind2 = dofList1[j];
                    localK[ind1, ind2] = planeStressPart[i, j];
                }
            }
            List<int> dofList2 = new List<int>();
            dofList2.Add(2);
            dofList2.Add(3);
            dofList2.Add(4);
            dofList2.Add(8);
            dofList2.Add(9);
            dofList2.Add(10);
            dofList2.Add(14);
            dofList2.Add(15);
            dofList2.Add(16);
            dofList2.Add(20);
            dofList2.Add(21);
            dofList2.Add(22);
            for (int i = 0; i < bendingPart.GetLength(0); i++)
            {
                for (int j = 0; j < bendingPart.GetLength(0); j++)
                {
                    int ind1 = dofList2[i];
                    int ind2 = dofList2[j];
                    localK[ind1, ind2] = bendingPart[i, j];
                }
            }
            transformationMat = MatrixOperations.CreateDiagonalMatrix(globalK.GetLength(0), tMat);
            globalK = MatrixOperations.MatrixProduct(MatrixOperations.Transpose(transformationMat),
                MatrixOperations.MatrixProduct(localK, transformationMat));
            //globalK = ArtificialStiffness(globalK, globalK[0, 0], xUpdated);
            return globalK;
        }
        private double[,] ArtificialStiffness(double[,] stifMat, double scalar, double[] xUpdated)
        {
            double[] positionVector1 = new double[] { xUpdated[0], xUpdated[1], xUpdated[2] };
            double[] positionVector2 = new double[] { xUpdated[6], xUpdated[7], xUpdated[8] };
            double[] positionVector4 = new double[] { xUpdated[18], xUpdated[19], xUpdated[20] };
            double[] v1 = VectorOperations.VectorVectorSubtraction(positionVector2, positionVector1);
            double[] v2 = VectorOperations.VectorVectorSubtraction(positionVector4, positionVector1);                           
            double[] v3 = VectorOperations.VectorCrossProduct(v1, v2);
            if(v3[0] == v3[1] && v3[1] == 0.0 && v3[2] != 0.0)
            {
                if(stifMat[5, 5]==0.0)stifMat[5,5] = scalar;
                if (stifMat[11, 11] == 0.0) stifMat[11, 11] = scalar;
                if (stifMat[17, 17] == 0.0) stifMat[17, 17] = scalar;
                if (stifMat[23, 23] == 0.0) stifMat[23, 23] = scalar;
            }
            return stifMat;
        }
        public double[,] CreateMassMatrix()
        {
            double[,] globalM = new double[24, 24];
            double[,] localM = new double[24, 24];
            double[,] transformationMat = new double[24, 24];
            double[] xUpdated = UpdateNodalCoordinates(DisplacementVector);
            double[,] tMat = transformationMatrix(xUpdated);
            double[,] planeStressPart = CreateLocalMassMatrixPlaneStressPart(xUpdated);
            double[,] bendingPart = CreateLocalMassMatrixBendingPart(xUpdated);
            List<int> dofList1 = new List<int>();
            dofList1.Add(0);
            dofList1.Add(1);
            dofList1.Add(6);
            dofList1.Add(7);
            dofList1.Add(12);
            dofList1.Add(13);
            dofList1.Add(18);
            dofList1.Add(19);
            for (int i = 0; i < planeStressPart.GetLength(0); i++)
            {
                for (int j = 0; j < planeStressPart.GetLength(0); j++)
                {
                    int ind1 = dofList1[i];
                    int ind2 = dofList1[j];
                    localM[ind1, ind2] = planeStressPart[i, j];
                }
            }
            List<int> dofList2 = new List<int>();
            dofList2.Add(2);
            dofList2.Add(3);
            dofList2.Add(4);
            dofList2.Add(8);
            dofList2.Add(9);
            dofList2.Add(10);
            dofList2.Add(14);
            dofList2.Add(15);
            dofList2.Add(16);
            dofList2.Add(20);
            dofList2.Add(21);
            dofList2.Add(22);
            for (int i = 0; i < bendingPart.GetLength(0); i++)
            {
                for (int j = 0; j < bendingPart.GetLength(0); j++)
                {
                    int ind1 = dofList2[i];
                    int ind2 = dofList2[j];
                    localM[ind1, ind2] = bendingPart[i, j];
                }
            }
            transformationMat = MatrixOperations.CreateDiagonalMatrix(globalM.GetLength(0), tMat);
            globalM = MatrixOperations.MatrixProduct(MatrixOperations.Transpose(transformationMat),
                MatrixOperations.MatrixProduct(localM, transformationMat));
            return globalM;
        }

        public double[,] CreateDampingMatrix()
        {
            return new double[24, 24];
        }

        public double[] CreateInternalGlobalForcesVector()
        {
            double[] fLocal = new double[24];
            double[] xUpdated = UpdateNodalCoordinates(DisplacementVector);
            double[,] tMat = transformationMatrix(xUpdated);
            double[,] stiffnessPlaneStressPart = CreateLocalStiffnessMatrixPlaneStressPart(xUpdated);
            double[,] stiffnessBendingPart = CreateLocalStiffnessMatrixBendingPart(xUpdated);
            double[,] transformationMat = MatrixOperations.CreateDiagonalMatrix(fLocal.GetLength(0), tMat);
            //double[] xUpdatedLocal = VectorOperations.MatrixVectorProduct(transformationMat, xUpdated);
            //double[] fullDisplacementVectorLocal = VectorOperations.MatrixVectorProduct(transformationMat, Properties.DU.ToArray());
            double[] fullDisplacementVectorLocal = VectorOperations.MatrixVectorProduct(transformationMat, DisplacementVector);
            double[] displacementVectorLocal1 = new double[8];
            double[] displacementVectorLocal2 = new double[12];
            List<int> dofList1 = new List<int>();
            dofList1.Add(0);
            dofList1.Add(1);
            dofList1.Add(6);
            dofList1.Add(7);
            dofList1.Add(12);
            dofList1.Add(13);
            dofList1.Add(18);
            dofList1.Add(19);
            for (int i = 0; i < stiffnessPlaneStressPart.GetLength(0); i++)
            {
                int index = dofList1[i];
                displacementVectorLocal1[i] = fullDisplacementVectorLocal[index];
            }
            List<int> dofList2 = new List<int>();
            dofList2.Add(2);
            dofList2.Add(3);
            dofList2.Add(4);
            dofList2.Add(8);
            dofList2.Add(9);
            dofList2.Add(10);
            dofList2.Add(14);
            dofList2.Add(15);
            dofList2.Add(16);
            dofList2.Add(20);
            dofList2.Add(21);
            dofList2.Add(22);
            for (int i = 0; i < stiffnessBendingPart.GetLength(0); i++)
            {
                int indx = dofList2[i];
                displacementVectorLocal2[i] = fullDisplacementVectorLocal[indx];
            }
            double[] planeStressPart = VectorOperations.MatrixVectorProduct(stiffnessPlaneStressPart, displacementVectorLocal1);
            double[] bendingPart = VectorOperations.MatrixVectorProduct(stiffnessBendingPart, displacementVectorLocal2);
            for (int i = 0; i < planeStressPart.GetLength(0); i++)
            {
                int indx = dofList1[i];
                fLocal[indx] = planeStressPart[i];
            }
            for (int i = 0; i < bendingPart.GetLength(0); i++)
            {
                int indx = dofList2[i];
                fLocal[indx] = bendingPart[i];
            }
            double[] fGlobal = VectorOperations.MatrixVectorProduct(MatrixOperations.Transpose(transformationMat), fLocal);
            //double[] fGlobal = VectorOperations.MatrixVectorProduct(
            //    CreateGlobalStiffnessMatrix(),
            //    Properties.DU.ToArray());
            //double[] fGlobal = VectorOperations.MatrixVectorProduct(
            //    CreateGlobalStiffnessMatrix(),
            //    DisplacementVector);
            return fGlobal;
        }
    }
}
