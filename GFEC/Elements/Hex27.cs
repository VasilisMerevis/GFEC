using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GFEC
{
    class Hex27 : IElement
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

        public Hex27(IElementProperties properties, Dictionary<int, INode> nodes)
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
            ElementFreedomSignature[9] = new bool[] { true, true, true, false, false, false };
            ElementFreedomSignature[10] = new bool[] { true, true, true, false, false, false };
            ElementFreedomSignature[11] = new bool[] { true, true, true, false, false, false };
            ElementFreedomSignature[12] = new bool[] { true, true, true, false, false, false };
            ElementFreedomSignature[13] = new bool[] { true, true, true, false, false, false };
            ElementFreedomSignature[14] = new bool[] { true, true, true, false, false, false };
            ElementFreedomSignature[15] = new bool[] { true, true, true, false, false, false };
            ElementFreedomSignature[16] = new bool[] { true, true, true, false, false, false };
            ElementFreedomSignature[17] = new bool[] { true, true, true, false, false, false };
            ElementFreedomSignature[18] = new bool[] { true, true, true, false, false, false };
            ElementFreedomSignature[19] = new bool[] { true, true, true, false, false, false };
            ElementFreedomSignature[20] = new bool[] { true, true, true, false, false, false };
            ElementFreedomSignature[21] = new bool[] { true, true, true, false, false, false };
            ElementFreedomSignature[22] = new bool[] { true, true, true, false, false, false };
            ElementFreedomSignature[23] = new bool[] { true, true, true, false, false, false };
            ElementFreedomSignature[24] = new bool[] { true, true, true, false, false, false };
            ElementFreedomSignature[25] = new bool[] { true, true, true, false, false, false };
            ElementFreedomSignature[26] = new bool[] { true, true, true, false, false, false };
            ElementFreedomSignature[27] = new bool[] { true, true, true, false, false, false };

            DisplacementVector = new double[81];
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

        public Dictionary<int, INode> NodesAtFinalState()
        {
            throw new Exception("Method not implemenented");
        }
        public void UpdateIncrementalDisplacements(double[] deltaU)
        {
            throw new Exception("Needs to be removed. Has beeb used only for testing purposes");
        }
        public List<double[]> GetStressVector()
        {
            List<double[]> GpointsStress = new List<double[]>();
            double[,] E = CalculateStressStrainMatrix(Properties.YoungMod, Properties.PoissonRatio);
            for (int i = 0; i <= 2; i++)
            {
                for (int j = 0; j <= 2; j++)
                {
                    for (int k = 0; k <= 2; k++)
                    {
                        double[] gP = GaussPoints(i, j, k).Item1;
                        //double[] gW = GaussPoints(i, j, k).Item2;
                        Dictionary<string, double[]> localdN = CalculateShapeFunctionsLocalDerivatives(gP);
                        double[,] J = CalculateJacobian(localdN);
                        double[,] invJ = CalculateInverseJacobian(J).Item1;
                        //double detJ = CalculateInverseJacobian(J).Item2;
                        Dictionary<int, double[]> globaldN = CalculateShapeFunctionsGlobalDerivatives(localdN, invJ);
                        double[,] B = CalculateBMatrix(globaldN);
                        double[] strainVector = CalculateStrainsVector(B);
                        double[] stressVector = CalculateStressVector(E, strainVector);
                        GpointsStress.Add(stressVector);
                    }
                }
            }
            return GpointsStress;
        }
        public List<double[]> GetStrainVector()
        {
            List<double[]> GpointsDeformation = new List<double[]>();
            for (int i = 0; i <= 2; i++)
            {
                for (int j = 0; j <= 2; j++)
                {
                    for (int k = 0; k <= 2; k++)
                    {
                        double[] gP = GaussPoints(i, j, k).Item1;
                        //double[] gW = GaussPoints(i, j, k).Item2;
                        Dictionary<string, double[]> localdN = CalculateShapeFunctionsLocalDerivatives(gP);
                        double[,] J = CalculateJacobian(localdN);
                        double[,] invJ = CalculateInverseJacobian(J).Item1;
                        //double detJ = CalculateInverseJacobian(J).Item2;
                        Dictionary<int, double[]> globaldN = CalculateShapeFunctionsGlobalDerivatives(localdN, invJ);
                        double[,] B = CalculateBMatrix(globaldN);
                        double[] strainVector = CalculateStrainsVector(B);
                        GpointsDeformation.Add(strainVector);
                    }
                }
            }
            return GpointsDeformation;
        }
        public List<double[]> GetGaussPointsInPhysicalSpace()
        {
            List<double[]> GpointsPhysicalCoordinates = new List<double[]>();
            double[] xUpdated = UpdateNodalCoordinates(DisplacementVector);
            for (int i = 0; i <= 2; i++)
            {
                for (int j = 0; j <= 2; j++)
                {
                    for (int k = 0; k <= 2; k++)
                    {
                        double[] gP = GaussPoints(i, j, k).Item1;
                        double[] gaussPoint = VectorOperations.MatrixVectorProduct(CalculateShapeFunctionMatrix(gP[0], gP[1], gP[2]), xUpdated);
                        GpointsPhysicalCoordinates.Add(gaussPoint);
                    }
                }
            }
            return GpointsPhysicalCoordinates;
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
        public List<double[]> GetStressFromElementsNodes()
        {
            List<double[]> NodalStessVector = new List<double[]>();
            double[,] E = CalculateStressStrainMatrix(Properties.YoungMod, Properties.PoissonRatio);
            int nodesCount = 27;
            List<double[]> nodalParametricCoordinates = new List<double[]>();
            nodalParametricCoordinates.Add(new double[] { -1.0, -1.0, -1.0 });
            nodalParametricCoordinates.Add(new double[] { 1.0, -1.0, -1.0 });
            nodalParametricCoordinates.Add(new double[] { 1.0, 1.0, -1.0 });
            nodalParametricCoordinates.Add(new double[] { -1.0, 1.0, -1.0 });
            nodalParametricCoordinates.Add(new double[] { -1.0, -1.0, 1.0 });
            nodalParametricCoordinates.Add(new double[] { 1.0, -1.0, 1.0 });
            nodalParametricCoordinates.Add(new double[] { 1.0, 1.0, 1.0 });
            nodalParametricCoordinates.Add(new double[] { -1.0, 1.0, 1.0 });
            nodalParametricCoordinates.Add(new double[] { 0.0, -1.0, -1.0 });

            nodalParametricCoordinates.Add(new double[] { 1.0, 0.0, -1.0 });
            nodalParametricCoordinates.Add(new double[] { 0.0, 1.0, -1.0 });
            nodalParametricCoordinates.Add(new double[] { -1.0, 0.0, -1.0 });
            nodalParametricCoordinates.Add(new double[] { 0.0, -1.0, 1.0 });
            nodalParametricCoordinates.Add(new double[] { 1.0, 0.0, 1.0 });
            nodalParametricCoordinates.Add(new double[] { 0.0, 1.0, 1.0 });
            nodalParametricCoordinates.Add(new double[] { -1.0, 0.0, 1.0 });
            nodalParametricCoordinates.Add(new double[] { -1.0, -1.0, 0.0 });
            nodalParametricCoordinates.Add(new double[] { 1.0, -1.0, 0.0 });

            nodalParametricCoordinates.Add(new double[] { 1.0, 1.0, 0.0 });
            nodalParametricCoordinates.Add(new double[] { -1.0, 1.0, 0.0 });
            nodalParametricCoordinates.Add(new double[] { 0.0, 0.0, -1.0 });
            nodalParametricCoordinates.Add(new double[] { 0.0, 0.0, 1.0 });
            nodalParametricCoordinates.Add(new double[] { 0.0, -1.0, 0.0 });
            nodalParametricCoordinates.Add(new double[] { 0.0, 1.0, 0.0 });
            nodalParametricCoordinates.Add(new double[] { -1.0, 0.0, 0.0 });
            nodalParametricCoordinates.Add(new double[] { 1.0, 0.0, 0.0 });
            nodalParametricCoordinates.Add(new double[] { 0.0, 0.0, 0.0 });

            for (int i = 0; i < nodesCount; i++)
            {
                double[] nodalParamCoord = nodalParametricCoordinates[i];
                Dictionary<string, double[]> localdN = CalculateShapeFunctionsLocalDerivatives(nodalParamCoord);
                double[,] J = CalculateJacobian(localdN);
                double[,] invJ = CalculateInverseJacobian(J).Item1;
                Dictionary<int, double[]> globaldN = CalculateShapeFunctionsGlobalDerivatives(localdN, invJ);
                double[,] B = CalculateBMatrix(globaldN);
                double[] strainVector = CalculateStrainsVector(B);
                double[] stressVector = CalculateStressVector(E, strainVector);
                NodalStessVector.Add(stressVector);
            }
            return NodalStessVector;
        }
        public List<double[]> GetStrainFromElementsNodes()
        {
            List<double[]> NodalDeformationsVector = new List<double[]>();
            int nodesCount = 27;
            List<double[]> nodalParametricCoordinates = new List<double[]>();
            nodalParametricCoordinates.Add(new double[] { -1.0, -1.0, -1.0 });
            nodalParametricCoordinates.Add(new double[] { 1.0, -1.0, -1.0 });
            nodalParametricCoordinates.Add(new double[] { 1.0, 1.0, -1.0 });
            nodalParametricCoordinates.Add(new double[] { -1.0, 1.0, -1.0 });
            nodalParametricCoordinates.Add(new double[] { -1.0, -1.0, 1.0 });
            nodalParametricCoordinates.Add(new double[] { 1.0, -1.0, 1.0 });
            nodalParametricCoordinates.Add(new double[] { 1.0, 1.0, 1.0 });
            nodalParametricCoordinates.Add(new double[] { -1.0, 1.0, 1.0 });
            nodalParametricCoordinates.Add(new double[] { 0.0, -1.0, -1.0 });

            nodalParametricCoordinates.Add(new double[] { 1.0, 0.0, -1.0 });
            nodalParametricCoordinates.Add(new double[] { 0.0, 1.0, -1.0 });
            nodalParametricCoordinates.Add(new double[] { -1.0, 0.0, -1.0 });
            nodalParametricCoordinates.Add(new double[] { 0.0, -1.0, 1.0 });
            nodalParametricCoordinates.Add(new double[] { 1.0, 0.0, 1.0 });
            nodalParametricCoordinates.Add(new double[] { 0.0, 1.0, 1.0 });
            nodalParametricCoordinates.Add(new double[] { -1.0, 0.0, 1.0 });
            nodalParametricCoordinates.Add(new double[] { -1.0, -1.0, 0.0 });
            nodalParametricCoordinates.Add(new double[] { 1.0, -1.0, 0.0 });

            nodalParametricCoordinates.Add(new double[] { 1.0, 1.0, 0.0 });
            nodalParametricCoordinates.Add(new double[] { -1.0, 1.0, 0.0 });
            nodalParametricCoordinates.Add(new double[] { 0.0, 0.0, -1.0 });
            nodalParametricCoordinates.Add(new double[] { 0.0, 0.0, 1.0 });
            nodalParametricCoordinates.Add(new double[] { 0.0, -1.0, 0.0 });
            nodalParametricCoordinates.Add(new double[] { 0.0, 1.0, 0.0 });
            nodalParametricCoordinates.Add(new double[] { -1.0, 0.0, 0.0 });
            nodalParametricCoordinates.Add(new double[] { 1.0, 0.0, 0.0 });
            nodalParametricCoordinates.Add(new double[] { 0.0, 0.0, 0.0 });

            for (int i = 0; i< nodesCount; i++)
            {
                double[] nodalParamCoord = nodalParametricCoordinates[i];
                //double[] gW = GaussPoints(i, j, k).Item2;
                Dictionary<string, double[]> localdN = CalculateShapeFunctionsLocalDerivatives(nodalParamCoord);
                double[,] J = CalculateJacobian(localdN);
                double[,] invJ = CalculateInverseJacobian(J).Item1;
                //double detJ = CalculateInverseJacobian(J).Item2;
                Dictionary<int, double[]> globaldN = CalculateShapeFunctionsGlobalDerivatives(localdN, invJ);
                double[,] B = CalculateBMatrix(globaldN);
                double[] strainVector = CalculateStrainsVector(B);
                NodalDeformationsVector.Add(strainVector);
            }
            return NodalDeformationsVector;
        }
        public List<double[]> GetStressFromElements(List<double[]> parametricCoordinatesVector)
        {
            List<double[]> StessVectorsList = new List<double[]>();
            double[,] E = CalculateStressStrainMatrix(Properties.YoungMod, Properties.PoissonRatio);
            int count = parametricCoordinatesVector.Count;
            for (int i = 0; i < count; i++)
            {
                double[] nodalParamCoord = parametricCoordinatesVector[i];
                Dictionary<string, double[]> localdN = CalculateShapeFunctionsLocalDerivatives(nodalParamCoord);
                double[,] J = CalculateJacobian(localdN);
                double[,] invJ = CalculateInverseJacobian(J).Item1;
                Dictionary<int, double[]> globaldN = CalculateShapeFunctionsGlobalDerivatives(localdN, invJ);
                double[,] B = CalculateBMatrix(globaldN);
                double[] strainVector = CalculateStrainsVector(B);
                double[] stressVector = CalculateStressVector(E, strainVector);
                StessVectorsList.Add(stressVector);
            }
            return StessVectorsList;
        }
        private double[] UpdateNodalCoordinates(double[] displacementVector)
        {
            double[] updatedCoor = new double[81];
            for (int i = 1; i <= 27; i++)
            {
                updatedCoor[3 * i - 3] = Nodes[i].XCoordinate + displacementVector[3 * i - 3];
                updatedCoor[3 * i - 2] = Nodes[i].YCoordinate + displacementVector[3 * i - 2];
                updatedCoor[3 * i - 1] = Nodes[i].ZCoordinate + displacementVector[3 * i - 1];
            }
            return updatedCoor;
        }

        private Dictionary<int, double> CalculateShapeFunctions(double ksi, double ihta, double zita)
        {
            Dictionary<int, double> shapeFunctions = new Dictionary<int, double>();
            double N1 = 1.0 / 8.0 * (Math.Pow(ksi, 2.0) - ksi) * (Math.Pow(ihta, 2.0) - ihta) * (Math.Pow(zita, 2.0) - zita); shapeFunctions.Add(1, N1);
            double N2 = 1.0 / 8.0 * (Math.Pow(ksi, 2.0) + ksi) * (Math.Pow(ihta, 2.0) - ihta) * (Math.Pow(zita, 2.0) - zita); shapeFunctions.Add(2, N2);
            double N3 = 1.0 / 8.0 * (Math.Pow(ksi, 2.0) + ksi) * (Math.Pow(ihta, 2.0) + ihta) * (Math.Pow(zita, 2.0) - zita); shapeFunctions.Add(3, N3);
            double N4 = 1.0 / 8.0 * (Math.Pow(ksi, 2.0) - ksi) * (Math.Pow(ihta, 2.0) + ihta) * (Math.Pow(zita, 2.0) - zita); shapeFunctions.Add(4, N4);
            double N5 = 1.0 / 8.0 * (Math.Pow(ksi, 2.0) - ksi) * (Math.Pow(ihta, 2.0) - ihta) * (Math.Pow(zita, 2.0) + zita); shapeFunctions.Add(5, N5);
            double N6 = 1.0 / 8.0 * (Math.Pow(ksi, 2.0) + ksi) * (Math.Pow(ihta, 2.0) - ihta) * (Math.Pow(zita, 2.0) + zita); shapeFunctions.Add(6, N6);
            double N7 = 1.0 / 8.0 * (Math.Pow(ksi, 2.0) + ksi) * (Math.Pow(ihta, 2.0) + ihta) * (Math.Pow(zita, 2.0) + zita); shapeFunctions.Add(7, N7);
            double N8 = 1.0 / 8.0 * (Math.Pow(ksi, 2.0) - ksi) * (Math.Pow(ihta, 2.0) + ihta) * (Math.Pow(zita, 2.0) + zita); shapeFunctions.Add(8, N8);
            double N9 = 1.0 / 4.0 * (1.0 - Math.Pow(ksi, 2.0)) * (Math.Pow(ihta, 2.0) - ihta) * (Math.Pow(zita, 2.0) - zita); shapeFunctions.Add(9, N9);
            double N10 = 1.0 / 4.0 * (Math.Pow(ksi, 2.0) + ksi) * (1.0 - Math.Pow(ihta, 2.0)) * (Math.Pow(zita, 2.0) - zita); shapeFunctions.Add(10, N10);
            double N11 = 1.0 / 4.0 * (1.0 - Math.Pow(ksi, 2.0)) * (Math.Pow(ihta, 2.0) + ihta) * (Math.Pow(zita, 2.0) - zita); shapeFunctions.Add(11, N11);
            double N12 = 1.0 / 4.0 * (Math.Pow(ksi, 2.0) - ksi) * (1.0 - Math.Pow(ihta, 2.0)) * (Math.Pow(zita, 2.0) - zita); shapeFunctions.Add(12, N12);
            double N13 = 1.0 / 4.0 * (1.0 - Math.Pow(ksi, 2.0)) * (Math.Pow(ihta, 2.0) - ihta) * (Math.Pow(zita, 2.0) + zita); shapeFunctions.Add(13, N13);
            double N14 = 1.0 / 4.0 * (Math.Pow(ksi, 2.0) + ksi) * (1.0 - Math.Pow(ihta, 2.0)) * (Math.Pow(zita, 2.0) + zita); shapeFunctions.Add(14, N14);
            double N15 = 1.0 / 4.0 * (1.0 - Math.Pow(ksi, 2.0)) * (Math.Pow(ihta, 2.0) + ihta) * (Math.Pow(zita, 2.0) + zita); shapeFunctions.Add(15, N15);
            double N16 = 1.0 / 4.0 * (Math.Pow(ksi, 2.0) - ksi) * (1.0 - Math.Pow(ihta, 2.0)) * (Math.Pow(zita, 2.0) + zita); shapeFunctions.Add(16, N16);
            double N17 = 1.0 / 4.0 * (Math.Pow(ksi, 2.0) - ksi) * (Math.Pow(ihta, 2.0) - ihta) * (1.0 - Math.Pow(zita, 2.0)); shapeFunctions.Add(17, N17);
            double N18 = 1.0 / 4.0 * (Math.Pow(ksi, 2.0) + ksi) * (Math.Pow(ihta, 2.0) - ihta) * (1.0 - Math.Pow(zita, 2.0)); shapeFunctions.Add(18, N18);
            double N19 = 1.0 / 4.0 * (Math.Pow(ksi, 2.0) + ksi) * (Math.Pow(ihta, 2.0) + ihta) * (1.0 - Math.Pow(zita, 2.0)); shapeFunctions.Add(19, N19);
            double N20 = 1.0 / 4.0 * (Math.Pow(ksi, 2.0) - ksi) * (Math.Pow(ihta, 2.0) + ihta) * (1.0 - Math.Pow(zita, 2.0)); shapeFunctions.Add(20, N20);
            double N21 = 1.0 / 2.0 * (1.0 - Math.Pow(ksi, 2.0)) * (1.0 - Math.Pow(ihta, 2.0)) * (Math.Pow(zita, 2.0) - zita); shapeFunctions.Add(21, N21);
            double N22 = 1.0 / 2.0 * (1.0 - Math.Pow(ksi, 2.0)) * (1.0 - Math.Pow(ihta, 2.0)) * (Math.Pow(zita, 2.0) + zita); shapeFunctions.Add(22, N22);
            double N23 = 1.0 / 2.0 * (1.0 - Math.Pow(ksi, 2.0)) * (Math.Pow(ihta, 2.0) - ihta) * (1.0 - Math.Pow(zita, 2.0)); shapeFunctions.Add(23, N23);
            double N24 = 1.0 / 2.0 * (1.0 - Math.Pow(ksi, 2.0)) * (Math.Pow(ihta, 2.0) + ihta) * (1.0 - Math.Pow(zita, 2.0)); shapeFunctions.Add(24, N24);
            double N25 = 1.0 / 2.0 * (Math.Pow(ksi, 2.0) - ksi) * (1.0 - Math.Pow(ihta, 2.0)) * (1.0 - Math.Pow(zita, 2.0)); shapeFunctions.Add(25, N25);
            double N26 = 1.0 / 2.0 * (Math.Pow(ksi, 2.0) + ksi) * (1.0 - Math.Pow(ihta, 2.0)) * (1.0 - Math.Pow(zita, 2.0)); shapeFunctions.Add(26, N26);
            double N27 = (1.0 - Math.Pow(ksi, 2.0)) * (1.0 - Math.Pow(ihta, 2.0)) * (1.0 - Math.Pow(zita, 2.0)); shapeFunctions.Add(27, N27);

            return shapeFunctions;
        }
        private double[,] CalculateShapeFunctionMatrix(double ksi, double ihta, double zita)
        {
            double N1 = 1.0 / 8.0 * (Math.Pow(ksi,2.0) - ksi) * (Math.Pow(ihta, 2.0) - ihta) * (Math.Pow(zita, 2.0) - zita);
            double N2 = 1.0 / 8.0 * (Math.Pow(ksi, 2.0) + ksi) * (Math.Pow(ihta, 2.0) - ihta) * (Math.Pow(zita, 2.0) - zita);
            double N3 = 1.0 / 8.0 * (Math.Pow(ksi, 2.0) + ksi) * (Math.Pow(ihta, 2.0) + ihta) * (Math.Pow(zita, 2.0) - zita);
            double N4 = 1.0 / 8.0 * (Math.Pow(ksi, 2.0) - ksi) * (Math.Pow(ihta, 2.0) + ihta) * (Math.Pow(zita, 2.0) - zita);
            double N5 = 1.0 / 8.0 * (Math.Pow(ksi, 2.0) - ksi) * (Math.Pow(ihta, 2.0) - ihta) * (Math.Pow(zita, 2.0) + zita);
            double N6 = 1.0 / 8.0 * (Math.Pow(ksi, 2.0) + ksi) * (Math.Pow(ihta, 2.0) - ihta) * (Math.Pow(zita, 2.0) + zita);
            double N7 = 1.0 / 8.0 * (Math.Pow(ksi, 2.0) + ksi) * (Math.Pow(ihta, 2.0) + ihta) * (Math.Pow(zita, 2.0) + zita);
            double N8 = 1.0 / 8.0 * (Math.Pow(ksi, 2.0) - ksi) * (Math.Pow(ihta, 2.0) + ihta) * (Math.Pow(zita, 2.0) + zita);

            double N9 = 1.0 / 4.0 * (1.0 - Math.Pow(ksi, 2.0)) * (Math.Pow(ihta, 2.0) - ihta) * (Math.Pow(zita, 2.0) - zita);
            double N10 = 1.0 / 4.0 * (Math.Pow(ksi, 2.0) + ksi) * (1.0 - Math.Pow(ihta, 2.0)) * (Math.Pow(zita, 2.0) - zita);
            double N11 = 1.0 / 4.0 * (1.0 - Math.Pow(ksi, 2.0)) * (Math.Pow(ihta, 2.0) + ihta) * (Math.Pow(zita, 2.0) - zita);
            double N12 = 1.0 / 4.0 * (Math.Pow(ksi, 2.0) - ksi) * (1.0 - Math.Pow(ihta, 2.0)) * (Math.Pow(zita, 2.0) - zita);
            double N13 = 1.0 / 4.0 * (1.0 - Math.Pow(ksi, 2.0)) * (Math.Pow(ihta, 2.0) - ihta) * (Math.Pow(zita, 2.0) + zita);
            double N14 = 1.0 / 4.0 * (Math.Pow(ksi, 2.0) + ksi) * (1.0 - Math.Pow(ihta, 2.0)) * (Math.Pow(zita, 2.0) + zita);
            double N15 = 1.0 / 4.0 * (1.0 - Math.Pow(ksi, 2.0)) * (Math.Pow(ihta, 2.0) + ihta) * (Math.Pow(zita, 2.0) + zita);
            double N16 = 1.0 / 4.0 * (Math.Pow(ksi, 2.0) - ksi) * (1.0 - Math.Pow(ihta, 2.0)) * (Math.Pow(zita, 2.0) + zita);
            double N17 = 1.0 / 4.0 * (Math.Pow(ksi, 2.0) - ksi) * (Math.Pow(ihta, 2.0) - ihta) * (1.0 - Math.Pow(zita, 2.0));
            double N18 = 1.0 / 4.0 * (Math.Pow(ksi, 2.0) + ksi) * (Math.Pow(ihta, 2.0) - ihta) * (1.0 - Math.Pow(zita, 2.0));
            double N19 = 1.0 / 4.0 * (Math.Pow(ksi, 2.0) + ksi) * (Math.Pow(ihta, 2.0) + ihta) * (1.0 - Math.Pow(zita, 2.0));
            double N20 = 1.0 / 4.0 * (Math.Pow(ksi, 2.0) - ksi) * (Math.Pow(ihta, 2.0) + ihta) * (1.0 - Math.Pow(zita, 2.0));

            double N21 = 1.0 / 2.0 * (1.0 - Math.Pow(ksi, 2.0)) * (1.0 - Math.Pow(ihta, 2.0)) * (Math.Pow(zita, 2.0) - zita);
            double N22 = 1.0 / 2.0 * (1.0 - Math.Pow(ksi, 2.0)) * (1.0 - Math.Pow(ihta, 2.0)) * (Math.Pow(zita, 2.0) + zita);
            double N23 = 1.0 / 2.0 * (1.0 - Math.Pow(ksi, 2.0)) * (Math.Pow(ihta, 2.0) - ihta) * (1.0 - Math.Pow(zita, 2.0));
            double N24 = 1.0 / 2.0 * (1.0 - Math.Pow(ksi, 2.0)) * (Math.Pow(ihta, 2.0) + ihta) * (1.0 - Math.Pow(zita, 2.0));
            double N25 = 1.0 / 2.0 * (Math.Pow(ksi, 2.0) - ksi) * (1.0 - Math.Pow(ihta, 2.0)) * (1.0 - Math.Pow(zita, 2.0));
            double N26 = 1.0 / 2.0 * (Math.Pow(ksi, 2.0) + ksi) * (1.0 - Math.Pow(ihta, 2.0)) * (1.0 - Math.Pow(zita, 2.0));
            double N27 = (1.0 - Math.Pow(ksi, 2.0)) * (1.0 - Math.Pow(ihta, 2.0)) * (1.0 - Math.Pow(zita, 2.0));

            double[,] shapeFunctionsMat = new double[,] {
                {N1, 0.0, 0.0, N2, 0.0, 0.0, N3, 0.0, 0.0, N4, 0.0, 0.0, N5, 0.0, 0.0, N6, 0.0, 0.0, N7, 0.0, 0.0, N8, 0.0, 0.0, N9, 0.0, 0.0, N10, 0.0, 0.0, N11, 0.0, 0.0, N12, 0.0, 0.0, N13, 0.0, 0.0, N14, 0.0, 0.0, N15, 0.0, 0.0, N16, 0.0, 0.0,N17, 0.0, 0.0, N18, 0.0, 0.0, N19, 0.0, 0.0, N20, 0.0, 0.0, N21, 0.0, 0.0, N22, 0.0, 0.0, N23, 0.0, 0.0, N24, 0.0, 0.0, N25, 0.0, 0.0, N26, 0.0, 0.0, N27, 0.0, 0.0 },
                {0.0, N1, 0.0, 0.0, N2, 0.0, 0.0, N3, 0.0, 0.0, N4, 0.0, 0.0, N5, 0.0, 0.0, N6, 0.0, 0.0, N7, 0.0, 0.0, N8, 0.0,0.0, N9, 0.0, 0.0, N10, 0.0, 0.0, N11, 0.0, 0.0, N12, 0.0, 0.0, N13, 0.0, 0.0, N14, 0.0, 0.0, N15, 0.0, 0.0, N16, 0.0, 0.0,N17, 0.0, 0.0, N18, 0.0, 0.0, N19, 0.0, 0.0, N20, 0.0, 0.0, N21, 0.0, 0.0, N22, 0.0, 0.0, N23, 0.0, 0.0, N24, 0.0, 0.0, N25, 0.0, 0.0, N26, 0.0, 0.0, N27, 0.0 },
                {0.0, 0.0, N1, 0.0, 0.0, N2, 0.0, 0.0, N3, 0.0, 0.0, N4, 0.0, 0.0, N5, 0.0, 0.0, N6, 0.0, 0.0, N7, 0.0, 0.0, N8,0.0,0.0, N9, 0.0, 0.0, N10, 0.0, 0.0, N11, 0.0, 0.0, N12, 0.0, 0.0, N13, 0.0, 0.0, N14, 0.0, 0.0, N15, 0.0, 0.0, N16, 0.0, 0.0,N17, 0.0, 0.0, N18, 0.0, 0.0, N19, 0.0, 0.0, N20, 0.0, 0.0, N21, 0.0, 0.0, N22, 0.0, 0.0, N23, 0.0, 0.0, N24, 0.0, 0.0, N25, 0.0, 0.0, N26, 0.0, 0.0, N27 }
            };
            return shapeFunctionsMat;
        }

        private Dictionary<string, double[]> CalculateShapeFunctionsLocalDerivatives(double[] naturalCoordinates)
        {
            double ksi = naturalCoordinates[0];
            double ihta = naturalCoordinates[1];
            double zita = naturalCoordinates[2];

            double[] dN_ksi = new double[]
            {
                (1.0 / 8.0 * (2 * ksi - 1.0) * (Math.Pow(ihta, 2.0) - ihta) * (Math.Pow(zita, 2.0) - zita)),
                (1.0 / 8.0 * (2 * ksi + 1.0) * (Math.Pow(ihta, 2.0) - ihta) * (Math.Pow(zita, 2.0) - zita)),
                (1.0 / 8.0 * (2 * ksi + 1.0) * (Math.Pow(ihta, 2.0) + ihta) * (Math.Pow(zita, 2.0) - zita)),
                (1.0 / 8.0 * (2 * ksi - 1.0) * (Math.Pow(ihta, 2.0) + ihta) * (Math.Pow(zita, 2.0) - zita)),
                (1.0 / 8.0 * (2 * ksi - 1.0) * (Math.Pow(ihta, 2.0) - ihta) * (Math.Pow(zita, 2.0) + zita)),
                (1.0 / 8.0 * (2 * ksi + 1.0) * (Math.Pow(ihta, 2.0) - ihta) * (Math.Pow(zita, 2.0) + zita)),
                (1.0 / 8.0 * (2 * ksi + 1.0) * (Math.Pow(ihta, 2.0) + ihta) * (Math.Pow(zita, 2.0) + zita)),
                (1.0 / 8.0 * (2 * ksi - 1.0) * (Math.Pow(ihta, 2.0) + ihta) * (Math.Pow(zita, 2.0) + zita)),
                (-1.0 / 2.0 * ksi * (Math.Pow(ihta, 2.0) - ihta) * (Math.Pow(zita, 2.0) - zita)),
                (1.0 / 4.0 * (2 * ksi + 1.0) * (1.0 - Math.Pow(ihta, 2.0)) * (Math.Pow(zita, 2.0) - zita)),
                (-1.0 / 2.0 * ksi * (Math.Pow(ihta, 2.0) + ihta) * (Math.Pow(zita, 2.0) - zita)),
                (1.0 / 4.0 * (2 * ksi - 1.0) * (1.0 - Math.Pow(ihta, 2.0)) * (Math.Pow(zita, 2.0) - zita)),
                (-1.0 / 2.0 * ksi * (Math.Pow(ihta, 2.0) - ihta) * (Math.Pow(zita, 2.0) + zita)),
                (1.0 / 4.0 * (2 * ksi + 1.0) * (1.0 - Math.Pow(ihta, 2.0)) * (Math.Pow(zita, 2.0) + zita)),
                (-1.0 / 2.0 * ksi * (Math.Pow(ihta, 2.0) + ihta) * (Math.Pow(zita, 2.0) + zita)),
                (1.0 / 4.0 * (2 * ksi - 1.0) * (1.0 - Math.Pow(ihta, 2.0)) * (Math.Pow(zita, 2.0) + zita)),
                (1.0 / 4.0 * (2 * ksi - 1.0) * (Math.Pow(ihta, 2.0) - ihta) * (1.0 - Math.Pow(zita, 2.0))),
                (1.0 / 4.0 * (2 * ksi + 1.0) * (Math.Pow(ihta, 2.0) - ihta) * (1.0 - Math.Pow(zita, 2.0))),
                (1.0 / 4.0 * (2 * ksi + 1.0) * (Math.Pow(ihta, 2.0) + ihta) * (1.0 - Math.Pow(zita, 2.0))),
                (1.0 / 4.0 * (2 * ksi - 1.0) * (Math.Pow(ihta, 2.0) + ihta) * (1.0 - Math.Pow(zita, 2.0))),
                (-ksi * (1.0 - Math.Pow(ihta, 2.0)) * (Math.Pow(zita, 2.0) - zita)),
                (-ksi * (1.0 - Math.Pow(ihta, 2.0)) * (Math.Pow(zita, 2.0) + zita)),
                (-ksi * (Math.Pow(ihta, 2.0) - ihta) * (1.0 - Math.Pow(zita, 2.0))),
                (- ksi * (Math.Pow(ihta, 2.0) + ihta) * (1.0 - Math.Pow(zita, 2.0))),
                (1.0 / 2.0 * (2 * ksi - 1.0) * (1.0 - Math.Pow(ihta, 2.0)) * (1.0 - Math.Pow(zita, 2.0))),
                (1.0 / 2.0 * (2 * ksi + 1.0) * (1.0 - Math.Pow(ihta, 2.0)) * (1.0 - Math.Pow(zita, 2.0))),
                (-2.0 * ksi * (1.0 - Math.Pow(ihta, 2.0)) * (1.0 - Math.Pow(zita, 2.0)))
            };

            double[] dN_ihta = new double[]
            {
                (1.0 / 8.0 * (Math.Pow(ksi, 2.0) - ksi) * (2 * ihta - 1.0) * (Math.Pow(zita, 2.0) - zita)),
                (1.0 / 8.0 * (Math.Pow(ksi, 2.0) + ksi) * (2 * ihta - 1.0) * (Math.Pow(zita, 2.0) - zita)),
                (1.0 / 8.0 * (Math.Pow(ksi, 2.0) + ksi) * (2 * ihta + 1.0) * (Math.Pow(zita, 2.0) - zita)),
                (1.0 / 8.0 * (Math.Pow(ksi, 2.0) - ksi) * (2 * ihta + 1.0) * (Math.Pow(zita, 2.0) - zita)),
                (1.0 / 8.0 * (Math.Pow(ksi, 2.0) - ksi) * (2 * ihta - 1.0) * (Math.Pow(zita, 2.0) + zita)),
                (1.0 / 8.0 * (Math.Pow(ksi, 2.0) + ksi) * (2 * ihta - 1.0) * (Math.Pow(zita, 2.0) + zita)),
                (1.0 / 8.0 * (Math.Pow(ksi, 2.0) + ksi) * (2 * ihta + 1.0) * (Math.Pow(zita, 2.0) + zita)),
                (1.0 / 8.0 * (Math.Pow(ksi, 2.0) - ksi) * (2 * ihta + 1.0) * (Math.Pow(zita, 2.0) + zita)),
                (1.0 / 4.0 * (1.0 - Math.Pow(ksi, 2.0)) * (2 * ihta - 1.0) * (Math.Pow(zita, 2.0) - zita)),
                (-1.0 / 2.0 * (Math.Pow(ksi, 2.0) + ksi) * ihta * (Math.Pow(zita, 2.0) - zita)),
                (1.0 / 4.0 * (1.0 - Math.Pow(ksi, 2.0)) * (2 * ihta + 1.0) * (Math.Pow(zita, 2.0) - zita)),
                (-1.0 / 2.0 * (Math.Pow(ksi, 2.0) - ksi) * ihta * (Math.Pow(zita, 2.0) - zita)),
                (1.0 / 4.0 * (1.0 - Math.Pow(ksi, 2.0)) * (2 * ihta - 1.0) * (Math.Pow(zita, 2.0) + zita)),
                (-1.0 / 2.0 * (Math.Pow(ksi, 2.0) + ksi) * ihta * (Math.Pow(zita, 2.0) + zita)),
                (1.0 / 4.0 * (1.0 - Math.Pow(ksi, 2.0)) * (2 * ihta + 1.0) * (Math.Pow(zita, 2.0) + zita)),
                (-1.0 / 2.0 * (Math.Pow(ksi, 2.0) - ksi) * ihta * (Math.Pow(zita, 2.0) + zita)),
                (1.0 / 4.0 * (Math.Pow(ksi, 2.0) - ksi) * (2 * ihta - 1.0) * (1.0 - Math.Pow(zita, 2.0))),
                (1.0 / 4.0 * (Math.Pow(ksi, 2.0) + ksi) * (2 * ihta - 1.0) * (1.0 - Math.Pow(zita, 2.0))),
                (1.0 / 4.0 * (Math.Pow(ksi, 2.0) + ksi) * (2 * ihta + 1.0) * (1.0 - Math.Pow(zita, 2.0))),
                (1.0 / 4.0 * (Math.Pow(ksi, 2.0) - ksi) * (2 * ihta + 1.0) * (1.0 - Math.Pow(zita, 2.0))),
                (-(1.0 - Math.Pow(ksi, 2.0)) * ihta * (Math.Pow(zita, 2.0) - zita)),
                (-(1.0 - Math.Pow(ksi, 2.0)) * ihta * (Math.Pow(zita, 2.0) + zita)),
                (1.0 / 2.0 * (1.0 - Math.Pow(ksi, 2.0)) * (2 * ihta - 1.0) * (1.0 - Math.Pow(zita, 2.0))),
                (1.0 / 2.0 * (1.0 - Math.Pow(ksi, 2.0)) * (2 * ihta + 1.0) * (1.0 - Math.Pow(zita, 2.0))),
                (- (Math.Pow(ksi, 2.0) - ksi) * ihta * (1.0 - Math.Pow(zita, 2.0))),
                (-(Math.Pow(ksi, 2.0) + ksi) * ihta * (1.0 - Math.Pow(zita, 2.0))),
                (-2.0 * (1.0 - Math.Pow(ksi, 2.0)) * ihta * (1.0 - Math.Pow(zita, 2.0)))
            };

            double[] dN_zita = new double[]
            {
                (1.0 / 8.0 * (Math.Pow(ksi, 2.0) - ksi) * (Math.Pow(ihta, 2.0) - ihta) * (2 * zita - 1.0)),
                (1.0 / 8.0 * (Math.Pow(ksi, 2.0) + ksi) * (Math.Pow(ihta, 2.0) - ihta) * (2 * zita - 1.0)),
                (1.0 / 8.0 * (Math.Pow(ksi, 2.0) + ksi) * (Math.Pow(ihta, 2.0) + ihta) * (2 * zita - 1.0)),
                (1.0 / 8.0 * (Math.Pow(ksi, 2.0) - ksi) * (Math.Pow(ihta, 2.0) + ihta) * (2 * zita - 1.0)),
                (1.0 / 8.0 * (Math.Pow(ksi, 2.0) - ksi) * (Math.Pow(ihta, 2.0) - ihta) * (2 * zita + 1.0)),
                (1.0 / 8.0 * (Math.Pow(ksi, 2.0) + ksi) * (Math.Pow(ihta, 2.0) - ihta) * (2 * zita + 1.0)),
                (1.0 / 8.0 * (Math.Pow(ksi, 2.0) + ksi) * (Math.Pow(ihta, 2.0) + ihta) * (2 * zita + 1.0)),
                (1.0 / 8.0 * (Math.Pow(ksi, 2.0) - ksi) * (Math.Pow(ihta, 2.0) + ihta) * (2 * zita + 1.0)),
                (1.0 / 4.0 * (1.0 - Math.Pow(ksi, 2.0)) * (Math.Pow(ihta, 2.0) - ihta) * (2 * zita - 1.0)),
                (1.0 / 4.0 * (Math.Pow(ksi, 2.0) + ksi) * (1.0 - Math.Pow(ihta, 2.0)) * (2 * zita - 1.0)),
                (1.0 / 4.0 * (1.0 - Math.Pow(ksi, 2.0)) * (Math.Pow(ihta, 2.0) + ihta) * (2 * zita - 1.0)),
                (1.0 / 4.0 * (Math.Pow(ksi, 2.0) - ksi) * (1.0 - Math.Pow(ihta, 2.0)) * (2 * zita - 1.0)),
                (1.0 / 4.0 * (1.0 - Math.Pow(ksi, 2.0)) * (Math.Pow(ihta, 2.0) - ihta) * (2 * zita + 1.0)),
                (1.0 / 4.0 * (Math.Pow(ksi, 2.0) + ksi) * (1.0 - Math.Pow(ihta, 2.0)) * (2 * zita + 1.0)),
                (1.0 / 4.0 * (1.0 - Math.Pow(ksi, 2.0)) * (Math.Pow(ihta, 2.0) + ihta) * (2 * zita + 1.0)),
                (1.0 / 4.0 * (Math.Pow(ksi, 2.0) - ksi) * (1.0 - Math.Pow(ihta, 2.0)) * (2 * zita + 1.0)),
                (-1.0 / 2.0 * (Math.Pow(ksi, 2.0) - ksi) * (Math.Pow(ihta, 2.0) - ihta) * zita),
                (-1.0 / 2.0 * (Math.Pow(ksi, 2.0) + ksi) * (Math.Pow(ihta, 2.0) - ihta) * zita),
                (-1.0 / 2.0 * (Math.Pow(ksi, 2.0) + ksi) * (Math.Pow(ihta, 2.0) + ihta) * zita),
                (-1.0 / 2.0 * (Math.Pow(ksi, 2.0) - ksi) * (Math.Pow(ihta, 2.0) + ihta) * zita),
                (1.0 / 2.0 * (1.0 - Math.Pow(ksi, 2.0)) * (1.0 - Math.Pow(ihta, 2.0)) * (2 * zita - 1.0)),
                (1.0 / 2.0 * (1.0 - Math.Pow(ksi, 2.0)) * (1.0 - Math.Pow(ihta, 2.0)) * (2 * zita + 1.0)),
                (-(1.0 - Math.Pow(ksi, 2.0)) * (Math.Pow(ihta, 2.0) - ihta) * zita),
                (- (1.0 - Math.Pow(ksi, 2.0)) * (Math.Pow(ihta, 2.0) + ihta) * zita),
                (- (Math.Pow(ksi, 2.0) - ksi) * (1.0 - Math.Pow(ihta, 2.0)) * zita),
                (-(Math.Pow(ksi, 2.0) + ksi) * (1.0 - Math.Pow(ihta, 2.0)) * zita),
                (-2.0 * (1.0 - Math.Pow(ksi, 2.0)) * (1.0 - Math.Pow(ihta, 2.0)) * zita)
            };

            Dictionary<string, double[]> dN = new Dictionary<string, double[]>();
            dN.Add("ksi", dN_ksi);
            dN.Add("ihta", dN_ihta);
            dN.Add("mhi", dN_zita);
            return dN;
        }

        private double[,] CalculateJacobian(Dictionary<string, double[]> dN)
        {
            double[,] jacobianMatrix = new double[3, 3];
            //DisplacementVector = new double[24];
            double[] xUpdated = UpdateNodalCoordinates(DisplacementVector);

            int k = 0;
            for (int i = 0; i < 27; i++)
            {
                jacobianMatrix[0, 0] = jacobianMatrix[0, 0] + xUpdated[k] * dN["ksi"][i];
                k += 3;
            }
            k = 1;
            for (int i = 0; i < 27; i++)
            {
                jacobianMatrix[0, 1] = jacobianMatrix[0, 1] + xUpdated[k] * dN["ksi"][i];
                k += 3;
            }
            k = 2;
            for (int i = 0; i < 27; i++)
            {
                jacobianMatrix[0, 2] = jacobianMatrix[0, 2] + xUpdated[k] * dN["ksi"][i];
                k += 3;
            }

            k = 0;
            for (int i = 0; i < 27; i++)
            {
                jacobianMatrix[1, 0] = jacobianMatrix[1, 0] + xUpdated[k] * dN["ihta"][i];
                k += 3;
            }
            k = 1;
            for (int i = 0; i < 27; i++)
            {
                jacobianMatrix[1, 1] = jacobianMatrix[1, 1] + xUpdated[k] * dN["ihta"][i];
                k += 3;
            }
            k = 2;
            for (int i = 0; i < 27; i++)
            {
                jacobianMatrix[1, 2] = jacobianMatrix[1, 2] + xUpdated[k] * dN["ihta"][i];
                k += 3;
            }

            k = 0;
            for (int i = 0; i < 27; i++)
            {
                jacobianMatrix[2, 0] = jacobianMatrix[2, 0] + xUpdated[k] * dN["mhi"][i];
                k += 3;
            }
            k = 1;
            for (int i = 0; i < 27; i++)
            {
                jacobianMatrix[2, 1] = jacobianMatrix[2, 1] + xUpdated[k] * dN["mhi"][i];
                k += 3;
            }
            k = 2;
            for (int i = 0; i < 27; i++)
            {
                jacobianMatrix[2, 2] = jacobianMatrix[2, 2] + xUpdated[k] * dN["mhi"][i];
                k += 3;
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

        private Dictionary<int, double[]> CalculateShapeFunctionsGlobalDerivatives(Dictionary<string, double[]> dN, double[,] Jinv)
        {
            Dictionary<int, double[]> dNg = new Dictionary<int, double[]>();

            for (int i = 0; i < 27; i++)
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

        private double[,] CalculateBMatrix(Dictionary<int, double[]> dNglobal)
        {
            double[,] Bmatrix = new double[6, 81];

            for (int i = 0; i < 27; i++)
            {
                Bmatrix[0, i * 3] = dNglobal[i][0];
                Bmatrix[1, i * 3 + 1] = dNglobal[i][1];
                Bmatrix[2, i * 3 + 2] = dNglobal[i][2];
                Bmatrix[3, i * 3] = dNglobal[i][1];
                Bmatrix[3, i * 3 + 1] = dNglobal[i][0];
                Bmatrix[4, i * 3 + 1] = dNglobal[i][2];
                Bmatrix[4, i * 3 + 2] = dNglobal[i][1];
                Bmatrix[5, i * 3] = dNglobal[i][2];
                Bmatrix[5, i * 3 + 2] = dNglobal[i][0];
            }
            //x
            //y
            //z
            //xy
            //yz
            //xz
            return Bmatrix;
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

        private Tuple<double[], double[]> GaussPoints(int i, int j, int k)
        {
            double[] gaussPoints = new double[] { -0.77459, 0.0, 0.77459 };
            double[] gaussWeights = new double[] { 0.55555, 0.88888, 0.55555 };

            double[] vectorWithPoints = new double[] { gaussPoints[i], gaussPoints[j], gaussPoints[k] };
            double[] vectorWithWeights = new double[] { gaussWeights[i], gaussWeights[j], gaussWeights[k] };
            return new Tuple<double[], double[]>(vectorWithPoints, vectorWithWeights);
        }

        public double[,] CreateGlobalStiffnessMatrix()
        {
            double[,] K = new double[81, 81];
            double[,] E = CalculateStressStrainMatrix(Properties.YoungMod, Properties.PoissonRatio);

            for (int i = 0; i <= 2; i++)
            {
                for (int j = 0; j <= 2; j++)
                {
                    for (int k = 0; k <= 2; k++)
                    {
                        double[] gP = GaussPoints(i, j, k).Item1;
                        double[] gW = GaussPoints(i, j, k).Item2;
                        Dictionary<string, double[]> localdN = CalculateShapeFunctionsLocalDerivatives(gP);
                        double[,] J = CalculateJacobian(localdN);
                        double[,] invJ = CalculateInverseJacobian(J).Item1;
                        double detJ = CalculateInverseJacobian(J).Item2;
                        Dictionary<int, double[]> globaldN = CalculateShapeFunctionsGlobalDerivatives(localdN, invJ);
                        double[,] B = CalculateBMatrix(globaldN);
                        K = MatrixOperations.MatrixAddition(K, MatrixOperations.ScalarMatrixProductNew(detJ * gW[0] * gW[1] * gW[2],
                            MatrixOperations.MatrixProduct(MatrixOperations.Transpose(B), MatrixOperations.MatrixProduct(E, B))));

                    }
                }
            }
            return K;
        }

        public double[,] CreateMassMatrix()
        {
            return new double[81, 81];
        }

        public double[,] CreateDampingMatrix()
        {
            throw new Exception("Not implemented");
        }

        private double[] CalculateStressVector(double[,] E, double[] strain)
        {
            double[] stressVector = VectorOperations.MatrixVectorProduct(E, strain);
            return stressVector;
        }

        public double[] CreateInternalGlobalForcesVector()
        {
            double[] fInt = new double[81];
            double[,] E = CalculateStressStrainMatrix(Properties.YoungMod, Properties.PoissonRatio);

            for (int i = 0; i <= 2; i++)
            {
                for (int j = 0; j <= 2; j++)
                {
                    for (int k = 0; k <= 2; k++)
                    {
                        double[] gP = GaussPoints(i, j, k).Item1;
                        double[] gW = GaussPoints(i, j, k).Item2;
                        Dictionary<string, double[]> localdN = CalculateShapeFunctionsLocalDerivatives(gP);
                        double[,] J = CalculateJacobian(localdN);
                        double[,] invJ = CalculateInverseJacobian(J).Item1;
                        double detJ = CalculateInverseJacobian(J).Item2;
                        Dictionary<int, double[]> globaldN = CalculateShapeFunctionsGlobalDerivatives(localdN, invJ);
                        double[,] B = CalculateBMatrix(globaldN);
                        double[] strainVector = CalculateStrainsVector(B);
                        double[] stressVector = CalculateStressVector(E, strainVector);
                        fInt = VectorOperations.VectorVectorAddition(fInt, VectorOperations.VectorScalarProductNew(
                            VectorOperations.MatrixVectorProduct(MatrixOperations.Transpose(B), stressVector), detJ * gW[0] * gW[1] * gW[2]));
                    }

                }
            }

            //double[,] Kstiff = CreateGlobalStiffnessMatrix();
            //double[] uDisp = DisplacementVector;
            //double[] fInt = VectorOperations.MatrixVectorProduct(Kstiff, uDisp);
            //fInt = VectorOperations.VectorScalarProductNew(fInt, 1.0);
            return fInt;
        }
    }
}