using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GFEC
{
    class ContactStS2D : IElement
    {
        public Dictionary<int, INode> Nodes { get; }
        public IElementProperties Properties { get; set; }
        public Dictionary<int, bool[]> ElementFreedomSignature { get; } = new Dictionary<int, bool[]>();
        public List<int> ElementFreedomList { get; set; }
        public double[] DisplacementVector { get; set; }
        public double[] AccelerationVector { get; set; }
        private double PenaltyFactor { get; set; }
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
        public ContactStS2D(IElementProperties properties, Dictionary<int, INode> nodes)
        {
            Properties = properties;
            this.Nodes = nodes;
            int amountOfNodes = Properties.SlaveSegmentPolynomialDegree + Properties.MasterSegmentPolynomialDegree + 2;
            for (int i = 1; i<= amountOfNodes; i++)
            {
                ElementFreedomSignature[i] = new bool[] { true, true, false, false, false, false };

            }
            DisplacementVector = new double[2* amountOfNodes];
            PenaltyFactor = properties.YoungMod * properties.PenaltyFactorRatio;
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
        public Dictionary<int, INode> NodesAtFinalState()
        {
            Dictionary<int, INode> finalNodes = new Dictionary<int, INode>();
            int amountOfNodes = Properties.SlaveSegmentPolynomialDegree + Properties.MasterSegmentPolynomialDegree + 2;
            int dofCount = new int();
            for (int i = 1; i<=amountOfNodes; i++)
            {
                finalNodes[i] = new Node(Nodes[i].XCoordinate + DisplacementVector[dofCount], Nodes[i].YCoordinate + DisplacementVector[dofCount + 1]);
                dofCount += 2;
            }
            return finalNodes;
        }
        private double[] NodalXUpdated()
        {
            int amountOfNodes = Properties.SlaveSegmentPolynomialDegree + Properties.MasterSegmentPolynomialDegree + 2;
            double[] x = new double[2 * amountOfNodes];
            for (int i = 1; i <= amountOfNodes; i++)
            {
                int count = (i - 1) * 2;
                x[count] = Nodes[i].XCoordinate + DisplacementVector[count];
                x[count + 1] = Nodes[i].YCoordinate + DisplacementVector[count + 1];
            }
            return x;
        }
        public List<double[]> GetStressVector()
        {
            List<double[]> l = new List<double[]>();
            l.Add(new double[] { 0.0, 0.0, 0.0 });
            return l;
        }
        public List<double[]> GetStrainVector()
        {
            List<double[]> l = new List<double[]>();
            l.Add(new double[] { 0.0, 0.0, 0.0 });
            return l;
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
            throw new Exception("Alternative method <Project> has been used for higher order elements");
        }

        private Tuple<double[,], double[,], double[,], double[,], double[,]> CalculatePositionMatrix(double ksi1, double ksi2)
        {
            if (Properties.SlaveSegmentPolynomialDegree == 1 && Properties.MasterSegmentPolynomialDegree == 1)
            {
                double N1 = 1.0 / 2.0 * (1.0 - ksi1);
                double N2 = 1.0 / 2.0 * (1.0 + ksi1);
                double N3 = 1.0 / 2.0 * (1.0 - ksi2);
                double N4 = 1.0 / 2.0 * (1.0 + ksi2);
                double dN11 = -1.0 / 2.0;
                double dN21 = 1.0 / 2.0;
                double dN32 = -1.0 / 2.0;
                double dN42 = 1.0 / 2.0;
                double[,] aMatrix = new double[,]
                    {
                    { -N1, 0.0, -N2, 0.0, N3, 0.0, N4, 0.0 },
                    { 0.0, -N1, 0.0, -N2, 0.0, N3, 0.0, N4 }
                    };

                double[,] da1Matrix = new double[,]
                    {
                    { -dN11, 0.0, -dN21, 0.0, 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, -dN11, 0.0, -dN21, 0.0, 0.0, 0.0, 0.0 }
                    };
                double[,] da11Matrix = new double[,]
                    {
                    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
                    };
                double[,] da2Matrix = new double[,]
                    {
                    { 0.0, 0.0, 0.0, 0.0, dN32, 0.0, dN42, 0.0 },
                    { 0.0, 0.0, 0.0, 0.0, 0.0, dN32, 0.0, dN42 }
                    };
                double[,] da22Matrix = new double[,]
                    {
                    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
                    };
                return new Tuple<double[,], double[,], double[,], double[,], double[,]>(aMatrix, da1Matrix, da11Matrix, da2Matrix, da22Matrix);
            }
            else if (Properties.MasterSegmentPolynomialDegree == 1 && Properties.SlaveSegmentPolynomialDegree == 2)
            {
                double N1 = 1.0 / 2.0 * (1.0 - ksi1);
                double N2 = 1.0 / 2.0 * (1.0 + ksi1);
                double N3 = 1.0 / 2.0 * (Math.Pow(ksi2, 2.0) - ksi2);
                double N4 = 1.0 - Math.Pow(ksi2, 2.0);
                double N5 = 1.0 / 2.0 * (Math.Pow(ksi2, 2.0) + ksi2);
                double dN11 = -1.0 / 2.0;
                double dN21 = 1.0 / 2.0;
                double dN32 = ksi2 - 0.5;
                double dN42 = -2.0 * ksi2;
                double dN52 = ksi2 + 0.5;
                double dN322 = 1.0;
                double dN422 = -2.0;
                double dN522 = 1.0;
                double[,] aMatrix = new double[,]
                    {
                    { -N1, 0.0, -N2, 0.0, N3, 0.0 , N4, 0.0, N5, 0.0 },
                    { 0.0, -N1, 0.0, -N2, 0.0, N3, 0.0, N4, 0.0, N5 }
                    };

                double[,] da1Matrix = new double[,]
                    {
                    { -dN11, 0.0, -dN21, 0.0, 0.0, 0.0 , 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, -dN11, 0.0, -dN21, 0.0, 0.0 , 0.0, 0.0, 0.0, 0.0 }
                    };
                double[,] da11Matrix = new double[,]
                    {
                    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 , 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 , 0.0, 0.0, 0.0, 0.0 }
                    };
                double[,] da2Matrix = new double[,]
                    {
                    { 0.0, 0.0, 0.0, 0.0, dN32, 0.0 , dN42, 0.0, dN52, 0.0 },
                    { 0.0, 0.0, 0.0, 0.0, 0.0, dN32 , 0.0, dN42, 0.0, dN52 }
                    };
                double[,] da22Matrix = new double[,]
                    {
                    { 0.0, 0.0, 0.0, 0.0, dN322, 0.0 , dN422, 0.0, dN522, 0.0 },
                    { 0.0, 0.0, 0.0, 0.0, 0.0, dN322 , 0.0, dN422, 0.0, dN522 }
                    };
                return new Tuple<double[,], double[,], double[,], double[,], double[,]>(aMatrix, da1Matrix, da11Matrix, da2Matrix, da22Matrix);
            }
            else if (Properties.MasterSegmentPolynomialDegree == 2 && Properties.SlaveSegmentPolynomialDegree == 1)
            {

                double N1 = 1.0 / 2.0 * (Math.Pow(ksi1, 2.0) - ksi1);
                double N2 = 1.0 - Math.Pow(ksi1, 2.0);
                double N3 = 1.0 / 2.0 * (Math.Pow(ksi1, 2.0) + ksi1);
                double N4 = 1.0 / 2.0 * (1.0 - ksi2);
                double N5 = 1.0 / 2.0 * (1.0 + ksi2);
                double dN11 = ksi1 - 0.5;
                double dN21 = -2.0 * ksi1;
                double dN31 = ksi1 + 0.5;
                double dN42 = -1.0 / 2.0;
                double dN52 = 1.0 / 2.0;
                double dN111 = 1.0;
                double dN211 = -2.0;
                double dN311 = 1.0;
                double[,] aMatrix = new double[,]
                    {
                    { -N1 ,0.0 ,-N2 ,0.0 ,-N3 , 0.0 , N4, 0.0, N5, 0.0 },
                    {0.0, -N1 , 0.0 ,-N2, 0.0, -N3, 0.0, N4, 0.0, N5 }
                    };
                double[,] da1Matrix = new double[,]
                    {
                    { -dN11 ,0.0 ,-dN21 ,0.0 ,-dN31 , 0.0 , 0.0, 0.0, 0.0, 0.0 },
                    {0.0, -dN11 , 0.0 ,-dN21, 0.0, -dN31, 0.0, 0.0, 0.0, 0.0 }
                    };
                double[,] da11Matrix = new double[,]
                    {
                    { -dN111 ,0.0 ,-dN211 ,0.0 ,-dN311 , 0.0 , 0.0, 0.0, 0.0, 0.0 },
                    {0.0, -dN111 , 0.0 ,-dN211, 0.0, -dN311, 0.0, 0.0, 0.0, 0.0 }
                    };
                double[,] da2Matrix = new double[,]
                    {
                    { 0.0, 0.0, 0.0 ,0.0 ,0.0 , 0.0 , dN42, 0.0, dN52, 0.0 },
                    {0.0, 0.0 , 0.0 ,0.0, 0.0, 0.0, 0.0, dN42, 0.0, dN52 }
                    };
                double[,] da22Matrix = new double[,]
                    {
                    { 0.0, 0.0 ,0.0 ,0.0 ,0.0 , 0.0 , 0.0, 0.0, 0.0, 0.0 },
                    {0.0, 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
                    };
                return new Tuple<double[,], double[,], double[,], double[,], double[,]>(aMatrix, da1Matrix, da11Matrix, da2Matrix, da22Matrix);
            }
            else
            {
                double N1 = 1.0 / 2.0 * (Math.Pow(ksi1, 2.0) - ksi1);
                double N2 = 1.0 - Math.Pow(ksi1, 2.0);
                double N3 = 1.0 / 2.0 * (Math.Pow(ksi1, 2.0) + ksi1);
                double N4 = 1.0 / 2.0 * (Math.Pow(ksi2, 2.0) - ksi2);
                double N5 = 1.0 - Math.Pow(ksi2, 2.0);
                double N6 = 1.0 / 2.0 * (Math.Pow(ksi2, 2.0) + ksi2);
                double dN11 = ksi1 - 0.5;
                double dN21 = -2.0 * ksi1;
                double dN31 = ksi1 + 0.5;
                double dN42 = ksi2 - 0.5;
                double dN52 = -2.0 * ksi2;
                double dN62 = ksi2 + 0.5;
                double dN111 = 1.0;
                double dN211 = -2.0;
                double dN311 = 1.0;
                double dN422 = 1.0;
                double dN522 = -2.0;
                double dN622 = 1.0;
                double[,] aMatrix = new double[,]
                    {
                    { -N1 ,0.0 ,-N2 ,0.0 ,-N3 , 0.0 , N4, 0.0, N5, 0.0, N6, 0.0 },
                    {0.0, -N1 , 0.0 ,-N2, 0.0, -N3, 0.0, N4, 0.0, N5, 0.0, N6 }
                    };

                double[,] da1Matrix = new double[,]
                    {
                    { -dN11 ,0.0 ,-dN21 ,0.0 ,-dN31, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                    {0.0, -dN11 , 0.0 ,-dN21, 0.0, -dN31, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
                    };
                double[,] da11Matrix = new double[,]
                    {
                    { -dN111 ,0.0 ,-dN211 ,0.0 ,-dN311, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                    {0.0, -dN111 , 0.0 ,-dN211, 0.0, -dN311, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
                    };
                double[,] da2Matrix = new double[,]
                    {
                    { 0.0, 0.0 ,0.0, 0.0, 0.0, 0.0, dN42, 0.0, dN52, 0.0, dN62, 0.0 },
                    { 0.0, 0.0 , 0.0, 0.0, 0.0, 0.0, 0.0, dN42, 0.0, dN52, 0.0, dN62 }
                    };
                double[,] da22Matrix = new double[,]
                    {
                    { 0.0, 0.0 ,0.0, 0.0, 0.0, 0.0, dN422, 0.0, dN522, 0.0, dN622, 0.0 },
                    { 0.0, 0.0 , 0.0, 0.0, 0.0, 0.0, 0.0, dN422, 0.0, dN522, 0.0, dN622 }
                    };
                return new Tuple<double[,], double[,], double[,], double[,], double[,]>(aMatrix, da1Matrix, da11Matrix, da2Matrix, da22Matrix);
            }
        }

        private Tuple<double[], double, double[], double[], double> MasterSegmentGeometry(double[,] daMatrix, double[,] da2Matrix)
        {
            double[] xupd = VectorOperations.VectorScalarProduct(NodalXUpdated(), -1);
            double[] surfaceVector = VectorOperations.MatrixVectorProduct(daMatrix, xupd);
            double[] surfaceVectorDerivative = VectorOperations.MatrixVectorProduct(da2Matrix, xupd);

            double detm = VectorOperations.VectorDotProduct(surfaceVector, surfaceVector);
            double m11 = 1.0 / detm;
            double[] vector = new double[2];
            double scalarCoef = new double();
            double scalarCoef2 = Math.Pow(m11, 2.0);
            double[] tangentVector = new double[2];
            double curvatureTensor = new double();

            if (Properties.MasterSegmentPolynomialDegree == 1)
            {
                double Xm1 = Nodes[1].XCoordinate + DisplacementVector[0];
                double Ym1 = Nodes[1].YCoordinate + DisplacementVector[1];
                double Xm2 = Nodes[2].XCoordinate + DisplacementVector[2];
                double Ym2 = Nodes[2].YCoordinate + DisplacementVector[3];
                vector[0] = Ym2 - Ym1;
                vector[1] = Xm1 - Xm2;
                scalarCoef = -1.0 / (2.0 * Math.Sqrt(detm));
            }
            else 
            {
                vector[0] = -surfaceVector[1];
                vector[1] = surfaceVector[0];
                scalarCoef = 1.0 / (Math.Sqrt(detm));
                tangentVector = VectorOperations.VectorScalarProductNew(surfaceVector, scalarCoef);
            }
            double[] normalUnitVec = VectorOperations.VectorScalarProductNew(vector, scalarCoef);
            curvatureTensor = scalarCoef2 * VectorOperations.VectorDotProduct(surfaceVectorDerivative, normalUnitVec);
            return new Tuple<double[], double, double[], double[], double>(surfaceVector, m11, normalUnitVec, tangentVector, curvatureTensor);
        }
        private Tuple<double[], double, double[], double[], double> SlaveSegmentGeometry(double[,] daMatrix, double[,] da2Matrix)
        {
            double[] xupd = NodalXUpdated();
            double[] surfaceVector = VectorOperations.MatrixVectorProduct(daMatrix, xupd);
            double[] surfaceVectorDerivative = VectorOperations.MatrixVectorProduct(da2Matrix, xupd);

            double detm = VectorOperations.VectorDotProduct(surfaceVector, surfaceVector);
            double m11 = 1.0 / detm;
            double[] vector = new double[] { -surfaceVector[1], surfaceVector[0] };
            double scalarCoef = 1.0 / (Math.Sqrt(detm));
            double[] normalUnitVec = VectorOperations.VectorScalarProductNew(vector, scalarCoef);
            double[] tangentVector = VectorOperations.VectorScalarProductNew(surfaceVector, scalarCoef);
            double scalarCoef2 = Math.Pow(m11, 2.0);
            double curvatureTensor = scalarCoef2 * VectorOperations.VectorDotProduct(surfaceVectorDerivative, normalUnitVec);

            return new Tuple<double[], double, double[], double[], double>(surfaceVector, detm, normalUnitVec, tangentVector, curvatureTensor);
        }
        private double CalculatePenetration(double[,] aMatrix, double[] n)
        {
            double[,] AT = MatrixOperations.Transpose(aMatrix);
            double[] AT_n = VectorOperations.MatrixVectorProduct(AT, n);
            double[] xupd = NodalXUpdated();
            double normalGap = VectorOperations.VectorDotProduct(xupd, AT_n);
            return normalGap;
        }
        private double CalculateDeltaKsi(double[] masterSlaveRelativeVector, double[] surfaceVector, double[] surfaceVectorDerivative)
        {
            double scalar1 = VectorOperations.VectorDotProduct(surfaceVector, masterSlaveRelativeVector);
            double scalar2 = VectorOperations.VectorDotProduct(surfaceVectorDerivative, masterSlaveRelativeVector) -
                             VectorOperations.VectorDotProduct(surfaceVector, surfaceVector);

            double deltaKsi = -scalar1 / scalar2;
            return deltaKsi;
        }
        private double Project(double ksi1Initial, double ksi2)
        {
            if (Properties.MasterSegmentPolynomialDegree == 1) 
            {
                double[,] aMatrix = CalculatePositionMatrix(ksi1Initial, ksi2).Item1;
                int m = Properties.SlaveSegmentPolynomialDegree + 1;
                double[,] slaveNMatrix = new double[2,2*m];
                double[] xUpdated = NodalXUpdated();
                List<double> list = new List<double>() ;
                for(int i = 4; i< list.Count; i++)
                {
                    list.Add(xUpdated[i]);
                }
                double[] x = list.ToArray();
                for (int i = 0; i<=1; i++)
                {
                    if(i == 0)
                    {
                        int countCols = 0;
                        for (int j = 4; j < aMatrix.GetLength(1) - 1; j += 2)
                        {
                            slaveNMatrix[i, countCols] = aMatrix[i, j];
                            countCols += 2;
                        }
                    }
                    else
                    {
                        int countCols = 1;
                        for (int j = 5; j < aMatrix.GetLength(1); j += 2)
                        {
                            slaveNMatrix[i, countCols] = aMatrix[i, j];
                            countCols += 2;
                        }
                    }
                }
                double[] slavePositionVector = VectorOperations.MatrixVectorProduct(slaveNMatrix, x);
                double xM1 = xUpdated[0];
                double yM1 = xUpdated[1];
                double xM2 = xUpdated[2];
                double yM2 = xUpdated[3];
                double xS = slavePositionVector[0];
                double yS = slavePositionVector[1];
                double ksi = (2*(xS*(xM2 - xM1) + yS * (yM2 - yM1)) - Math.Pow(xM2,2) - Math.Pow(yM2, 2) + Math.Pow(xM1, 2) + Math.Pow(yM1, 2)) /(Math.Pow(xM2 - xM1, 2) + Math.Pow(yM2 - yM1, 2));
                return ksi;
            }
            else
            {
                int maxIterations = 1000;
                double tol = Math.Pow(10.0, -6.0);
                double deltaKsi = 0.0;
                double ksi = ksi1Initial;
                double[] xUpdated = NodalXUpdated();
                for (int i = 1; i <= maxIterations; i++)
                {
                    Tuple<double[,], double[,], double[,], double[,], double[,]> aMatrices = CalculatePositionMatrix(ksi, ksi2);
                    double[] masterSlaveRelativeVector = VectorOperations.MatrixVectorProduct(aMatrices.Item1, xUpdated);
                    double[] surfaceVector = VectorOperations.VectorScalarProduct(VectorOperations.MatrixVectorProduct(aMatrices.Item2, xUpdated), -1);
                    double[] surfaceVectorDerivative = VectorOperations.VectorScalarProduct(VectorOperations.MatrixVectorProduct(aMatrices.Item3, xUpdated), -1);

                    deltaKsi = CalculateDeltaKsi(masterSlaveRelativeVector, surfaceVector, surfaceVectorDerivative);
                    ksi += deltaKsi;
                    if (Math.Abs(deltaKsi) <= tol)
                    {
                        break;
                    }
                }
                if (Math.Abs(deltaKsi) > tol)
                {
                    throw new Exception("CPP not found in current iterations");
                }
                else
                {
                    return ksi;
                }
            }
        }
        private Tuple<double, double> GaussPoints(int i)
        {
            int iP = Properties.IntegrationPoints;
            double[] gaussPoints = new double[iP];
            double[] gaussWeights = new double[iP];
            if(iP == 1)
            {
                gaussPoints[0] = 0.0;
                gaussWeights[0] = 2.0;
            }
            else if(iP == 2)
            {
                gaussPoints[0] = -1.0 / Math.Sqrt(3);
                gaussPoints[1] = 1.0 / Math.Sqrt(3);
                gaussWeights[0] = 1.0;
                gaussWeights[1] = 1.0;
            }
            else if(iP == 3)
            {
                gaussPoints[0] = -0.77459;
                gaussPoints[1] = 0.0;
                gaussPoints[2] = 0.77459;
                gaussWeights[0] = 0.55555;
                gaussWeights[1] = 0.88888;
                gaussWeights[2] = 0.55555;
            }
            else if (iP == 4)
            {
                gaussPoints[0] = -0.86113;
                gaussPoints[1] = -0.33998;
                gaussPoints[2] = 0.33998;
                gaussPoints[3] = 0.86113;
                gaussWeights[0] = 0.34785;
                gaussWeights[1] = 0.65214;
                gaussWeights[2] = 0.65214;
                gaussWeights[3] = 0.34785;
            }
            else if(iP == 5)
            {
                gaussPoints[0] = -0.90617;
                gaussPoints[1] = -0.53846;
                gaussPoints[2] = 0.0;
                gaussPoints[3] = 0.53846;
                gaussPoints[4] = 0.90617;
                gaussWeights[0] = 0.23692;
                gaussWeights[1] = 0.47862;
                gaussWeights[2] = 0.56888;
                gaussWeights[3] = 0.47862;
                gaussWeights[4] = 0.23692;
            }
            else if(iP == 6)
            {
                gaussPoints[0] = -0.93246;
                gaussPoints[1] = -0.66120;
                gaussPoints[2] = -0.23861;
                gaussPoints[3] = 0.23861;
                gaussPoints[4] = 0.66120;
                gaussPoints[5] = 0.93246;
                gaussWeights[0] = 0.17132;
                gaussWeights[1] = 0.36076;
                gaussWeights[2] = 0.46791;
                gaussWeights[3] = 0.46791;
                gaussWeights[4] = 0.36076;
                gaussWeights[5] = 0.17132;
            }
            else if (iP == 7)
            {
                gaussPoints[0] = -0.94910;
                gaussPoints[1] = -0.74153;
                gaussPoints[2] = -0.40584;
                gaussPoints[3] = 0.0;
                gaussPoints[4] = 0.40584;
                gaussPoints[5] = 0.74153;
                gaussPoints[6] = 0.94910;
                gaussWeights[0] = 0.12948;
                gaussWeights[1] = 0.27970;
                gaussWeights[2] = 0.38183;
                gaussWeights[3] = 0.41795;
                gaussWeights[4] = 0.38183;
                gaussWeights[5] = 0.27970;
                gaussWeights[6] = 0.12948;
            }
            else if (iP == 8)
            {
                gaussPoints[0] = -0.96028;
                gaussPoints[1] = -0.79666;
                gaussPoints[2] = -0.52553;
                gaussPoints[3] = -0.18343;
                gaussPoints[4] = 0.18343;
                gaussPoints[5] = 0.52553;
                gaussPoints[6] = 0.79666;
                gaussPoints[7] = 0.96028;
                gaussWeights[0] = 0.10122;
                gaussWeights[1] = 0.22238;
                gaussWeights[2] = 0.31370;
                gaussWeights[3] = 0.36268;
                gaussWeights[4] = 0.36268;
                gaussWeights[5] = 0.31370;
                gaussWeights[6] = 0.22238;
                gaussWeights[7] = 0.10122;
            }
            else if (iP == 9)
            {
                gaussPoints[0] = -0.96816;
                gaussPoints[1] = -0.83603;
                gaussPoints[2] = -0.61337;
                gaussPoints[3] = -0.32425;
                gaussPoints[4] = 0.0;
                gaussPoints[5] = 0.32425;
                gaussPoints[6] = 0.61337;
                gaussPoints[7] = 0.83603;
                gaussPoints[8] = 0.96816;
                gaussWeights[0] = 0.08127;
                gaussWeights[1] = 0.18064;
                gaussWeights[2] = 0.26061;
                gaussWeights[3] = 0.31234;
                gaussWeights[4] = 0.33023;
                gaussWeights[5] = 0.31234;
                gaussWeights[6] = 0.26061;
                gaussWeights[7] = 0.18064;
                gaussWeights[8] = 0.08127;
            }
            else if (iP == 10)
            {
                gaussPoints[0] = -0.97390;
                gaussPoints[1] = -0.86506;
                gaussPoints[2] = -0.67940;
                gaussPoints[3] = -0.43339;
                gaussPoints[4] = -0.14887;
                gaussPoints[5] = 0.14887;
                gaussPoints[6] = 0.43339;
                gaussPoints[7] = 0.67940;
                gaussPoints[8] = 0.86506;
                gaussPoints[9] = 0.97390;
                gaussWeights[0] = 0.06667;
                gaussWeights[1] = 0.14945;
                gaussWeights[2] = 0.21908;
                gaussWeights[3] = 0.26926;
                gaussWeights[4] = 0.29552;
                gaussWeights[5] = 0.29552;
                gaussWeights[6] = 0.26926;
                gaussWeights[7] = 0.21908;
                gaussWeights[8] = 0.14945;
                gaussWeights[9] = 0.06667;
            }
            else if (iP == 11)
            {
                gaussPoints[0] = -0.97823;
                gaussPoints[1] = -0.88706;
                gaussPoints[2] = -0.73015;
                gaussPoints[3] = -0.51909;
                gaussPoints[4] = -0.26954;
                gaussPoints[5] = 0.0;
                gaussPoints[6] = 0.26954;
                gaussPoints[7] = 0.51909;
                gaussPoints[8] = 0.73015;
                gaussPoints[9] = 0.88706;
                gaussPoints[10] = 0.97823;
                gaussWeights[0] = 0.05567;
                gaussWeights[1] = 0.12558;
                gaussWeights[2] = 0.18629;
                gaussWeights[3] = 0.23319;
                gaussWeights[4] = 0.26280;
                gaussWeights[5] = 0.27292;
                gaussWeights[6] = 0.26280;
                gaussWeights[7] = 0.23319;
                gaussWeights[8] = 0.18629;
                gaussWeights[9] = 0.12558;
                gaussWeights[10] = 0.05567;
            }
            else if (iP == 12)
            {
                gaussPoints[0] = -0.98156;
                gaussPoints[1] = -0.90411;
                gaussPoints[2] = -0.76990;
                gaussPoints[3] = -0.58732;
                gaussPoints[4] = -0.36783;
                gaussPoints[5] = -0.12523;
                gaussPoints[6] = 0.12523;
                gaussPoints[7] = 0.36783;
                gaussPoints[8] = 0.58732;
                gaussPoints[9] = 0.76990;
                gaussPoints[10] = 0.90411;
                gaussPoints[11] = 0.98156;
                gaussWeights[0] = 0.047175;
                gaussWeights[1] = 0.106939;
                gaussWeights[2] = 0.16008;
                gaussWeights[3] = 0.20317;
                gaussWeights[4] = 0.23349;
                gaussWeights[5] = 0.249147;
                gaussWeights[6] = 0.249147;
                gaussWeights[7] = 0.23349;
                gaussWeights[8] = 0.20317;
                gaussWeights[9] = 0.16008;
                gaussWeights[10] = 0.106939;
                gaussWeights[11] = 0.047175;
            }
            else if (iP == 13)
            {
                gaussPoints[0] = -0.98418;
                gaussPoints[1] = -0.917598;
                gaussPoints[2] = -0.801578;
                gaussPoints[3] = -0.642349;
                gaussPoints[4] = -0.448493;
                gaussPoints[5] = -0.230458;
                gaussPoints[6] = 0.0;
                gaussPoints[7] = 0.230458;
                gaussPoints[8] = 0.448493;
                gaussPoints[9] = 0.642349;
                gaussPoints[10] = 0.801578;
                gaussPoints[11] = 0.917598;
                gaussPoints[12] = 0.98418;
                gaussWeights[0] = 0.040484;
                gaussWeights[1] = 0.0921215;
                gaussWeights[2] = 0.138874;
                gaussWeights[3] = 0.178146;
                gaussWeights[4] = 0.207816;
                gaussWeights[5] = 0.226283;
                gaussWeights[6] = 0.232552;
                gaussWeights[7] = 0.226283;
                gaussWeights[8] = 0.207816;
                gaussWeights[9] = 0.178146;
                gaussWeights[10] = 0.138874;
                gaussWeights[11] = 0.0921215;
                gaussWeights[12] = 0.040484;
            }
            else if (iP == 14)
            {
                gaussPoints[0] = -0.986284;
                gaussPoints[1] = -0.928435;
                gaussPoints[2] = -0.827201;
                gaussPoints[3] = -0.687293;
                gaussPoints[4] = -0.515249;
                gaussPoints[5] = -0.319112;
                gaussPoints[6] = -0.108055;
                gaussPoints[7] = 0.108055;
                gaussPoints[8] = 0.319112;
                gaussPoints[9] = 0.515249;
                gaussPoints[10] = 0.687293;
                gaussPoints[11] = 0.827201;
                gaussPoints[12] = 0.928435;
                gaussPoints[13] = 0.986284;
                gaussWeights[0] = 0.0351195;
                gaussWeights[1] = 0.0801581;
                gaussWeights[2] = 0.121519;
                gaussWeights[3] = 0.157203;
                gaussWeights[4] = 0.185538;
                gaussWeights[5] = 0.205198;
                gaussWeights[6] = 0.215264;
                gaussWeights[7] = 0.215264;
                gaussWeights[8] = 0.205198;
                gaussWeights[9] = 0.185538;
                gaussWeights[10] = 0.157203;
                gaussWeights[11] = 0.121519;
                gaussWeights[12] = 0.0801581;
                gaussWeights[13] = 0.0351195;
            }
            else
            {
                gaussPoints[0] = -0.987993;
                gaussPoints[1] = -0.937273;
                gaussPoints[2] = -0.848207;
                gaussPoints[3] = -0.724418;
                gaussPoints[4] = -0.570972;
                gaussPoints[5] = -0.394151;
                gaussPoints[6] = -0.201194;
                gaussPoints[7] = 0.0;
                gaussPoints[8] = 0.201194;
                gaussPoints[9] = 0.394151;
                gaussPoints[10] = 0.570972;
                gaussPoints[11] = 0.724418;
                gaussPoints[12] = 0.848207;
                gaussPoints[13] = 0.937273;
                gaussPoints[14] = 0.987993;
                gaussWeights[0] = 0.0307532;
                gaussWeights[1] = 0.070366;
                gaussWeights[2] = 0.107159;
                gaussWeights[3] = 0.139571;
                gaussWeights[4] = 0.166269;
                gaussWeights[5] = 0.186161;
                gaussWeights[6] = 0.198431;
                gaussWeights[7] = 0.202578;
                gaussWeights[8] = 0.198431;
                gaussWeights[9] = 0.186161;
                gaussWeights[10] = 0.166269;
                gaussWeights[11] = 0.139571;
                gaussWeights[12] = 0.107159;
                gaussWeights[13] = 0.070366;
                gaussWeights[14] = 0.0307532;
            }
            double GaussPoint = gaussPoints[i];
            double Weight = gaussWeights[i];
            return new Tuple<double, double>(GaussPoint, Weight);
        }
        private double[,] CalculateMainStiffnessPart(double ksi1, double ksi2, double[] n)
        {
            int numberOfNodes = Properties.MasterSegmentPolynomialDegree + Properties.SlaveSegmentPolynomialDegree + 2;
            double[,] mainStiffnessMatrix = new double[2 * numberOfNodes, 2* numberOfNodes];
            Tuple<double[,], double[,], double[,], double[,], double[,]> positionMatrices = CalculatePositionMatrix(ksi1, ksi2);
            double[,] A = positionMatrices.Item1;
            double[,] nxn = VectorOperations.VectorVectorTensorProduct(n, n);
            double[,] nxn_A = MatrixOperations.MatrixProduct(nxn, A);
            double[,] AT_nxn_A = MatrixOperations.MatrixProduct(MatrixOperations.Transpose(A), nxn_A);
            mainStiffnessMatrix = MatrixOperations.ScalarMatrixProductNew(PenaltyFactor, AT_nxn_A);
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
            return rotationalPart;
        }
        private double[,] CalculateCurvatureStiffnessPart(double[,] A, double ksi3, double m11, double[] dRho, double h11)
        {
            double coef = PenaltyFactor * ksi3 * m11 * h11;
            double[,] curvaturePart;
            double[,] dRho_x_dRho = VectorOperations.VectorVectorTensorProduct(dRho, dRho);
            double[,] Matrix = MatrixOperations.MatrixProduct(
                                                                    MatrixOperations.Transpose(A),
                                                                    MatrixOperations.MatrixProduct(dRho_x_dRho, A)
                                                                    );
            curvaturePart = MatrixOperations.ScalarMatrixProductNew(
                                                                        coef,
                                                                        Matrix
                                                                        );
            return curvaturePart;
        }
        public double[,] CreateGlobalStiffnessMatrix()
        {
            int nodesNumber = Properties.MasterSegmentPolynomialDegree + Properties.SlaveSegmentPolynomialDegree + 2;
            double[,] globalStifnessMatrix = new double[2* nodesNumber, 2* nodesNumber];
            for (int i = 0; i< Properties.IntegrationPoints; i++)
            {
                double ksi2 = GaussPoints(i).Item1;
                double gW = GaussPoints(i).Item2;
                double ksi1 = Project(0.0, ksi2);
                if (Math.Abs(ksi1) <= 1.05)
                {
                    Tuple<double[,], double[,], double[,], double[,], double[,]> positionMatrices = CalculatePositionMatrix(ksi1, ksi2);
                    double[,] aMatrix = positionMatrices.Item1;
                    double[,] daMatrix = positionMatrices.Item2;
                    double[,] da2Matrix = positionMatrices.Item3;

                    Tuple<double[], double, double[], double[], double> masterSurfaceCharacteristics = MasterSegmentGeometry(daMatrix, da2Matrix);
                    double m11 = masterSurfaceCharacteristics.Item2;
                    double[] dRho = masterSurfaceCharacteristics.Item1;
                    double[] n = masterSurfaceCharacteristics.Item3;
                    double h11 = masterSurfaceCharacteristics.Item5;
                    double ksi3 = CalculatePenetration(aMatrix, n);
                    if (ksi3 <= 0)
                    {
                        double slaveMetricTensor = SlaveSegmentGeometry(positionMatrices.Item4, positionMatrices.Item5).Item2;
                        double[,] mainPart = CalculateMainStiffnessPart(ksi1, ksi2, n);
                        double[,] rotationalPart = CalculateRotationalStiffnessPart(aMatrix, daMatrix, n, ksi3, m11, dRho);
                        double[,] curvaturePart = CalculateCurvatureStiffnessPart(aMatrix, ksi3, m11, dRho, h11);
                        double scalar = Math.Pow(slaveMetricTensor,0.5) * gW;
                        double[,] StifnessMatrix = MatrixOperations.ScalarMatrixProductNew(scalar, MatrixOperations.MatrixAddition(MatrixOperations.MatrixAddition(mainPart, rotationalPart),
                            curvaturePart));
                        globalStifnessMatrix = MatrixOperations.MatrixAddition(globalStifnessMatrix, StifnessMatrix);
                    }
                }
            }
            return globalStifnessMatrix;
        }

        public double[] CreateInternalGlobalForcesVector()
        {
            int nodesNumber = Properties.MasterSegmentPolynomialDegree + Properties.SlaveSegmentPolynomialDegree + 2;
            double[] internalGlobalForcesVector = new double[2* nodesNumber];
            for (int i = 0; i < Properties.IntegrationPoints; i++)
            {
                double ksi2 = GaussPoints(i).Item1;
                double gW = GaussPoints(i).Item2;
                double ksi1 = Project(0.0, ksi2);
                if (Math.Abs(ksi1) <= 1.05)
                {
                    Tuple<double[,], double[,], double[,], double[,], double[,]> positionMatrices = CalculatePositionMatrix(ksi1, ksi2);
                    double[,] aMatrix = positionMatrices.Item1;
                    double[,] daMatrix = positionMatrices.Item2;
                    double[,] da2Matrix = positionMatrices.Item3;

                    Tuple<double[], double, double[], double[], double> surfaceCharacteristics = MasterSegmentGeometry(daMatrix, da2Matrix);
                    double[] n = surfaceCharacteristics.Item3;
                    double ksi3 = CalculatePenetration(aMatrix, n);
                    if (ksi3 <= 0)
                    {
                        double slaveMetricTensor = SlaveSegmentGeometry(positionMatrices.Item4, positionMatrices.Item5).Item2;
                        double scalar = Math.Pow(slaveMetricTensor, 0.5) * gW;
                        double[,] AT = MatrixOperations.Transpose(aMatrix);
                        double[] AT_n = VectorOperations.MatrixVectorProduct(AT, n);
                        double[] internalLocalForcesVector = VectorOperations.VectorScalarProductNew(VectorOperations.VectorScalarProductNew(AT_n, PenaltyFactor * ksi3),scalar);
                        internalGlobalForcesVector = VectorOperations.VectorVectorAddition(internalGlobalForcesVector, internalLocalForcesVector);
                    }
                }
            }
            return internalGlobalForcesVector;
        }

        public double[,] CreateMassMatrix()
        {
            int nodesNumber = Properties.MasterSegmentPolynomialDegree + Properties.SlaveSegmentPolynomialDegree + 2;
            return new double[2* nodesNumber, 2* nodesNumber];
        }

        public double[,] CreateDampingMatrix()
        {
            int nodesNumber = Properties.MasterSegmentPolynomialDegree + Properties.SlaveSegmentPolynomialDegree + 2;
            return new double[2* nodesNumber, 2* nodesNumber];
        }
    }
}

