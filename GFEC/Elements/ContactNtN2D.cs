﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GFEC
{
    class ContactNtN2D : IElement
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
        public ContactNtN2D(IElementProperties properties, Dictionary<int, INode> nodes)
        {
            Properties = properties;
            this.Nodes = nodes;
            ElementFreedomSignature[1] = new bool[] { true, true, false, false, false, false };
            ElementFreedomSignature[2] = new bool[] { true, true, false, false, false, false };
            DisplacementVector = new double[4];
            //PenaltyFactor = properties.YoungMod * 1000.0 * properties.SectionArea;
            PenaltyFactor = properties.YoungMod * 10.0;
            //PenaltyFactor = properties.YoungMod * 3000.0 * properties.SectionArea;
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
        public Dictionary<int, INode> NodesAtFinalState()
        {
            Dictionary<int, INode> finalNodes = new Dictionary<int, INode>();
            finalNodes[1] = new Node(Nodes[1].XCoordinate + DisplacementVector[0], Nodes[1].YCoordinate + DisplacementVector[1]);
            finalNodes[2] = new Node(Nodes[2].XCoordinate + DisplacementVector[2], Nodes[2].YCoordinate + DisplacementVector[3]);
            return finalNodes;
        }

        private double[] CalculateNormalUnitVector()
        {
            double X1 = Nodes[1].XCoordinate;
            double Y1 = Nodes[1].YCoordinate;
            double X2 = Nodes[2].XCoordinate;
            double Y2 = Nodes[2].YCoordinate;
            double[] normalVector = new double[] { X2 - X1, Y2 - Y1 };
            double normalVectorLength = VectorOperations.VectorNorm2(normalVector);
            double[] normalUnitVec = new double[] { normalVector[0] / normalVectorLength, normalVector[1] / normalVectorLength };
            return normalUnitVec;
        }

        private double[,] CalculatePositionMatrix()
        {
            double[,] aMatrix = new double[,]
                {
                    { -1,0,1,0},
                    {0,-1,0,1 }
                };
            return aMatrix;
        }

        private double CalculateNormalGap()
        {
            double[,] A = CalculatePositionMatrix();
            double[,] AT = MatrixOperations.Transpose(A);
            double[] n = CalculateNormalUnitVector();
            double[] AT_n = VectorOperations.MatrixVectorProduct(AT, n);
            double[] xupd = new double[] {
                Nodes[1].XCoordinate + DisplacementVector[0],
                Nodes[1].YCoordinate + DisplacementVector[1],
                Nodes[2].XCoordinate + DisplacementVector[2],
                Nodes[2].YCoordinate + DisplacementVector[3]
            };
            double normalGap = VectorOperations.VectorDotProduct(xupd, AT_n);
            return normalGap;
        }

        public double[,] CreateGlobalStiffnessMatrix()
        {
            double penetration = CalculateNormalGap();
            if (penetration <= 0)
            {
                double[] n = CalculateNormalUnitVector();
                double[,] A = CalculatePositionMatrix();
                double[,] AT = MatrixOperations.Transpose(A);
                double[,] nxn = VectorOperations.VectorVectorTensorProduct(n, n);
                double[,] nxn_A = MatrixOperations.MatrixProduct(nxn, A);
                double[,] AT_nxn_A = MatrixOperations.MatrixProduct(AT, nxn_A);
                double[,] globalStiffnessMatrix = MatrixOperations.ScalarMatrixProductNew(PenaltyFactor, AT_nxn_A);
                return globalStiffnessMatrix;
            }
            else
            {
                double[,] globalStifnessMatrix = new double[4, 4];
                return globalStifnessMatrix;
            }

        }

        public double[] CreateInternalGlobalForcesVector()
        {
            double penetration = CalculateNormalGap();
            if (penetration <= 0)
            {
                double[,] A = CalculatePositionMatrix();
                double[,] AT = MatrixOperations.Transpose(A);
                double[] n = CalculateNormalUnitVector();
                double[] AT_n = VectorOperations.MatrixVectorProduct(AT, n);
                double ksi = CalculateNormalGap();
                double[] ksi_AT_n = VectorOperations.VectorScalarProductNew(AT_n, ksi);
                double[] e_ksi_AT_n = VectorOperations.VectorScalarProductNew(ksi_AT_n, PenaltyFactor);
                return e_ksi_AT_n;
            }
            else
            {
                double[] internalGlobalForcesVector = new double[4];
                return internalGlobalForcesVector;
            }
        }

        public double[,] CreateMassMatrix()
        {
            return new double[4, 4];
        }

        public double[,] CreateDampingMatrix()
        {
            return new double[4, 4];
        }
    }
}

