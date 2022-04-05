﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GFEC
{
    class ContactNtS2DTh : IElement
    {
        public Dictionary<int, INode> Nodes { get; }
        public IElementProperties Properties { get; set; }
        public Dictionary<int, bool[]> ElementFreedomSignature { get; } = new Dictionary<int, bool[]>();
        public List<int> ElementFreedomList { get; set; }
        public double[] DisplacementVector { get; set; }
        public double[] AccelerationVector { get; set; }
        private double PenaltyFactor { get; set; }
        private double ContactArea { get; set; }
        private double ContactPressure { get; set; }
        private double ContactThermalConductivity { get; set; }
        private double SurfaceRoughness { get; set; }
        private double YieldStrength { get; set; }
        private double Dx1 { get; set; }
        private double Dx2 { get; set; }
        private double Dx { get; set; }
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
        public ContactNtS2DTh(IElementProperties properties, Dictionary<int, INode> nodes)
        {
            Properties = properties;
            this.Nodes = nodes;
            ElementFreedomSignature[1] = new bool[] { true, false, false, false, false, false };
            ElementFreedomSignature[2] = new bool[] { true, false, false, false, false, false };
            ElementFreedomSignature[3] = new bool[] { true, false, false, false, false, false };
            DisplacementVector = new double[3];
            ContactArea = properties.SectionArea;
            ContactPressure = properties.ContactForceValue / properties.SectionArea;
            SurfaceRoughness = properties.SurfaceRoughness;
            ContactThermalConductivity = properties.ContactThermalConductivity;
            YieldStrength = properties.YieldStrength;
            Dx = properties.Dx;
            Dx1 = properties.Dx1;
            Dx2 = properties.Dx2;
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
        public List<double[]> GetStressVector()
        {
            throw new Exception("Needs to be removed. Has beeb used only for testing purposes");
        }
        public List<double[]> GetStrainVector()
        {
            throw new Exception("Needs to be removed. Has beeb used only for testing purposes");
        }
        public List<double[]> GetStressFromElementsNodes()
        {
            throw new Exception("Method not implemenented");

        }
        public List<double[]> GetStrainFromElementsNodes()
        {
            throw new Exception("Method not implemenented");

        }
        public List<double[]> GetGaussPointsInPhysicalSpace()
        {
            List<double[]> l = new List<double[]>();
            l.Add(new double[] { 0.0, 0.0 });
            return l;
        }
        public double CalculateConductivity()
        {
            //double k = 19.2;
            double m = 0.1259 * Math.Pow(SurfaceRoughness * Math.Pow(10, 6), 0.402);
            //double c1 = 6271.0 * Math.Pow(10, 6);
            //double c2 = -0.229;
            //double sigma = 0.478 * Math.Pow(10, -6);
            //double cc = (1.25 * k * m / sigma) * Math.Pow((ContactPressure / c1) * Math.Pow(1.6177 * 1000000 * sigma / m, -c2), 0.95 / (1 + 0.0711 * c2));
            double cc = 1.25 * ContactThermalConductivity * (m / (SurfaceRoughness)) * Math.Pow(ContactPressure / (3.0 * YieldStrength), 0.95);// * 2.582260191 * Math.Pow(10, -4);
            double cH = cc * ContactArea;
            return cH;
        }

        public double CalculateCoefficient1()
        {
            double W1;
            W1 = 1.0 - Dx1;  //(Dx - Dx1) / Dx;
            return W1;
        }
        public double CalculateCoefficient2()
        {
            double W2;
            double W1 = CalculateCoefficient1();
            W2 = 1.0 - W1;
            //W2 = 1.0 - Dx2 / Dx;
            return W2;
        }
        private double CalculateTemperatureJump()
        {
            double W1 = CalculateCoefficient1();
            double W2 = CalculateCoefficient2();
            double theta1 = 0.0;
            double theta2 = 0.0;
            double theta3 = 0.0;
            double gH = (theta3 + DisplacementVector[2]) - (W1 * (theta1 + DisplacementVector[0]) + W2 * (theta2 + DisplacementVector[1]));
            return gH;
        }
        private double CalculateTemperatureJump1()
        {
            //double W1 = CalculateCoefficient1();
            double theta1 = 0.0;
            //double theta2 = 0.0;
            double theta3 = 0.0;
            double gH1 = (theta3 + DisplacementVector[2]) - (theta1 + DisplacementVector[0]);
            return gH1;
        }
        private double CalculateTemperatureJump2()
        {
            //double W2 = CalculateCoefficient2();
            //double theta1 = 0.0;
            double theta2 = 0.0;
            double theta3 = 0.0;
            double gH2 = (theta3 + DisplacementVector[2]) - (theta2 + DisplacementVector[1]);
            return gH2;
        }
        public double[,] CreateGlobalStiffnessMatrix()
        {
            double cH = CalculateConductivity();
            double W1 = CalculateCoefficient1();
            double W2 = CalculateCoefficient2();
            double[,] stiffMatrix = new double[,]
            {
                { W1 * cH, 0.0, -W1 * cH},
                { 0.0, W2 * cH, -W2 * cH},
                { -W1 * cH, -W2 * cH, 1.0 * cH},
            };
            return stiffMatrix;
            //{
            //    {0.50 * W1 * cH, 0.50 * W2 * cH, -0.50 * cH},
            //    {0.50 * W1 * cH, 0.50 * W2 * cH, -0.50 * cH},
            //    {-W1 * cH, -W2 * cH, 1.0 * cH},
            //};
        }

        public double[] CreateInternalGlobalForcesVector()
        {
            double cH = CalculateConductivity();
            double gH = CalculateTemperatureJump();
            double gH1 = CalculateTemperatureJump1();
            double gH2 = CalculateTemperatureJump2();
            double W1 = CalculateCoefficient1();
            double W2 = CalculateCoefficient2();
            double[] heatFlux = new double[]
            {
                -W1 * cH * gH1 ,
                -W2 * cH * gH2 ,
                1.0 * cH * gH
            };
            return heatFlux;
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

