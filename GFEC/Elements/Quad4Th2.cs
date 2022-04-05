using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GFEC
{
    class Quad4Th2 : IElement
    {
        public Dictionary<int, INode> Nodes { get; }
        public IElementProperties Properties { get; set; }
        public Dictionary<int, bool[]> ElementFreedomSignature { get; } = new Dictionary<int, bool[]>();
        public List<int> ElementFreedomList { get; set; }
        public double[] DisplacementVector { get; set; }
        public double[] AccelerationVector { get; set; }
        public double kc;
        private double A { get; set; }
        private double B { get; set; }
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
        //private double thickness = 1.0; //To be included in Element Properties
        //private double density = 1.0; //To be included in Element Properties

        public Quad4Th2(IElementProperties properties, Dictionary<int, INode> nodes)
        {
            Properties = properties;
            this.Nodes = nodes;
            ElementFreedomSignature[1] = new bool[] { true, false, false, false, false, false };
            ElementFreedomSignature[2] = new bool[] { true, false, false, false, false, false };
            ElementFreedomSignature[3] = new bool[] { true, false, false, false, false, false };
            ElementFreedomSignature[4] = new bool[] { true, false, false, false, false, false };
            DisplacementVector = new double[4];
            A = properties.A;
            B = properties.B;
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
            //    double[] positionVector = VectorOperations.MatrixVectorProduct(CalculateShapeFunctionMatrix(parametricCoordinatesVec[0], parametricCoordinatesVec[1]), xUpdated);
            //    PositionVectorsList.Add(positionVector);
            //}
            return PositionVectorsList;
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
        public List<double[]> GetGaussPointsInPhysicalSpace()
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
        public double[,] CreateGlobalStiffnessMatrix()
        {
            kc = Properties.ThermalConductivity;
            double[,] K = new double[4, 4];

            K[0, 0] = 2 * (Math.Pow(A, 2) + Math.Pow(B, 2)) * kc / (6.0 * A * B);
            K[0, 1] = (Math.Pow(A, 2) - 2 * Math.Pow(B, 2)) * kc / (6.0 * A * B);
            K[0, 2] = -(Math.Pow(A, 2) + Math.Pow(B, 2)) * kc / (6.0 * A * B);
            K[0, 3] = (Math.Pow(B, 2) - 2 * Math.Pow(A, 2)) * kc / (6.0 * A * B);

            K[1, 0] = K[0, 1];
            K[1, 1] = 2 * (Math.Pow(A, 2) + Math.Pow(B, 2)) * kc / (6.0 * A * B);
            K[1, 2] = (Math.Pow(B, 2) - 2 * Math.Pow(A, 2)) * kc / (6.0 * A * B);
            K[1, 3] = -(Math.Pow(A, 2) + Math.Pow(B, 2)) * kc / (6.0 * A * B);

            K[2, 0] = K[0, 2];
            K[2, 1] = K[1, 2];
            K[2, 2] = 2 * (Math.Pow(A, 2) + Math.Pow(B, 2)) * kc / (6.0 * A * B);
            K[2, 3] = (Math.Pow(A, 2) - 2 * Math.Pow(B, 2)) * kc / (6.0 * A * B);

            K[3, 0] = K[0, 3];
            K[3, 1] = K[1, 3];
            K[3, 2] = K[2, 3];
            K[3, 3] = 2 * (Math.Pow(A, 2) + Math.Pow(B, 2)) * kc / (6.0 * A * B);

            return K;
        }

        public double[,] CreateMassMatrix()
        {
            throw new Exception("Mass matrix not implemented for Quad4Th element");
        }

        public double[,] CreateDampingMatrix()
        {
            throw new Exception("Damping matrix not implemented for Quad4Th element");
        }

        public double[] CreateInternalGlobalForcesVector()
        {
            double[] intForces;
            double[,] stiff = CreateGlobalStiffnessMatrix();

            intForces = VectorOperations.MatrixVectorProduct(stiff, DisplacementVector);
            intForces = VectorOperations.VectorScalarProductNew(intForces, 1.0);
            return intForces;
        }
    }
}

