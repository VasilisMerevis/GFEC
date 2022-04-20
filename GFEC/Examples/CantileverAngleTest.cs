using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading;

namespace GFEC
{
    public static class CantileverAngleTest
    {
        private const int totalNodes = 648;

        private const int totalContactElements = 40;//42;//20;//8;
        private static int actualTotalContactElements;

        private const int totalElements = 480;
        private const int nodesInXCoor = 81;
        private const int nodesInYCoor = 4;
        private const double scaleFactor = 1.0;
        private const double xIntervals = 0.375;
        private const double yIntervals = 0.41;

        public static ISolver structuralSolution;
        public static ISolver thermalSolution;

        private const double angleInDegrees = 75;

        private const double angle = (Math.PI / 180) * angleInDegrees; //tested: 2.2, 2.57, 2.40, 2.12
        private static double offset = (nodesInXCoor - 1) * xIntervals * Math.Sin(angle) - totalContactElements * xIntervals + 0.1; //6.95;//8.1;//9.3// tested: 7.0 - 0.05, 7.0 -0.45, 7.0 -0.25, 7.0 -0.0
        private static double offsetInY = (nodesInXCoor - 1) * xIntervals * Math.Cos(angle);
        private static double gap = offsetInY + 0.05; //tested: 1.14, 2.75, 2.10, 0.75
        private static int loadStepsNumber = 40;

       
       
        //Model2
        static int[] structuralBoundaryConditions; // = new int[] { 1, 203, 505, 707, 909, 1012, 1014, 1016, 1018, 1020, 1022, 1024, 1026, 1211, 1413, 1615, 1817, 2019 };
        




        //External loads
        const double externalStructuralLoad = -2.85 * 4;//-2.6 * 4 - totalContactElements;//-2.6, -2,85
        //const double externalHeatLoad = 3300 * 100.0/2;//2500.0 * 1e-9;
        const double T0 = 100.0;
        






        static List<int> loadedStructuralDOFs; // = new List<int>(new int[] { 995, 997, 999, 1001, 1003, 1005, 1007, 1009 });
        static double[] externalForcesStructuralVector; // = new double[2020];

        const double YoungMod = 1.45e6;
        //const double density = 8000.0;
        //const double area = 0.01;
        const double thickness = 0.38;
        const double solidThermalCond = 3300.0;
        const double roughness = 0.0075;
        const double contactCond = 3300.0;
        const double yieldStrength = 60.0e3;



        private static void CreateStructuralBoundaryConditions()
        {
            List<int> boundedDofs = new List<int>();
            //for (int i = 0; i < nodesInYCoor; i++)
            //{
            //    boundedDofs.Add(i * 2 * nodesInXCoor + 1); //upper beam left side support
            //}

            //for (int i = 0; i < nodesInYCoor; i++)
            //{
            //    boundedDofs.Add(i * 2 * nodesInXCoor + 2 * nodesInXCoor - 1); //upper beam right side support
            //}

            for (int i = 0; i < nodesInYCoor; i++) //upper beam left side support
            {
                //boundedDofs.Add(i * nodesInXCoor * 2 + 1);
                boundedDofs.Add(i * nodesInXCoor * 2 + 2);
            }

            //for (int i = 1; i <= totalContactElements; i++)
            //{
            //    boundedDofs.Add(nodesInXCoor * nodesInYCoor * 2 + 2 * i); //lower beam lower side support
            //}
            //for (int i = 0; i < nodesInYCoor; i++)
            //{
            //    boundedDofs.Add(nodesInYCoor * nodesInXCoor * 2 + nodesInXCoor * 2 * (i+1) - 1); //lower beam right side support
            //}

            for (int i = 0; i < totalNodes; i++)
            {
                boundedDofs.Add(i * 2 + 1); //support for all nodes at X direction
            }

            //for (int i = totalNodes / 2 + 1; i < totalNodes; i++)
            //{
            //    boundedDofs.Add(i * 2 + 0);
            //    boundedDofs.Add(i * 2 + 1); //lower beam support for all nodes
            //}

            for (int i = 0; i < nodesInXCoor; i++)
            {
                boundedDofs.Add((totalNodes / 2 + i) * 2 + 2); //lower beam lower side y support
            }

            //for (int i = 0; i < nodesInYCoor; i++)
            //{
            //    boundedDofs.Add((totalNodes / 2 + (i + 1) * nodesInXCoor) * 2 - 1); //lower beam right side x support
            //    boundedDofs.Add((totalNodes / 2 + i * nodesInXCoor) * 2 + 1); //lower beam left side x support
            //}
            structuralBoundaryConditions = boundedDofs.ToArray<int>();
        }

        

        private static void CreateStructuralLoadVector()
        {
            loadedStructuralDOFs = new List<int>();
            //for (int i = 0; i < totalContactElements; i++)
            //{
            //    loadedStructuralDOFs.Add(nodesInXCoor * nodesInYCoor * 2 - 2 * i); //load in contact area
            //}

            for (int i = 0; i < nodesInXCoor - 1; i++)
            {
                loadedStructuralDOFs.Add((nodesInXCoor * (nodesInYCoor - 1) + i) * 2 + 2); //load at upper surface of upper beam
            }
            externalForcesStructuralVector = new double[totalNodes * 2];
        }

        private static Dictionary<int, INode> CreateNodes()
        {
            Dictionary<int, INode> nodes = new Dictionary<int, INode>();
            //Upper cantilever
            int k;
            k = 1;
            for (int j = 0; j < nodesInYCoor; j++)
            {
                for (int i = 0; i < nodesInXCoor; i++)
                {
                    nodes[k] = new Node(i * xIntervals * scaleFactor * Math.Sin(angle) + j * yIntervals * scaleFactor * Math.Cos(angle), -Math.Cos(angle) * i * xIntervals * scaleFactor + Math.Sin(angle) * j * yIntervals * scaleFactor);//j * yIntervals * scaleFactor - (j * yIntervals * scaleFactor * Math.Cos(angle)));
                    k += 1;
                }
            }

            //Lower cantilever
            for (int j = 0; j < nodesInYCoor; j++)
            {
                for (int i = 0; i < nodesInXCoor; i++)
                {
                    nodes[k] = new Node(i * xIntervals * scaleFactor + offset, j * yIntervals * scaleFactor - ((nodesInYCoor - 1) * yIntervals + gap));
                    k += 1;
                }
            }
            return nodes;
        }

        private static Dictionary<int, Dictionary<int, int>> CreateConnectivity()
        {

            Dictionary<int, Dictionary<int, int>> connectivity = new Dictionary<int, Dictionary<int, int>>();
            int k = 1;
            for (int j = 0; j <= nodesInYCoor - 2; j++)
            {
                for (int i = 1; i <= nodesInXCoor - 1; i++)
                {
                    connectivity[k] = new Dictionary<int, int>() { { 1, i + j * nodesInXCoor }, { 2, i + 1 + j * nodesInXCoor }, { 3, i + 1 + nodesInXCoor + j * nodesInXCoor }, { 4, i + nodesInXCoor + j * nodesInXCoor } };
                    k += 1;
                }

            }

            for (int j = 0 + nodesInYCoor; j <= nodesInYCoor - 2 + nodesInYCoor; j++)
            {
                for (int i = 1; i <= nodesInXCoor - 1; i++)
                {
                    connectivity[k] = new Dictionary<int, int>() { { 1, i + j * nodesInXCoor }, { 2, i + 1 + j * nodesInXCoor }, { 3, i + 1 + nodesInXCoor + j * nodesInXCoor }, { 4, i + nodesInXCoor + j * nodesInXCoor } };
                    k += 1;
                }

            }

            //Contact elements
            for (int i = 1; i <= totalContactElements - 1; i++)
            {
                int lowerMiddleNode = 2 * nodesInXCoor * nodesInYCoor - nodesInXCoor + i + 1;
                int lowerLeftNode = lowerMiddleNode - 1;
                int lowerRightNode = lowerMiddleNode + 1;
                int upperNode = nodesInXCoor - totalContactElements + i + 1;
                connectivity[totalElements + i] = new Dictionary<int, int>() { { 1, lowerLeftNode }, { 2, lowerRightNode }, { 3, upperNode } };
            }

            int count = connectivity.Count;
            for (int i = 1; i <= totalContactElements - 1 - 1; i++)
            {
                int lowerMiddleNode = 2 * nodesInXCoor * nodesInYCoor - nodesInXCoor + i + 1;
                int lowerLeftNode = lowerMiddleNode - 1;
                int lowerRightNode = lowerMiddleNode + 1;
                int upperNode = nodesInXCoor - totalContactElements + i + 1 + 1;
                connectivity[count + i] = new Dictionary<int, int>() { { 1, lowerLeftNode }, { 2, lowerRightNode }, { 3, upperNode } };
            }

            count = connectivity.Count;
            for (int i = 1; i <= totalContactElements - 1 - 1 - 1; i++)
            {
                int lowerMiddleNode = 2 * nodesInXCoor * nodesInYCoor - nodesInXCoor + i + 1;
                int lowerLeftNode = lowerMiddleNode - 1;
                int lowerRightNode = lowerMiddleNode + 1;
                int upperNode = nodesInXCoor - totalContactElements + i + 1 + 1 + 1;
                connectivity[count + i] = new Dictionary<int, int>() { { 1, lowerLeftNode }, { 2, lowerRightNode }, { 3, upperNode } };
            }

            count = connectivity.Count;
            for (int i = 1; i <= totalContactElements - 1 - 1 - 1 - 1; i++)
            {
                int lowerMiddleNode = 2 * nodesInXCoor * nodesInYCoor - nodesInXCoor + i + 1;
                int lowerLeftNode = lowerMiddleNode - 1;
                int lowerRightNode = lowerMiddleNode + 1;
                int upperNode = nodesInXCoor - totalContactElements + i + 1 + 1 + 1 + 1;
                connectivity[count + i] = new Dictionary<int, int>() { { 1, lowerLeftNode }, { 2, lowerRightNode }, { 3, upperNode } };
            }
            actualTotalContactElements = totalContactElements * 4 - 1 - 2 - 3 - 4;
            return connectivity;
        }

        private static Dictionary<int, bool[]> CreateNodeFAT()
        {
            Dictionary<int, bool[]> nodeFAT = new Dictionary<int, bool[]>();
            for (int i = 1; i <= totalNodes; i++)
            {
                nodeFAT[i] = new bool[] { true, true, false, false, false, false };
            }
            return nodeFAT;
        }

        
        private static Dictionary<int, IElementProperties> CreateElementProperties()
        {
            double E = YoungMod;
            //double A = area;
            string type = "Quad4";
            string type2 = "ContactNtS2D";

            Dictionary<int, IElementProperties> elementProperties = new Dictionary<int, IElementProperties>();
            for (int i = 1; i <= totalElements; i++)
            {
                elementProperties[i] = new ElementProperties(E, type);
            }

            for (int i = 1; i <= totalElements; i++)
            {
                //elementProperties[i].Density = density;
                elementProperties[i].Thickness = thickness;
            }

            for (int i = totalElements + 1; i <= totalElements + totalContactElements - 1 + totalContactElements - 2 + totalContactElements - 3 + totalContactElements - 4; i++)
            {
                elementProperties[i] = new ElementProperties(E, type2);
                //elementProperties[i].Density = density;
                elementProperties[i].Thickness = thickness;
            }
            return elementProperties;
        }

        private static IAssembly CreateAssembly()
        {
            IAssembly assembly = new Assembly();
            assembly.Nodes = CreateNodes();
            assembly.ElementsConnectivity = CreateConnectivity();
            assembly.ElementsProperties = CreateElementProperties();
            assembly.NodeFreedomAllocationList = CreateNodeFAT();
            CreateStructuralBoundaryConditions();
            CreateStructuralLoadVector();
            assembly.BoundedDOFsVector = structuralBoundaryConditions;
            return assembly;
        }

        public static Results RunStaticExample()
        {
            #region Structural
            IAssembly elementsAssembly = CreateAssembly();
            elementsAssembly.CreateElementsAssembly();
            elementsAssembly.ActivateBoundaryConditions = true;
            double[,] globalStiffnessMatrix = elementsAssembly.CreateTotalStiffnessMatrix();
            int countContactElements = elementsAssembly.CountElementsOfSameType(typeof(ContactNtS2D));
            elementsAssembly.SeperateContactDoF();
            //Gnuplot graphs
            ShowToGUI.PlotInitialGeometry(elementsAssembly);
            //ExportToFile.ExportMatlabInitialGeometry(elementsAssembly);
            Dictionary<int, INode> initialNodes = elementsAssembly.Nodes;
            double[] initialXCoord = Assembly.NodalCoordinatesToVectors(initialNodes).Item1;
            double[] initialYCoord = Assembly.NodalCoordinatesToVectors(initialNodes).Item2;

            double[] Xvec1Initial = new double[totalNodes / 2];
            double[] Yvec1Initial = new double[totalNodes / 2];
            double[] Xvec2Initial = new double[totalNodes / 2];
            double[] Yvec2Initial = new double[totalNodes / 2];
            double[] Ζvec1Initial = Enumerable.Repeat(1.0, totalNodes / 2).ToArray();
            double[] Ζvec2Initial = Enumerable.Repeat(1.0, totalNodes / 2).ToArray();

            Array.Copy(initialXCoord, 0, Xvec1Initial, 0, totalNodes / 2);
            Array.Copy(initialYCoord, 0, Yvec1Initial, 0, totalNodes / 2);

            Array.Copy(initialXCoord, totalNodes / 2, Xvec2Initial, 0, totalNodes / 2);
            Array.Copy(initialYCoord, totalNodes / 2, Yvec2Initial, 0, totalNodes / 2);
            string pathForContour1 = @"C:\Users\Public\Documents\Total\1";
            string pathForContour2 = @"C:\Users\Public\Documents\Total\2";
            //ExportToFile.CreateContourDataForMatlab(Xvec1Initial, Yvec1Initial, Ζvec1Initial, nodesInYCoor, nodesInXCoor, pathForContour1);
            //ExportToFile.CreateContourDataForMatlab(Xvec2Initial, Yvec2Initial, Ζvec2Initial, nodesInYCoor, nodesInXCoor, pathForContour2);




            ///structuralSolution = new StaticSolver();
            structuralSolution.LinearScheme = new MMCPCGSolver();
            //structuralSolution.NonLinearScheme = new LoadControlledNewtonRaphson();
            structuralSolution.NonLinearScheme.Tolerance = 1e-5;
            structuralSolution.ActivateNonLinearSolver = true;
            structuralSolution.NonLinearScheme.numberOfLoadSteps = loadStepsNumber;

            double[] externalForces3 = externalForcesStructuralVector;
            foreach (var dof in loadedStructuralDOFs)
            {
                externalForces3[dof - 1] = externalStructuralLoad;
            }



            double[] reducedExternalForces3 = elementsAssembly.MMCPGCreateReducedFromFullVector(externalForces3);
            structuralSolution.AssemblyData = elementsAssembly;
            structuralSolution.Solve(reducedExternalForces3);
            double[] solvector3 = structuralSolution.GetSolution();
            elementsAssembly.UpdateDisplacements(solvector3);
            ShowToGUI.PlotFinalGeometry(elementsAssembly);
            double[] fullSolVector3 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(solvector3, elementsAssembly.BoundedDOFsVector);
            Dictionary<int, INode> finalNodes = Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullSolVector3);
            double[] xFinalNodalCoor = Assembly.NodalCoordinatesToVectors(finalNodes).Item1;
            double[] yFinalNodalCoor = Assembly.NodalCoordinatesToVectors(finalNodes).Item2;
            Dictionary<int, double[]> allStepsSolutions = structuralSolution.GetAllStepsSolutions();

            Dictionary<int, Dictionary<int, double[]>> allStepsContactForces = new Dictionary<int, Dictionary<int, double[]>>();
            Dictionary<int, Dictionary<int, double>> allStepsProjectionPoints = new Dictionary<int, Dictionary<int, double>>();
            Dictionary<int, double[]> elementsInternalContactForcesVector;
            Dictionary<int, double> projectionPointForEachElement;
            for (int i = 1; i <= allStepsSolutions.Count; i++)
            {
                elementsInternalContactForcesVector = new Dictionary<int, double[]>();
                projectionPointForEachElement = new Dictionary<int, double>();
                elementsAssembly.UpdateDisplacements(allStepsSolutions[i]);
                for (int j = totalElements + 1; j <= totalElements + countContactElements; j++)
                {
                    elementsInternalContactForcesVector[j] = elementsAssembly.ElementsAssembly[j].CreateInternalGlobalForcesVector();
                    projectionPointForEachElement[j] = elementsAssembly.ElementsAssembly[j].ClosestPointProjection();
                }
                allStepsContactForces[i] = elementsInternalContactForcesVector;
                allStepsProjectionPoints[i] = projectionPointForEachElement;
            }



            List<double[]> structuralSolutions = new List<double[]>();

            //ExportToFile.ExportMatlabInitialGeometry(elementsAssembly);
            #endregion






            return new Results() { NonlinearSolution = structuralSolutions, SelectedDOF = 2, SolutionType = "Nonlinear" };
        }

        public static void RunDynamicExample()
        {
            IAssembly elementsAssembly = CreateAssembly();
            elementsAssembly.CreateElementsAssembly();
            elementsAssembly.ActivateBoundaryConditions = true;

            InitialConditions initialValues = new InitialConditions();
            initialValues.InitialAccelerationVector = new double[6];
            initialValues.InitialDisplacementVector = new double[6];
            //initialValues.InitialDisplacementVector[7] = -0.02146;
            initialValues.InitialVelocityVector = new double[6];
            initialValues.InitialTime = 0.0;

            ExplicitSolver newSolver = new ExplicitSolver(1.0, 10000);
            newSolver.Assembler = elementsAssembly;

            newSolver.InitialValues = initialValues;
            newSolver.ExternalForcesVector = new double[] { 0, 0, 0, 0, -50000, -50000 };
            newSolver.LinearSolver = new CholeskyFactorization();
            newSolver.ActivateNonLinearSolution = true;
            newSolver.SolveNewmark();
            newSolver.PrintExplicitSolution();//
        }

    }
}

