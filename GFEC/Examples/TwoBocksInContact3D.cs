using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GFEC
{
    public static class TwoBlocksInContact3D
    {
        private const double BlockLength = 1.0;
        private const double ElementSize  = 0.1;
        private static int ElementsNumber = 9;
        private const double Gap  = 0.01;
        private const int nodesPerSide = 3;
        private static Dictionary<int, INode> nodes;
        public static ISolver newSolu;

        private static Dictionary<int, INode> CreateNodes()
        {
            //nodesPerSide = (int)(BlockLength / ElementSize + 1.0);
            //ElementsNumber = (int)Math.Pow(BlockLength / ElementSize, 2.0);
            nodes = new Dictionary<int, INode>();

            int nodeNumber = 0;
            for (int k = 0; k < 2; k++)
            {
                for (int i = 0; i < nodesPerSide; i++)
                {
                    for (int j = 0; j < nodesPerSide; j++)
                    {
                        nodeNumber = nodeNumber + 1;
                        double x = i * ElementSize;
                        double y = j * ElementSize;
                        double z = k * ElementSize;
                        nodes[nodeNumber] = new Node(x, y, z);
                    }
                }
            }

            nodes[19] = new Node(ElementSize / 2.0, ElementSize / 2.0, ElementSize + Gap);
            nodes[20] = new Node(3.0 * ElementSize / 2.0, ElementSize / 2.0, ElementSize + Gap);
            nodes[21] = new Node(ElementSize / 2.0, 3.0 * ElementSize / 2.0, ElementSize + Gap);
            nodes[22] = new Node(3.0 * ElementSize / 2.0, 3.0 * ElementSize / 2.0, ElementSize + Gap);

            nodes[23] = new Node(ElementSize / 2.0, ElementSize / 2.0, 2.0 * ElementSize + Gap);
            nodes[24] = new Node(3.0 * ElementSize / 2.0, ElementSize / 2.0, 2.0 * ElementSize + Gap);
            nodes[25] = new Node(ElementSize / 2.0, 3.0 * ElementSize / 2.0, 2.0 * ElementSize + Gap);
            nodes[26] = new Node(3.0 * ElementSize / 2.0, 3.0 * ElementSize / 2.0, 2.0 * ElementSize + Gap);
            //for (int k = 0; k < nodesPerSide; k++)
            //{
            //    for (int i = 0; i < nodesPerSide; i++)
            //    {
            //        for (int j = 0; j < nodesPerSide; j++)
            //        {
            //            nodeNumber = nodeNumber + 1;
            //            double x = i * ElementSize;
            //            double y = j * ElementSize;
            //            double z = k * ElementSize + BlockLength + Gap;
            //            nodes[nodeNumber] = new Node(x, y, z);
            //        }
            //    }
            //}
            return nodes;
        }

        private static Dictionary<int, Dictionary<int, int>> CreateConnectivity()
        {
            Dictionary<int, Dictionary<int, int>> connectivity = new Dictionary<int, Dictionary<int, int>>();
            //int element = 1;
            //for (int i = 1; i < nodesPerSide; i++)
            //{
            //    for (int j = 1; j < nodesPerSide; j++)
            //    {
            //        int localNode1 = (i - 1) * nodesPerSide + j;
            //        int localNode2 = (i - 1) * nodesPerSide + j + 1;
            //        int localNode3 = i * nodesPerSide + j;
            //        int localNode4 = i * nodesPerSide + j + 1;
            //        int localNode5 = i * nodesPerSide + j + 1;
            //        int localNode6 = i * nodesPerSide + j + 1;
            //        int localNode7 = i * nodesPerSide + j + 1;
            //        int localNode8 = i * nodesPerSide + j + 1;
            //        connectivity[element] = new Dictionary<int, int>() { { 1, localNode1 }, { 2, localNode2 }, { 3, localNode3 }, { 4, localNode4 } };
            //        element = element + 1;
            //    }
            //}
            //int firstBodyTotalNodes = nodesPerSide * nodesPerSide;
            //for (int i = 1; i < nodesPerSide; i++)
            //{
            //    for (int j = 1; j < nodesPerSide; j++)
            //    {
            //        int localNode1 = (i - 1) * nodesPerSide + j + firstBodyTotalNodes;
            //        int localNode2 = (i - 1) * nodesPerSide + j + 1 + firstBodyTotalNodes;
            //        int localNode3 = i * nodesPerSide + j + firstBodyTotalNodes;
            //        int localNode4 = i * nodesPerSide + j + 1 + firstBodyTotalNodes;
            //        connectivity[element] = new Dictionary<int, int>() { { 1, localNode1 }, { 2, localNode2 }, { 3, localNode3 }, { 4, localNode4 } };
            //        element = element + 1;
            //    }
            //}

            //for (int i = 1; i <= nodesPerSide; i++)
            //{
            //    int localNode1 = i + firstBodyTotalNodes - nodesPerSide;
            //    int localNode2 = i + firstBodyTotalNodes;
            //    connectivity[element] = new Dictionary<int, int>() { { 1, localNode1 }, { 2, localNode2 } };
            //    element = element + 1;
            //}
            connectivity[1] = new Dictionary<int, int>() { { 1, 1 }, { 2, 2 }, { 3, 4 }, { 4, 5 }, { 5, 10 }, { 6, 11 }, { 7, 13 }, { 8, 14 } };
            connectivity[2] = new Dictionary<int, int>() { { 1, 2 }, { 2, 3 }, { 3, 5 }, { 4, 6 }, { 5, 11 }, { 6, 12 }, { 7, 14 }, { 8, 15 } };
            connectivity[3] = new Dictionary<int, int>() { { 1, 4 }, { 2, 5 }, { 3, 7 }, { 4, 8 }, { 5, 13 }, { 6, 14 }, { 7, 16 }, { 8, 17 } };
            connectivity[4] = new Dictionary<int, int>() { { 1, 5 }, { 2, 6 }, { 3, 8 }, { 4, 9 }, { 5, 14 }, { 6, 15 }, { 7, 17 }, { 8, 18 } };

            connectivity[5] = new Dictionary<int, int>() { { 1, 19 }, { 2, 20 }, { 3, 21 }, { 4, 22 }, { 5, 23 }, { 6, 24}, { 7, 25 }, { 8, 26 } };

            connectivity[6] = new Dictionary<int, int>() { { 1, 10 }, { 2, 11 }, { 3, 13 }, { 4, 14 }, { 5, 19 } };
            connectivity[7] = new Dictionary<int, int>() { { 1, 11 }, { 2, 12 }, { 3, 14 }, { 4, 15 }, { 5, 20 } };
            connectivity[8] = new Dictionary<int, int>() { { 1, 13 }, { 2, 14 }, { 3, 16 }, { 4, 17 }, { 5, 21 } };
            connectivity[9] = new Dictionary<int, int>() { { 1, 14 }, { 2, 15 }, { 3, 17 }, { 4, 18 }, { 5, 22 } };


            ElementsNumber = connectivity.Count;
            return connectivity;
        }

        private static Dictionary<int, bool[]> CreateNodeFAT()
        {
            Dictionary<int, bool[]> nodeFAT = new Dictionary<int, bool[]>();
            for (int i = 1; i <= nodes.Count; i++)
            {
                nodeFAT[i] = new bool[] { true, true, true, false, false, false };
            }
            return nodeFAT;
        }

        private static Dictionary<int, IElementProperties> CreateElementProperties()
        {
            double E = 200.0e9;
            double A = 1.0;
            string type = "Hex8";
            string type2 = "ContactNtS3D";

            Dictionary<int, IElementProperties> elementProperties = new Dictionary<int, IElementProperties>();
            for (int i = 1; i <= 5; i++)
            {
                elementProperties[i] = new ElementProperties(E, A, type);
                elementProperties[i].Density = 8000.0;
                //elementProperties[i].Thickness = 0.1;
            }
            for (int i = 6; i <= 9; i++)
            {
                elementProperties[i] = new ElementProperties(E, A, type2);
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

            assembly.BoundedDOFsVector = new int[27];
            for (int i = 1; i <= 27; i++)
            {
                assembly.BoundedDOFsVector[i - 1] = i;
            }
            return assembly;
        }

        public static Results RunStaticExample()
        {
            IAssembly elementsAssembly = CreateAssembly();
            elementsAssembly.CreateElementsAssembly();
            elementsAssembly.ActivateBoundaryConditions = true;
            double[,] globalStiffnessMatrix = elementsAssembly.CreateTotalStiffnessMatrix();

            //ISolver newSolu = new StaticSolver();
            newSolu.LinearScheme = new PCGSolver();
            //newSolu.NonLinearScheme = new LoadControlledNewtonRaphson();
            newSolu.ActivateNonLinearSolver = true;
            newSolu.NonLinearScheme.numberOfLoadSteps = 20;

            double[] externalForces = new double[78];
            externalForces[67] = -100.0;
            externalForces[70] = -100.0;
            externalForces[73] = -100.0;
            externalForces[76] = -100.0;

            double[] reducedExternalFVector = BoundaryConditionsImposition.ReducedVector(externalForces, elementsAssembly.BoundedDOFsVector);

            newSolu.AssemblyData = elementsAssembly;
            newSolu.Solve(reducedExternalFVector);
            newSolu.PrintSolution();

            return new Results() { NonlinearSolution = new List<double[]>(), SelectedDOF = 2, SolutionType = "Nonlinear" };
        }

        public static void RunDynamicExample()
        {
            IAssembly elementsAssembly = CreateAssembly();
            elementsAssembly.CreateElementsAssembly();
            elementsAssembly.ActivateBoundaryConditions = true;



            InitialConditions initialValues = new InitialConditions();
            initialValues.InitialAccelerationVector = new double[462];
            initialValues.InitialDisplacementVector = new double[462];
            //initialValues.InitialDisplacementVector[7] = -0.02146;
            initialValues.InitialVelocityVector = new double[462];
            initialValues.InitialTime = 0.0;

            ExplicitSolver newSolver = new ExplicitSolver(1.0, 1000000);
            newSolver.Assembler = elementsAssembly;

            newSolver.InitialValues = initialValues;
            newSolver.ExternalForcesVector = new double[462];
            for (int i = 441; i <= 462; i += 2)
            {
                newSolver.ExternalForcesVector[i] = -10000.0;
            }
            newSolver.LinearSolver = new CholeskyFactorization();
            newSolver.ActivateNonLinearSolution = true;
            newSolver.SolveExplicit();
            //newSolver.PrintExplicitSolution();
        }

    }
}

