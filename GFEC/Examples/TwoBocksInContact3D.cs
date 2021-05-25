using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GFEC
{
    public static class TwoBlocksInContact3D
    {
        private const double BlockLength = 1.0;
        private const double ElementSize  = 1.0;
        private static int ElementsNumber = 9;
        private const double Gap  = 0.01;
        private const int nodesPerSide = 3;
        private static Dictionary<int, INode> nodes;
        public static ISolver newSolu;

        private static Dictionary<int, INode> CreateNodes()
        {
            
            nodes = new Dictionary<int, INode>();

            nodes[1] = new Node(0.0, 0.0, 0.0);
            nodes[2] = new Node(ElementSize, 0.0, 0.0);
            nodes[3] = new Node(2.0*ElementSize, 0.0, 0.0);

            nodes[4] = new Node(0.0, 0.0, ElementSize);
            nodes[5] = new Node(ElementSize, 0.0, ElementSize);
            nodes[6] = new Node(2.0*ElementSize, 0.0, ElementSize);

            nodes[7] = new Node(0.0, 0.0, 2.0 * ElementSize);
            nodes[8] = new Node(ElementSize, 0.0, 2.0*ElementSize);
            nodes[9] = new Node(2.0 * ElementSize, 0.0, 2.0*ElementSize);

            nodes[10] = new Node(0.0, ElementSize, 0.0);
            nodes[11] = new Node(ElementSize, ElementSize, 0.0);
            nodes[12] = new Node(2.0 * ElementSize, ElementSize, 0.0);

            nodes[13] = new Node(0.0, ElementSize, ElementSize);
            nodes[14] = new Node(ElementSize, ElementSize, ElementSize);
            nodes[15] = new Node(2.0 * ElementSize, ElementSize, ElementSize);

            nodes[16] = new Node(0.0, ElementSize, 2.0 * ElementSize);
            nodes[17] = new Node(ElementSize, ElementSize, 2.0 * ElementSize);
            nodes[18] = new Node(2.0 * ElementSize, ElementSize, 2.0 * ElementSize);

            nodes[19] = new Node(ElementSize / 2.0, ElementSize+Gap, ElementSize/2.0);
            nodes[20] = new Node(1.5 * ElementSize, ElementSize+Gap, ElementSize/2.0);
            nodes[21] = new Node(ElementSize / 2.0, ElementSize+Gap, 1.5*ElementSize);
            nodes[22] = new Node(1.5 * ElementSize, ElementSize+Gap, 1.5*ElementSize);

            nodes[23] = new Node(ElementSize / 2.0,2.0* ElementSize + Gap, ElementSize / 2.0);
            nodes[24] = new Node(1.5 * ElementSize,2.0* ElementSize + Gap, ElementSize / 2.0);
            nodes[25] = new Node(ElementSize / 2.0,2.0* ElementSize + Gap, 1.5 * ElementSize);
            nodes[26] = new Node(1.5 * ElementSize,2.0* ElementSize + Gap, 1.5 * ElementSize);

            return nodes;
        }

        private static Dictionary<int, Dictionary<int, int>> CreateConnectivity()
        {
            Dictionary<int, Dictionary<int, int>> connectivity = new Dictionary<int, Dictionary<int, int>>();
           
            connectivity[1] = new Dictionary<int, int>() { { 1, 7 }, { 2, 8 }, { 3, 5 }, { 4, 4 }, { 5, 16 }, { 6, 17 }, { 7, 14 }, { 8, 13 } };
            connectivity[2] = new Dictionary<int, int>() { { 1, 8 }, { 2, 9 }, { 3, 6 }, { 4, 5 }, { 5, 17 }, { 6, 18 }, { 7, 15 }, { 8, 14 } };
            connectivity[3] = new Dictionary<int, int>() { { 1, 4 }, { 2, 5 }, { 3, 2 }, { 4, 1 }, { 5, 13 }, { 6, 14 }, { 7, 11 }, { 8, 10 } };
            connectivity[4] = new Dictionary<int, int>() { { 1, 5 }, { 2, 6 }, { 3, 3 }, { 4, 2 }, { 5, 14 }, { 6, 15 }, { 7, 12 }, { 8, 11 } };

            connectivity[5] = new Dictionary<int, int>() { { 1, 21 }, { 2, 22 }, { 3, 20 }, { 4, 19 }, { 5, 25 }, { 6, 26}, { 7, 24 }, { 8, 23 } };

            connectivity[6] = new Dictionary<int, int>() { { 1, 16 }, { 2, 17 }, { 3, 14 }, { 4, 13 }, { 5, 21 } };
            connectivity[7] = new Dictionary<int, int>() { { 1, 17 }, { 2, 18 }, { 3, 15 }, { 4, 14 }, { 5, 22 } };
            connectivity[8] = new Dictionary<int, int>() { { 1, 13 }, { 2, 14 }, { 3, 11 }, { 4, 10 }, { 5, 19 } };
            connectivity[9] = new Dictionary<int, int>() { { 1, 14 }, { 2, 15 }, { 3, 12 }, { 4, 11 }, { 5, 20 } };


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
            newSolu.LinearScheme = new LUFactorization();
            //newSolu.NonLinearScheme = new LoadControlledNewtonRaphson();
            newSolu.ActivateNonLinearSolver = true;
            newSolu.NonLinearScheme.numberOfLoadSteps = 20;

            double[] externalForces = new double[78];
            externalForces[76] = -1000000000.0;


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

