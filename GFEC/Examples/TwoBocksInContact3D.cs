using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GFEC
{
    public static class TwoBlocksInContact3D2
    {
        private const double BlockLength = 1.0;
        private const double ElementSize  = 1.0;
        private const double ElementSizeb = 1.0;
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
            nodes[3] = new Node(ElementSize, 0.0, ElementSize);
            nodes[4] = new Node(0.0, 0.0, ElementSize);

            nodes[5] = new Node(0.0, ElementSize, 0.0);
            nodes[6] = new Node(ElementSize, ElementSize, 0.0);
            nodes[7] = new Node(ElementSize, ElementSize, ElementSize);
            nodes[8] = new Node(0.0, ElementSize, ElementSize);

            nodes[9] = new Node(0.0+0.2, ElementSize+Gap, 0.0);
            nodes[10] = new Node(ElementSize-0.2, ElementSize + Gap, 0.0);
            nodes[11] = new Node(ElementSize-0.2, ElementSize + Gap, ElementSize);
            nodes[12] = new Node(0.0+0.2, ElementSize + Gap, ElementSize);

            nodes[13] = new Node(0.0+0.2, 2*ElementSize + Gap, 0.0);
            nodes[14] = new Node(ElementSize-0.2, 2*ElementSize + Gap, 0.0);
            nodes[15] = new Node(ElementSize-0.2, 2*ElementSize + Gap, ElementSize);
            nodes[16] = new Node(0.0+0.2, 2*ElementSize + Gap, ElementSize);

            return nodes;
        }

        private static Dictionary<int, Dictionary<int, int>> CreateConnectivity()
        {
            Dictionary<int, Dictionary<int, int>> connectivity = new Dictionary<int, Dictionary<int, int>>();
           
            connectivity[1] = new Dictionary<int, int>() { { 1, 1 }, { 2, 2 }, { 3, 3 }, { 4, 4 }, { 5, 5 }, { 6, 6 }, { 7, 7 }, { 8, 8 } };
            connectivity[2] = new Dictionary<int, int>() { { 1, 9 }, { 2, 10 }, { 3, 11 }, { 4, 12 }, { 5, 13 }, { 6, 14 }, { 7, 15 }, { 8, 16 } };
            
            connectivity[3] = new Dictionary<int, int>() { { 1, 5 }, { 2, 6 }, { 3, 7 }, { 4, 8 }, { 5, 9 } };
            connectivity[4] = new Dictionary<int, int>() { { 1, 5 }, { 2, 6 }, { 3, 7 }, { 4, 8 }, { 5, 10 } };
            connectivity[5] = new Dictionary<int, int>() { { 1, 5 }, { 2, 6 }, { 3, 7 }, { 4, 8 }, { 5, 11 } };
            connectivity[6] = new Dictionary<int, int>() { { 1, 5 }, { 2, 6 }, { 3, 7 }, { 4, 8 }, { 5, 12 } };

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
            string type2 = "ContactNtS3Df";

            Dictionary<int, IElementProperties> elementProperties = new Dictionary<int, IElementProperties>();
            for (int i = 1; i <= 2; i++)
            {
                elementProperties[i] = new ElementProperties(E, A, type);
                elementProperties[i].Density = 8000.0;
                //elementProperties[i].Thickness = 0.1;
            }
            for (int i = 3; i <= 6; i++)
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

            assembly.BoundedDOFsVector = new int[12];
            for (int i = 1; i <= 12; i++)
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
            newSolu.LinearScheme = new BiCGSTABSolver();
            //newSolu.NonLinearScheme = new LoadControlledNewtonRaphson();
            newSolu.ActivateNonLinearSolver = true;
            newSolu.NonLinearScheme.numberOfLoadSteps = 10;

            double[] externalForces = new double[48];
            externalForces[46] = -1000.0;


            double[] reducedExternalFVector = BoundaryConditionsImposition.ReducedVector(externalForces, elementsAssembly.BoundedDOFsVector);

            newSolu.AssemblyData = elementsAssembly;
            newSolu.Solve(reducedExternalFVector);
            double[] solutionVector = newSolu.GetSolution();
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

