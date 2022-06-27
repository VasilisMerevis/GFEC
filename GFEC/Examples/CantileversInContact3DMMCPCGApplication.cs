using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GFEC
{
    public static class CantileversInContact3DMMCPCGApplication
    {

        public static ISolver structuralSolution;
        static int[] structuralBoundaryConditions;
        const double thickness = 0.1;
        const int nodesInX = 8;
        const int nodesInY = 2;
        const int nodesInZ = 2;
        const int nodesNumber = 64;
        const int elementsNumber = 14;
        const int contactElements = 4;

        //const double xInterv1 = 0.20;
        const double xInterv2 = 0.10;
        const double yInterv = 0.10;
        const double zInterv = 0.10;
        //private const double angleInDegrees = 88.5;
        //private const double angle = (Math.PI / 180) * angleInDegrees;
        const double gap = 0.02;
        const double offset = 0.50;

        //External loads
        const double externalStructuralLoad = -0.60;

        static List<int> loadedStructuralDOFs;
        static double[] externalForcesStructuralVector;

        const double YoungMod = 1.0 * 1e5;

        const double poissonRatio = 0.25;
        const double density = 8000.0;
        const double area = 1.0;
        const double contactArea = thickness * xInterv2;

        //Friction coefficients
        //const double miS = 0.20;
        //const double miD = 0.20;


        private static void CreateStructuralBoundaryConditions()
        {
            List<int> boundedDofs = new List<int>();
            boundedDofs.Add(1);
            boundedDofs.Add(2);
            boundedDofs.Add(3);
            boundedDofs.Add(4);
            boundedDofs.Add(5);
            boundedDofs.Add(6);
            boundedDofs.Add(7);
            boundedDofs.Add(8);
            boundedDofs.Add(9);
            boundedDofs.Add(10);
            boundedDofs.Add(11);
            boundedDofs.Add(12);
            for(int i = 1; i<= 12; i++)
            {
                boundedDofs.Add(180 + i);
            }
            structuralBoundaryConditions = boundedDofs.ToArray<int>();
        }

        private static void CreateStructuralLoadVector()
        {
            loadedStructuralDOFs = new List<int>();
            loadedStructuralDOFs.Add(90);
            loadedStructuralDOFs.Add(96);
            externalForcesStructuralVector = new double[nodesNumber * 3];
        }

        private static Dictionary<int, INode> CreateNodes()
        {

            Dictionary<int, INode> nodes = new Dictionary<int, INode>();
            int l;
            l = 1;
            //First beam
            for (int i = 0; i < nodesInX; i++)
            {
                for (int j = 0; j < nodesInY; j++)
                {
                    for (int k = 0; k < nodesInZ; k++)
                    {
                        nodes[l] = new Node(i * xInterv2, j * yInterv, k * zInterv);
                        l += 1;
                    }
                }
            }
            //Second beam
            for (int i = 0; i < nodesInX; i++)
            {
                for (int j = 0; j < nodesInY; j++)
                {
                    for (int k = 0; k < nodesInZ; k++)
                    {
                        nodes[l] = new Node(i * xInterv2 + offset, j * yInterv, (k - 1) * zInterv - gap);
                        l += 1;
                    }
                }
            }
            return nodes;
        }

        private static Dictionary<int, Dictionary<int, int>> CreateConnectivity()
        {

            Dictionary<int, Dictionary<int, int>> connectivity = new Dictionary<int, Dictionary<int, int>>();
            int l = 1;
            for (int i = 1; i <= nodesInX - 1; i++)
            {
                for (int j = 1; j <= nodesInY - 1; j++)
                {
                    for (int k = 1; k <= nodesInZ - 1; k++)
                    {
                        int first = (i - 1) * nodesInY * nodesInZ + (j - 1) * nodesInZ + k;
                        connectivity[l] = new Dictionary<int, int>() { { 1, first }, { 2, first + nodesInY * nodesInZ },
                            { 3, first + nodesInY * nodesInZ + nodesInZ }, { 4, first + nodesInZ },
                            { 5, first + 1 },{ 6, first + nodesInY * nodesInZ + 1 },{ 7, first + (nodesInY + 1) * nodesInZ + 1 },
                            { 8, first + 1 + nodesInZ } };
                        l += 1;
                    }
                }
            }
            for (int i = 1; i <= nodesInX - 1; i++)
            {
                for (int j = 1; j <= nodesInY - 1; j++)
                {
                    for (int k = 1; k <= nodesInZ - 1; k++)
                    {
                        int first = (i - 1) * nodesInY * nodesInZ + (j - 1) * nodesInZ + k +
                            nodesInX * nodesInY * nodesInZ;
                        connectivity[l] = new Dictionary<int, int>() { { 1, first }, { 2, first + nodesInY * nodesInZ },
                            { 3, first + nodesInY * nodesInZ + nodesInZ }, { 4, first + nodesInZ },
                            { 5, first + 1 },{ 6, first + nodesInY * nodesInZ + 1 },{ 7, first + (nodesInY + 1) * nodesInZ + 1 },
                            { 8, first + 1 + nodesInZ } };
                        l += 1;
                    }
                }
            }
            //Contact elements
            connectivity[l] = new Dictionary<int, int>() { { 1, 32 + 2 }, { 2, 32 + 6 }, { 3, 32 + 8 },
                                                             {4, 32 + 4 },
                                                             { 5, 20 + 1 }, { 6, 20 + 5 }, { 7, 20 + 7 },
                                                             {8, 20 + 3 } };
            l += 1;
            connectivity[l] = new Dictionary<int, int>() { { 1, 36 + 2 }, { 2, 36 + 6 }, { 3, 36 + 8 },
                                                             {4, 36 + 4 },
                                                             { 5, 20 + 1 }, { 6, 20 + 5 }, { 7, 20 + 7 },
                                                             {8, 20 + 3 } };
            l += 1;
            connectivity[l] = new Dictionary<int, int>() { { 1, 32 + 2 }, { 2, 32 + 6 }, { 3, 32 + 8 },
                                                             {4, 32 + 4 },
                                                             { 5, 24 + 1 }, { 6, 24 + 5 }, { 7, 24 + 7 },
                                                             {8, 24 + 3 } };
            l += 1;
            connectivity[l] = new Dictionary<int, int>() { { 1, 36 + 2 }, { 2, 36 + 6 }, { 3, 36 + 8 },
                                                             {4, 36 + 4 },
                                                             { 5, 24 + 1 }, { 6, 24 + 5 }, { 7, 24 + 7 },
                                                             {8, 24 + 3 } };
            return connectivity;
        }

        private static Dictionary<int, bool[]> CreateNodeFAT()
        {
            Dictionary<int, bool[]> nodeFAT = new Dictionary<int, bool[]>();
            for (int i = 1; i <= nodesNumber; i++)
            {
                nodeFAT[i] = new bool[] { true, true, true, false, false, false };
            }
            return nodeFAT;
        }
        private static Dictionary<int, IElementProperties> CreateElementProperties()
        {
            double E = YoungMod;

            double A = area;
            //double CA = contactArea;

            string type = "Hex8";
            string type3 = "ContactStS3D";
            //string type3 = "ContactNtS3D";
            Dictionary<int, IElementProperties> elementProperties = new Dictionary<int, IElementProperties>();
            for (int i = 1; i <= elementsNumber; i++)
            {
                elementProperties[i] = new ElementProperties(E, poissonRatio, A, thickness, density, type);
            }
            for (int i = elementsNumber + 1; i <= elementsNumber + contactElements; i++)
            {
                //elementProperties[i] = new ElementProperties(E, type3, 5.0, 9, 1, 1, 5.0, miS, miD);
                elementProperties[i] = new ElementProperties(E, A, type3, 5.0, 9, 1, 1);
                //elementProperties[i] = new ElementProperties(E, A, type3);
                elementProperties[i].Density = density;
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
            //
            elementsAssembly.SeperateContactDoF();
            //
            double[,] globalStiffnessMatrix = elementsAssembly.CreateTotalStiffnessMatrix();
            //int countContactElements = elementsAssembly.CountElementsOfSameType(typeof(ContactStS2D));
            //ShowToGUI.PlotInitialGeometry(elementsAssembly);
            //
            //structuralSolution.LinearScheme = new LUFactorization();
            //
            structuralSolution.LinearScheme = new MMCPCGSolver();
            //
            structuralSolution.NonLinearScheme.Tolerance = 1e-5;
            structuralSolution.ActivateNonLinearSolver = true;
            structuralSolution.NonLinearScheme.numberOfLoadSteps = 200;

            double[] externalForces3 = externalForcesStructuralVector;
            foreach (var dof in loadedStructuralDOFs)
            {
                externalForces3[dof - 1] = externalStructuralLoad;
            }
            double[] reducedExternalForces3 = elementsAssembly.MMCPGCreateReducedFromFullVector(externalForces3);
            structuralSolution.AssemblyData = elementsAssembly;
            structuralSolution.Solve(reducedExternalForces3);
            double[] solvector3 = structuralSolution.GetSolution();
            Dictionary<int, double[]> allStepsSolutions = structuralSolution.GetAllStepsSolutions();
            elementsAssembly.UpdateDisplacements(solvector3);
            double[] fullSolVector3 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(solvector3, elementsAssembly.BoundedDOFsVector);

            Dictionary<int, double[]> allStepsFullSolutions = new Dictionary<int, double[]>();
            for (int i = 0; i < allStepsSolutions.Count; i++)
            {
                allStepsFullSolutions.Add(i + 1, BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions.Single(m => m.Key == i + 1).Value, elementsAssembly.BoundedDOFsVector));
                int j = i + 1;
                string name = "solution" + j.ToString() + ".dat";
                VectorOperations.PrintVectorToFile(allStepsFullSolutions.Single(m => m.Key == i + 1).Value, @"C:\Users\Public\Documents\" + name);
            }
            List<double[]> structuralSolutions = new List<double[]>();

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
