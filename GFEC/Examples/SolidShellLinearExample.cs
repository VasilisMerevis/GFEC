using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GFEC
{
    class SolidShellLinearExample
    {
        public static ISolver structuralSolution;
        static int[] structuralBoundaryConditions;
        const double shellThickness = 0.006;
        const int nodesInX = 21;
        const int nodesInY = 21;
        const int nodesInZ = 2;
        const int nodesNumber = 882;
        const int elementsNumber = 400;

        //const double xInterv1 = 0.20;
        const double xInterv = 0.05;
        const double yInterv = 0.05;

        //External loads
        const double externalStructuralLoad = -30000.0;

        static List<int> loadedStructuralDOFs;
        static double[] externalForcesStructuralVector;

        const double YoungMod = 200.0 * 1e9;

        const double poissonRatio = 0.25;
        const double density = 8000.0;
        const double area = 1.0;
        private static void CreateStructuralBoundaryConditions()
        {
            List<int> boundedDofs = new List<int>();
            boundedDofs.Add(1);
            boundedDofs.Add(2);
            for (int i = 0; i < nodesInY - 1; i++)
            {
                boundedDofs.Add(i * 3 + 3);
            }
            boundedDofs.Add(61);
            boundedDofs.Add(62);
            boundedDofs.Add(63);
            for (int i = 1; i < nodesInX - 1; i++)
            {
                boundedDofs.Add(i * nodesInY * 3 + 3);
                boundedDofs.Add((i * nodesInY + nodesInY - 1) * 3 + 3);
            }
            boundedDofs.Add(1261);
            boundedDofs.Add(1262);
            for (int i = 0; i < nodesInY; i++)
            {
                boundedDofs.Add((nodesInX - 1) * nodesInY * 3 + i * 3 + 3);
            }
            boundedDofs.Add(1321);
            boundedDofs.Add(1322);
            //--------------------------------------------------
            int add = nodesInX * nodesInY * 3;
            boundedDofs.Add(add + 1);
            boundedDofs.Add(add + 2);
            for (int i = 0; i < nodesInY - 1; i++)
            {
                boundedDofs.Add(add + i * 3 + 3);
            }
            boundedDofs.Add(add + 61);
            boundedDofs.Add(add + 62);
            boundedDofs.Add(add + 63);
            for (int i = 1; i < nodesInX - 1; i++)
            {
                boundedDofs.Add(add + i * nodesInY * 3 + 3);
                boundedDofs.Add(add + (i * nodesInY + nodesInY - 1) * 3 + 3);
            }
            boundedDofs.Add(add + 1261);
            boundedDofs.Add(add + 1262);
            for (int i = 0; i < nodesInY; i++)
            {
                boundedDofs.Add(add + (nodesInX - 1) * nodesInY * 3 + i * 3 + 3);
            }
            boundedDofs.Add(add + 1321);
            boundedDofs.Add(add + 1322);
            structuralBoundaryConditions = boundedDofs.ToArray<int>();
        }
        private static void CreateStructuralLoadVector()
        {
            int add = nodesInX * nodesInY * 3;
            loadedStructuralDOFs = new List<int>();
            for (int j = 1; j < 20; j++)
            {
                for (int i = 1; i < 20; i++)
                {
                    loadedStructuralDOFs.Add(add + j * 21 * 3 + i * 3 + 3);
                }
            }
            externalForcesStructuralVector = new double[nodesNumber * 3];
        }

        private static Dictionary<int, INode> CreateNodes()
        {

            Dictionary<int, INode> nodes = new Dictionary<int, INode>();
            int l;
            l = 1;
            for (int i = 0; i < nodesInX; i++)
            {
                for (int j = 0; j < nodesInY; j++)
                {
                    nodes[l] = new Node(i * xInterv, j * yInterv, -shellThickness / 2.0);
                    l += 1;
                }
            }
            for (int i = 0; i < nodesInX; i++)
            {
                for (int j = 0; j < nodesInY; j++)
                {
                    nodes[l] = new Node(i * xInterv, j * yInterv, shellThickness / 2.0);
                    l += 1;
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
                    int first = (i - 1) * nodesInY + j;
                    int first2 = first + nodesNumber / 2;
                    connectivity[l] = new Dictionary<int, int>()
                    {
                            { 1,  first2 },
                            { 2,  first2 + nodesInY },
                            { 3,  first2 + nodesInY + 1},
                            { 4,  first2 + 1},
                            { 5,  first },
                            { 6,  first + nodesInY},
                            { 7,  first + nodesInY + 1},
                            { 8,  first + 1}
                    };
                    l += 1;
                }
            }
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
            string type = "ANSSolidShell8EAS";

            Dictionary<int, IElementProperties> elementProperties = new Dictionary<int, IElementProperties>();
            for (int i = 1; i <= elementsNumber; i++)
            {
                elementProperties[i] = new ElementProperties(E, poissonRatio, A, shellThickness, density, type);
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
            ExportToFile.ExportMatlabInitialGeometry(elementsAssembly);
            double[,] globalStiffnessMatrix = elementsAssembly.CreateTotalStiffnessMatrix();
            structuralSolution.LinearScheme = new Skyline();
            structuralSolution.ActivateNonLinearSolver = false;
            //structuralSolution.NonLinearScheme.Tolerance = 1e-5;//
            //structuralSolution.NonLinearScheme.MaxIterations = 100;//
            //structuralSolution.NonLinearScheme.numberOfLoadSteps = 6;//

            double[] externalForces3 = externalForcesStructuralVector;
            foreach (var dof in loadedStructuralDOFs)
            {
                externalForces3[dof - 1] = externalStructuralLoad * xInterv * yInterv;

            }
            double[] reducedExternalForces3 = BoundaryConditionsImposition.ReducedVector(externalForces3, elementsAssembly.BoundedDOFsVector);
            structuralSolution.AssemblyData = elementsAssembly;
            structuralSolution.Solve(reducedExternalForces3);
            //---------------------------------------------------
            //Case Linear Solver
            double[] fullSolVector = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(structuralSolution.GetSolution(),
                elementsAssembly.BoundedDOFsVector);
            string name = "SolidShellLsolution" + ".dat";
            VectorOperations.PrintVectorToFile(fullSolVector, @"C:\Users\Public\Documents\" + name);
            //---------------------------------------------------
            ////Case NonLinear Solver
            //Dictionary<int, double[]> solvectors = structuralSolution.GetAllStepsSolutions();
            //for (int i = solvectors.Keys.Min(); i <= solvectors.Keys.Max(); i++)
            //{
            //    double[] fullSolVector = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(solvectors.Single(s => s.Key == i).Value,
            //    elementsAssembly.BoundedDOFsVector);
            //    string name = "SolidShellLsolution" + i.ToString() + ".dat";
            //    VectorOperations.PrintVectorToFile(fullSolVector, @"C:\Users\Public\Documents\" + name);
            //}
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