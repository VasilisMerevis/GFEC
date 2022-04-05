using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GFEC
{
    class DegenerateShellElementsLinearExample
    {
        public static ISolver structuralSolution;
        static int[] structuralBoundaryConditions;
        const double thickness = 0.006;
        const int nodesInX = 21;
        const int nodesInY = 21;
        const int nodesInZ = 2;
        const int nodesNumber = 882;
        const int elementsNumber = 100;

        //const double xInterv1 = 0.20;
        const double xInterv = 0.05;
        const double yInterv = 0.05;

        //External loads
        const double externalStructuralLoad = -20000.0;

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
                boundedDofs.Add(i * 5 + 3);
            }
            boundedDofs.Add(101);
            boundedDofs.Add(102);
            boundedDofs.Add(103);
            for (int i = 1; i < nodesInX - 1; i++)
            {
                boundedDofs.Add(i * nodesInY * 5 + 3);
                boundedDofs.Add((i * nodesInY + nodesInY - 1) * 5 + 3);
            }
            boundedDofs.Add(2101);
            boundedDofs.Add(2102);
            for (int i = 0; i < nodesInY; i++)
            {
                boundedDofs.Add((nodesInX - 1) * nodesInY * 5 + i * 5 + 3);
            }
            boundedDofs.Add(2201);
            boundedDofs.Add(2202);
            //boundedDofs.Add(2053);
            structuralBoundaryConditions = boundedDofs.ToArray<int>();
        }
        private static void CreateStructuralLoadVector()
        {
            loadedStructuralDOFs = new List<int>();
            for (int j = 1; j < 20; j++)
            {
                for (int i = 1; i < 20; i++)
                {
                    loadedStructuralDOFs.Add(j * 21 * 5 + i * 5 + 3);
                }
            }
            externalForcesStructuralVector = new double[nodesNumber/2 * 5];
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
                    nodes[l] = new Node(i * xInterv, j * yInterv, -thickness / 2.0);
                    l += 1;
                }
            }
            for (int i = 0; i < nodesInX; i++)
            {
                for (int j = 0; j < nodesInY; j++)
                {
                    nodes[l] = new Node(i * xInterv, j * yInterv, thickness / 2.0);
                    l += 1;
                }
            }
            return nodes;
        }
        private static Dictionary<int, Dictionary<int, int>> CreateConnectivity()
        {

            Dictionary<int, Dictionary<int, int>> connectivity = new Dictionary<int, Dictionary<int, int>>();
            int l = 1;
            for (int i = 1; i <= nodesInX - 2; i+=2)
            {
                for (int j = 1; j <= nodesInY - 2; j+=2)
                {
                    int first = (i - 1) * nodesInY + j;
                    int first2 = first + nodesNumber / 2;
                    connectivity[l] = new Dictionary<int, int>() { { 1, first }, { 2, first + nodesInY },
                            { 3,  first + nodesInY + nodesInY}, { 4,  first + nodesInY + nodesInY + 1},
                            { 5,  first + nodesInY + nodesInY + 2}, { 6,  first + nodesInY + 2},
                            { 7,  first + 2}, { 8,  first + 1},
                            { 9,  first + nodesInY + 1},
                            { 10, first2 }, { 11, first2 + nodesInY },
                            { 12,  first2 + nodesInY + nodesInY}, { 13,  first2 + nodesInY + nodesInY + 1},
                            { 14,  first2 + nodesInY + nodesInY + 2}, { 15,  first2 + nodesInY + 2},
                            { 16,  first2 + 2}, { 17,  first2 + 1},
                            { 18,  first2 + nodesInY + 1}};
                    l += 1;
                }
            }
            return connectivity;
        }

        private static Dictionary<int, bool[]> CreateNodeFAT()
        {
            Dictionary<int, bool[]> nodeFAT = new Dictionary<int, bool[]>();
            for (int i = 1; i <= nodesNumber/2; i++)
            {
                nodeFAT[i] = new bool[] { true, true, true, true, true, false };
            }
            return nodeFAT;
        }
        private static Dictionary<int, IElementProperties> CreateElementProperties()
        {
            double E = YoungMod;
            double A = area;
            string type = "IsoparamShell18";

            Dictionary<int, IElementProperties> elementProperties = new Dictionary<int, IElementProperties>();
            for (int i = 1; i <= elementsNumber; i++)
            {
                elementProperties[i] = new ElementProperties(E, poissonRatio, A, thickness, density, type);
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
            structuralSolution.LinearScheme = new LUFactorization();
            structuralSolution.ActivateNonLinearSolver = true;
            structuralSolution.NonLinearScheme.Tolerance = 1e-5;//
            structuralSolution.NonLinearScheme.numberOfLoadSteps = 10;//

            double[] externalForces3 = externalForcesStructuralVector;
            //foreach (var dof in loadedStructuralDOFs)
            //{
            //    externalForces3[dof - 1] = externalStructuralLoad;
            //}
            int count = 1;
            int count2 = 1;
            foreach (var dof in loadedStructuralDOFs)
            {
                if (count2 % 2 != 0)
                {
                    if (count % 2 != 0)
                    {
                        externalForces3[dof - 1] = 16.0 * externalStructuralLoad * 4.0 *  xInterv * yInterv / 36.0;
                        count += 1;
                    }
                    else
                    {
                        externalForces3[dof - 1] = 8.0 * externalStructuralLoad * 4.0 * xInterv * yInterv / 36.0;
                        count += 1;
                    }
                    if (count == 20)
                    {
                        count = 1;
                        count2 += 1;
                    }
                }
                else
                {
                    if (count % 2 != 0)
                    {
                        externalForces3[dof - 1] = 8.0 * externalStructuralLoad * 4.0 * xInterv * yInterv / 36.0;
                        count += 1;
                    }
                    else
                    {
                        externalForces3[dof - 1] = 4.0 * externalStructuralLoad * 4.0 * xInterv * yInterv / 36.0;
                        count += 1;
                    }
                    if (count == 20)
                    {
                        count = 1;
                        count2 += 1;
                    }
                }
            }
            double[] reducedExternalForces3 = BoundaryConditionsImposition.ReducedVector(externalForces3, elementsAssembly.BoundedDOFsVector);
            structuralSolution.AssemblyData = elementsAssembly;
            structuralSolution.Solve(reducedExternalForces3);
            Dictionary<int, double[]> solvectors = structuralSolution.GetAllStepsSolutions();

            //double[]  = structuralSolution.GetSolution();
            for(int i = solvectors.Keys.Min(); i<=solvectors.Keys.Max(); i++)
            {
                double[] fullSolVector = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(solvectors.Single(s=>s.Key == i).Value,
                elementsAssembly.BoundedDOFsVector);
                string name = "DegenElementsNLsolution" + i.ToString() + ".dat";
                VectorOperations.PrintVectorToFile(fullSolVector, @"C:\Users\Public\Documents\" + name);
            }
            //double[] fullSolVector = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(solvector,
            //    elementsAssembly.BoundedDOFsVector);

            //string name = "DegenElementssolution" + ".dat";
            //VectorOperations.PrintVectorToFile(fullSolVector, @"C:\Users\Public\Documents\" + name);
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