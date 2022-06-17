using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading;

namespace GFEC
{
    public static class SolidShellThinCylinderConnectivity2
    {
        public static ISolver structuralSolution;
        static int[] structuralBoundaryConditions;
        static double[] externalForcesStructuralVector;
        const double externalStructuralLoad = -250000.0;
        const int nodesNumber = 3066;
        const int elmntsNumber = 1440;
        const double radius1 = 0.5;
        const double radius2 = 0.49;
        const double initialTheta = 0.0;
        const double deltaTheta = 2.5 * Math.PI / 180.0;
        const double zInterv = 0.05;
        private static Dictionary<int, INode> CreateNodes()
        {
            Dictionary<int, INode> nodes = new Dictionary<int, INode>();
            int l;
            l = 1;
            for (int k = 0; k < 21; k++)
            {
                for (int i = 0; i < 73; i++)
                {
                    nodes[l] = new Node(radius1 * Math.Cos(initialTheta + deltaTheta * i),
                        radius1 * Math.Sin(initialTheta + deltaTheta * i),
                        k * zInterv);
                    l += 1;
                }
            }
            for (int k = 0; k < 21; k++)
            {
                for (int i = 0; i < 73; i++)
                {
                    nodes[l] = new Node(radius2 * Math.Cos(initialTheta + deltaTheta * i),
                        radius2 * Math.Sin(initialTheta + deltaTheta * i),
                        k * zInterv);
                    l += 1;
                }
            }
            return nodes;
        }
        private static void CreateStructuralBoundaryConditions()
        {
            int addedNodes = 1533;
            List<int> boundedDofs = new List<int>();
            for (int node = 1; node <= 73; node++)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            for (int node = 74; node <= 1461; node+=73)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            for (int node = 146; node <= 1533; node+=73)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            for (int node = 1462; node < 1533; node++)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            //--------------------------------------
            for (int node = addedNodes + 1; node <= addedNodes + 73; node++)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            for (int node = addedNodes + 74; node <= addedNodes + 1461; node += 73)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            for (int node = addedNodes + 146; node <= addedNodes + 1533; node += 73)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            for (int node = addedNodes + 1462; node < addedNodes + 1533; node++)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            boundedDofs = boundedDofs.Distinct().ToList().OrderBy(x => x).ToList();
            structuralBoundaryConditions = boundedDofs.ToArray<int>();
        }
        private static Dictionary<int, Dictionary<int, int>> CreateConnectivity()
        {
            Dictionary<int, Dictionary<int, int>> connectivity = new Dictionary<int, Dictionary<int, int>>();
            int addedNodes = 1533;
            int l = 1;
            for (int i = 73; i >= 2; i-=1)
            {
                for (int j = 1; j <= 21 - 1; j++)
                {
                    int firstUpp = i + (j - 1) * 73;
                    connectivity[l] = new Dictionary<int, int>() {
                            { 1, firstUpp + 73 },
                            { 2, firstUpp + 72 },
                            { 3,  firstUpp - 1},
                            { 4,  firstUpp},
                            { 5, firstUpp + addedNodes + 73 },
                            { 6, firstUpp + 72 + addedNodes },
                            { 7,  firstUpp + addedNodes - 1},
                            { 8,  firstUpp + addedNodes}
                        };
                    l += 1;
                }
            }
            return connectivity;
        }
        private static void CreateStructuralLoadVector()
        {
            externalForcesStructuralVector = new double[nodesNumber * 3];
        }
        private static Dictionary<int, bool[]> CreateNodeFAT(Dictionary<int, INode> nodes)
        {
            int totalNodes = nodes.Count;
            Dictionary<int, bool[]> nodeFAT = new Dictionary<int, bool[]>();
            for (int i = 1; i <= totalNodes; i++)
            {
                nodeFAT[i] = new bool[] { true, true, true, false, false, false };
            }
            return nodeFAT;
        }

        private static Dictionary<int, IElementProperties> CreateElementProperties(Dictionary<int, Dictionary<int, int>> elementsConnectivity)
        {
            double E = 200.0 * 1e9;
            string type = "ANSSolidShell8LEAS1RI";
            Dictionary<int, IElementProperties> elementProperties = new Dictionary<int, IElementProperties>();
            int totalElements = elementsConnectivity.Count;
            for (int i = 1; i <= totalElements; i++)
            {
                elementProperties[i] = new ElementProperties(E, 0.28, 1.0, 0.005, 100, type);
            }
            return elementProperties;
        }


        private static IAssembly CreateAssembly()
        {
            IAssembly assembly = new Assembly();
            assembly.Nodes = CreateNodes();
            assembly.ElementsConnectivity = CreateConnectivity();
            assembly.ElementsProperties = CreateElementProperties(assembly.ElementsConnectivity);
            assembly.NodeFreedomAllocationList = CreateNodeFAT(assembly.Nodes);
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
            ExportToFile.ExportMatlabInitialGeometry(elementsAssembly);
            elementsAssembly.ActivateBoundaryConditions = true;
            double[,] globalStiffnessMatrix = elementsAssembly.CreateTotalStiffnessMatrix();
            if (!globalStiffnessMatrix.Cast<double>().Any(d => double.IsNaN(d) || double.IsInfinity(d)))
            {
                bool noInfiniteValues = true;
            }
            else
            {
                bool noInfiniteValues = false;
            }
            int count = 0;
            List<int> Indices1 = new List<int>();
            List<int> Indices2 = new List<int>();

            for (int i = 0; i < globalStiffnessMatrix.GetLength(0); i++)
            {
                for (int j = 0; j < globalStiffnessMatrix.GetLength(1); j++)
                {
                    if (double.IsNaN(globalStiffnessMatrix[i, j]) || double.IsInfinity(globalStiffnessMatrix[i, j]))
                    {
                        count += 1;
                        Indices1.Add(i);
                        Indices2.Add(j);

                    }
                }
            }
            List<int> noDupes1 = Indices1.Distinct().ToList().OrderBy(x => x).ToList();
            List<int> noDupes2 = Indices2.Distinct().ToList().OrderBy(x => x).ToList();
            for (int i = 0; i < noDupes1.Count; i++)
            {
                noDupes1[i] = (noDupes1[i]) / 3;
            }
            List<int> noDupes = noDupes1.Distinct().ToList();
            int[] indices = noDupes.ToArray();
            int[] indices12 = noDupes2.ToArray();
            structuralSolution.LinearScheme = new LUFactorization();
            structuralSolution.ActivateNonLinearSolver = false;
            double[] externalForcesVector = externalForcesStructuralVector;
            externalForcesVector[767 * 3 - 2] = externalStructuralLoad;
            double[] reducedExternalForces3 = BoundaryConditionsImposition.ReducedVector(externalForcesVector, elementsAssembly.BoundedDOFsVector);
            structuralSolution.AssemblyData = elementsAssembly;
            structuralSolution.Solve(reducedExternalForces3);
            double[] solvector = structuralSolution.GetSolution();
            elementsAssembly.UpdateDisplacements(solvector);
            //ShowToGUI.PlotFinalGeometry(elementsAssembly);
            double[] fullSolVector = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(solvector, elementsAssembly.BoundedDOFsVector);
            VectorOperations.PrintVectorToFile(fullSolVector, @"C:\Users\Public\Documents\" + "SolidShellElementsThinCylinderLinearSolution.dat");
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

