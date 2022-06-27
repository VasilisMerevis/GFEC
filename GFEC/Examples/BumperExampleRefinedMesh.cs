using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading;

namespace GFEC
{
    public static class BumperExampleRefinedMesh
    {
        public static ISolver structuralSolution;
        static int[] structuralBoundaryConditions;
        //public static  Dictionary<int, INode> nodes;
        //public static Dictionary<int, Dictionary<int, int>> elementsConnectivity;

        static List<int> loadedStructuralDOFs;
        static double[] externalForcesStructuralVector;
        const double externalStructuralLoad = 100.0;
        const int nodesNumber = 9315;
        const int elmntsNumber = 7040;
        private static void CreateStructuralBoundaryConditions(List<int> fixedNodes)
        {
            List<int> boundedDofs = new List<int>();
            foreach(var node in fixedNodes)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            bool duplicates = new bool();
            if (boundedDofs.Count != boundedDofs.Distinct().Count())
            {
                duplicates = true;
                boundedDofs = boundedDofs.Distinct().ToList();
            }
            structuralBoundaryConditions = boundedDofs.ToArray<int>();
        }

        private static void CreateStructuralLoadVector(List<int> loadedNodes)
        {
            loadedStructuralDOFs = new List<int>();
            foreach (var node in loadedNodes)
            {
                loadedStructuralDOFs.Add(3 * node);
            }
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


        private static IAssembly CreateAssembly(Dictionary<int, INode> nodes, Dictionary<int, Dictionary<int, int>> elementsConnectivity,
            List<int> fixedNodes, List<int> loadedNodes)
        {
            IAssembly assembly = new Assembly();
            assembly.Nodes = nodes;
            assembly.ElementsConnectivity = elementsConnectivity;
            assembly.ElementsProperties = CreateElementProperties(elementsConnectivity);
            assembly.NodeFreedomAllocationList = CreateNodeFAT(nodes);
            CreateStructuralBoundaryConditions(fixedNodes);
            CreateStructuralLoadVector(loadedNodes);
            assembly.BoundedDOFsVector = structuralBoundaryConditions;
            return assembly;
        }
        private static IAssembly CreateAssembly()
        {
            IAssembly assembly = new Assembly();
            //assembly.ElementsProperties = CreateElementProperties(elementsConnectivity);
            //assembly.NodeFreedomAllocationList = CreateNodeFAT();
            //CreateStructuralBoundaryConditions();
            //CreateStructuralLoadVector();
            //assembly.BoundedDOFsVector = structuralBoundaryConditions;
            return assembly;
        }

        public static Results RunStaticExample(Dictionary<int, INode> nodes, Dictionary<int, Dictionary<int, int>> elementsConnectivity,
            List<int> fixedNodes, List<int> loadedNodes)
        {
            #region Structural
            IAssembly elementsAssembly = CreateAssembly(nodes, elementsConnectivity, fixedNodes, loadedNodes);
            elementsAssembly.CreateElementsAssembly();
            ExportToFile.ExportMatlabInitialGeometry(elementsAssembly);
            elementsAssembly.ActivateBoundaryConditions = true;
            //double[,] globalStiffnessMatrix = elementsAssembly.CreateTotalStiffnessMatrix();
            //if (!globalStiffnessMatrix.Cast<double>().Any(d => double.IsNaN(d) || double.IsInfinity(d)))
            //{
            //    bool noInfiniteValues = true;
            //}
            //else
            //{
            //    bool noInfiniteValues = false;
            //}
            //int count = 0;
            //List<int> Indices1 = new List<int>();
            //List<int> Indices2 = new List<int>();

            //for (int i = 0; i < globalStiffnessMatrix.GetLength(0); i++)
            //{
            //    for (int j = 0; j < globalStiffnessMatrix.GetLength(1); j++)
            //    {
            //        if (double.IsNaN(globalStiffnessMatrix[i, j]) || double.IsInfinity(globalStiffnessMatrix[i, j]))
            //        {
            //            count += 1;
            //            Indices1.Add(i);
            //            Indices2.Add(j);

            //        }
            //    }
            //}
            //List<int> noDupes1 = Indices1.Distinct().ToList().OrderBy(x => x).ToList();
            //List<int> noDupes2 = Indices2.Distinct().ToList().OrderBy(x => x).ToList();
            //for (int i = 0; i < noDupes1.Count; i++)
            //{
            //    noDupes1[i] = (noDupes1[i]) / 3;
            //}
            //List<int> noDupes = noDupes1.Distinct().ToList();
            //int[] indices = noDupes.ToArray();
            //int[] indices12 = noDupes2.ToArray();
            //VectorOperations.PrintIntVectorToFile(indices, @"C:\Users\Public\Documents\" + "indices.dat");
            //VectorOperations.PrintIntVectorToFile(indices12, @"C:\Users\Public\Documents\" + "indices2.dat");

            structuralSolution.LinearScheme = new CholeskyFactorization();
            //structuralSolution.NonLinearScheme.Tolerance = 1e-4;
            structuralSolution.ActivateNonLinearSolver = false;
            //structuralSolution.NonLinearScheme.numberOfLoadSteps = 40;
            double[] externalForces3 = externalForcesStructuralVector;
            foreach (var dof in loadedStructuralDOFs)
            {
                externalForces3[dof - 1] = externalStructuralLoad;
            }
            double[] reducedExternalForces3 = BoundaryConditionsImposition.ReducedVector(externalForces3, elementsAssembly.BoundedDOFsVector);
            structuralSolution.AssemblyData = elementsAssembly;
            structuralSolution.Solve(reducedExternalForces3);
            double[] solvector = structuralSolution.GetSolution();
            elementsAssembly.UpdateDisplacements(solvector);
            //ShowToGUI.PlotFinalGeometry(elementsAssembly);
            double[] fullSolVector = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(solvector, elementsAssembly.BoundedDOFsVector);
            VectorOperations.PrintVectorToFile(fullSolVector, @"C:\Users\Public\Documents\" + "BumperLinearSolution.dat");
            //Dictionary<int, INode> finalNodes = Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullSolVector);
            //double[] xFinalNodalCoor = Assembly.NodalCoordinatesToVectors(finalNodes).Item1;
            //double[] yFinalNodalCoor = Assembly.NodalCoordinatesToVectors(finalNodes).Item2;
            //Dictionary<int, double[]> allStepsSolutions = structuralSolution.GetAllStepsSolutions();
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

