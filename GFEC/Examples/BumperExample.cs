using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading;

namespace GFEC
{
    public static class BumperExample
    {
        public static ISolver structuralSolution;
        static int[] structuralBoundaryConditions;
        //public static  Dictionary<int, INode> nodes;
        //public static Dictionary<int, Dictionary<int, int>> elementsConnectivity;

        static List<int> loadedStructuralDOFs;
        static double[] externalForcesStructuralVector;
        const double externalStructuralLoad = 100.0;
        const int nodesNumber = 3726;
        const int elmntsNumber = 1760;
        private static void CreateStructuralBoundaryConditions()
        {
            List<int> boundedDofs = new List<int>();
            int initialNode = 1169;
            for (int node = initialNode; node <= initialNode + 2; node++)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            initialNode = 1199;
            for (int node = initialNode; node <= initialNode + 2; node++)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            initialNode = 1249;
            for (int node = initialNode; node <= initialNode + 2; node++)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            initialNode = 1279;
            for (int node = initialNode; node <= initialNode + 2; node++)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            initialNode = 1328;
            for (int node = initialNode; node <= initialNode + 2; node++)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            initialNode = 1358;
            for (int node = initialNode; node <= initialNode + 2; node++)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            initialNode = 1807;
            for (int node = initialNode; node <= initialNode + 2; node++)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            initialNode = 1837;
            for (int node = initialNode; node <= initialNode + 2; node++)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            initialNode = 1917;
            for (int node = initialNode; node <= initialNode + 2; node++)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            initialNode = 1947;
            for (int node = initialNode; node <= initialNode + 2; node++)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            initialNode = 1997;
            for (int node = initialNode; node <= initialNode + 2; node++)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            initialNode = 2027;
            for (int node = initialNode; node <= initialNode + 2; node++)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            initialNode = 2076;
            for (int node = initialNode; node <= initialNode + 2; node++)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            initialNode = 2106;
            for (int node = initialNode; node <= initialNode + 2; node++)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            initialNode = 2156;
            for (int node = initialNode; node <= initialNode + 2; node++)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            initialNode = 2186;
            for (int node = initialNode; node <= initialNode + 2; node++)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            initialNode = 2236;
            for (int node = initialNode; node <= initialNode + 2; node++)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            initialNode = 2266;
            for (int node = initialNode; node <= initialNode + 2; node++)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            initialNode = 2316;
            for (int node = initialNode; node <= initialNode + 2; node++)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            initialNode = 2346;
            for (int node = initialNode; node <= initialNode + 2; node++)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            initialNode = 2395;
            for (int node = initialNode; node <= initialNode + 2; node++)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            initialNode = 2425;
            for (int node = initialNode; node <= initialNode + 2; node++)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            initialNode = 2475;
            for (int node = initialNode; node <= initialNode + 2; node++)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            initialNode = 2505;
            for (int node = initialNode; node <= initialNode + 2; node++)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            initialNode = 2555;
            for (int node = initialNode; node <= initialNode + 2; node++)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            initialNode = 2585;
            for (int node = initialNode; node <= initialNode + 2; node++)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            initialNode = 2635;
            for (int node = initialNode; node <= initialNode + 2; node++)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            initialNode = 2665;
            for (int node = initialNode; node <= initialNode + 2; node++)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            initialNode = 2715;
            for (int node = initialNode; node <= initialNode + 2; node++)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            initialNode = 2745;
            for (int node = initialNode; node <= initialNode + 2; node++)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            initialNode = 2797;
            for (int node = initialNode; node <= initialNode + 2; node++)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            initialNode = 2827;
            for (int node = initialNode; node <= initialNode + 2; node++)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            initialNode = 2877;
            for (int node = initialNode; node <= initialNode + 2; node++)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            initialNode = 2907;
            for (int node = initialNode; node <= initialNode + 2; node++)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            initialNode = 2960;
            for (int node = initialNode; node <= initialNode + 2; node++)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            initialNode = 2990;
            for (int node = initialNode; node <= initialNode + 2; node++)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            initialNode = 3039;
            for (int node = initialNode; node <= initialNode + 2; node++)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            initialNode = 3069;
            for (int node = initialNode; node <= initialNode + 2; node++)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            initialNode = 3197;
            for (int node = initialNode; node <= initialNode + 2; node++)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            initialNode = 3227;
            for (int node = initialNode; node <= initialNode + 2; node++)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            initialNode = 3276;
            for (int node = initialNode; node <= initialNode + 2; node++)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            initialNode = 3306;
            for (int node = initialNode; node <= initialNode + 2; node++)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            initialNode = 3513;
            for (int node = initialNode; node <= initialNode + 2; node++)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            initialNode = 3543;
            for (int node = initialNode; node <= initialNode + 2; node++)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            initialNode = 3592;
            for (int node = initialNode; node <= initialNode + 2; node++)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            initialNode = 3622;
            for (int node = initialNode; node <= initialNode + 2; node++)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            initialNode = 3671;
            for (int node = initialNode; node <= initialNode + 2; node++)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            initialNode = 3701;
            for (int node = initialNode; node <= initialNode + 2; node++)
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

        private static void CreateStructuralLoadVector()
        {
            loadedStructuralDOFs = new List<int>();
            for (int i = 744; i <= 986; i ++)
            {
                loadedStructuralDOFs.Add(3 * i);
            }
            externalForcesStructuralVector = new double[3726 * 3];
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
            string type = "ANSSolidShell8EAS";
            Dictionary<int, IElementProperties> elementProperties = new Dictionary<int, IElementProperties>();
            int totalElements = elementsConnectivity.Count;
            for (int i = 1; i <= totalElements; i++)
            {
                elementProperties[i] = new ElementProperties(E, 0.28, 1.0, 0.005, 100, type);
            }
            return elementProperties;
        }
        

        private static IAssembly CreateAssembly(Dictionary<int, INode> nodes, Dictionary<int, Dictionary<int, int>> elementsConnectivity)
        {
            IAssembly assembly = new Assembly();
            assembly.Nodes = nodes;
            assembly.ElementsConnectivity = elementsConnectivity;
            assembly.ElementsProperties = CreateElementProperties(elementsConnectivity);
            assembly.NodeFreedomAllocationList = CreateNodeFAT(nodes);
            CreateStructuralBoundaryConditions();
            CreateStructuralLoadVector();
            assembly.BoundedDOFsVector = structuralBoundaryConditions;
            return assembly;
        }
        private static IAssembly CreateAssembly()
        {
            IAssembly assembly = new Assembly();
            //assembly.ElementsProperties = CreateElementProperties(elementsConnectivity);
            //assembly.NodeFreedomAllocationList = CreateNodeFAT();
            CreateStructuralBoundaryConditions();
            CreateStructuralLoadVector();
            //assembly.BoundedDOFsVector = structuralBoundaryConditions;
            return assembly;
        }

        public static Results RunStaticExample(Dictionary<int, INode> nodes, Dictionary<int, Dictionary<int, int>> elementsConnectivity)
        {
            #region Structural
            IAssembly elementsAssembly = CreateAssembly(nodes, elementsConnectivity);
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

            structuralSolution.LinearScheme = new Skyline();
            //structuralSolution.NonLinearScheme.Tolerance = 1e-4;
            structuralSolution.ActivateNonLinearSolver = false;
            //structuralSolution.NonLinearScheme.numberOfLoadSteps = 40;
            double[] externalForces3 = externalForcesStructuralVector;
            foreach (var dof in loadedStructuralDOFs)
            {
                externalForces3[dof - 1] =  externalStructuralLoad;
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