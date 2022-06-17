using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading;

namespace GFEC
{
    public static class SolidShellElementsThinCylinder
    {
        public static ISolver structuralSolution;
        static int[] structuralBoundaryConditions;
        static double[] externalForcesStructuralVector;
        const double externalStructuralLoad = -250000;
        const int nodesNumber = 374;
        const int elmntsNumber = 160;
        private static void CreateStructuralBoundaryConditions()
        {
            List<int> boundedDofs = new List<int>();
            for (int node = 1; node <= 11; node++)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            for (int node = 30; node <= 45; node++)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            for (int node = 190; node <= 191; node++)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            for (int node = 201; node <= 208; node++)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            for (int node = 216; node <= 223; node++)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            for (int node = 231; node <= 239; node++)
            {
                boundedDofs.Add(3 * node - 2);
                boundedDofs.Add(3 * node - 1);
                boundedDofs.Add(3 * node);
            }
            structuralBoundaryConditions = boundedDofs.ToArray<int>();
            //VectorOperations.PrintIntVectorToFile(structuralBoundaryConditions, @"C:\Users\Public\Documents\" + "BDoF.dat");
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
            string type = "ANSSolidShell8LEAS7";
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
            //VectorOperations.PrintIntVectorToFile(indices, @"C:\Users\Public\Documents\" + "indices.dat");
            //VectorOperations.PrintIntVectorToFile(indices12, @"C:\Users\Public\Documents\" + "indices2.dat");

            structuralSolution.LinearScheme = new LUFactorization();
            structuralSolution.ActivateNonLinearSolver = false;
            double[] externalForcesVector = externalForcesStructuralVector;
            externalForcesVector[166 * 3 - 3] = externalStructuralLoad;
            //externalForcesVector[265 * 3 - 1] = externalStructuralLoad;
            //externalForcesVector[249 * 3 - 3] = -externalStructuralLoad;
            //externalForcesVector[510 * 3 - 3] = externalStructuralLoad;
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

