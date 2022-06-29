using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading;

namespace GFEC
{
    public static class TruckBumperImpact
    {
        public static ISolver structuralSolution;
        static int[] structuralBoundaryConditions;
        //public static  Dictionary<int, INode> nodes;
        //public static Dictionary<int, Dictionary<int, int>> elementsConnectivity;

        static List<int> loadedStructuralDOFs;
        static double[] externalForcesStructuralVector;
        const double externalStructuralLoad = 100.0;
        const int nodesNumber = 5134;
        const int nodesNumberBumper = 3726;
        const int nodesNumberRigidObject = 1408;

        const int elmntsNumber = 2820;
        const int elmntsNumberSolidShell = 1760;
        const int elmntsNumberSolid = 1060;

        private static void CreateStructuralBoundaryConditions(List<int> fixedNodes)
        {
            List<int> boundedDofs = new List<int>();
            foreach (var node in fixedNodes)
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
            double E1 = 200.0 * 1e9;
            double E2 = 1.0 * 1e12;
            double density1 = 7850.0;
            double density2 = 7850.0;
            double poissonRatio1 = 0.28;
            double poissonRatio2 = 0.20;
            string type = "ANSSolidShell8EAS";
            string type2 = "Hex8";
            string type3 = "ContactStS3D";

            Dictionary<int, IElementProperties> elementProperties = new Dictionary<int, IElementProperties>();
            int totalElements = elementsConnectivity.Count;
            for (int i = 1; i <= elmntsNumberSolidShell; i++)
            {
                elementProperties[i] = new ElementProperties(E1, poissonRatio1, 1.0, 0.005, density1, type);
            }
            for (int i = elmntsNumberSolidShell + 1; i <= elmntsNumberSolidShell + elmntsNumberSolid; i++)
            {
                elementProperties[i] = new ElementProperties(E2, poissonRatio2, 1.0, 0.005, density2, type2);
            }
            for (int i = elmntsNumberSolidShell + elmntsNumberSolid + 1; i <= totalElements; i++)
            {
                elementProperties[i] = new ElementProperties(E1, 1.0, type3, 4d, 3, 1, 1);
            }
            return elementProperties;
        }
        private static Dictionary<int, Dictionary<int, int>> CreateConnectivity(Dictionary<int, Dictionary<int, int>> elementsConnectivity,
           Dictionary<int, Dictionary<int, int>> masterElementsConnectivity,
            Dictionary<int, Dictionary<int, int>> slaveElementsConnectivity)
        {

            Dictionary<int, Dictionary<int, int>> connectivity = elementsConnectivity;
            //Contact elements
            int connectCount = connectivity.Count + 1;
            foreach (var nodeList in masterElementsConnectivity)
            {
                if (nodeList.Key <= 8)
                {
                    for (int i = 65; i <= 80; i++)
                    {
                        connectivity[connectCount] = new Dictionary<int, int>()
                        {
                        { 1, slaveElementsConnectivity[i][1] },
                        { 2, slaveElementsConnectivity[i][2] },
                        { 3, slaveElementsConnectivity[i][3]},
                        { 4, slaveElementsConnectivity[i][4]},
                        { 5, nodeList.Value[1] },
                        { 6, nodeList.Value[2] },
                        { 7, nodeList.Value[3] },
                        { 8, nodeList.Value[4] }

                        };
                        connectCount += 1;
                    }
                }
                if ((nodeList.Key >= 9 && nodeList.Key <= 32) ||
                    (nodeList.Key >= 65 && nodeList.Key <= 72) ||
                   (nodeList.Key >= 97 && nodeList.Key <= 104))
                {
                    for (int i = 57; i <= 80; i++)
                    {
                        connectivity[connectCount] = new Dictionary<int, int>()
                        {
                        { 1, slaveElementsConnectivity[i][1] },
                        { 2, slaveElementsConnectivity[i][2] },
                        { 3, slaveElementsConnectivity[i][3]},
                        { 4, slaveElementsConnectivity[i][4]},
                        { 5, nodeList.Value[1] },
                        { 6, nodeList.Value[2] },
                        { 7, nodeList.Value[3] },
                        { 8, nodeList.Value[4] }
                        };
                        connectCount += 1;
                    }
                }
                if (nodeList.Key >= 65 && nodeList.Key <= 112)
                {
                    for (int i = 25; i <= 56; i++)
                    {
                        connectivity[connectCount] = new Dictionary<int, int>()
                        {
                        { 1, slaveElementsConnectivity[i][1] },
                        { 2, slaveElementsConnectivity[i][2] },
                        { 3, slaveElementsConnectivity[i][3]},
                        { 4, slaveElementsConnectivity[i][4]},
                        { 5, nodeList.Value[1] },
                        { 6, nodeList.Value[2] },
                        { 7, nodeList.Value[3] },
                        { 8, nodeList.Value[4] }
                        };
                        connectCount += 1;
                    }
                }
                if ((nodeList.Key >= 33 && nodeList.Key <= 40))
                {
                    for (int i = 1; i <= 16; i++)
                    {
                        connectivity[connectCount] = new Dictionary<int, int>()
                        {
                        { 1, slaveElementsConnectivity[i][1] },
                        { 2, slaveElementsConnectivity[i][2] },
                        { 3, slaveElementsConnectivity[i][3]},
                        { 4, slaveElementsConnectivity[i][4]},
                        { 5, nodeList.Value[1] },
                        { 6, nodeList.Value[2] },
                        { 7, nodeList.Value[3] },
                        { 8, nodeList.Value[4] }                
                        };
                        connectCount += 1;
                    }
                }
                if ((nodeList.Key >= 41 && nodeList.Key <= 64) ||
                   (nodeList.Key >= 73 && nodeList.Key <= 80) ||
                   (nodeList.Key >= 105 && nodeList.Key <= 112))
                {
                    for (int i = 1; i <= 24; i++)
                    {
                        connectivity[connectCount] = new Dictionary<int, int>()
                        {
                        { 1, slaveElementsConnectivity[i][1] },
                        { 2, slaveElementsConnectivity[i][2] },
                        { 3, slaveElementsConnectivity[i][3]},
                        { 4, slaveElementsConnectivity[i][4]},
                        { 5, nodeList.Value[1] },
                        { 6, nodeList.Value[2] },
                        { 7, nodeList.Value[3] },
                        { 8, nodeList.Value[4] }
                        };
                        connectCount += 1;
                    }
                }
            }
            return connectivity;
        }

        private static IAssembly CreateAssembly(Dictionary<int, INode> nodes, Dictionary<int, Dictionary<int, int>> elementsConnectivity,
            List<int> fixedNodes,
            Dictionary<int, Dictionary<int, int>> masterElementsConnectivity,
            Dictionary<int, Dictionary<int, int>> slaveElementsConnectivity)
        {
            IAssembly assembly = new Assembly();
            assembly.Nodes = nodes;
            assembly.ElementsConnectivity = CreateConnectivity(elementsConnectivity, masterElementsConnectivity, slaveElementsConnectivity);
            assembly.ElementsProperties = CreateElementProperties(assembly.ElementsConnectivity);
            assembly.NodeFreedomAllocationList = CreateNodeFAT(nodes);
            CreateStructuralBoundaryConditions(fixedNodes);
            CreateStructuralLoadVector();
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

        public static Results RunStaticExample(Dictionary<int, INode> nodes, Dictionary<int, Dictionary<int, int>> elementsConnectivity, List<int> fixedNodes,
            Dictionary<int, Dictionary<int, int>> masterElementsConnectivity,
            Dictionary<int, Dictionary<int, int>> slaveElementsConnectivity)
        {
            #region Structural
            IAssembly elementsAssembly = CreateAssembly(nodes, elementsConnectivity, fixedNodes, masterElementsConnectivity, slaveElementsConnectivity);
            elementsAssembly.CreateElementsAssembly();
            ExportToFile.ExportMatlabInitialGeometry(elementsAssembly);
            elementsAssembly.ActivateBoundaryConditions = true;

           
            double[,] globalStiffnessMatrix = elementsAssembly.CreateTotalStiffnessMatrix();
           

           

            double[] stiffnessVector = new double[globalStiffnessMatrix.GetLength(1)];
            for (int i = 0; i < globalStiffnessMatrix.GetLength(1); i++)
            {
                stiffnessVector[i] = globalStiffnessMatrix[i, i];
            }
            VectorOperations.PrintVectorToFile(stiffnessVector, @"C:\Users\Public\Documents\StiffDiagonal.dat");
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
            VectorOperations.PrintIntVectorToFile(indices, @"C:\Users\Public\Documents\" + "indices.dat");
            VectorOperations.PrintIntVectorToFile(indices12, @"C:\Users\Public\Documents\" + "indices2.dat");

            structuralSolution.LinearScheme = new Skyline();
            //structuralSolution.NonLinearScheme.Tolerance = 1e-4;
            structuralSolution.ActivateNonLinearSolver = true;
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

        public static Results RunDynamicExample(Dictionary<int, INode> nodes, Dictionary<int, Dictionary<int, int>> elementsConnectivity, List<int> fixedNodes,
            Dictionary<int, Dictionary<int, int>> masterElementsConnectivity,
            Dictionary<int, Dictionary<int, int>> slaveElementsConnectivity)
        {
            IAssembly elementsAssembly = CreateAssembly(nodes, elementsConnectivity, fixedNodes, masterElementsConnectivity, slaveElementsConnectivity);
            elementsAssembly.CreateElementsAssembly();
            ExportToFile.ExportMatlabInitialGeometry(elementsAssembly);

            
            Stopwatch watch1 = new Stopwatch();
            watch1.Start();
            double[,] globalStiffnessMatrix = elementsAssembly.CreateTotalStiffnessMatrix();
            watch1.Stop();
            long watch1Time = watch1.ElapsedMilliseconds;

            elementsAssembly.ActivateParallelCalculations = true;
            Stopwatch watch2 = new Stopwatch();
            watch2.Start();
            double[,] globalStiffnessMatrixParallel = elementsAssembly.CreateTotalStiffnessMatrix();
            watch2.Stop();
            long watch2Time = watch2.ElapsedMilliseconds;
            //double[,] massMatrix = elementsAssembly.CreateTotalMassMatrix();
            //double[] internalForces = elementsAssembly.CreateTotalInternalForcesVector();
            //VectorOperations.PrintVectorToFile(internalForces, @"C:\Users\Public\Documents\" + "internalForces.dat");

            //double[] massVector = new double[massMatrix.GetLength(1)];
            //for (int i = 0; i < massMatrix.GetLength(1); i++)
            //{
            //    massVector[i] = massMatrix[i, i];
            //}
            //VectorOperations.PrintVectorToFile(massVector, @"C:\Users\Public\Documents\Mass.dat");
            //if (!massMatrix.Cast<double>().Any(d => double.IsNaN(d) || double.IsInfinity(d)))
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

            //for (int i = 0; i < massMatrix.GetLength(0); i++)
            //{
            //    for (int j = 0; j < massMatrix.GetLength(1); j++)
            //    {
            //        if (double.IsNaN(massMatrix[i, j]) || double.IsInfinity(massMatrix[i, j]))
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
            elementsAssembly.ActivateBoundaryConditions = true;
            var AccelerationVector = new double[nodesNumber * 3];
            var DisplacementVector = new double[nodesNumber * 3];
            var VelocityVector = new double[nodesNumber * 3];
            for (int i = nodesNumberBumper * 3 + 2; i <= nodesNumber * 3 - 1; i += 3)
            {
                VelocityVector[i] =  25.0;
            }
            InitialConditions initialValues = new InitialConditions();
            initialValues.InitialAccelerationVector = BoundaryConditionsImposition.ReducedVector(AccelerationVector, elementsAssembly.BoundedDOFsVector);
            initialValues.InitialDisplacementVector = BoundaryConditionsImposition.ReducedVector(DisplacementVector, elementsAssembly.BoundedDOFsVector);
            initialValues.InitialVelocityVector = BoundaryConditionsImposition.ReducedVector(VelocityVector, elementsAssembly.BoundedDOFsVector);
            initialValues.InitialTime = 0.0;
            ExplicitSolver newSolver = new ExplicitSolver(0.005, 12);
            newSolver.Assembler = elementsAssembly;
            double[] externalForces = externalForcesStructuralVector;
            newSolver.InitialValues = initialValues;
            newSolver.ExternalForcesVector = BoundaryConditionsImposition.ReducedVector(externalForces, elementsAssembly.BoundedDOFsVector);
            //newSolver.LinearSolver = new LUFactorization();
            newSolver.LinearSolver = new Skyline();
            newSolver.ActivateNonLinearSolution = true;
            newSolver.SolveNewmark();
            //newSolver.SolveExplicit();
            Tuple<Dictionary<int, double[]>, Dictionary<int, double>> solvectors = newSolver.GetResults();
            Dictionary<int, double[]> allStepsSolutions = solvectors.Item1;

            //for (int i = 0; i< solutionsIndices.Count; i++)
            //{
            //    int index = solutionsIndices[i];
            //    double[] fullDynamicSol = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[index], elementsAssembly.BoundedDOFsVector);
            //    var k = index + 1;
            //    VectorOperations.PrintVectorToFile(fullDynamicSol, @"C:\Users\Public\Documents\Results" + k.ToString() + ".dat");

            //}

            for (int i = 0; i <= allStepsSolutions.Keys.Max(); i++)
            {
                double[] fullDynamicSol = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[i], elementsAssembly.BoundedDOFsVector);
                var k = i + 1;
                VectorOperations.PrintVectorToFile(fullDynamicSol, @"C:\Users\Public\Documents\Results" + k.ToString() + ".dat");

            }
            Results finalResults = new Results() { DynamicSolution = newSolver.explicitSolution, TimeSteps = newSolver.TimeAtEachStep, SelectedDOF = 1, SelectedInterval = 1, SolutionType = "Dynamic" };
            return finalResults;
        }

    }
}

