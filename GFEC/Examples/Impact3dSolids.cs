using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GFEC
{
    class Impact3dSolids
    {
        //1.0 X 1.0 /0.05-->21^2 elms -> 21^2 * 6 dof
        static int[] structuralBoundaryConditions;
        const double thickness = 0.006;
        const double uniformGap = 0.0005;
        const double offsetX = 0.25;
        const double offsetY = 0.075;
        const double offsetZ = 0.575;

        const int nodesInX = 3;
        const int nodesInY = 3;
        const int nodesInZ = 11;

        const int nodesNumber1 = 99;

        const int solidElementsNumber1 = 40;
        //
        const int nodesInX2 = 2;
        const int nodesInY2 = 2;
        const int nodesInZ2 = 2;

        const int nodesNumber2 = 8;
        const int solidElementsNumber2 = 1;
        const int contactElements = 4;
        const double xInterv = 0.125;
        const double yInterv = 0.125;
        const double zInterv = 0.125;
        //
        const double xInterv2 = 0.10;
        const double yInterv2 = 0.10;
        const double zInterv2 = 0.10;

        //External loads
        //const double externalStructuralLoad = -10000.0;

        static List<int> loadedStructuralDOFs;
        static double[] externalForcesStructuralVector;

        const double YoungMod = 4.0 * 1e9;
        const double YoungMod2 = 200.0 * 1e9;

        const double poissonRatio = 0.25;
        const double density = 1000.0;
        const double density2 = 8000.0;
        const double area = 1.0;
        private static void CreateStructuralBoundaryConditions()
        {
            List<int> boundedDofs = new List<int>();
            boundedDofs.Add(1);
            boundedDofs.Add(2);
            boundedDofs.Add(3);
            boundedDofs.Add(34);
            boundedDofs.Add(35);
            boundedDofs.Add(36);
            boundedDofs.Add(67);
            boundedDofs.Add(68);
            boundedDofs.Add(69);
            boundedDofs.Add(100);
            boundedDofs.Add(101);
            boundedDofs.Add(102);
            boundedDofs.Add(133);
            boundedDofs.Add(134);
            boundedDofs.Add(135);
            boundedDofs.Add(166);
            boundedDofs.Add(167);
            boundedDofs.Add(168);
            boundedDofs.Add(199);
            boundedDofs.Add(200);
            boundedDofs.Add(201);
            boundedDofs.Add(232);
            boundedDofs.Add(233);
            boundedDofs.Add(234);
            boundedDofs.Add(265);
            boundedDofs.Add(266);
            boundedDofs.Add(267);
            structuralBoundaryConditions = boundedDofs.ToArray<int>();
        }
        private static void CreateStructuralLoadVector()
        {
            loadedStructuralDOFs = new List<int>();
            externalForcesStructuralVector = new double[nodesNumber1 * 3 +
                nodesNumber2 * 3];
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
                    for (int k = 0; k < nodesInZ; k++)
                    {
                        nodes[l] = new Node(i * xInterv, j * yInterv, k * zInterv);
                        l += 1;
                    }
                }
            }
            for (int i = 0; i < nodesInX2; i++)
            {
                for (int j = 0; j < nodesInY2; j++)
                {
                    for (int k = 0; k < nodesInZ2; k++)
                    {
                        nodes[l] = new Node(offsetX + uniformGap + i * xInterv2,
                            offsetY + j * yInterv2,
                            offsetZ + k * zInterv2);
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
                        int firstDown = (i - 1) * nodesInY * nodesInZ +
                            (j - 1) * nodesInZ +
                            k;
                        int firstUpp = firstDown + 1;
                        connectivity[l] = new Dictionary<int, int>() {
                            { 1, firstDown },
                            { 2, firstDown + nodesInY * nodesInZ },
                            { 3,  firstDown + nodesInY * nodesInZ + nodesInZ},
                            { 4,  firstDown + nodesInZ},
                            { 5, firstUpp },
                            { 6, firstUpp + nodesInY * nodesInZ },
                            { 7,  firstUpp + + nodesInY * nodesInZ + nodesInZ},
                            { 8,  firstUpp + nodesInZ}
                        };
                        l += 1;
                    }
                }
            }
            for (int i = 1; i <= nodesInX2 - 1; i++)
            {
                for (int j = 1; j <= nodesInY2 - 1; j++)
                {
                    for (int k = 1; k <= nodesInZ2 - 1; k++)
                    {
                        int firstDown = nodesNumber1 +
                            (i - 1) * nodesInY2 * nodesInZ2 + (j - 1) * nodesInZ2 +
                            k;
                        int firstUpp = firstDown + 1;
                        connectivity[l] = new Dictionary<int, int>() {
                            { 1, firstDown },
                            { 2, firstDown + nodesInY2 * nodesInZ2 },
                            { 3,  firstDown + nodesInY2 * nodesInZ2 + nodesInZ2},
                            { 4,  firstDown + nodesInZ2},
                            { 5, firstUpp },
                            { 6, firstUpp + nodesInY2 * nodesInZ2 },
                            { 7,  firstUpp + + nodesInY2 * nodesInZ2 + nodesInZ2},
                            { 8,  firstUpp + nodesInZ2}
                        };
                        l += 1;
                    }
                }
            }
            //contact elements
            int slaveNode1 = nodesNumber1 + 1;
            int slaveNode2 = nodesNumber1 + 3;
            int slaveNode3 = nodesNumber1 + 4;
            int slaveNode4 = nodesNumber1 + 2;
            int masterNode1 = 71;
            int masterNode2 = 82;
            int masterNode3 = 83;
            int masterNode4 = 72;
            connectivity[l] = new Dictionary<int, int>()
            {
                { 1, masterNode1 },
                { 2, masterNode2 },
                { 3, masterNode3 },
                { 4, masterNode4 },
                { 5, slaveNode1 },
                { 6, slaveNode2 },
                { 7, slaveNode3 },
                { 8, slaveNode4 }
            };
            l += 1;
            masterNode1 = 72;
            masterNode2 = 83;
            masterNode3 = 84;
            masterNode4 = 73;
            connectivity[l] = new Dictionary<int, int>()
            {
                { 1, masterNode1 },
                { 2, masterNode2 },
                { 3, masterNode3 },
                { 4, masterNode4 },
                { 5, slaveNode1 },
                { 6, slaveNode2 },
                { 7, slaveNode3 },
                { 8, slaveNode4 }
            };
            l += 1;
            masterNode1 = 82;
            masterNode2 = 93;
            masterNode3 = 94;
            masterNode4 = 83;
            connectivity[l] = new Dictionary<int, int>()
            {
                { 1, masterNode1 },
                { 2, masterNode2 },
                { 3, masterNode3 },
                { 4, masterNode4 },
                { 5, slaveNode1 },
                { 6, slaveNode2 },
                { 7, slaveNode3 },
                { 8, slaveNode4 }
            };
            l += 1;
            masterNode1 = 83;
            masterNode2 = 94;
            masterNode3 = 95;
            masterNode4 = 84;
            connectivity[l] = new Dictionary<int, int>()
            {
                { 1, masterNode1 },
                { 2, masterNode2 },
                { 3, masterNode3 },
                { 4, masterNode4 },
                { 5, slaveNode1 },
                { 6, slaveNode2 },
                { 7, slaveNode3 },
                { 8, slaveNode4 }
            };
            return connectivity;
        }

        private static Dictionary<int, bool[]> CreateNodeFAT()
        {
            Dictionary<int, bool[]> nodeFAT = new Dictionary<int, bool[]>();
            for (int i = 1; i <= nodesNumber1 + nodesNumber2; i++)
            {
                nodeFAT[i] = new bool[] { true, true, true, false, false, false };
            }
            return nodeFAT;
        }
        private static Dictionary<int, IElementProperties> CreateElementProperties()
        {
            double E1 = YoungMod;
            double E2 = YoungMod2;
            double A = area;
            string type = "Hex8";
            string type2 = "ContactStS3D";
            Dictionary<int, IElementProperties> elementProperties = new Dictionary<int, IElementProperties>();
            for (int i = 1; i <= solidElementsNumber1; i++)
            {
                elementProperties[i] = new ElementProperties(E1, poissonRatio, A, thickness, density, type);
            }
            for (int i = solidElementsNumber1 + 1; i <= solidElementsNumber1 + solidElementsNumber2; i++)
            {
                elementProperties[i] = new ElementProperties(E2, poissonRatio,
                    A, thickness, density2, type);
            }
            for (int i = solidElementsNumber1 + solidElementsNumber2 + 1; i <= solidElementsNumber1 + solidElementsNumber2 + contactElements; i++)
            {
                elementProperties[i] = new ElementProperties(E1, A, type2, 10.0, 5, 1, 1);
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
        public static Results RunDynamicExample()
        {
            IAssembly elementsAssembly = CreateAssembly();
            elementsAssembly.CreateElementsAssembly();
            ExportToFile.ExportMatlabInitialGeometry(elementsAssembly);
            elementsAssembly.ActivateBoundaryConditions = true;
            var AccelerationVector = new double[nodesNumber1 * 3 +
                nodesNumber2 * 3];
            var DisplacementVector = new double[nodesNumber1 * 3 +
                nodesNumber2 * 3];
            var VelocityVector = new double[nodesNumber1 * 3 +
                nodesNumber2 * 3];
            for (int i = nodesNumber1 * 3; i <= nodesNumber1 * 3 +
                nodesNumber2 * 3 - 3; i += 3)
            {
                VelocityVector[i] = -50.0;
            }
            InitialConditions initialValues = new InitialConditions();
            initialValues.InitialAccelerationVector = BoundaryConditionsImposition.ReducedVector(AccelerationVector, elementsAssembly.BoundedDOFsVector);
            initialValues.InitialDisplacementVector = BoundaryConditionsImposition.ReducedVector(DisplacementVector, elementsAssembly.BoundedDOFsVector);
            initialValues.InitialVelocityVector = BoundaryConditionsImposition.ReducedVector(VelocityVector, elementsAssembly.BoundedDOFsVector);
            initialValues.InitialTime = 0.0;
            ExplicitSolver newSolver = new ExplicitSolver(0.001, 100);
            newSolver.Assembler = elementsAssembly;
            double[] externalForces = externalForcesStructuralVector;
            newSolver.InitialValues = initialValues;
            newSolver.ExternalForcesVector = BoundaryConditionsImposition.ReducedVector(externalForces, elementsAssembly.BoundedDOFsVector);
            newSolver.LinearSolver = new LUFactorization();
            newSolver.ActivateNonLinearSolution = true;
            newSolver.SolveNewmark();
            //newSolver.SolveExplicit();
            Tuple<Dictionary<int, double[]>, Dictionary<int, double>> solvectors = newSolver.GetResults();
            Dictionary<int, double[]> allStepsSolutions = solvectors.Item1;
            //List<int> solutionsIndices = new List<int>();
            //solutionsIndices.Add(0);
            //solutionsIndices.Add(99);
            //solutionsIndices.Add(199);
            //solutionsIndices.Add(299);
            //solutionsIndices.Add(399);
            //solutionsIndices.Add(499);
            //solutionsIndices.Add(599);
            //solutionsIndices.Add(699);
            //solutionsIndices.Add(799);
            //solutionsIndices.Add(899);
            //solutionsIndices.Add(999);

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
            //get stress
            //List<double[]> parametricCoordinatesVectors = new List<double[]>();
            //for (int i = 0; i < 19; i++)
            //{
            //    for (int j = 0; j < 19; j++)
            //    {
            //        double[] paramCoordinates = new double[]
            //        {
            //            1.0, -0.90 + 0.1 * j, -0.90 + 0.1 * i
            //        };
            //        parametricCoordinatesVectors.Add(paramCoordinates);
            //    }
            //}
            //for (int j = 21; j <= 40; j++)
            //{
            //    string name2 = "Element" + j.ToString() + "Stress.dat";
            //    string name3 = "Element" + j.ToString() + "Coordinates.dat";

            //    List<double[]> stress = elementsAssembly.ElementsAssembly[j].GetStressFromElements(parametricCoordinatesVectors);
            //    VectorOperations.PrintListofVectorsToFile(stress, @"C:\Users\Public\Documents\Results\StressResults\timeStep" + 61.ToString() + name2);
            //    List<double[]> physicalSpaceCoordinates = elementsAssembly.ElementsAssembly[j].GetphysicalCoordinatesFromElements(parametricCoordinatesVectors);
            //    VectorOperations.PrintListofVectorsToFile(physicalSpaceCoordinates, @"C:\Users\Public\Documents\Results\StressResults\timeStep" + 61.ToString() + name3);
            //}
            //get stress
            //double[] fullDynamicSol1 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[1], elementsAssembly.BoundedDOFsVector);
            //double[] fullDynamicSol2 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[2], elementsAssembly.BoundedDOFsVector);
            //double[] fullDynamicSol3 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[3], elementsAssembly.BoundedDOFsVector);
            //double[] fullDynamicSol40 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[4], elementsAssembly.BoundedDOFsVector);
            //double[] fullDynamicSol50 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[49], elementsAssembly.BoundedDOFsVector);
            //double[] fullDynamicSol75 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[74], elementsAssembly.BoundedDOFsVector);
            //double[] fullDynamicSol100 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[99], elementsAssembly.BoundedDOFsVector);
            //double[] fullDynamicSol999 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[999], elementsAssembly.BoundedDOFsVector);
            //VectorOperations.PrintVectorToFile(fullDynamicSol1, @"C:\Users\Public\Documents\Results1.dat");
            //VectorOperations.PrintVectorToFile(fullDynamicSol2, @"C:\Users\Public\Documents\Results2.dat");
            //VectorOperations.PrintVectorToFile(fullDynamicSol3, @"C:\Users\Public\Documents\Results3.dat");
            //VectorOperations.PrintVectorToFile(fullDynamicSol40, @"C:\Users\Public\Documents\Results40.dat");
            //VectorOperations.PrintVectorToFile(fullDynamicSol50, @"C:\Users\Public\Documents\Results50.dat");
            //VectorOperations.PrintVectorToFile(fullDynamicSol75, @"C:\Users\Public\Documents\Results75.dat");
            //VectorOperations.PrintVectorToFile(fullDynamicSol100, @"C:\Users\Public\Documents\Results100.dat");

            //VectorOperations.PrintVectorToFile(fullDynamicSol999, @"C:\Users\Public\Documents\Results999.dat");
            Results finalResults = new Results() { DynamicSolution = newSolver.explicitSolution, TimeSteps = newSolver.TimeAtEachStep, SelectedDOF = 1, SelectedInterval = 1, SolutionType = "Dynamic" };
            return finalResults;
        }
    }
}