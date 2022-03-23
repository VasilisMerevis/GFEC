using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GFEC
{
    class DegenerateShellElementsImpactExample
    {
        static int[] structuralBoundaryConditions;
        const double thickness = 0.006;
        const double uniformGap = 0.0001;
        const double offsetX = 0.25;
        const double offsetY = 0.125;
        const double offsetZ = 0.5779;

        const int nodesInX = 3;
        const int nodesInY = 3;
        const int nodesInZ = 11;

        const int SolidElementsNodesNumber = 99;

        const int solidElementsNumber1 = 40;
        //
        const int nodesInXY = 21;

        //const int nodesInY2 = 2;
        const int nodesInZ2 = 7;

        const int shellElementsNodesNumber = 147;
        const int shellElementsNumber = 30;
        const int contactElements = 120;
        const double xInterv = 0.125;
        const double yInterv = 0.125;
        const double zInterv = 0.125;
        //
        const double radius = 0.20;
        const double initialTheta = 225.0 * Math.PI / 180.0;
        const double deltaTheta = -4.5 * Math.PI / 180.0;
        const double zInterv2 = 0.0157;

        //External loads
        //const double externalStructuralLoad = -10000.0;

        static List<int> loadedStructuralDOFs;
        static double[] externalForcesStructuralVector;

        const double YoungMod = 30.0 * 1e6;
        const double YoungMod2 = 200.0 * 1e6;

        const double poissonRatio = 0.20;
        const double poissonRatio2 = 0.25;
        //const double density = 2400.0;
        //const double density2 = 8000.0;
        const double density = 2.40;
        const double density2 = 8.00;
        const double area = 1.0;
        private static void CreateStructuralBoundaryConditions()
        {
            List<int> boundedDofs = new List<int>();
            //boundedDofs.Add(shellElementsNodesNumber * 5 + 1);
            //boundedDofs.Add(shellElementsNodesNumber * 5 + 2);
            //boundedDofs.Add(shellElementsNodesNumber * 5 + 3);
            //boundedDofs.Add(shellElementsNodesNumber * 5 + 34);
            //boundedDofs.Add(shellElementsNodesNumber * 5 + 35);
            //boundedDofs.Add(shellElementsNodesNumber * 5 + 36);
            //boundedDofs.Add(shellElementsNodesNumber * 5 + 67);
            //boundedDofs.Add(shellElementsNodesNumber * 5 + 68);
            //boundedDofs.Add(shellElementsNodesNumber * 5 + 69);
            //boundedDofs.Add(shellElementsNodesNumber * 5 + 100);
            //boundedDofs.Add(shellElementsNodesNumber * 5 + 101);
            //boundedDofs.Add(shellElementsNodesNumber * 5 + 102);
            //boundedDofs.Add(shellElementsNodesNumber * 5 + 133);
            //boundedDofs.Add(shellElementsNodesNumber * 5 + 134);
            //boundedDofs.Add(shellElementsNodesNumber * 5 + 135);
            //boundedDofs.Add(shellElementsNodesNumber * 5 + 166);
            //boundedDofs.Add(shellElementsNodesNumber * 5 + 167);
            //boundedDofs.Add(shellElementsNodesNumber * 5 + 168);
            //boundedDofs.Add(shellElementsNodesNumber * 5 + 199);
            //boundedDofs.Add(shellElementsNodesNumber * 5 + 200);
            //boundedDofs.Add(shellElementsNodesNumber * 5 + 201);
            //boundedDofs.Add(shellElementsNodesNumber * 5 + 232);
            //boundedDofs.Add(shellElementsNodesNumber * 5 + 233);
            //boundedDofs.Add(shellElementsNodesNumber * 5 + 234);
            //boundedDofs.Add(shellElementsNodesNumber * 5 + 265);
            //boundedDofs.Add(shellElementsNodesNumber * 5 + 266);
            //boundedDofs.Add(shellElementsNodesNumber * 5 + 267);
            for(int i = 0; i < nodesInXY; i++)
            {
                if (i == 0 || i == nodesInXY - 1)
                {
                    for (int j = 1; j <= nodesInZ2; j++)
                    {
                        boundedDofs.Add(i * nodesInZ2 * 5 + (j - 1) * 5 + 1);
                        boundedDofs.Add(i * nodesInZ2 * 5 + (j - 1) * 5 + 2);
                        boundedDofs.Add(i * nodesInZ2 * 5 + (j - 1) * 5 + 3);
                        boundedDofs.Add(i * nodesInZ2 * 5 + (j - 1) * 5 + 4);
                        boundedDofs.Add(i * nodesInZ2 * 5 + (j - 1) * 5 + 5);
                    }
                }
                //else
                //{
                //    boundedDofs.Add(i * nodesInZ2 * 5 + 1);
                //    boundedDofs.Add(i * nodesInZ2 * 5 + 2);
                //    boundedDofs.Add(i * nodesInZ2 * 5 + 3);
                //    boundedDofs.Add(i * nodesInZ2 * 5 + 4);
                //    boundedDofs.Add(i * nodesInZ2 * 5 + 5);
                //    boundedDofs.Add(i * nodesInZ2 * 5 + (nodesInZ2 - 1) * 5 + 1);
                //    boundedDofs.Add(i * nodesInZ2 * 5 + (nodesInZ2 - 1) * 5 + 2);
                //    boundedDofs.Add(i * nodesInZ2 * 5 + (nodesInZ2 - 1) * 5 + 3);
                //    boundedDofs.Add(i * nodesInZ2 * 5 + (nodesInZ2 - 1) * 5 + 4);
                //    boundedDofs.Add(i * nodesInZ2 * 5 + (nodesInZ2 - 1) * 5 + 5);
                //}
            }
            structuralBoundaryConditions = boundedDofs.ToArray<int>();
        }
        private static void CreateStructuralLoadVector()
        {
            loadedStructuralDOFs = new List<int>();
            externalForcesStructuralVector = new double[SolidElementsNodesNumber * 3 +
                shellElementsNodesNumber * 5];
        }

        private static Dictionary<int, INode> CreateNodes()
        {

            Dictionary<int, INode> nodes = new Dictionary<int, INode>();
            int l;
            l = 1;
            for (int i = 0; i < nodesInXY; i++)
            {
                for (int k = 0; k < nodesInZ2; k++)
                {
                    nodes[l] = new Node(
                        offsetX + uniformGap + thickness / 2.0 + radius
                        + (radius + thickness / 2.0) * Math.Cos(initialTheta + deltaTheta * i),
                        offsetY + (radius + thickness / 2.0) *
                        Math.Sin(initialTheta + deltaTheta * i),
                        offsetZ + k * zInterv2);
                    l += 1;
                }
            }
            for (int i = 0; i < nodesInXY; i++)
            {
                for (int k = 0; k < nodesInZ2; k++)
                {
                    nodes[l] = new Node(
                        offsetX + uniformGap + thickness / 2.0 + radius
                        + (radius - thickness / 2.0) * Math.Cos(initialTheta + deltaTheta * i),
                        offsetY + (radius - thickness / 2.0) *
                        Math.Sin(initialTheta + deltaTheta * i),
                        offsetZ + k * zInterv2);
                    l += 1;
                }
            }
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
                        int firstDown = 2 * shellElementsNodesNumber + (i - 1) * nodesInY * nodesInZ +
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
            for (int i = 1; i <= nodesInXY - 2; i+=2)
            {
                for (int k = 1; k <= nodesInZ2 - 2; k+=2)
                {
                    int firstNode = (i - 1) * nodesInZ2 + k;
                    int firstNode2 = firstNode + shellElementsNodesNumber;
                    connectivity[l] = new Dictionary<int, int>() 
                    {
                            { 1, firstNode },
                            { 2, firstNode + nodesInZ2 },
                            { 3,  firstNode + nodesInZ2 + nodesInZ2},
                            { 4,  firstNode + nodesInZ2 + nodesInZ2 + 1},
                            { 5, firstNode + nodesInZ2 + nodesInZ2 + 2 },
                            { 6, firstNode + nodesInZ2 + 2 },
                            { 7,  firstNode + 2},
                            { 8,  firstNode + 1},
                            { 9,  firstNode + nodesInZ2 + 1},
                            { 10, firstNode2 },
                            { 11, firstNode2 + nodesInZ2 },
                            { 12,  firstNode2 + nodesInZ2 + nodesInZ2},
                            { 13,  firstNode2 + nodesInZ2 + nodesInZ2 + 1},
                            { 14, firstNode2 + nodesInZ2 + nodesInZ2 + 2 },
                            { 15, firstNode2 + nodesInZ2 + 2 },
                            { 16,  firstNode2 + 2},
                            { 17,  firstNode2 + 1},
                            { 18,  firstNode2 + nodesInZ2 + 1}
                    };
                    l += 1;
                }
            }
            //contact elements
            for (int i = 1; i <= nodesInXY - 2; i += 2)
            {
                for (int k = 1; k <= nodesInZ2 - 2; k += 2)
                {
                    int masterNode1 = 2 * shellElementsNodesNumber + 71;
                    int masterNode2 = 2 * shellElementsNodesNumber + 82;
                    int masterNode3 = 2 * shellElementsNodesNumber + 83;
                    int masterNode4 = 2 * shellElementsNodesNumber + 72;
                    int firstNode = (i - 1) * nodesInZ2 + k;
                    int firstNode2 = firstNode + shellElementsNodesNumber;
                    connectivity[l] = new Dictionary<int, int>()
                    {
                            { 1, masterNode1 },
                            { 2, masterNode2 },
                            { 3, masterNode3 },
                            { 4, masterNode4 },
                            { 5, firstNode },
                            { 6, firstNode + nodesInZ2 },
                            { 7,  firstNode + nodesInZ2 + nodesInZ2},
                            { 8,  firstNode + nodesInZ2 + nodesInZ2 + 1},
                            { 9, firstNode + nodesInZ2 + nodesInZ2 + 2 },
                            { 10, firstNode + nodesInZ2 + 2 },
                            { 11,  firstNode + 2},
                            { 12,  firstNode + 1},
                            { 13,  firstNode + nodesInZ2 + 1}
                    };
                    l += 1;
                    masterNode1 = 2 * shellElementsNodesNumber + 72;
                    masterNode2 = 2 * shellElementsNodesNumber + 83;
                    masterNode3 = 2 * shellElementsNodesNumber + 84;
                    masterNode4 = 2 * shellElementsNodesNumber + 73;
                    connectivity[l] = new Dictionary<int, int>()
                    {
                            { 1, masterNode1 },
                            { 2, masterNode2 },
                            { 3, masterNode3 },
                            { 4, masterNode4 },
                            { 5, firstNode },
                            { 6, firstNode + nodesInZ2 },
                            { 7,  firstNode + nodesInZ2 + nodesInZ2},
                            { 8,  firstNode + nodesInZ2 + nodesInZ2 + 1},
                            { 9, firstNode + nodesInZ2 + nodesInZ2 + 2 },
                            { 10, firstNode + nodesInZ2 + 2 },
                            { 11,  firstNode + 2},
                            { 12,  firstNode + 1},
                            { 13,  firstNode + nodesInZ2 + 1}
                    };
                    l += 1;
                    masterNode1 = 2 * shellElementsNodesNumber + 82;
                    masterNode2 = 2 * shellElementsNodesNumber + 93;
                    masterNode3 = 2 * shellElementsNodesNumber + 94;
                    masterNode4 = 2 * shellElementsNodesNumber + 83;
                    connectivity[l] = new Dictionary<int, int>()
                    {
                            { 1, masterNode1 },
                            { 2, masterNode2 },
                            { 3, masterNode3 },
                            { 4, masterNode4 },
                            { 5, firstNode },
                            { 6, firstNode + nodesInZ2 },
                            { 7,  firstNode + nodesInZ2 + nodesInZ2},
                            { 8,  firstNode + nodesInZ2 + nodesInZ2 + 1},
                            { 9, firstNode + nodesInZ2 + nodesInZ2 + 2 },
                            { 10, firstNode + nodesInZ2 + 2 },
                            { 11,  firstNode + 2},
                            { 12,  firstNode + 1},
                            { 13,  firstNode + nodesInZ2 + 1}
                    };
                    l += 1;
                    masterNode1 = 2 * shellElementsNodesNumber + 83;
                    masterNode2 = 2 * shellElementsNodesNumber + 94;
                    masterNode3 = 2 * shellElementsNodesNumber + 95;
                    masterNode4 = 2 * shellElementsNodesNumber + 84;
                    connectivity[l] = new Dictionary<int, int>()
                    {
                            { 1, masterNode1 },
                            { 2, masterNode2 },
                            { 3, masterNode3 },
                            { 4, masterNode4 },
                            { 5, firstNode },
                            { 6, firstNode + nodesInZ2 },
                            { 7,  firstNode + nodesInZ2 + nodesInZ2},
                            { 8,  firstNode + nodesInZ2 + nodesInZ2 + 1},
                            { 9, firstNode + nodesInZ2 + nodesInZ2 + 2 },
                            { 10, firstNode + nodesInZ2 + 2 },
                            { 11,  firstNode + 2},
                            { 12,  firstNode + 1},
                            { 13,  firstNode + nodesInZ2 + 1}
                    };
                    l += 1;
                }
            }
            return connectivity;
        }

        private static Dictionary<int, bool[]> CreateNodeFAT()
        {
            Dictionary<int, bool[]> nodeFAT = new Dictionary<int, bool[]>();
            for (int i = 1; i <= shellElementsNodesNumber; i++)
            {
                nodeFAT[i] = new bool[] { true, true, true, true, true, false };
            }
            for (int i = shellElementsNodesNumber + 1; i <= SolidElementsNodesNumber + shellElementsNodesNumber; i++)
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
            string type2 = "IsoparamShell18";
            string type3 = "ContactStS3D";
            Dictionary<int, IElementProperties> elementProperties = new Dictionary<int, IElementProperties>();
            for (int i = 1; i <= solidElementsNumber1; i++)
            {
                elementProperties[i] = new ElementProperties(E1, poissonRatio, A, thickness, density, type);
            }
            for (int i = solidElementsNumber1 + 1; i <= solidElementsNumber1 + shellElementsNumber; i++)
            {
                elementProperties[i] = new ElementProperties(E2, poissonRatio2, A, thickness, density2, type2);

            }
            for (int i = solidElementsNumber1 + shellElementsNumber + 1; i <= solidElementsNumber1 + shellElementsNumber + contactElements; i++)
            {
                elementProperties[i] = new ElementProperties(E2, A, type3, 2.0, 5, 2, 1);
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
            var AccelerationVector = new double[SolidElementsNodesNumber * 3 +
                shellElementsNodesNumber * 5];
            var DisplacementVector = new double[SolidElementsNodesNumber * 3 +
                shellElementsNodesNumber * 5];
            var VelocityVector = new double[SolidElementsNodesNumber * 3 +
                shellElementsNodesNumber * 5];
            //for (int i = 0; i <= shellElementsNodesNumber * 5 - 5; i += 5)
            //{
            //    VelocityVector[i] = -50.0;
            //}
            for (int i = shellElementsNodesNumber * 5;
                i <= shellElementsNodesNumber * 5 + SolidElementsNodesNumber * 3 - 3; i += 3)
            {
                VelocityVector[i] = 5.0;
            }
            InitialConditions initialValues = new InitialConditions();
            initialValues.InitialAccelerationVector = BoundaryConditionsImposition.ReducedVector(AccelerationVector, elementsAssembly.BoundedDOFsVector);
            initialValues.InitialDisplacementVector = BoundaryConditionsImposition.ReducedVector(DisplacementVector, elementsAssembly.BoundedDOFsVector);
            initialValues.InitialVelocityVector = BoundaryConditionsImposition.ReducedVector(VelocityVector, elementsAssembly.BoundedDOFsVector);
            initialValues.InitialTime = 0.0;
            ExplicitSolver newSolver = new ExplicitSolver(0.0015, 500);
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

            for (int i = 0; i <= allStepsSolutions.Keys.Max(); i++)
            {
                double[] fullDynamicSol = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[i], elementsAssembly.BoundedDOFsVector);
                var k = i + 1;
                VectorOperations.PrintVectorToFile(fullDynamicSol, @"C:\Users\Public\Documents\DegenerateShellElementsImpactExampleResults" + k.ToString() + ".dat");
            }
            //for (int i = 9000; i <= 10000; i++)
            //{
            //    double[] fullDynamicSol = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[i], elementsAssembly.BoundedDOFsVector);
            //    var k = i + 1;
            //    VectorOperations.PrintVectorToFile(fullDynamicSol, @"C:\Users\Public\Documents\DegenerateShellElementsImpactExampleResults" + k.ToString() + ".dat");
            //}
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
            Results finalResults = new Results() { DynamicSolution = newSolver.explicitSolution, TimeSteps = newSolver.TimeAtEachStep, SelectedDOF = 1, SelectedInterval = 1, SolutionType = "Dynamic" };
            return finalResults;
        }
    }
}