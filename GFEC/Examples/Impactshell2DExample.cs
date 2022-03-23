using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GFEC
{
    class Impactshell2DExample
    {
        //1.0 X 1.0 /0.05-->21^2 elms -> 21^2 * 6 dof
        static int[] structuralBoundaryConditions;
        const double thickness = 0.006;
        const double uniformGap = 0.00001;
        const double offset = 0.40;
        const int nodesInX = 21;
        const int nodesInY = 21;
        const int nodesNumber = 441;
        const int shellElementsNumber = 400;
        //
        const int nodesInX2 = 3;
        const int nodesInY2 = 3;
        const int nodesInZ2 = 4;

        const int nodesNumber2 = 36;
        const int solidElementsNumber = 12;
        const int contactElements = 64;
        const double xInterv = 0.05;
        const double yInterv = 0.05;
        //
        const double xInterv2 = 0.10;
        const double yInterv2 = 0.10;
        const double zInterv = 0.125;

        //External loads
        //const double externalStructuralLoad = -10000.0;

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
                boundedDofs.Add(i * 6 + 3);
                //boundedDofs.Add(i * 6 + 6);
            }
            boundedDofs.Add(121);
            boundedDofs.Add(122);
            boundedDofs.Add(123);
            //boundedDofs.Add(126);
            for (int i = 1; i < nodesInX - 1; i++)
            {
                boundedDofs.Add(i * nodesInY * 6 + 3);
                //boundedDofs.Add(i * nodesInY * 6 + 6);
                boundedDofs.Add((i * nodesInY + nodesInY - 1) * 6 + 3);
                //boundedDofs.Add((i * nodesInY + nodesInY - 1) * 6 + 6);
            }
            boundedDofs.Add(2521);
            boundedDofs.Add(2522);
            for (int i = 0; i < nodesInY; i++)
            {
                boundedDofs.Add((nodesInX - 1) * nodesInY * 6 + i * 6 + 3);
                //boundedDofs.Add((nodesInX - 1) * nodesInY * 6 + i * 6 + 6);
            }
            boundedDofs.Add(2641);
            boundedDofs.Add(2642);
            for (int i = 0; i < nodesNumber; i++)
            {
                boundedDofs.Add(i * 6 + 6);
            }
            structuralBoundaryConditions = boundedDofs.ToArray<int>();
        }
        private static void CreateStructuralLoadVector()
        {
            loadedStructuralDOFs = new List<int>();
            //loadedStructuralDOFs.Add(4);
            //loadedStructuralDOFs.Add(5);
            //loadedStructuralDOFs.Add(124);
            //loadedStructuralDOFs.Add(125);
            //loadedStructuralDOFs.Add(2524);
            //loadedStructuralDOFs.Add(2525);
            //loadedStructuralDOFs.Add(2644);
            //loadedStructuralDOFs.Add(2645);
            //for (int i = 1; i < 20; i++)
            //{
            //    loadedStructuralDOFs.Add(i * 6 + 4);
            //    loadedStructuralDOFs.Add(i * 6 + 5);
            //}
            //for (int i = 1; i < 20; i++)
            //{
            //    loadedStructuralDOFs.Add(21 * i * 6 + 4);
            //    loadedStructuralDOFs.Add(21 * i * 6 + 5);
            //    loadedStructuralDOFs.Add(21 * i * 6 + 124);
            //    loadedStructuralDOFs.Add(21 * i * 6 + 125);
            //}
            //for (int i = 1; i < 20; i++)
            //{
            //    loadedStructuralDOFs.Add(2520 + i * 6 + 4);
            //    loadedStructuralDOFs.Add(2520 + i * 6 + 5);
            //}
            //for (int j = 1; j < 20; j++)
            //{
            //    for (int i = 1; i < 20; i++)
            //    {
            //        loadedStructuralDOFs.Add(j * 21 * 6 + i * 6 + 3);
            //        loadedStructuralDOFs.Add(j * 21 * 6 + i * 6 + 4);
            //        loadedStructuralDOFs.Add(j * 21 * 6 + i * 6 + 5);
            //    }
            //}
            externalForcesStructuralVector = new double[nodesNumber * 6 +
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
                    nodes[l] = new Node(i * xInterv, j * yInterv, thickness, 0.0, 0.0, 0.0);
                    l += 1;
                }
            }
            for (int i = 0; i < nodesInX2; i++)
            {
                for (int j = 0; j < nodesInY2; j++)
                {
                    for (int k = 0; k < nodesInZ2; k++)
                    {
                        nodes[l] = new Node(offset + i * xInterv2,
                            offset + j * yInterv2,
                            thickness + uniformGap + k * zInterv);
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
                    int first = (i - 1) * nodesInY + j;
                    connectivity[l] = new Dictionary<int, int>() { { 1, first }, { 2, first + nodesInY },
                            { 3,  first + nodesInY + 1}, { 4,  first + 1}};
                    l += 1;
                }
            }
            for (int i = 1; i <= nodesInX2 - 1; i++)
            {
                for (int j = 1; j <= nodesInY2 - 1; j++)
                {
                    for (int k = 1; k <= nodesInZ2 - 1; k++)
                    {
                        int firstDown = nodesNumber + 
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
            int slaveNode1 = nodesNumber + 1;
            for(int i = 8; i < 12; i++)
            {
                for (int j = 8; j < 12; j++)
                {
                    int firstMaster = (i - 1) * nodesInY + j;
                    connectivity[l] = new Dictionary<int, int>() { { 1, firstMaster },
                                                { 2, firstMaster + nodesInY },
                                                { 3, firstMaster + nodesInY + 1 },
                                                {4, firstMaster + 1 },
                                                { 5, slaveNode1 },
                                                { 6, slaveNode1 + nodesInY2 * nodesInZ2 },
                                                { 7, slaveNode1 + nodesInY2 * nodesInZ2 + nodesInZ2 },
                                                {8, slaveNode1 + nodesInZ2 } };
                    l += 1;
                }
            }
            int slaveNode2 = nodesNumber + nodesInZ2 + 1;
            for (int i = 8; i < 12; i++)
            {
                for (int j = 10; j < 14; j++)
                {
                    int firstMaster = (i - 1) * nodesInY + j;
                    connectivity[l] = new Dictionary<int, int>() { { 1, firstMaster },
                                                { 2, firstMaster + nodesInY },
                                                { 3, firstMaster + nodesInY + 1 },
                                                {4, firstMaster + 1 },
                                                { 5, slaveNode2 },
                                                { 6, slaveNode2 + nodesInY2 * nodesInZ2 },
                                                { 7, slaveNode2 + nodesInY2 * nodesInZ2 + nodesInZ2 },
                                                {8, slaveNode2 + nodesInZ2 } };
                    l += 1;
                }
            }
            int slaveNode3 = nodesNumber +nodesInY2 * nodesInZ2 + 1;
            for (int i = 10; i < 14; i++)
            {
                for (int j = 8; j < 12; j++)
                {
                    int firstMaster = (i - 1) * nodesInY + j;
                    connectivity[l] = new Dictionary<int, int>() { { 1, firstMaster },
                                                { 2, firstMaster + nodesInY },
                                                { 3, firstMaster + nodesInY + 1 },
                                                {4, firstMaster + 1 },
                                                { 5, slaveNode3 },
                                                { 6, slaveNode3 + nodesInY2 * nodesInZ2 },
                                                { 7, slaveNode3 + nodesInY2 * nodesInZ2 + nodesInZ2 },
                                                {8, slaveNode3 + nodesInZ2 } };
                    l += 1;
                }
            }
            int slaveNode4 = nodesNumber + nodesInY2 * nodesInZ2 + nodesInZ2 + 1;
            for (int i = 10; i < 14; i++)
            {
                for (int j = 10; j < 14; j++)
                {
                    int firstMaster = (i - 1) * nodesInY + j;
                    connectivity[l] = new Dictionary<int, int>() { { 1, firstMaster },
                                                { 2, firstMaster + nodesInY },
                                                { 3, firstMaster + nodesInY + 1 },
                                                {4, firstMaster + 1 },
                                                { 5, slaveNode4 },
                                                { 6, slaveNode4 + nodesInY2 * nodesInZ2 },
                                                { 7, slaveNode4 + nodesInY2 * nodesInZ2 + nodesInZ2 },
                                                {8, slaveNode4 + nodesInZ2 } };
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
                nodeFAT[i] = new bool[] { true, true, true, true, true, true };
            }
            for (int i = nodesNumber + 1; i <= nodesNumber + nodesNumber2; i++)
            {
                nodeFAT[i] = new bool[] { true, true, true, false, false, false };
            }
            return nodeFAT;
        }
        private static Dictionary<int, IElementProperties> CreateElementProperties()
        {
            double E = YoungMod;
            double A = area;
            string type = "Shell2DQuadratic4";
            string type2 = "Hex8";
            string type3 = "ContactStS3D";
            Dictionary<int, IElementProperties> elementProperties = new Dictionary<int, IElementProperties>();
            for (int i = 1; i <= shellElementsNumber; i++)
            {
                elementProperties[i] = new ElementProperties(E, poissonRatio, A, thickness, density, type);
            }
            for (int i = shellElementsNumber + 1; i <= shellElementsNumber + solidElementsNumber; i++)
            {
                elementProperties[i] = new ElementProperties(100.0 * E, poissonRatio,
                    A, thickness, density, type2);
            }
            for (int i = shellElementsNumber + solidElementsNumber + 1; i <= shellElementsNumber + solidElementsNumber + contactElements; i++)
            {
                elementProperties[i] = new ElementProperties(E, A, type3, 5.0, 10, 1, 1);
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
            elementsAssembly.ActivateBoundaryConditions = true;
            var AccelerationVector = new double[nodesNumber * 6 +
                nodesNumber2 * 3];
            var DisplacementVector = new double[nodesNumber * 6 +
                nodesNumber2 * 3];
            var VelocityVector = new double[nodesNumber * 6 +
                nodesNumber2 * 3];
            for (int i = nodesNumber * 6 + 2; i <= nodesNumber * 6 +
                nodesNumber2 * 3 - 1; i += 3)
            {
                VelocityVector[i] = -10.0;
            }
            InitialConditions initialValues = new InitialConditions();
            initialValues.InitialAccelerationVector = BoundaryConditionsImposition.ReducedVector(AccelerationVector, elementsAssembly.BoundedDOFsVector);
            initialValues.InitialDisplacementVector = BoundaryConditionsImposition.ReducedVector(DisplacementVector, elementsAssembly.BoundedDOFsVector);
            initialValues.InitialVelocityVector = BoundaryConditionsImposition.ReducedVector(VelocityVector, elementsAssembly.BoundedDOFsVector);
            initialValues.InitialTime = 0.0;
            ExplicitSolver newSolver = new ExplicitSolver(0.0001, 100);
            newSolver.Assembler = elementsAssembly;
            double[] externalForces = externalForcesStructuralVector;
            //foreach (var dof in loadedStructuralDOFs)
            //{
            //    externalForces[dof - 1] = externalForce;

            //}
            newSolver.InitialValues = initialValues;
            newSolver.ExternalForcesVector = BoundaryConditionsImposition.ReducedVector(externalForces, elementsAssembly.BoundedDOFsVector);
            newSolver.LinearSolver = new LUFactorization();
            newSolver.ActivateNonLinearSolution = true;
            //newSolver.SolveNewmark();
            newSolver.SolveExplicit();
            Tuple<Dictionary<int, double[]>, Dictionary<int, double>> solvectors = newSolver.GetResults();
            Dictionary<int, double[]> allStepsSolutions = solvectors.Item1;
            double[] fullDynamicSol10 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[9], elementsAssembly.BoundedDOFsVector);
            double[] fullDynamicSol20 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[19], elementsAssembly.BoundedDOFsVector);
            double[] fullDynamicSol30 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[29], elementsAssembly.BoundedDOFsVector);
            double[] fullDynamicSol40 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[39], elementsAssembly.BoundedDOFsVector);
            double[] fullDynamicSol50 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[49], elementsAssembly.BoundedDOFsVector);
            double[] fullDynamicSol75 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[74], elementsAssembly.BoundedDOFsVector);
            double[] fullDynamicSol100 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[99], elementsAssembly.BoundedDOFsVector);
            //double[] fullDynamicSol999 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[999], elementsAssembly.BoundedDOFsVector);
            VectorOperations.PrintVectorToFile(fullDynamicSol10, @"C:\Users\Public\Documents\Results10.dat");
            VectorOperations.PrintVectorToFile(fullDynamicSol20, @"C:\Users\Public\Documents\Results20.dat");
            VectorOperations.PrintVectorToFile(fullDynamicSol30, @"C:\Users\Public\Documents\Results30.dat");
            VectorOperations.PrintVectorToFile(fullDynamicSol40, @"C:\Users\Public\Documents\Results40.dat");
            VectorOperations.PrintVectorToFile(fullDynamicSol50, @"C:\Users\Public\Documents\Results50.dat");
            VectorOperations.PrintVectorToFile(fullDynamicSol75, @"C:\Users\Public\Documents\Results75.dat");
            VectorOperations.PrintVectorToFile(fullDynamicSol100, @"C:\Users\Public\Documents\Results100.dat");

            //VectorOperations.PrintVectorToFile(fullDynamicSol999, @"C:\Users\Public\Documents\Results999.dat");
            Results finalResults = new Results() { DynamicSolution = newSolver.explicitSolution, TimeSteps = newSolver.TimeAtEachStep, SelectedDOF = 1, SelectedInterval = 1, SolutionType = "Dynamic" };
            return finalResults;
        }
    }
}