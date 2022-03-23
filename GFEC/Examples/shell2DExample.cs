using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GFEC
{
    class shell2DExample
    {
        public static ISolver structuralSolution;
        static int[] structuralBoundaryConditions;
        const double thickness = 0.006;
        const int nodesInX = 21;
        const int nodesInY = 21;
        const int nodesNumber = 441;
        const int elementsNumber = 400;

        //const double xInterv1 = 0.20;
        const double xInterv = 0.05;
        const double yInterv = 0.05;

        //External loads
        const double externalStructuralLoad = -1000.0;

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
            loadedStructuralDOFs.Add(4);
            loadedStructuralDOFs.Add(5);
            loadedStructuralDOFs.Add(124);
            loadedStructuralDOFs.Add(125);
            loadedStructuralDOFs.Add(2524);
            loadedStructuralDOFs.Add(2525);
            loadedStructuralDOFs.Add(2644);
            loadedStructuralDOFs.Add(2645);
            for (int i = 1; i < 20; i++)
            {
                loadedStructuralDOFs.Add(i * 6 + 4);
                loadedStructuralDOFs.Add(i * 6 + 5);
            }
            for (int i = 1; i < 20; i++)
            {
                loadedStructuralDOFs.Add(21 * i * 6 + 4);
                loadedStructuralDOFs.Add(21 * i * 6 + 5);
                loadedStructuralDOFs.Add(21 * i * 6 + 124);
                loadedStructuralDOFs.Add(21 * i * 6 + 125);
            }
            for (int i = 1; i < 20; i++)
            {
                loadedStructuralDOFs.Add(2520 + i * 6 + 4);
                loadedStructuralDOFs.Add(2520 + i * 6 + 5);
            }
            for (int j = 1; j < 20; j++)
            {
                for (int i = 1; i < 20; i++)
                {
                    loadedStructuralDOFs.Add(j * 21 * 6 + i * 6 + 3);
                    loadedStructuralDOFs.Add(j * 21 * 6 + i * 6 + 4);
                    loadedStructuralDOFs.Add(j * 21 * 6 + i * 6 + 5);
                }
            }
            externalForcesStructuralVector = new double[nodesNumber * 6];
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
                    nodes[l] = new Node(i * xInterv, j * yInterv, thickness, 0.0, 0.0,0.0);
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
                    connectivity[l] = new Dictionary<int, int>() { { 1, first }, { 2, first + nodesInY },
                            { 3,  first + nodesInY + 1}, { 4,  first + 1}};
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
            return nodeFAT;
        }
        private static Dictionary<int, IElementProperties> CreateElementProperties()
        {
            double E = YoungMod;
            double A = area;
            string type = "Shell2DQuadratic4";

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
            structuralSolution.LinearScheme = new PCGSolver();
            //structuralSolution.NonLinearScheme.Tolerance = 1e-5;//
            //structuralSolution.ActivateNonLinearSolver = true;
            //structuralSolution.NonLinearScheme.numberOfLoadSteps = 10;//

            double[] externalForces3 = externalForcesStructuralVector;
            for (int i = 0; i < loadedStructuralDOFs.Count; i++)
            {
                var dof = loadedStructuralDOFs[i];
                //if (i == 0)
                //{
                //    externalForces3[dof - 1] = -externalStructuralLoad * xInterv * yInterv *
                //        yInterv / 3.0;
                //}
                //if (i == 1)
                //{
                //    externalForces3[dof - 1] = externalStructuralLoad * xInterv * yInterv *
                //        xInterv / 3.0;
                //}
                //if (i == 2)
                //{
                //    externalForces3[dof - 1] = externalStructuralLoad * xInterv * yInterv *
                //        yInterv / 3.0;
                //}
                //if (i == 3)
                //{
                //    externalForces3[dof - 1] = externalStructuralLoad * xInterv * yInterv *
                //        xInterv / 3.0;
                //}
                //if (i == 4)
                //{
                //    externalForces3[dof - 1] = -externalStructuralLoad * xInterv * yInterv *
                //        yInterv / 3.0;
                //}
                //if (i == 5)
                //{
                //    externalForces3[dof - 1] = -externalStructuralLoad * xInterv * yInterv *
                //        xInterv / 3.0;
                //}
                //if (i == 6)
                //{
                //    externalForces3[dof - 1] = externalStructuralLoad * xInterv * yInterv *
                //        yInterv / 3.0;
                //}
                //if (i == 7)
                //{
                //    externalForces3[dof - 1] = -externalStructuralLoad * xInterv * yInterv *
                //        xInterv / 3.0;
                //}
                //if (i == 8 || i == 10 || i == 12 || i == 14 || i == 16 || i == 18 ||
                //    i == 20 || i == 22 || i == 24 || i == 26 || i == 28 || i == 30 ||
                //    i == 32 || i == 34 || i == 36 || i == 38 || i == 40 || i == 42 ||
                //    i == 44)
                //{
                //    externalForces3[dof - 1] = 0.0;
                //}
                //if (i == 9 || i == 11 || i == 13 || i == 15 || i == 17 || i == 19 ||
                //    i == 21 || i == 23 || i == 25 || i == 27 || i == 29 || i == 31 ||
                //    i == 33 || i == 35 || i == 37 || i == 39 || i == 41 || i == 43 ||
                //    i == 45)
                //{
                //    externalForces3[dof - 1] = 0.0;
                //}
                //if (i == 46 || i == 50 || i == 54 || i == 58 || i == 62 || i == 66 ||
                //    i == 70 || i == 74 || i == 78 || i == 82 || i == 86 || i == 90 ||
                //    i == 94 || i == 98 || i == 102 || i == 106 || i == 110 || i == 114 ||
                //    i == 118)
                //{
                //    externalForces3[dof - 1] = 0.0;
                //}
                //if (i == 47 || i == 51 || i == 55 || i == 59 || i == 63 || i == 67 ||
                //    i == 71 || i == 75 || i == 79 || i == 83 || i == 87 || i == 91 ||
                //    i == 95 || i == 99 || i == 103 || i == 107 || i == 111 || i == 115 ||
                //    i == 119)
                //{
                //    externalForces3[dof - 1] = externalForces3[dof - 1] = externalStructuralLoad * xInterv * yInterv *
                //        2.0 * xInterv / 3.0;
                //}
                //if (i == 48 || i == 52 || i == 56 || i == 60 || i == 64 || i == 68 ||
                //    i == 72 || i == 76 || i == 80 || i == 84 || i == 88 || i == 92 ||
                //    i == 96 || i == 100 || i == 104 || i == 108 || i == 112 || i == 116 ||
                //    i == 120)
                //{
                //    externalForces3[dof - 1] = 0.0;
                //}
                //if (i == 49 || i == 53 || i == 57 || i == 61 || i == 65 || i == 69 ||
                //    i == 73 || i == 77 || i == 81 || i == 85 || i == 89 || i == 93 ||
                //    i == 97 || i == 101 || i == 105 || i == 109 || i == 113 || i == 117 ||
                //    i == 121)
                //{
                //    externalForces3[dof - 1] = -externalStructuralLoad * xInterv * yInterv *
                //        2.0 * xInterv / 3.0;
                //}
                //if (i == 122 || i == 124 || i == 126 || i == 128 || i == 130 || i == 132 ||
                //    i == 134 || i == 136 || i == 138 || i == 140 || i == 142 || i == 144 ||
                //    i == 146 || i == 148 || i == 150 || i == 152 || i == 154 || i == 156 ||
                //    i == 158)
                //{
                //    externalForces3[dof - 1] = 0.0;
                //}
                //if (i == 123 || i == 125 || i == 127 || i == 129 || i == 131 || i == 133 ||
                //    i == 135 || i == 137 || i == 139 || i == 141 || i == 143 || i == 145 ||
                //    i == 147 || i == 149 || i == 151 || i == 153 || i == 155 || i == 157 ||
                //    i == 159)
                //{
                //    externalForces3[dof - 1] = 0.0;
                //}
                if (i >= 160)
                {
                    if(i == 160 ||
                        (i - 160) % 3  == 0)
                    {
                        externalForces3[dof - 1] = externalStructuralLoad * xInterv * yInterv;
                    }
                }
            }
            double[] reducedExternalForces3 = BoundaryConditionsImposition.ReducedVector(externalForces3, elementsAssembly.BoundedDOFsVector);
            structuralSolution.AssemblyData = elementsAssembly;
            structuralSolution.Solve(reducedExternalForces3);
            double[] solvector = structuralSolution.GetSolution();
            double[] fullSolVector = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(solvector, 
                elementsAssembly.BoundedDOFsVector);

            string name = "solution" + ".dat";
            VectorOperations.PrintVectorToFile(fullSolVector, @"C:\Users\Public\Documents\" + name);

            //Dictionary<int, double[]> allStepsSolutions = structuralSolution.GetAllStepsSolutions();

            //elementsAssembly.UpdateDisplacements(solvector3);
            //double[] fullSolVector3 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(solvector3, elementsAssembly.BoundedDOFsVector);
            //Dictionary<int, double[]> allStepsFullSolutions = new Dictionary<int, double[]>();
            //for (int i = 0; i < allStepsSolutions.Count; i++)
            //{
            //    allStepsFullSolutions.Add(i + 1, BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions.Single(m => m.Key == i + 1).Value, elementsAssembly.BoundedDOFsVector));
            //    int j = i + 1;
            //    string name = "solution" + j.ToString() + ".dat";
            //    VectorOperations.PrintVectorToFile(allStepsFullSolutions.Single(m => m.Key == i + 1).Value, @"C:\Users\Public\Documents\" + name);
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