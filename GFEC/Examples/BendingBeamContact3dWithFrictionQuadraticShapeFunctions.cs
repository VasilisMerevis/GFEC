using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GFEC
{
    public static class BendingBeamContact3dWithFrictionQuadraticShapeFunctions
    {

        public static ISolver structuralSolution;
        static int[] structuralBoundaryConditions;
        const double gap = 0.01;
        const double thickness = 0.1;
        const int nodesInX = 21;
        const int nodesInY = 3;
        const int nodesInZ = 3;
        const int nodesNumber = 378;
        //const int nodesNumber = 189;

        //const int elementsNumber = 10;

        const int elementsNumber = 20;
        const int contactElements = 28;

        //const double xInterv1 = 0.20;
        const double xInterv2 = 0.20;
        const double yInterv = 0.15;
        const double zInterv = 0.20;

        const double offset = 0.0;

        //External loads
        const double externalStructuralLoad = -30.0;

        static List<int> loadedStructuralDOFs;
        static double[] externalForcesStructuralVector;

        const double YoungMod = 1.0 * 1e5;

        const double poissonRatio = 0.25;
        const double density = 8000.0;
        const double area = 1.0;
        const double contactArea = thickness * xInterv2;

        //Friction coefficients
        //const double miS = 0.20;
        //const double miD = 0.20;


        private static void CreateStructuralBoundaryConditions()
        {
            List<int> boundedDofs = new List<int>();
            boundedDofs.Add(1);
            boundedDofs.Add(2);
            boundedDofs.Add(3);
            boundedDofs.Add(4);
            boundedDofs.Add(5);
            boundedDofs.Add(6);
            boundedDofs.Add(7);
            boundedDofs.Add(8);
            boundedDofs.Add(9);
            //boundedDofs.Add(10);
            //boundedDofs.Add(11);
            //boundedDofs.Add(12);
            //boundedDofs.Add(13);
            //boundedDofs.Add(14);
            //boundedDofs.Add(15);
            //boundedDofs.Add(16);
            //boundedDofs.Add(17);
            //boundedDofs.Add(18);
            boundedDofs.Add(542);
            boundedDofs.Add(543);
            boundedDofs.Add(545);
            boundedDofs.Add(546);
            boundedDofs.Add(548);
            boundedDofs.Add(549);

            boundedDofs.Add(568);
            boundedDofs.Add(569);
            boundedDofs.Add(570);
            boundedDofs.Add(571);
            boundedDofs.Add(572);
            boundedDofs.Add(573);
            boundedDofs.Add(574);
            boundedDofs.Add(575);
            boundedDofs.Add(576);
            boundedDofs.Add(577);
            boundedDofs.Add(578);
            boundedDofs.Add(579);
            boundedDofs.Add(580);
            boundedDofs.Add(581);
            boundedDofs.Add(582);
            boundedDofs.Add(583);
            boundedDofs.Add(584);
            boundedDofs.Add(585);
            boundedDofs.Add(586);
            boundedDofs.Add(587);
            boundedDofs.Add(588);
            boundedDofs.Add(589);
            boundedDofs.Add(590);
            boundedDofs.Add(591);
            boundedDofs.Add(592);
            boundedDofs.Add(593);
            boundedDofs.Add(594);
            boundedDofs.Add(1108);
            boundedDofs.Add(1109);
            boundedDofs.Add(1110);
            boundedDofs.Add(1111);
            boundedDofs.Add(1112);
            boundedDofs.Add(1113);
            boundedDofs.Add(1114);
            boundedDofs.Add(1115);
            boundedDofs.Add(1116);
            boundedDofs.Add(1117);
            boundedDofs.Add(1118);
            boundedDofs.Add(1119);
            boundedDofs.Add(1120);
            boundedDofs.Add(1121);
            boundedDofs.Add(1122);
            boundedDofs.Add(1123);
            boundedDofs.Add(1124);
            boundedDofs.Add(1125);
            boundedDofs.Add(1126);
            boundedDofs.Add(1127);
            boundedDofs.Add(1128);
            boundedDofs.Add(1129);
            boundedDofs.Add(1130);
            boundedDofs.Add(1131);
            boundedDofs.Add(1132);
            boundedDofs.Add(1133);
            boundedDofs.Add(1134);
            structuralBoundaryConditions = boundedDofs.ToArray<int>();
        }

        private static void CreateStructuralLoadVector()
        {
            loadedStructuralDOFs = new List<int>();
            loadedStructuralDOFs.Add(293);
            externalForcesStructuralVector = new double[nodesNumber * 3];
        }

        private static Dictionary<int, INode> CreateNodes()
        {

            Dictionary<int, INode> nodes = new Dictionary<int, INode>();
            int l;
            l = 1;
            //First beam
            for (int i = 0; i < nodesInX; i++)
            {
                for (int j = 0; j < nodesInY; j++)
                {
                    for (int k = 0; k < nodesInZ; k++)
                    {
                        nodes[l] = new Node(i * xInterv2, j * yInterv, k * zInterv);
                        l += 1;
                    }
                }
            }
            //Second beam
            for (int i = 0; i < nodesInX; i++)
            {
                for (int j = 0; j < nodesInY; j++)
                {
                    for (int k = 0; k < nodesInZ; k++)
                    {
                        nodes[l] = new Node(offset + i * xInterv2, j * yInterv - gap - (nodesInY - 1) * yInterv, k * zInterv);
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
            for (int i = 1; i <= 10; i++)
            {
                int first = (i - 1) * 2 * nodesInY * nodesInZ + 1;
                connectivity[l] = new Dictionary<int, int>() { { 1, first }, { 2, first + 2 * nodesInY * nodesInZ }, { 3, first + 2 * nodesInY * nodesInZ + 6 },
                    { 4, first + 6 },
                    { 5, first + 2 },{ 6, first + 2 * nodesInY * nodesInZ + 2 },
                    { 7, first + 2 * nodesInY * nodesInZ + 8 }, { 8, first + 8 },
                    { 9, first + nodesInY * nodesInZ }, { 10, first + 2 * nodesInY * nodesInZ + 3 },
                    { 11, first + 6 + nodesInY * nodesInZ }, { 12, first + 3 }, { 13, first + nodesInY * nodesInZ + 2 },
                    { 14, first + 2 * nodesInY * nodesInZ + 5 }, { 15, first + 8 + nodesInY * nodesInZ }, { 16, first + 5 },
                    { 17, first + 1 }, { 18, first + 1 + 2 * nodesInY * nodesInZ }, { 19, first + 2 * nodesInY * nodesInZ + 7 },
                    { 20, first + 7 }, { 21, first + 3 + nodesInY * nodesInZ }, { 22, first + nodesInY * nodesInZ + 5 },
                    { 23, first + nodesInY * nodesInZ + 1 }, { 24, first + nodesInY * nodesInZ + 7 }, { 25, first + 4 },
                    { 26, first + 2 * nodesInY * nodesInZ + 4 }, { 27, first + 4 + nodesInY * nodesInZ }};
                l += 1;
            }
            for (int i = 1; i <= 10; i++)
            {
                int first = nodesInX * nodesInY * nodesInZ + (i - 1) * 2 * nodesInY * nodesInZ + 1;
                connectivity[l] = new Dictionary<int, int>() { { 1, first }, { 2, first + 2 * nodesInY * nodesInZ }, { 3, first + 2 * nodesInY * nodesInZ + 6 },
                    { 4, first + 6 },
                    { 5, first + 2 },{ 6, first + 2 * nodesInY * nodesInZ + 2 },
                    { 7, first + 2 * nodesInY * nodesInZ + 8 }, { 8, first + 8 },
                    { 9, first + nodesInY * nodesInZ }, { 10, first + 2 * nodesInY * nodesInZ + 3 },
                    { 11, first + 6 + nodesInY * nodesInZ }, { 12, first + 3 }, { 13, first + nodesInY * nodesInZ + 2 },
                    { 14, first + 2 * nodesInY * nodesInZ + 5 }, { 15, first + 8 + nodesInY * nodesInZ }, { 16, first + 5 },
                    { 17, first + 1 }, { 18, first + 1 + 2 * nodesInY * nodesInZ }, { 19, first + 2 * nodesInY * nodesInZ + 7 },
                    { 20, first + 7 }, { 21, first + 3 + nodesInY * nodesInZ }, { 22, first + nodesInY * nodesInZ + 5 },
                    { 23, first + nodesInY * nodesInZ + 1 }, { 24, first + nodesInY * nodesInZ + 7 }, { 25, first + 4 },
                    { 26, first + 2 * nodesInY * nodesInZ + 4 }, { 27, first + 4 + nodesInY * nodesInZ }};
                l += 1;
            }
            //Contact elements
            for (int j = 0; j < 10; j++)
            {
                int firstSlaveNode = j * 2 * nodesInY * nodesInZ + 1;
                int firstMasterNode = firstSlaveNode + nodesInX * nodesInY * nodesInZ + 6;
                connectivity[l] = new Dictionary<int, int>()
                    {
                        { 1, firstMasterNode }, { 2, firstMasterNode + 1 }, { 3, firstMasterNode + 2 },
                        {4, firstMasterNode + 11 }, {5, firstMasterNode + 20 }, {6, firstMasterNode + 19 },
                        {7, firstMasterNode + 18 }, {8, firstMasterNode + 9 }, {9, firstMasterNode + 10 },
                        { 10, firstSlaveNode }, { 11, firstSlaveNode + 1 }, { 12, firstSlaveNode + 2 },
                        {13, firstSlaveNode + 11 }, {14, firstSlaveNode + 20 }, {15, firstSlaveNode + 19 },
                        {16, firstSlaveNode + 18 }, {17, firstSlaveNode + 9 }, {18, firstSlaveNode + 10 }
                    };
                l += 1;
            }
            for (int j = 0; j < 9; j++)
            {
                int firstSlaveNode = j * 2 * nodesInY * nodesInZ + 1;
                int firstMasterNode = firstSlaveNode + nodesInX * nodesInY * nodesInZ + nodesInY * nodesInZ + 6;
                connectivity[l] = new Dictionary<int, int>()
                    {
                        { 1, firstMasterNode }, { 2, firstMasterNode + 1 }, { 3, firstMasterNode + 2 },
                        {4, firstMasterNode + 11 }, {5, firstMasterNode + 20 }, {6, firstMasterNode + 19 },
                        {7, firstMasterNode + 18 }, {8, firstMasterNode + 9 }, {9, firstMasterNode + 10 },
                        { 10, firstSlaveNode }, { 11, firstSlaveNode + 1 }, { 12, firstSlaveNode + 2 },
                        {13, firstSlaveNode + 11 }, {14, firstSlaveNode + 20 }, {15, firstSlaveNode + 19 },
                        {16, firstSlaveNode + 18 }, {17, firstSlaveNode + 9 }, {18, firstSlaveNode + 10 }
                    };
                l += 1;
            }
            for (int j = 1; j < 10; j++)
            {
                int firstSlaveNode = j * 2 * nodesInY * nodesInZ + 1;
                int firstMasterNode = firstSlaveNode + nodesInX * nodesInY * nodesInZ - nodesInY * nodesInZ + 6;
                connectivity[l] = new Dictionary<int, int>()
                    {
                        { 1, firstMasterNode }, { 2, firstMasterNode + 1 }, { 3, firstMasterNode + 2 },
                        {4, firstMasterNode + 11 }, {5, firstMasterNode + 20 }, {6, firstMasterNode + 19 },
                        {7, firstMasterNode + 18 }, {8, firstMasterNode + 9 }, {9, firstMasterNode + 10 },
                        { 10, firstSlaveNode }, { 11, firstSlaveNode + 1 }, { 12, firstSlaveNode + 2 },
                        {13, firstSlaveNode + 11 }, {14, firstSlaveNode + 20 }, {15, firstSlaveNode + 19 },
                        {16, firstSlaveNode + 18 }, {17, firstSlaveNode + 9 }, {18, firstSlaveNode + 10 }
                    };
                l += 1;
            }
            return connectivity;
        }

        private static Dictionary<int, bool[]> CreateNodeFAT()
        {
            Dictionary<int, bool[]> nodeFAT = new Dictionary<int, bool[]>();
            for (int i = 1; i <= nodesNumber; i++)
            {
                nodeFAT[i] = new bool[] { true, true, true, false, false, false };
            }
            return nodeFAT;
        }
        private static Dictionary<int, IElementProperties> CreateElementProperties()
        {
            double E = YoungMod;

            double A = area;
            //double CA = contactArea;

            string type = "Hex27";
            //string type3 = "ContactStS2D";
            //string type3 = "ContactStS3Df";
            string type3 = "ContactStS3D";
            Dictionary<int, IElementProperties> elementProperties = new Dictionary<int, IElementProperties>();
            for (int i = 1; i <= elementsNumber; i++)
            {
                elementProperties[i] = new ElementProperties(E, poissonRatio, A, thickness, density, type);
            }
            for (int i = elementsNumber + 1; i <= elementsNumber + contactElements; i++)
            {
                //elementProperties[i] = new ElementProperties(E, type3, 5.0, 9, 1, 1, 5.0, miS, miD);
                elementProperties[i] = new ElementProperties(E, A, type3, 5.0, 15, 2, 2);
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
        public static Results RunStaticExample()
        {
            #region Structural
            IAssembly elementsAssembly = CreateAssembly();
            elementsAssembly.CreateElementsAssembly();
            elementsAssembly.ActivateBoundaryConditions = true;
            double[,] globalStiffnessMatrix = elementsAssembly.CreateTotalStiffnessMatrix();
            structuralSolution.LinearScheme = new LUFactorization();
            structuralSolution.NonLinearScheme.Tolerance = 1e-5;
            structuralSolution.ActivateNonLinearSolver = true;
            structuralSolution.NonLinearScheme.numberOfLoadSteps = 15;

            double[] externalForces3 = externalForcesStructuralVector;
            foreach (var dof in loadedStructuralDOFs)
            {
                externalForces3[dof - 1] = externalStructuralLoad;
            }
            double[] reducedExternalForces3 = BoundaryConditionsImposition.ReducedVector(externalForces3, elementsAssembly.BoundedDOFsVector);
            structuralSolution.AssemblyData = elementsAssembly;
            structuralSolution.Solve(reducedExternalForces3);
            double[] solvector3 = structuralSolution.GetSolution();
            Dictionary<int, double[]> allStepsSolutions = structuralSolution.GetAllStepsSolutions();
            //Dictionary<int, List<double[]>> gPointsStress = new Dictionary<int, List<double[]>>();
            //Dictionary<int, List<double[]>> gPointsStrain = new Dictionary<int, List<double[]>>();
            //Dictionary<int, List<double[]>> gPoints = new Dictionary<int, List<double[]>>();
            //Dictionary<int, List<double[]>> nodalStress = new Dictionary<int, List<double[]>>();
            //Dictionary<int, List<double[]>> nodalStrain = new Dictionary<int, List<double[]>>();
            //for (int i = 1; i <= allStepsSolutions.Count; i++)
            //{
            //    string name = "NodalCoordinates" + i.ToString() + ".dat";
            //    gPointsStress = elementsAssembly.GetElementsStresses(allStepsSolutions[i]);
            //    gPointsStrain = elementsAssembly.GetElementsStains(allStepsSolutions[i]);
            //    gPoints = elementsAssembly.GetElementsGaussPoints(allStepsSolutions[i]);
            //    nodalStress = elementsAssembly.GetElementsNodesStresses(allStepsSolutions[i]);
            //    nodalStrain = elementsAssembly.GetElementsNodesStains(allStepsSolutions[i]);
            //    string name1 = "GPointsStress" + i.ToString() + ".dat";
            //    string name2 = "GPointsStrain" + i.ToString() + ".dat";
            //    string name3 = "GPointsCoordinates" + i.ToString() + ".dat";
            //    string name4 = "NodalStress" + i.ToString() + ".dat";
            //    string name5 = "NodalStrain" + i.ToString() + ".dat";

            //    //VectorOperations.PrintDictionaryofListsofVectorsToFile(gPointsStress, @"C:\Users\Public\Documents\" + name1);
            //    //VectorOperations.PrintDictionaryofListsofVectorsToFile(gPointsStrain, @"C:\Users\Public\Documents\" + name2);
            //    //VectorOperations.PrintDictionaryofListsofVectorsToFile(gPoints, @"C:\Users\Public\Documents\" + name3);
            //    //VectorOperations.PrintDictionaryofListsofVectorsToFile(nodalStress, @"C:\Users\Public\Documents\" + name4);
            //    //VectorOperations.PrintDictionaryofListsofVectorsToFile(nodalStrain, @"C:\Users\Public\Documents\" + name5);
            //    if (i == allStepsSolutions.Count)
            //    {
            //        ExportToFile.ExportUpdatedNodalCoordinates(elementsAssembly, BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions.Single(m => m.Key == i).Value, elementsAssembly.BoundedDOFsVector), name);
            //        VectorOperations.PrintDictionaryofListsofVectorsToFile(gPointsStress, @"C:\Users\Public\Documents\" + name1);
            //        VectorOperations.PrintDictionaryofListsofVectorsToFile(gPointsStrain, @"C:\Users\Public\Documents\" + name2);
            //        VectorOperations.PrintDictionaryofListsofVectorsToFile(gPoints, @"C:\Users\Public\Documents\" + name3);
            //        VectorOperations.PrintDictionaryofListsofVectorsToFile(nodalStress, @"C:\Users\Public\Documents\" + name4);
            //        VectorOperations.PrintDictionaryofListsofVectorsToFile(nodalStrain, @"C:\Users\Public\Documents\" + name5);
            //    }
            //}
            elementsAssembly.UpdateDisplacements(solvector3);
            List<double[]> parametricCoordinatesVectors = new List<double[]>();
            for(int i = 0;  i < 11; i++)
            {
                for (int j = 0; j < 11; j++)
                {
                    double[] paramCoordinates = new double[]
                    {
                        -1.0 + 0.2 * i, 1.0, -1.0 + 0.2 * j
                    };
                    parametricCoordinatesVectors.Add(paramCoordinates);
                }
            }
            //ShowToGUI.PlotFinalGeometry(elementsAssembly);
            double[] fullSolVector3 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(solvector3, elementsAssembly.BoundedDOFsVector);
            //Dictionary<int, INode> finalNodes = Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullSolVector3);
            //double[] xFinalNodalCoor = Assembly.NodalCoordinatesToVectors(finalNodes).Item1;
            //double[] yFinalNodalCoor = Assembly.NodalCoordinatesToVectors(finalNodes).Item2;
            Dictionary<int, double[]> allStepsFullSolutions = new Dictionary<int, double[]>();
            //Dictionary<int, Dictionary<int, double[]>> allStepsContactForces = new Dictionary<int, Dictionary<int, double[]>>();
            //Dictionary<int, double[]> elementsContactForcesVector;
            for (int i = 1; i <= allStepsSolutions.Count; i++)
            {
                //elementsContactForcesVector = new Dictionary<int, double[]>();
                //elementsAssembly.UpdateDisplacements(allStepsSolutions[i]);
                //for (int j = 1; j <= contactElements; j++)
                //{
                //    elementsContactForcesVector[elementsNumber + j] = elementsAssembly.ElementsAssembly[elementsNumber + j].CreateInternalGlobalForcesVector();
                //}
                //allStepsContactForces[i] = elementsContactForcesVector;
                //string name = "ContactForces" + i.ToString() + ".dat";
                //double[] Vector = new double[contactElements * 27];
                //int count = 0;
                //for (int j = 1; j <= contactElements; j++)
                //{
                //    Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == elementsNumber + j).Value[0];
                //    count += 1;
                //    Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == elementsNumber + j).Value[1];
                //    count += 1;
                //    Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == elementsNumber + j).Value[2];
                //    count += 1;
                //    Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == elementsNumber + j).Value[3];
                //    count += 1;
                //    Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == elementsNumber + j).Value[4];
                //    count += 1;
                //    Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == elementsNumber + j).Value[5];
                //    count += 1;
                //    Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == elementsNumber + j).Value[6];
                //    count += 1;
                //    Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == elementsNumber + j).Value[7];
                //    count += 1;
                //    Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == elementsNumber + j).Value[8];
                //    count += 1;
                //    Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == elementsNumber + j).Value[9];
                //    count += 1;
                //    Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == elementsNumber + j).Value[10];
                //    count += 1;
                //    Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == elementsNumber + j).Value[11];
                //    count += 1;
                //    Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == elementsNumber + j).Value[12];
                //    count += 1;
                //    Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == elementsNumber + j).Value[13];
                //    count += 1;
                //    Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == elementsNumber + j).Value[14];
                //    count += 1;
                //    Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == elementsNumber + j).Value[15];
                //    count += 1;
                //    Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == elementsNumber + j).Value[16];
                //    count += 1;
                //    Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == elementsNumber + j).Value[17];
                //    count += 1;
                //    Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == elementsNumber + j).Value[18];
                //    count += 1;
                //    Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == elementsNumber + j).Value[19];
                //    count += 1;
                //    Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == elementsNumber + j).Value[20];
                //    count += 1;
                //    Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == elementsNumber + j).Value[21];
                //    count += 1;
                //    Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == elementsNumber + j).Value[22];
                //    count += 1;
                //    Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == elementsNumber + j).Value[23];
                //    count += 1;
                //    Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == elementsNumber + j).Value[24];
                //    count += 1;
                //    Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == elementsNumber + j).Value[25];
                //    count += 1;
                //    Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == elementsNumber + j).Value[26];
                //    count += 1;
                //}
                //VectorOperations.PrintVectorToFile(Vector, @"C:\Users\Public\Documents\" + name);
                if (i == allStepsSolutions.Count)
                {
                    for (int j = 11; j <= 20; j++)
                    {
                        //string name2 = "Element" + j.ToString() + "Stress.dat";
                        //List<double[]> stress = elementsAssembly.ElementsAssembly[j].GetStressFromElementsNodes();
                        //VectorOperations.PrintListofVectorsToFile(stress, @"C:\Users\Public\Documents\" + name2);
                        string name2 = "Element" + j.ToString() + "Stress.dat";
                        string name3 = "Element" + j.ToString() + "Coordinates.dat";

                        List<double[]> stress = elementsAssembly.ElementsAssembly[j].GetStressFromElements(parametricCoordinatesVectors);
                        VectorOperations.PrintListofVectorsToFile(stress, @"C:\Users\Public\Documents\" + name2);
                        List<double[]> physicalSpaceCoordinates = elementsAssembly.ElementsAssembly[j].GetphysicalCoordinatesFromElements(parametricCoordinatesVectors);
                        VectorOperations.PrintListofVectorsToFile(physicalSpaceCoordinates, @"C:\Users\Public\Documents\" + name3);
                    }
                }
            }

            for (int i = 0; i < allStepsSolutions.Count; i++)
            {
                allStepsFullSolutions.Add(i + 1, BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions.Single(m => m.Key == i + 1).Value, elementsAssembly.BoundedDOFsVector));
                int j = i + 1;
                string name = "solution" + j.ToString() + ".dat";
                //VectorOperations.PrintVectorToFile(allStepsFullSolutions.Single(m => m.Key == i + 1).Value, @"C:\Users\Public\Documents\" + name);
                VectorOperations.PrintVectorToFile(allStepsFullSolutions.Single(m => m.Key == i + 1).Value, @"C:\Users\Public\Documents\" + name);
                //if (i == allStepsSolutions.Count - 1)
                //{
                //    VectorOperations.PrintVectorToFile(allStepsFullSolutions.Single(m => m.Key == i + 1).Value, @"C:\Users\Public\Documents\" + name);

                //}
            }
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
