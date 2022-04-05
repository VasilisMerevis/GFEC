using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GFEC
{
    public static class Blocks3dContactSlidingQuadratic
    {

        public static ISolver structuralSolution;
        static int[] structuralBoundaryConditions;
        const double gap = 0.00001;
        const double thickness = 0.1;
        const int nodesInX = 13;
        const int nodesInY = 3;
        const int nodesInZ = 3;
        const int nodesNumber = 234;
        const int elementsNumber = 12;
        const int contactElements = 16;

        //const double xInterv1 = 0.20;
        const double xInterv2 = 0.25;
        const double yInterv = 0.25;
        const double zInterv = 0.25;

        const double offset = 0.0;

        //External loads
        const double externalStructuralLoad = 3000.0;

        static List<int> loadedStructuralDOFs;
        static double[] externalForcesStructuralVector;

        const double YoungMod = 1.0 * 1e6;

        const double poissonRatio = 0.25;
        const double density = 8000.0;
        const double area = 1.0;
        const double contactArea = thickness * xInterv2;

        //Friction coefficients
        const double miS = 1.00;
        const double miD = 1.00;


        private static void CreateStructuralBoundaryConditions()
        {
            List<int> boundedDofs = new List<int>();
            //boundedDofs.Add(1);
            boundedDofs.Add(2);
            boundedDofs.Add(3);
            //boundedDofs.Add(4);
            boundedDofs.Add(5);
            boundedDofs.Add(6);
            //boundedDofs.Add(7);
            boundedDofs.Add(8);
            boundedDofs.Add(9);
            //boundedDofs.Add(10);
            //boundedDofs.Add(11);
            //boundedDofs.Add(12);
            ////boundedDofs.Add(13);
            ////boundedDofs.Add(14);
            ////boundedDofs.Add(15);
            //boundedDofs.Add(16);
            //boundedDofs.Add(17);
            //boundedDofs.Add(18);
            //boundedDofs.Add(19);
            //boundedDofs.Add(20);
            //boundedDofs.Add(21);
            //boundedDofs.Add(22);
            //boundedDofs.Add(23);
            //boundedDofs.Add(24);
            //boundedDofs.Add(25);
            //boundedDofs.Add(26);
            //boundedDofs.Add(27);
            //boundedDofs.Add(325);
            //

            //boundedDofs.Add(326);
            //boundedDofs.Add(327);
            ////boundedDofs.Add(328);
            //boundedDofs.Add(329);
            //boundedDofs.Add(330);
            ////boundedDofs.Add(331);
            //boundedDofs.Add(332);
            //boundedDofs.Add(333);
            for (int i = 1; i <= 27; i++)
            {
                boundedDofs.Add(324 + i);
            }
            for (int i = 1; i <= 27; i++)
            {
                boundedDofs.Add(351 + i);
            }
            for (int i = 1; i <= 9; i++)
            {
                boundedDofs.Add(162 + (i - 1) * 3 + 3);
            }
            for (int i = 1; i <= 27; i++)
            {
                boundedDofs.Add(675 + i);
            }
            structuralBoundaryConditions = boundedDofs.ToArray<int>();
        }

        private static void CreateStructuralLoadVector()
        {
            loadedStructuralDOFs = new List<int>();
            for (int i = 0; i < 12; i++)
            {
                loadedStructuralDOFs.Add(27 * i + 20);
                loadedStructuralDOFs.Add(27 * i + 23);
                loadedStructuralDOFs.Add(27 * i + 26);
            }
            //for (int i = 0; i < 9; i++)
            //{
            //    loadedStructuralDOFs.Add(27 * 12 + i * 3 + 1);
            //}
            for (int i = 0; i < 9; i++)
            {
                loadedStructuralDOFs.Add(i * 3 + 1);
            }
            externalForcesStructuralVector = new double[nodesNumber * 3];
        }

        private static Dictionary<int, INode> CreateNodes()
        {

            Dictionary<int, INode> nodes = new Dictionary<int, INode>();
            int l;
            l = 1;
            //First block
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
            //Second block
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
            for (int i = 1; i <= 6; i++)
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
            for (int i = 1; i <= 6; i++)
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
            for (int j = 0; j < 6; j++)
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
            for (int j = 0; j < 5; j++)
            {
                int firstSlaveNode = j * 2 * nodesInY * nodesInZ + 1;
                int firstMasterNode = firstSlaveNode + nodesInX * nodesInY * nodesInZ + 2 * nodesInY * nodesInZ + 6;
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
            for (int j = 1; j < 6; j++)
            {
                int firstSlaveNode = j * 2 * nodesInY * nodesInZ + 1;
                int firstMasterNode = firstSlaveNode + nodesInX * nodesInY * nodesInZ -2 * nodesInY * nodesInZ + 6;
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
            string type3 = "ContactStS3Df";
            //string type3 = "ContactStS3D";
            Dictionary<int, IElementProperties> elementProperties = new Dictionary<int, IElementProperties>();
            for (int i = 1; i <= elementsNumber; i++)
            {
                elementProperties[i] = new ElementProperties(E, poissonRatio, A, thickness, density, type);
            }
            for (int i = elementsNumber + 1; i <= elementsNumber + contactElements; i++)
            {
                elementProperties[i] = new ElementProperties(E, type3, 5.0, 15, 2, 2, 1.20, miS, miD);
                //elementProperties[i] = new ElementProperties(E, A, type3, 5.0, 15, 2, 2);
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
            int countContactElements = elementsAssembly.CountElementsOfSameType(typeof(ContactStS2D));
            //ShowToGUI.PlotInitialGeometry(elementsAssembly);
            structuralSolution.LinearScheme = new LUFactorization();
            structuralSolution.NonLinearScheme.Tolerance = 1e-5;
            structuralSolution.ActivateNonLinearSolver = true;
            structuralSolution.NonLinearScheme.numberOfLoadSteps = 30;

            double[] externalForces3 = externalForcesStructuralVector;
            foreach (var dof in loadedStructuralDOFs)
            {
                //if (dof == 20 || dof == 26 || dof == 344 || dof == 350)
                //{
                //    externalForces3[dof - 1] = -externalStructuralLoad * xInterv2 * zInterv / 4.0;
                //}
                //else if (dof == 23 || dof == 347 ||
                //    dof % 27 == 20 || dof % 27 == 26)
                //{
                //    externalForces3[dof - 1] = -externalStructuralLoad * xInterv2 * zInterv / 2.0;
                //}
                if (dof == 20 || dof == 26)
                {
                    externalForces3[dof - 1] = -750.0 / 36.0;
                }
                else if (dof == 23)
                {
                    externalForces3[dof - 1] = -750.0 * 4.0 / 36.0;
                }
                else if (dof == 74 || dof == 128 || dof == 182 || dof == 236 || dof == 290 ||
                         dof == 80 || dof == 134 || dof == 188 || dof == 242 || dof == 296
                         || dof == 317 || dof == 323)
                {
                    externalForces3[dof - 1] = -750.0 * 2.0 / 36.0;
                }
                else if (dof == 47 || dof == 101 || dof == 155 || dof == 209 || dof == 263 ||
                         dof == 53 || dof == 107 || dof == 161 || dof == 215 || dof == 269)
                {
                    externalForces3[dof - 1] = -750.0 * 4.0 / 36.0;
                }
                else if (dof == 77 || dof == 131 || dof == 185 || dof == 239 || dof == 293 || dof == 320)
                {
                    externalForces3[dof - 1] = -750.0 * 8.0 / 36.0;
                }
                else if (dof == 50 || dof == 104 || dof == 158 || dof == 212 || dof == 266)
                {
                    externalForces3[dof - 1] = -750.0 * 16.0 / 36.0;
                }
                else if (dof == 1 || dof == 7 ||
                        dof == 19 || dof == 25)
                {
                    externalForces3[dof - 1] = -208.3125 - 1.0 / 48.0;
                }
                else if (dof == 4 || dof == 10 ||
                        dof == 16 || dof == 22)
                {
                    externalForces3[dof - 1] = -833.25 - 1.0 / 12.0;
                }
                else if (dof == 13)
                {
                    externalForces3[dof - 1] = -3333 - 1.0 / 3.0;
                }
                //else if (dof < 324)
                //{
                //    externalForces3[dof - 1] = -externalStructuralLoad * xInterv2 * zInterv;
                //}
                //else if (dof == 325 || dof == 331 ||
                //    dof == 343 || dof == 349)
                //{
                //    externalForces3[dof - 1] = 10 * externalStructuralLoad * yInterv * zInterv / 4.0;
                //}
                //else if (dof == 328 || dof == 334 ||
                //    dof == 340 || dof == 346)
                //{
                //    externalForces3[dof - 1] = 10 * externalStructuralLoad * yInterv * zInterv / 2.0;
                //}
                //else if (dof == 337)
                //{
                //    externalForces3[dof - 1] = 10 * externalStructuralLoad * yInterv * zInterv;
                //}
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
            for (int i = 0; i < 11; i++)
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
            //Dictionary<int, double[]> elementsInternalContactForcesVector;
            //for (int i = 1; i <= allStepsSolutions.Count; i++)
            //{
            //    elementsInternalContactForcesVector = new Dictionary<int, double[]>();
            //    elementsAssembly.UpdateDisplacements(allStepsSolutions[i]);
            //    for (int j = 1; j <= contactElements; j++)
            //    {
            //        elementsInternalContactForcesVector[elementsNumber + j] = elementsAssembly.ElementsAssembly[elementsNumber + j].CreateInternalGlobalForcesVector();
            //    }
            //    allStepsContactForces[i] = elementsInternalContactForcesVector;
            //    string name = "ContactForces" + i.ToString() + ".dat";
            //    double[] Vector = new double[contactElements * 12];
            //    int count = 0;
            //    for (int j = 1; j <= contactElements; j++)
            //    {
            //        Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == elementsNumber + j).Value[0];
            //        count += 1;
            //        Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == elementsNumber + j).Value[1];
            //        count += 1;
            //        Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == elementsNumber + j).Value[2];
            //        count += 1;
            //        Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == elementsNumber + j).Value[3];
            //        count += 1;
            //        Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == elementsNumber + j).Value[4];
            //        count += 1;
            //        Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == elementsNumber + j).Value[5];
            //        count += 1;
            //        Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == elementsNumber + j).Value[6];
            //        count += 1;
            //        Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == elementsNumber + j).Value[7];
            //        count += 1;
            //        Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == elementsNumber + j).Value[8];
            //        count += 1;
            //        Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == elementsNumber + j).Value[9];
            //        count += 1;
            //        Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == elementsNumber + j).Value[10];
            //        count += 1;
            //        Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == elementsNumber + j).Value[11];
            //        count += 1;
            //    }
            //    //VectorOperations.PrintVectorToFile(Vector, @"C:\Users\Public\Documents\" + name);
            //    if (i == allStepsSolutions.Count)
            //    {
            //        VectorOperations.PrintVectorToFile(Vector, @"C:\Users\Public\Documents\" + name);

            //    }
            //}

            for (int i = 1; i <= allStepsSolutions.Count; i++)
            {
                allStepsFullSolutions.Add(i, BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions.Single(m => m.Key == i).Value, elementsAssembly.BoundedDOFsVector));
                //int j = i + 1;
                string name = "solution" + i.ToString() + ".dat";
                //VectorOperations.PrintVectorToFile(allStepsFullSolutions.Single(m => m.Key == i + 1).Value, @"C:\Users\Public\Documents\" + name);
                VectorOperations.PrintVectorToFile(allStepsFullSolutions.Single(m => m.Key == i).Value, @"C:\Users\Public\Documents\" + name);
                //if (i == allStepsSolutions.Count - 1)
                //{
                //    VectorOperations.PrintVectorToFile(allStepsFullSolutions.Single(m => m.Key == i + 1).Value, @"C:\Users\Public\Documents\" + name);

                //}
                if (i == allStepsSolutions.Count)
                {
                    for (int k = 7; k <= 12; k++)
                    {
                        //string name2 = "Element" + j.ToString() + "Stress.dat";
                        //List<double[]> stress = elementsAssembly.ElementsAssembly[j].GetStressFromElementsNodes();
                        //VectorOperations.PrintListofVectorsToFile(stress, @"C:\Users\Public\Documents\" + name2);
                        string name2 = "Element" + k.ToString() + "Stress.dat";
                        string name3 = "Element" + k.ToString() + "Coordinates.dat";

                        List<double[]> stress = elementsAssembly.ElementsAssembly[k].GetStressFromElements(parametricCoordinatesVectors);
                        VectorOperations.PrintListofVectorsToFile(stress, @"C:\Users\Public\Documents\" + name2);
                        List<double[]> physicalSpaceCoordinates = elementsAssembly.ElementsAssembly[k].GetphysicalCoordinatesFromElements(parametricCoordinatesVectors);
                        VectorOperations.PrintListofVectorsToFile(physicalSpaceCoordinates, @"C:\Users\Public\Documents\" + name3);
                    }
                }
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
