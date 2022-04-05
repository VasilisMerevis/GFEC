using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GFEC
{
    public static class BendingBeamContact3dWithFrictionRefinedMesh
    {

        public static ISolver structuralSolution;
        static int[] structuralBoundaryConditions;
        const double gap = 0.01;
        const double thickness = 0.1;
        const int nodesInX = 41;
        const int nodesInY = 3;
        const int nodesInZ = 3;
        const int nodesNumber = 738;
        const int elementsNumber = 320;
        const int contactElements = 472;

        //const double xInterv1 = 0.20;
        const double xInterv = 0.10;
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
        //const double contactArea = thickness * xInterv;

        //Friction coefficients
        const double miS = 0.20;
        const double miD = 0.20;


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
            //boundedDofs.Add(1081);
            boundedDofs.Add(1082);
            boundedDofs.Add(1083);
            //boundedDofs.Add(1084);
            boundedDofs.Add(1085);
            boundedDofs.Add(1086);
            //boundedDofs.Add(1087);
            boundedDofs.Add(1088);
            boundedDofs.Add(1089);
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
            boundedDofs.Add(2188);
            boundedDofs.Add(2189);
            boundedDofs.Add(2190);
            boundedDofs.Add(2191);
            boundedDofs.Add(2192);
            boundedDofs.Add(2193);
            boundedDofs.Add(2194);
            boundedDofs.Add(2195);
            boundedDofs.Add(2196);
            boundedDofs.Add(2197);
            boundedDofs.Add(2198);
            boundedDofs.Add(2199);
            boundedDofs.Add(2200);
            boundedDofs.Add(2201);
            boundedDofs.Add(2202);
            boundedDofs.Add(2203);
            boundedDofs.Add(2204);
            boundedDofs.Add(2205);
            boundedDofs.Add(2206);
            boundedDofs.Add(2207);
            boundedDofs.Add(2208);
            boundedDofs.Add(2209);
            boundedDofs.Add(2210);
            boundedDofs.Add(2211);
            boundedDofs.Add(2212);
            boundedDofs.Add(2213);
            boundedDofs.Add(2214);
            structuralBoundaryConditions = boundedDofs.ToArray<int>();
        }

        private static void CreateStructuralLoadVector()
        {
            loadedStructuralDOFs = new List<int>();
            loadedStructuralDOFs.Add(563);
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
                        nodes[l] = new Node(i * xInterv, j * yInterv, k * zInterv);
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
                        nodes[l] = new Node(offset + i * xInterv, j * yInterv - gap - (nodesInY - 1) * yInterv, k * zInterv);
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
                        int first = (i - 1) * nodesInY * nodesInZ + (j - 1) * nodesInZ + k;
                        connectivity[l] = new Dictionary<int, int>() { { 1, first }, { 2, first + nodesInY * nodesInZ },
                            { 3, first + nodesInY * nodesInZ + nodesInZ }, { 4, first + nodesInZ },
                            { 5, first + 1 },{ 6, first + nodesInY * nodesInZ + 1 },{ 7, first + (nodesInY + 1) * nodesInZ + 1 },
                            { 8, first + 1 + nodesInZ } };
                        l += 1;
                    }
                }
            }
            for (int i = 1; i <= nodesInX - 1; i++)
            {
                for (int j = 1; j <= nodesInY - 1; j++)
                {
                    for (int k = 1; k <= nodesInZ - 1; k++)
                    {
                        int first = (i - 1) * nodesInY * nodesInZ + (j - 1) * nodesInZ + k +
                            nodesInX * nodesInY * nodesInZ;
                        connectivity[l] = new Dictionary<int, int>() { { 1, first }, { 2, first + nodesInY * nodesInZ },
                            { 3, first + nodesInY * nodesInZ + nodesInZ }, { 4, first + nodesInZ },
                            { 5, first + 1 },{ 6, first + nodesInY * nodesInZ + 1 },{ 7, first + (nodesInY + 1) * nodesInZ + 1 },
                            { 8, first + 1 + nodesInZ } };
                        l += 1;
                    }
                }
            }
            //Contact elements
            for (int j = 0; j < 40; j++)
            {
                int firstSlaveNode = j * nodesInY * nodesInZ + 3;
                int firstMasterNode = firstSlaveNode + nodesInX * nodesInY * nodesInZ + 6;
                for (int i = 0; i <= 1; i++)
                {
                    connectivity[l] = new Dictionary<int, int>() { { 1, firstMasterNode }, { 2, firstMasterNode + nodesInY * nodesInZ }, { 3, firstMasterNode + nodesInY * nodesInZ - 1 },
                                                             {4, firstMasterNode - 1 },
                                                             { 5, firstSlaveNode }, { 6, firstSlaveNode + nodesInY * nodesInZ }, { 7, firstSlaveNode + nodesInY * nodesInZ - 1 },
                                                             {8, firstSlaveNode - 1 } };
                    l += 1;
                    firstSlaveNode += -1;
                    firstMasterNode = firstSlaveNode + nodesInX * nodesInY * nodesInZ + 6;
                }
            }

            for (int j = 0; j < 40; j++)
            {
                int firstSlaveNode = j * nodesInY * nodesInZ + 3;
                int firstMasterNode = firstSlaveNode + nodesInX * nodesInY * nodesInZ + 6 - 1;
                for (int i = 0; i <= 1; i++)
                {
                    connectivity[l] = new Dictionary<int, int>() { { 1, firstMasterNode }, { 2, firstMasterNode + nodesInY * nodesInZ }, { 3, firstMasterNode + nodesInY * nodesInZ - 1 },
                                                             {4, firstMasterNode - 1 },
                                                             { 5, firstSlaveNode }, { 6, firstSlaveNode + nodesInY * nodesInZ }, { 7, firstSlaveNode + nodesInY * nodesInZ - 1 },
                                                             {8, firstSlaveNode - 1 } };
                    l += 1;
                    firstSlaveNode += -1;
                    firstMasterNode = firstSlaveNode + nodesInX * nodesInY * nodesInZ + 6 + 1;
                }
            }

            for (int j = 0; j < 39; j++)
            {
                int firstSlaveNode = j * nodesInY * nodesInZ + 3;
                int firstMasterNode = firstSlaveNode + nodesInX * nodesInY * nodesInZ + 6 + nodesInY * nodesInZ;
                for (int i = 0; i <= 1; i++)
                {
                    connectivity[l] = new Dictionary<int, int>() { { 1, firstMasterNode }, { 2, firstMasterNode + nodesInY * nodesInZ }, { 3, firstMasterNode + nodesInY * nodesInZ - 1 },
                                                             {4, firstMasterNode - 1 },
                                                             { 5, firstSlaveNode }, { 6, firstSlaveNode + nodesInY * nodesInZ }, { 7, firstSlaveNode + nodesInY * nodesInZ - 1 },
                                                             {8, firstSlaveNode - 1 } };
                    l += 1;
                    firstSlaveNode += -1;
                    firstMasterNode = firstSlaveNode + nodesInX * nodesInY * nodesInZ + 6 + nodesInY * nodesInZ;
                }
            }

            for (int j = 0; j < 39; j++)
            {
                int firstSlaveNode = j * nodesInY * nodesInZ + 3;
                int firstMasterNode = firstSlaveNode + nodesInX * nodesInY * nodesInZ + 6 + nodesInY * nodesInZ - 1;
                for (int i = 0; i <= 1; i++)
                {
                    connectivity[l] = new Dictionary<int, int>() { { 1, firstMasterNode }, { 2, firstMasterNode + nodesInY * nodesInZ }, { 3, firstMasterNode + nodesInY * nodesInZ - 1 },
                                                             {4, firstMasterNode - 1 },
                                                             { 5, firstSlaveNode }, { 6, firstSlaveNode + nodesInY * nodesInZ }, { 7, firstSlaveNode + nodesInY * nodesInZ - 1 },
                                                             {8, firstSlaveNode - 1 } };
                    l += 1;
                    firstSlaveNode += -1;
                    firstMasterNode = firstSlaveNode + nodesInX * nodesInY * nodesInZ + 6 + nodesInY * nodesInZ + 1;
                }
            }

            for (int j = 1; j < 40; j++)
            {
                int firstSlaveNode = j * nodesInY * nodesInZ + 3;
                int firstMasterNode = firstSlaveNode + nodesInX * nodesInY * nodesInZ + 6 - nodesInY * nodesInZ;
                for (int i = 0; i <= 1; i++)
                {
                    connectivity[l] = new Dictionary<int, int>() { { 1, firstMasterNode }, { 2, firstMasterNode + nodesInY * nodesInZ }, { 3, firstMasterNode + nodesInY * nodesInZ - 1 },
                                                             {4, firstMasterNode - 1 },
                                                             { 5, firstSlaveNode }, { 6, firstSlaveNode + nodesInY * nodesInZ }, { 7, firstSlaveNode + nodesInY * nodesInZ - 1 },
                                                             {8, firstSlaveNode - 1 } };
                    l += 1;
                    firstSlaveNode += -1;
                    firstMasterNode = firstSlaveNode + nodesInX * nodesInY * nodesInZ + 6 - nodesInY * nodesInZ;
                }
            }

            for (int j = 1; j < 40; j++)
            {
                int firstSlaveNode = j * nodesInY * nodesInZ + 3;
                int firstMasterNode = firstSlaveNode + nodesInX * nodesInY * nodesInZ + 6 - nodesInY * nodesInZ - 1;
                for (int i = 0; i <= 1; i++)
                {
                    connectivity[l] = new Dictionary<int, int>() { { 1, firstMasterNode }, { 2, firstMasterNode + nodesInY * nodesInZ }, { 3, firstMasterNode + nodesInY * nodesInZ - 1 },
                                                             {4, firstMasterNode - 1 },
                                                             { 5, firstSlaveNode }, { 6, firstSlaveNode + nodesInY * nodesInZ }, { 7, firstSlaveNode + nodesInY * nodesInZ - 1 },
                                                             {8, firstSlaveNode - 1 } };
                    l += 1;
                    firstSlaveNode += -1;
                    firstMasterNode = firstSlaveNode + nodesInX * nodesInY * nodesInZ + 6 - nodesInY * nodesInZ + 1;
                }
            }

            //for (int j = 0; j < 38; j++)
            //{
            //    int firstSlaveNode = j * nodesInY * nodesInZ + 3;
            //    int firstMasterNode = firstSlaveNode + nodesInX * nodesInY * nodesInZ + 6 + 2 * nodesInY * nodesInZ;
            //    for (int i = 0; i <= 1; i++)
            //    {
            //        connectivity[l] = new Dictionary<int, int>() { { 1, firstMasterNode }, { 2, firstMasterNode + nodesInY * nodesInZ }, { 3, firstMasterNode + nodesInY * nodesInZ - 1 },
            //                                                 {4, firstMasterNode - 1 },
            //                                                 { 5, firstSlaveNode }, { 6, firstSlaveNode + nodesInY * nodesInZ }, { 7, firstSlaveNode + nodesInY * nodesInZ - 1 },
            //                                                 {8, firstSlaveNode - 1 } };
            //        l += 1;
            //        firstSlaveNode += -1;
            //        firstMasterNode = firstSlaveNode + nodesInX * nodesInY * nodesInZ + 6 + 2 * nodesInY * nodesInZ;
            //    }
            //}

            //for (int j = 0; j < 38; j++)
            //{
            //    int firstSlaveNode = j * nodesInY * nodesInZ + 3;
            //    int firstMasterNode = firstSlaveNode + nodesInX * nodesInY * nodesInZ + 6 + 2 * nodesInY * nodesInZ - 1;
            //    for (int i = 0; i <= 1; i++)
            //    {
            //        connectivity[l] = new Dictionary<int, int>() { { 1, firstMasterNode }, { 2, firstMasterNode + nodesInY * nodesInZ }, { 3, firstMasterNode + nodesInY * nodesInZ - 1 },
            //                                                 {4, firstMasterNode - 1 },
            //                                                 { 5, firstSlaveNode }, { 6, firstSlaveNode + nodesInY * nodesInZ }, { 7, firstSlaveNode + nodesInY * nodesInZ - 1 },
            //                                                 {8, firstSlaveNode - 1 } };
            //        l += 1;
            //        firstSlaveNode += -1;
            //        firstMasterNode = firstSlaveNode + nodesInX * nodesInY * nodesInZ + 6 + 2 * nodesInY * nodesInZ + 1;
            //    }
            //}

            //for (int j = 2; j < 40; j++)
            //{
            //    int firstSlaveNode = j * nodesInY * nodesInZ + 3;
            //    int firstMasterNode = firstSlaveNode + nodesInX * nodesInY * nodesInZ + 6 - 2 * nodesInY * nodesInZ;
            //    for (int i = 0; i <= 1; i++)
            //    {
            //        connectivity[l] = new Dictionary<int, int>() { { 1, firstMasterNode }, { 2, firstMasterNode + nodesInY * nodesInZ }, { 3, firstMasterNode + nodesInY * nodesInZ - 1 },
            //                                                 {4, firstMasterNode - 1 },
            //                                                 { 5, firstSlaveNode }, { 6, firstSlaveNode + nodesInY * nodesInZ }, { 7, firstSlaveNode + nodesInY * nodesInZ - 1 },
            //                                                 {8, firstSlaveNode - 1 } };
            //        l += 1;
            //        firstSlaveNode += -1;
            //        firstMasterNode = firstSlaveNode + nodesInX * nodesInY * nodesInZ + 6 - 2 * nodesInY * nodesInZ;
            //    }
            //}

            //for (int j = 2; j < 40; j++)
            //{
            //    int firstSlaveNode = j * nodesInY * nodesInZ + 3;
            //    int firstMasterNode = firstSlaveNode + nodesInX * nodesInY * nodesInZ + 6 - 2 * nodesInY * nodesInZ - 1;
            //    for (int i = 0; i <= 1; i++)
            //    {
            //        connectivity[l] = new Dictionary<int, int>() { { 1, firstMasterNode }, { 2, firstMasterNode + nodesInY * nodesInZ }, { 3, firstMasterNode + nodesInY * nodesInZ - 1 },
            //                                                 {4, firstMasterNode - 1 },
            //                                                 { 5, firstSlaveNode }, { 6, firstSlaveNode + nodesInY * nodesInZ }, { 7, firstSlaveNode + nodesInY * nodesInZ - 1 },
            //                                                 {8, firstSlaveNode - 1 } };
            //        l += 1;
            //        firstSlaveNode += -1;
            //        firstMasterNode = firstSlaveNode + nodesInX * nodesInY * nodesInZ + 6 - 2 * nodesInY * nodesInZ + 1;
            //    }
            //}
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

            string type = "Hex8";
            //string type2 = "ContactStS3Df";
            string type2 = "ContactStS3D";

            Dictionary<int, IElementProperties> elementProperties = new Dictionary<int, IElementProperties>();
            for (int i = 1; i <= elementsNumber; i++)
            {
                elementProperties[i] = new ElementProperties(E, poissonRatio, A, thickness, density, type);
            }
            for (int i = elementsNumber + 1; i <= elementsNumber + contactElements; i++)
            {
                //elementProperties[i] = new ElementProperties(E, type2, 2.0, 9, 1, 1, 2.0, miS, miD);
                elementProperties[i] = new ElementProperties(E, A, type2, 5.0, 9, 1, 1);
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
            //
            //elementsAssembly.SeperateContactDoF();
            //
            //ShowToGUI.PlotInitialGeometry(elementsAssembly);
            structuralSolution.LinearScheme = new PCGSolver();
            structuralSolution.NonLinearScheme.Tolerance = 1e-5;
            structuralSolution.ActivateNonLinearSolver = true;
            structuralSolution.NonLinearScheme.numberOfLoadSteps = 15;

            double[] externalForces3 = externalForcesStructuralVector;
            foreach (var dof in loadedStructuralDOFs)
            {
                externalForces3[dof - 1] = externalStructuralLoad;
            }
            double[] reducedExternalForces3 = BoundaryConditionsImposition.ReducedVector(externalForces3, elementsAssembly.BoundedDOFsVector);
            //double[] reducedExternalForces3 = elementsAssembly.MMCPGCreateReducedFromFullVector(externalForces3);

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
