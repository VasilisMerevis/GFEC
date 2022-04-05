using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GFEC
{
    class SolidShellElementsContactExample
    {
        public static ISolver structuralSolution;
        static int[] structuralBoundaryConditions;
        const double thickness = 0.006;
        const double uniformGap = 0.014;
        const double offsetX = 0.45;
        const double offsetY = -0.50;
        const int nodesInX = 21;
        const int nodesInY = 21;
        const int nodesNumberShellEmements = 882;
        const int ShellelementsNumber = 400;
        const int ContactElementsNumber = 180;
        //
        const int nodesInX2 = 3;
        const int nodesInY2 = 21;
        const int nodesInZ2 = 3;
        const int nodesNumberSolidEmements = 189;
        const int solidElementsNumber = 80;
        //const double xInterv1 = 0.20;
        const double xInterv = 0.05;
        const double yInterv = 0.05;
        const double zInterv = 0.05;

        //External loads
        const double externalStructuralLoad = -40000.0;

        static List<int> loadedStructuralDOFs;
        static double[] externalForcesStructuralVector;

        const double YoungMod = 200.0 * 1e9;
        const double YoungMod2 = 30.0 * 1e9;

        const double poissonRatio = 0.25;
        const double poissonRatio2 = 0.20;
        const double density = 8000.0;
        const double density2 = 2400.0;
        const double area = 1.0;
        private static void CreateStructuralBoundaryConditions()
        {
            List<int> boundedDofs = new List<int>();
            boundedDofs.Add(1);
            boundedDofs.Add(2);
            boundedDofs.Add(3);

            boundedDofs.Add(61);
            boundedDofs.Add(62);
            boundedDofs.Add(63);

            boundedDofs.Add(1261);
            boundedDofs.Add(1262);
            boundedDofs.Add(1263);

            boundedDofs.Add(1321);
            boundedDofs.Add(1322);
            boundedDofs.Add(1323);

            int add = nodesInX * nodesInY * 3;

            boundedDofs.Add(add + 1);
            boundedDofs.Add(add + 2);
            boundedDofs.Add(add + 3);

            boundedDofs.Add(add + 61);
            boundedDofs.Add(add + 62);
            boundedDofs.Add(add + 63);

            boundedDofs.Add(add + 1261);
            boundedDofs.Add(add + 1262);
            boundedDofs.Add(add + 1263);

            boundedDofs.Add(add + 1321);
            boundedDofs.Add(add + 1322);
            boundedDofs.Add(add + 1323);
            for (int i = 1; i <= 9; i++)
            {
                boundedDofs.Add(nodesNumberShellEmements * 3 + i);

            }
            for (int i = 1; i <= 9; i++)
            {
                boundedDofs.Add(nodesNumberShellEmements * 3 + 63 * 3 + i);

            }
            for (int i = 1; i <= 9; i++)
            {
                boundedDofs.Add(nodesNumberShellEmements * 3 + 126 * 3 + i);

            }
            structuralBoundaryConditions = boundedDofs.ToArray<int>();
        }
        private static void CreateStructuralLoadVector()
        {
            loadedStructuralDOFs = new List<int>();
            for (int i = 2; i <= 21; i++)
            {
                loadedStructuralDOFs.Add(nodesNumberShellEmements * 3 + i * 9);
            }
            for (int i = 2; i <= 21; i++)
            {
                loadedStructuralDOFs.Add(nodesNumberShellEmements * 3 + 189 + i * 9);
            }
            for (int i = 2; i <= 21; i++)
            {
                loadedStructuralDOFs.Add(nodesNumberShellEmements * 3 + 2 * 189 + i * 9);
            }
            externalForcesStructuralVector = new double[nodesNumberShellEmements * 3 +
                3 * nodesNumberSolidEmements];
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
                    nodes[l] = new Node(i * xInterv, j * yInterv, -thickness / 2.0);
                    l += 1;
                }
            }
            for (int i = 0; i < nodesInX; i++)
            {
                for (int j = 0; j < nodesInY; j++)
                {
                    nodes[l] = new Node(i * xInterv, j * yInterv, thickness / 2.0);
                    l += 1;
                }
            }
            for (int i = 0; i < nodesInX2; i++)
            {
                for (int j = 0; j < nodesInY2; j++)
                {
                    for (int k = 0; k < nodesInZ2; k++)
                    {
                        nodes[l] = new Node(offsetX + i * xInterv,
                            offsetY + j * yInterv,
                            thickness / 2.0 + uniformGap + k * zInterv);
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
            //shel elements
            for (int i = 1; i <= nodesInX - 1; i ++)
            {
                for (int j = 1; j <= nodesInY - 1; j ++)
                {
                    int first = (i - 1) * nodesInY + j;
                    int first2 = first + nodesNumberShellEmements / 2;
                    connectivity[l] = new Dictionary<int, int>()
                    {
                        { 1, first },
                        { 2, first + nodesInY },
                        { 3,  first + nodesInY + 1},
                        { 4,  first + 1},
                        { 5, first2 },
                        { 6, first2 + nodesInY },
                        { 7,  first2 + nodesInY + 1},
                        { 8,  first2 + 1}
                    };
                    l += 1;
                }
            }
            //solid elements
            for (int i = 1; i <= nodesInX2 - 1; i++)
            {
                for (int j = 1; j <= nodesInY2 - 1; j++)
                {
                    for (int k = 1; k <= nodesInZ2 - 1; k++)
                    {
                        int first = nodesNumberShellEmements +
                            (i - 1) * nodesInY2 * nodesInZ2 + (j - 1) * nodesInZ2 + k;
                        connectivity[l] = new Dictionary<int, int>() {
                            { 1, first },
                            { 2, first + nodesInY2 * nodesInZ2 },
                            { 3,  first + nodesInY2 * nodesInZ2 + nodesInZ2},
                            { 4,  first + nodesInZ2},
                            { 5, first + 1 },
                            { 6, first + nodesInY2 * nodesInZ2 + 1 },
                            { 7,  first + nodesInY2 * nodesInZ2 + nodesInZ2 + 1},
                            { 8,  first + nodesInZ2 + 1}};
                        l += 1;
                    }
                }
            }
            //contact elements
            int slaveNode1 = nodesNumberShellEmements + 28;
            for (int masterNode1 = nodesNumberShellEmements / 2 + 8 * 21 + 1;
                masterNode1 <= nodesNumberShellEmements / 2 + 8 * 21 + 11; masterNode1++)
            {
                connectivity[l] = new Dictionary<int, int>()
                    {
                        { 1, masterNode1 },
                        { 2, masterNode1 + nodesInY },
                        { 3, masterNode1 + nodesInY + 1 },
                        { 4, masterNode1 + 1 },
                        { 5, slaveNode1 },
                        { 6, slaveNode1 + 63 },
                        { 7, slaveNode1 + 66 },
                        { 8, slaveNode1 + 3}
                    };
                l += 1;
                slaveNode1 += 3;
            }
            slaveNode1 = nodesNumberShellEmements + 31;
            for (int masterNode1 = nodesNumberShellEmements / 2 + 8 * 21 + 1;
                masterNode1 <= nodesNumberShellEmements / 2 + 8 * 21 + 10; masterNode1++)
            {
                connectivity[l] = new Dictionary<int, int>()
                    {
                        { 1, masterNode1 },
                        { 2, masterNode1 + nodesInY },
                        { 3, masterNode1 + nodesInY + 1 },
                        { 4, masterNode1 + 1 },
                        { 5, slaveNode1 },
                        { 6, slaveNode1 + 63 },
                        { 7, slaveNode1 + 66 },
                        { 8, slaveNode1 + 3}
                    };
                l += 1;
                slaveNode1 += 3;
            }
            slaveNode1 = nodesNumberShellEmements + 34;
            for (int masterNode1 = nodesNumberShellEmements / 2 + 8 * 21 + 1;
                masterNode1 <= nodesNumberShellEmements / 2 + 8 * 21 + 9; masterNode1++)
            {
                connectivity[l] = new Dictionary<int, int>()
                    {
                        { 1, masterNode1 },
                        { 2, masterNode1 + nodesInY },
                        { 3, masterNode1 + nodesInY + 1 },
                        { 4, masterNode1 + 1 },
                        { 5, slaveNode1 },
                        { 6, slaveNode1 + 63 },
                        { 7, slaveNode1 + 66 },
                        { 8, slaveNode1 + 3}
                    };
                l += 1;
                slaveNode1 += 3;
            }
            //----------------------------------------------------------------------------------------------------
            slaveNode1 = nodesNumberShellEmements + 91;
            for (int masterNode1 = nodesNumberShellEmements / 2 + 11 * 21 + 1;
                masterNode1 <= nodesNumberShellEmements / 2 + 11 * 21 + 11; masterNode1++)
            {
                connectivity[l] = new Dictionary<int, int>()
                    {
                        { 1, masterNode1 },
                        { 2, masterNode1 + nodesInY },
                        { 3, masterNode1 + nodesInY + 1 },
                        { 4, masterNode1 + 1 },
                        { 5, slaveNode1 },
                        { 6, slaveNode1 + 63 },
                        { 7, slaveNode1 + 66 },
                        { 8, slaveNode1 + 3}
                    };
                l += 1;
                slaveNode1 += 3;
            }
            slaveNode1 = nodesNumberShellEmements + 94;
            for (int masterNode1 = nodesNumberShellEmements / 2 + 11 * 21 + 1;
                masterNode1 <= nodesNumberShellEmements / 2 + 11 * 21 + 10; masterNode1++)
            {
                connectivity[l] = new Dictionary<int, int>()
                    {
                        { 1, masterNode1 },
                        { 2, masterNode1 + nodesInY },
                        { 3, masterNode1 + nodesInY + 1 },
                        { 4, masterNode1 + 1 },
                        { 5, slaveNode1 },
                        { 6, slaveNode1 + 63 },
                        { 7, slaveNode1 + 66 },
                        { 8, slaveNode1 + 3}
                    };
                l += 1;
                slaveNode1 += 3;
            }
            slaveNode1 = nodesNumberShellEmements + 97;
            for (int masterNode1 = nodesNumberShellEmements / 2 + 11 * 21 + 1;
                masterNode1 <= nodesNumberShellEmements / 2 + 11 * 21 + 9; masterNode1++)
            {
                connectivity[l] = new Dictionary<int, int>()
                    {
                        { 1, masterNode1 },
                        { 2, masterNode1 + nodesInY },
                        { 3, masterNode1 + nodesInY + 1 },
                        { 4, masterNode1 + 1 },
                        { 5, slaveNode1 },
                        { 6, slaveNode1 + 63 },
                        { 7, slaveNode1 + 66 },
                        { 8, slaveNode1 + 3}
                    };
                l += 1;
                slaveNode1 += 3;
            }
            //----------------------------------------------------------------------------------------------------
            slaveNode1 = nodesNumberShellEmements + 28;
            for (int masterNode1 = nodesNumberShellEmements / 2 + 9 * 21 + 1;
                masterNode1 <= nodesNumberShellEmements / 2 + 9 * 21 + 11; masterNode1++)
            {
                connectivity[l] = new Dictionary<int, int>()
                    {
                        { 1, masterNode1 },
                        { 2, masterNode1 + nodesInY },
                        { 3, masterNode1 + nodesInY + 1 },
                        { 4, masterNode1 + 1 },
                        { 5, slaveNode1 },
                        { 6, slaveNode1 + 63 },
                        { 7, slaveNode1 + 66 },
                        { 8, slaveNode1 + 3}
                    };
                l += 1;
                slaveNode1 += 3;
            }
            slaveNode1 = nodesNumberShellEmements + 31;
            for (int masterNode1 = nodesNumberShellEmements / 2 + 9 * 21 + 1;
                masterNode1 <= nodesNumberShellEmements / 2 + 9 * 21 + 10; masterNode1++)
            {
                connectivity[l] = new Dictionary<int, int>()
                    {
                        { 1, masterNode1 },
                        { 2, masterNode1 + nodesInY },
                        { 3, masterNode1 + nodesInY + 1 },
                        { 4, masterNode1 + 1 },
                        { 5, slaveNode1 },
                        { 6, slaveNode1 + 63 },
                        { 7, slaveNode1 + 66 },
                        { 8, slaveNode1 + 3}
                    };
                l += 1;
                slaveNode1 += 3;
            }
            slaveNode1 = nodesNumberShellEmements + 34;
            for (int masterNode1 = nodesNumberShellEmements / 2 + 9 * 21 + 1;
                masterNode1 <= nodesNumberShellEmements / 2 + 9 * 21 + 9; masterNode1++)
            {
                connectivity[l] = new Dictionary<int, int>()
                    {
                        { 1, masterNode1 },
                        { 2, masterNode1 + nodesInY },
                        { 3, masterNode1 + nodesInY + 1 },
                        { 4, masterNode1 + 1 },
                        { 5, slaveNode1 },
                        { 6, slaveNode1 + 63 },
                        { 7, slaveNode1 + 66 },
                        { 8, slaveNode1 + 3}
                    };
                l += 1;
                slaveNode1 += 3;
            }
            //-----
            slaveNode1 = nodesNumberShellEmements + 91;
            for (int masterNode1 = nodesNumberShellEmements / 2 + 9 * 21 + 1;
                masterNode1 <= nodesNumberShellEmements / 2 + 9 * 21 + 11; masterNode1++)
            {
                connectivity[l] = new Dictionary<int, int>()
                    {
                        { 1, masterNode1 },
                        { 2, masterNode1 + nodesInY },
                        { 3, masterNode1 + nodesInY + 1 },
                        { 4, masterNode1 + 1 },
                        { 5, slaveNode1 },
                        { 6, slaveNode1 + 63 },
                        { 7, slaveNode1 + 66 },
                        { 8, slaveNode1 + 3}
                    };
                l += 1;
                slaveNode1 += 3;
            }
            slaveNode1 = nodesNumberShellEmements + 94;
            for (int masterNode1 = nodesNumberShellEmements / 2 + 9 * 21 + 1;
                masterNode1 <= nodesNumberShellEmements / 2 + 9 * 21 + 10; masterNode1++)
            {
                connectivity[l] = new Dictionary<int, int>()
                    {
                        { 1, masterNode1 },
                        { 2, masterNode1 + nodesInY },
                        { 3, masterNode1 + nodesInY + 1 },
                        { 4, masterNode1 + 1 },
                        { 5, slaveNode1 },
                        { 6, slaveNode1 + 63 },
                        { 7, slaveNode1 + 66 },
                        { 8, slaveNode1 + 3}
                    };
                l += 1;
                slaveNode1 += 3;
            }
            slaveNode1 = nodesNumberShellEmements + 97;
            for (int masterNode1 = nodesNumberShellEmements / 2 + 9 * 21 + 1;
                masterNode1 <= nodesNumberShellEmements / 2 + 9 * 21 + 9; masterNode1++)
            {
                connectivity[l] = new Dictionary<int, int>()
                    {
                        { 1, masterNode1 },
                        { 2, masterNode1 + nodesInY },
                        { 3, masterNode1 + nodesInY + 1 },
                        { 4, masterNode1 + 1 },
                        { 5, slaveNode1 },
                        { 6, slaveNode1 + 63 },
                        { 7, slaveNode1 + 66 },
                        { 8, slaveNode1 + 3}
                    };
                l += 1;
                slaveNode1 += 3;
            }
            //----------------------------------------------------------------------------------------------------
            slaveNode1 = nodesNumberShellEmements + 28;
            for (int masterNode1 = nodesNumberShellEmements / 2 + 10 * 21 + 1;
                masterNode1 <= nodesNumberShellEmements / 2 + 10 * 21 + 11; masterNode1++)
            {
                connectivity[l] = new Dictionary<int, int>()
                    {
                        { 1, masterNode1 },
                        { 2, masterNode1 + nodesInY },
                        { 3, masterNode1 + nodesInY + 1 },
                        { 4, masterNode1 + 1 },
                        { 5, slaveNode1 },
                        { 6, slaveNode1 + 63 },
                        { 7, slaveNode1 + 66 },
                        { 8, slaveNode1 + 3}
                    };
                l += 1;
                slaveNode1 += 3;
            }
            slaveNode1 = nodesNumberShellEmements + 31;
            for (int masterNode1 = nodesNumberShellEmements / 2 + 10 * 21 + 1;
                masterNode1 <= nodesNumberShellEmements / 2 + 10 * 21 + 10; masterNode1++)
            {
                connectivity[l] = new Dictionary<int, int>()
                    {
                        { 1, masterNode1 },
                        { 2, masterNode1 + nodesInY },
                        { 3, masterNode1 + nodesInY + 1 },
                        { 4, masterNode1 + 1 },
                        { 5, slaveNode1 },
                        { 6, slaveNode1 + 63 },
                        { 7, slaveNode1 + 66 },
                        { 8, slaveNode1 + 3}
                    };
                l += 1;
                slaveNode1 += 3;
            }
            slaveNode1 = nodesNumberShellEmements + 34;
            for (int masterNode1 = nodesNumberShellEmements / 2 + 10 * 21 + 1;
                masterNode1 <= nodesNumberShellEmements / 2 + 10 * 21 + 9; masterNode1++)
            {
                connectivity[l] = new Dictionary<int, int>()
                    {
                        { 1, masterNode1 },
                        { 2, masterNode1 + nodesInY },
                        { 3, masterNode1 + nodesInY + 1 },
                        { 4, masterNode1 + 1 },
                        { 5, slaveNode1 },
                        { 6, slaveNode1 + 63 },
                        { 7, slaveNode1 + 66 },
                        { 8, slaveNode1 + 3}
                    };
                l += 1;
                slaveNode1 += 3;
            }
            //-----
            slaveNode1 = nodesNumberShellEmements + 91;
            for (int masterNode1 = nodesNumberShellEmements / 2 + 10 * 21 + 1;
                masterNode1 <= nodesNumberShellEmements / 2 + 10 * 21 + 11; masterNode1++)
            {
                connectivity[l] = new Dictionary<int, int>()
                    {
                        { 1, masterNode1 },
                        { 2, masterNode1 + nodesInY },
                        { 3, masterNode1 + nodesInY + 1 },
                        { 4, masterNode1 + 1 },
                        { 5, slaveNode1 },
                        { 6, slaveNode1 + 63 },
                        { 7, slaveNode1 + 66 },
                        { 8, slaveNode1 + 3}
                    };
                l += 1;
                slaveNode1 += 3;
            }
            slaveNode1 = nodesNumberShellEmements + 94;
            for (int masterNode1 = nodesNumberShellEmements / 2 + 10 * 21 + 1;
                masterNode1 <= nodesNumberShellEmements / 2 + 10 * 21 + 10; masterNode1++)
            {
                connectivity[l] = new Dictionary<int, int>()
                    {
                        { 1, masterNode1 },
                        { 2, masterNode1 + nodesInY },
                        { 3, masterNode1 + nodesInY + 1 },
                        { 4, masterNode1 + 1 },
                        { 5, slaveNode1 },
                        { 6, slaveNode1 + 63 },
                        { 7, slaveNode1 + 66 },
                        { 8, slaveNode1 + 3}
                    };
                l += 1;
                slaveNode1 += 3;
            }
            slaveNode1 = nodesNumberShellEmements + 97;
            for (int masterNode1 = nodesNumberShellEmements / 2 + 10 * 21 + 1;
                masterNode1 <= nodesNumberShellEmements / 2 + 10 * 21 + 9; masterNode1++)
            {
                connectivity[l] = new Dictionary<int, int>()
                    {
                        { 1, masterNode1 },
                        { 2, masterNode1 + nodesInY },
                        { 3, masterNode1 + nodesInY + 1 },
                        { 4, masterNode1 + 1 },
                        { 5, slaveNode1 },
                        { 6, slaveNode1 + 63 },
                        { 7, slaveNode1 + 66 },
                        { 8, slaveNode1 + 3}
                    };
                l += 1;
                slaveNode1 += 3;
            }
            return connectivity;
        }

        private static Dictionary<int, bool[]> CreateNodeFAT()
        {
            Dictionary<int, bool[]> nodeFAT = new Dictionary<int, bool[]>();
            for (int i = 1; i <= nodesNumberShellEmements; i++)
            {
                nodeFAT[i] = new bool[] { true, true, true, false, false, false };
            }
            for (int i = nodesNumberShellEmements + 1; i <= nodesNumberShellEmements + nodesNumberSolidEmements; i++)
            {
                nodeFAT[i] = new bool[] { true, true, true, false, false, false };
            }
            return nodeFAT;
        }
        private static Dictionary<int, IElementProperties> CreateElementProperties()
        {
            double E = YoungMod;
            double E2 = YoungMod2;

            double A = area;
            string type = "ANSSolidShell8EAS";
            string type2 = "Hex8";
            string type3 = "ContactStS3D";

            Dictionary<int, IElementProperties> elementProperties = new Dictionary<int, IElementProperties>();
            for (int i = 1; i <= ShellelementsNumber; i++)
            {
                elementProperties[i] = new ElementProperties(E, poissonRatio, A, thickness, density, type);
            }
            for (int i = ShellelementsNumber + 1; i <= ShellelementsNumber + solidElementsNumber; i++)
            {
                elementProperties[i] = new ElementProperties(E2, poissonRatio2, A, thickness, density2, type2);
            }
            for (int i = ShellelementsNumber + solidElementsNumber + 1; i <= ShellelementsNumber + solidElementsNumber +
                ContactElementsNumber; i++)
            {
                elementProperties[i] = new ElementProperties(E, A, type3, 10, 5, 1, 1);
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
            ExportToFile.ExportMatlabInitialGeometry(elementsAssembly);
            elementsAssembly.ActivateBoundaryConditions = true;
            double[,] globalStiffnessMatrix = elementsAssembly.CreateTotalStiffnessMatrix();
            structuralSolution.LinearScheme = new LUFactorization();
            structuralSolution.ActivateNonLinearSolver = true;
            structuralSolution.NonLinearScheme.Tolerance = 1e-5;//
            structuralSolution.NonLinearScheme.numberOfLoadSteps = 20;//

            double[] externalForces3 = externalForcesStructuralVector;
            int count = 1;
            foreach (var dof in loadedStructuralDOFs)
            {
                if (count < 20)
                {
                    externalForces3[dof - 1] = 2.0 * 0.01 * externalStructuralLoad;
                }
                else if (count == 20 ||
                    count == 60)
                {
                    externalForces3[dof - 1] = 0.01 * externalStructuralLoad;
                }
                else if (count == 40)
                {
                    externalForces3[dof - 1] = 2.0 * 0.01 * externalStructuralLoad;
                }
                else if (count > 20 && count < 40)
                {
                    externalForces3[dof - 1] = 4.0 * 0.01 * externalStructuralLoad;
                }
                else if (count > 40 && count < 60)
                {
                    externalForces3[dof - 1] = 2.0 * 0.01 * externalStructuralLoad;
                }
                count += 1;
            }
            //double r = 0;
            //for(int i = 0; i< externalForces3.GetLength(0); i++)
            //{
            //    r += externalForces3[i];
            //}
            double[] reducedExternalForces3 = BoundaryConditionsImposition.ReducedVector(externalForces3, elementsAssembly.BoundedDOFsVector);
            structuralSolution.AssemblyData = elementsAssembly;
            structuralSolution.Solve(reducedExternalForces3);
            Dictionary<int, double[]> solvectors = structuralSolution.GetAllStepsSolutions();

            //double[]  = structuralSolution.GetSolution();
            for (int i = solvectors.Keys.Min(); i <= solvectors.Keys.Max(); i++)
            {
                double[] fullSolVector = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(solvectors.Single(s => s.Key == i).Value,
                elementsAssembly.BoundedDOFsVector);
                string name = "SolidElementsContactsolution" + i.ToString() + ".dat";
                VectorOperations.PrintVectorToFile(fullSolVector, @"C:\Users\Public\Documents\" + name);
            }
            //double[] fullSolVector = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(solvector,
            //    elementsAssembly.BoundedDOFsVector);

            //string name = "DegenElementssolution" + ".dat";
            //VectorOperations.PrintVectorToFile(fullSolVector, @"C:\Users\Public\Documents\" + name);
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