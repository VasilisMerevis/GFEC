using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GFEC
{
    public static class Blocks3dContactSlidingMeshRefined
    {

        public static ISolver structuralSolution;
        static int[] structuralBoundaryConditions;
        const double gap = 0.00001;
        const double thickness = 0.1;
        const int nodesInX = 25;
        const int nodesInY = 5;
        const int nodesInZ = 5;
        const int nodesNumber = 1250;
        const int elementsNumber = 768;
        const int contactElements = 700;

        //const double xInterv1 = 0.20;
        const double xInterv2 = 0.125;
        const double yInterv = 0.125;
        const double zInterv = 0.125;

        const double offset = 0.0;

        //External loads
        const double externalStructuralLoad = 1000.0;

        static List<int> loadedStructuralDOFs;
        static double[] externalForcesStructuralVector;

        const double YoungMod = 1.0 * 1e6;

        const double poissonRatio = 0.25;
        const double density = 8000.0;
        const double area = 1.0;
        const double contactArea = thickness * xInterv2;

        //Friction coefficients
        //const double miS = 1.50;
        //const double miD = 1.50;


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
            boundedDofs.Add(10);
            boundedDofs.Add(11);
            boundedDofs.Add(12);
            boundedDofs.Add(13);
            boundedDofs.Add(14);
            boundedDofs.Add(15);
            boundedDofs.Add(16);
            boundedDofs.Add(17);
            boundedDofs.Add(18);
            boundedDofs.Add(28);
            boundedDofs.Add(29);
            boundedDofs.Add(30);
            boundedDofs.Add(31);
            boundedDofs.Add(32);
            boundedDofs.Add(33);
            boundedDofs.Add(43);
            boundedDofs.Add(44);
            boundedDofs.Add(45);
            boundedDofs.Add(46);
            boundedDofs.Add(47);
            boundedDofs.Add(48);
            boundedDofs.Add(58);
            boundedDofs.Add(59);
            boundedDofs.Add(60);
            boundedDofs.Add(61);
            boundedDofs.Add(62); 
            boundedDofs.Add(63);
            boundedDofs.Add(64);
            boundedDofs.Add(65);
            boundedDofs.Add(66);
            boundedDofs.Add(67);
            boundedDofs.Add(68);
            boundedDofs.Add(69);
            boundedDofs.Add(70);
            boundedDofs.Add(71);
            boundedDofs.Add(72);
            boundedDofs.Add(73);
            boundedDofs.Add(74);
            boundedDofs.Add(75);

            //boundedDofs.Add(1801);
            boundedDofs.Add(1802);
            boundedDofs.Add(1803);
            //boundedDofs.Add(1804);
            boundedDofs.Add(1805);
            boundedDofs.Add(1806);
            //boundedDofs.Add(1807);
            boundedDofs.Add(1808);
            boundedDofs.Add(1809);
            //boundedDofs.Add(1810);
            boundedDofs.Add(1811);
            boundedDofs.Add(1812);
            //boundedDofs.Add(1813);
            boundedDofs.Add(1814);
            boundedDofs.Add(1815);
            for (int i = 1; i <= 75; i++)
            {
                boundedDofs.Add(1875 + i);
            }
            for (int i = 1; i <= 75; i++)
            {
                boundedDofs.Add(3675 + i);
            }
            structuralBoundaryConditions = boundedDofs.ToArray<int>();
        }

        private static void CreateStructuralLoadVector()
        {
            loadedStructuralDOFs = new List<int>();
            for (int i = 1; i < 25; i++)
            {
                loadedStructuralDOFs.Add(75 * i + 62);
                loadedStructuralDOFs.Add(75 * i + 65);
                loadedStructuralDOFs.Add(75 * i + 68);
                loadedStructuralDOFs.Add(75 * i + 71);
                loadedStructuralDOFs.Add(75 * i + 74);
            }
            for (int i = 0; i < 25; i++)
            {
                loadedStructuralDOFs.Add(75 * 24 + i * 3 + 1);
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
            for (int j = 0; j < 24; j++)
            {
                int firstSlaveNode = j * nodesInY * nodesInZ + 3;
                int firstMasterNode = firstSlaveNode + nodesInX * nodesInY * nodesInZ + 20;
                for (int i = 0; i <= 1; i++)
                {
                    connectivity[l] = new Dictionary<int, int>() { { 1, firstMasterNode - 1 }, { 2, firstMasterNode }, { 3, firstMasterNode + nodesInY * nodesInZ },
                                                             {4, firstMasterNode + nodesInY * nodesInZ - 1 },
                                                             { 5, firstSlaveNode - 1 }, { 6, firstSlaveNode }, { 7, firstSlaveNode + nodesInY * nodesInZ},
                                                             {8, firstSlaveNode + nodesInY * nodesInZ - 1 } };
                    l += 1;
                    firstSlaveNode += -1;
                    firstMasterNode = firstSlaveNode + nodesInX * nodesInY * nodesInZ + 20;
                }
            }

            for (int j = 0; j < 24; j++)
            {
                int firstSlaveNode = j * nodesInY * nodesInZ + 3;
                int firstMasterNode = firstSlaveNode + nodesInX * nodesInY * nodesInZ + 20 - 1;
                for (int i = 0; i <= 1; i++)
                {
                    connectivity[l] = new Dictionary<int, int>() { { 1, firstMasterNode - 1 }, { 2, firstMasterNode }, { 3, firstMasterNode + nodesInY * nodesInZ },
                                                             {4, firstMasterNode + nodesInY * nodesInZ - 1 },
                                                             { 5, firstSlaveNode - 1 }, { 6, firstSlaveNode }, { 7, firstSlaveNode + nodesInY * nodesInZ},
                                                             {8, firstSlaveNode + nodesInY * nodesInZ - 1 } };
                    l += 1;
                    firstSlaveNode += -1;
                    firstMasterNode = firstSlaveNode + nodesInX * nodesInY * nodesInZ + 20 + 1;
                }
            }

            for (int j = 0; j < 23; j++)
            {
                int firstSlaveNode = j * nodesInY * nodesInZ + 3;
                int firstMasterNode = firstSlaveNode + nodesInX * nodesInY * nodesInZ + 20 + nodesInY * nodesInZ;
                for (int i = 0; i <= 1; i++)
                {
                    connectivity[l] = new Dictionary<int, int>() { { 1, firstMasterNode - 1 }, { 2, firstMasterNode }, { 3, firstMasterNode + nodesInY * nodesInZ },
                                                             {4, firstMasterNode + nodesInY * nodesInZ - 1 },
                                                             { 5, firstSlaveNode - 1 }, { 6, firstSlaveNode }, { 7, firstSlaveNode + nodesInY * nodesInZ},
                                                             {8, firstSlaveNode + nodesInY * nodesInZ - 1 } };
                    l += 1;
                    firstSlaveNode += -1;
                    firstMasterNode = firstSlaveNode + nodesInX * nodesInY * nodesInZ + 20 + nodesInY * nodesInZ;
                }
            }

            for (int j = 0; j < 23; j++)
            {
                int firstSlaveNode = j * nodesInY * nodesInZ + 3;
                int firstMasterNode = firstSlaveNode + nodesInX * nodesInY * nodesInZ + 20 + nodesInY * nodesInZ - 1;
                for (int i = 0; i <= 1; i++)
                {
                    connectivity[l] = new Dictionary<int, int>() { { 1, firstMasterNode - 1 }, { 2, firstMasterNode }, { 3, firstMasterNode + nodesInY * nodesInZ },
                                                             {4, firstMasterNode + nodesInY * nodesInZ - 1 },
                                                             { 5, firstSlaveNode - 1 }, { 6, firstSlaveNode }, { 7, firstSlaveNode + nodesInY * nodesInZ},
                                                             {8, firstSlaveNode + nodesInY * nodesInZ - 1 } };
                    l += 1;
                    firstSlaveNode += -1;
                    firstMasterNode = firstSlaveNode + nodesInX * nodesInY * nodesInZ + 20 + nodesInY * nodesInZ + 1;
                }
            }

            for (int j = 1; j < 24; j++)
            {
                int firstSlaveNode = j * nodesInY * nodesInZ + 3;
                int firstMasterNode = firstSlaveNode + nodesInX * nodesInY * nodesInZ + 20 - nodesInY * nodesInZ;
                for (int i = 0; i <= 1; i++)
                {
                    connectivity[l] = new Dictionary<int, int>() { { 1, firstMasterNode - 1 }, { 2, firstMasterNode }, { 3, firstMasterNode + nodesInY * nodesInZ },
                                                             {4, firstMasterNode + nodesInY * nodesInZ - 1 },
                                                             { 5, firstSlaveNode - 1 }, { 6, firstSlaveNode }, { 7, firstSlaveNode + nodesInY * nodesInZ},
                                                             {8, firstSlaveNode + nodesInY * nodesInZ - 1 } };
                    l += 1;
                    firstSlaveNode += -1;
                    firstMasterNode = firstSlaveNode + nodesInX * nodesInY * nodesInZ + 20 - nodesInY * nodesInZ;
                }
            }

            for (int j = 1; j < 24; j++)
            {
                int firstSlaveNode = j * nodesInY * nodesInZ + 3;
                int firstMasterNode = firstSlaveNode + nodesInX * nodesInY * nodesInZ + 20 - nodesInY * nodesInZ - 1;
                for (int i = 0; i <= 1; i++)
                {
                    connectivity[l] = new Dictionary<int, int>() { { 1, firstMasterNode - 1 }, { 2, firstMasterNode }, { 3, firstMasterNode + nodesInY * nodesInZ },
                                                             {4, firstMasterNode + nodesInY * nodesInZ - 1 },
                                                             { 5, firstSlaveNode - 1 }, { 6, firstSlaveNode }, { 7, firstSlaveNode + nodesInY * nodesInZ},
                                                             {8, firstSlaveNode + nodesInY * nodesInZ - 1 } };
                    l += 1;
                    firstSlaveNode += -1;
                    firstMasterNode = firstSlaveNode + nodesInX * nodesInY * nodesInZ + 20 - nodesInY * nodesInZ + 1;
                }
            }

            for (int j = 0; j < 24; j++)
            {
                int firstSlaveNode = j * nodesInY * nodesInZ + 5;
                int firstMasterNode = firstSlaveNode + nodesInX * nodesInY * nodesInZ + 20;
                for (int i = 0; i <= 1; i++)
                {
                    connectivity[l] = new Dictionary<int, int>() { { 1, firstMasterNode - 1 }, { 2, firstMasterNode }, { 3, firstMasterNode + nodesInY * nodesInZ },
                                                             {4, firstMasterNode + nodesInY * nodesInZ - 1 },
                                                             { 5, firstSlaveNode - 1 }, { 6, firstSlaveNode }, { 7, firstSlaveNode + nodesInY * nodesInZ},
                                                             {8, firstSlaveNode + nodesInY * nodesInZ - 1 } };
                    l += 1;
                    firstSlaveNode += -1;
                    firstMasterNode = firstSlaveNode + nodesInX * nodesInY * nodesInZ + 20;
                }
            }

            for (int j = 0; j < 24; j++)
            {
                int firstSlaveNode = j * nodesInY * nodesInZ + 5;
                int firstMasterNode = firstSlaveNode + nodesInX * nodesInY * nodesInZ + 20 - 1;
                for (int i = 0; i <= 1; i++)
                {
                    connectivity[l] = new Dictionary<int, int>() { { 1, firstMasterNode - 1 }, { 2, firstMasterNode }, { 3, firstMasterNode + nodesInY * nodesInZ },
                                                             {4, firstMasterNode + nodesInY * nodesInZ - 1 },
                                                             { 5, firstSlaveNode - 1 }, { 6, firstSlaveNode }, { 7, firstSlaveNode + nodesInY * nodesInZ},
                                                             {8, firstSlaveNode + nodesInY * nodesInZ - 1 } };
                    l += 1;
                    firstSlaveNode += -1;
                    firstMasterNode = firstSlaveNode + nodesInX * nodesInY * nodesInZ + 20 + 1;
                }
            }

            for (int j = 0; j < 23; j++)
            {
                int firstSlaveNode = j * nodesInY * nodesInZ + 5;
                int firstMasterNode = firstSlaveNode + nodesInX * nodesInY * nodesInZ + 20 + nodesInY * nodesInZ;
                for (int i = 0; i <= 1; i++)
                {
                    connectivity[l] = new Dictionary<int, int>() { { 1, firstMasterNode - 1 }, { 2, firstMasterNode }, { 3, firstMasterNode + nodesInY * nodesInZ },
                                                             {4, firstMasterNode + nodesInY * nodesInZ - 1 },
                                                             { 5, firstSlaveNode - 1 }, { 6, firstSlaveNode }, { 7, firstSlaveNode + nodesInY * nodesInZ},
                                                             {8, firstSlaveNode + nodesInY * nodesInZ - 1 } };
                    l += 1;
                    firstSlaveNode += -1;
                    firstMasterNode = firstSlaveNode + nodesInX * nodesInY * nodesInZ + 20 + nodesInY * nodesInZ;
                }
            }

            for (int j = 0; j < 23; j++)
            {
                int firstSlaveNode = j * nodesInY * nodesInZ + 5;
                int firstMasterNode = firstSlaveNode + nodesInX * nodesInY * nodesInZ + 20 + nodesInY * nodesInZ - 1;
                for (int i = 0; i <= 1; i++)
                {
                    connectivity[l] = new Dictionary<int, int>() { { 1, firstMasterNode - 1 }, { 2, firstMasterNode }, { 3, firstMasterNode + nodesInY * nodesInZ },
                                                             {4, firstMasterNode + nodesInY * nodesInZ - 1 },
                                                             { 5, firstSlaveNode - 1 }, { 6, firstSlaveNode }, { 7, firstSlaveNode + nodesInY * nodesInZ},
                                                             {8, firstSlaveNode + nodesInY * nodesInZ - 1 } };
                    l += 1;
                    firstSlaveNode += -1;
                    firstMasterNode = firstSlaveNode + nodesInX * nodesInY * nodesInZ + 20 + nodesInY * nodesInZ + 1;
                }
            }

            for (int j = 1; j < 24; j++)
            {
                int firstSlaveNode = j * nodesInY * nodesInZ + 5;
                int firstMasterNode = firstSlaveNode + nodesInX * nodesInY * nodesInZ + 20 - nodesInY * nodesInZ;
                for (int i = 0; i <= 1; i++)
                {
                    connectivity[l] = new Dictionary<int, int>() { { 1, firstMasterNode - 1 }, { 2, firstMasterNode }, { 3, firstMasterNode + nodesInY * nodesInZ },
                                                             {4, firstMasterNode + nodesInY * nodesInZ - 1 },
                                                             { 5, firstSlaveNode - 1 }, { 6, firstSlaveNode }, { 7, firstSlaveNode + nodesInY * nodesInZ},
                                                             {8, firstSlaveNode + nodesInY * nodesInZ - 1 } };
                    l += 1;
                    firstSlaveNode += -1;
                    firstMasterNode = firstSlaveNode + nodesInX * nodesInY * nodesInZ + 20 - nodesInY * nodesInZ;
                }
            }

            for (int j = 1; j < 24; j++)
            {
                int firstSlaveNode = j * nodesInY * nodesInZ + 5;
                int firstMasterNode = firstSlaveNode + nodesInX * nodesInY * nodesInZ + 20 - nodesInY * nodesInZ - 1;
                for (int i = 0; i <= 1; i++)
                {
                    connectivity[l] = new Dictionary<int, int>() { { 1, firstMasterNode - 1 }, { 2, firstMasterNode }, { 3, firstMasterNode + nodesInY * nodesInZ },
                                                             {4, firstMasterNode + nodesInY * nodesInZ - 1 },
                                                             { 5, firstSlaveNode - 1 }, { 6, firstSlaveNode }, { 7, firstSlaveNode + nodesInY * nodesInZ},
                                                             {8, firstSlaveNode + nodesInY * nodesInZ - 1 } };
                    l += 1;
                    firstSlaveNode += -1;
                    firstMasterNode = firstSlaveNode + nodesInX * nodesInY * nodesInZ + 20 - nodesInY * nodesInZ + 1;
                }
            }

            for (int j = 0; j < 24; j++)
            {
                int firstSlaveNode = j * nodesInY * nodesInZ + 4;
                int firstMasterNode = firstSlaveNode + nodesInX * nodesInY * nodesInZ + 20 - 1;
                for (int i = 0; i <= 1; i++)
                {
                    connectivity[l] = new Dictionary<int, int>() { { 1, firstMasterNode - 1 }, { 2, firstMasterNode }, { 3, firstMasterNode + nodesInY * nodesInZ },
                                                             {4, firstMasterNode + nodesInY * nodesInZ - 1 },
                                                             { 5, firstSlaveNode - 1 }, { 6, firstSlaveNode }, { 7, firstSlaveNode + nodesInY * nodesInZ},
                                                             {8, firstSlaveNode + nodesInY * nodesInZ - 1 } };
                    l += 1;
                    firstSlaveNode += -1;
                    firstMasterNode = firstSlaveNode + nodesInX * nodesInY * nodesInZ + 20 + 1;
                }
            }

            for (int j = 0; j < 23; j++)
            {
                int firstSlaveNode = j * nodesInY * nodesInZ + 4;
                int firstMasterNode = firstSlaveNode + nodesInX * nodesInY * nodesInZ + 20 + nodesInY * nodesInZ - 1;
                for (int i = 0; i <= 1; i++)
                {
                    connectivity[l] = new Dictionary<int, int>() { { 1, firstMasterNode - 1 }, { 2, firstMasterNode }, { 3, firstMasterNode + nodesInY * nodesInZ },
                                                             {4, firstMasterNode + nodesInY * nodesInZ - 1 },
                                                             { 5, firstSlaveNode - 1 }, { 6, firstSlaveNode }, { 7, firstSlaveNode + nodesInY * nodesInZ},
                                                             {8, firstSlaveNode + nodesInY * nodesInZ - 1 } };
                    l += 1;
                    firstSlaveNode += -1;
                    firstMasterNode = firstSlaveNode + nodesInX * nodesInY * nodesInZ + 20 + nodesInY * nodesInZ + 1;
                }
            }

            for (int j = 1; j < 24; j++)
            {
                int firstSlaveNode = j * nodesInY * nodesInZ + 4;
                int firstMasterNode = firstSlaveNode + nodesInX * nodesInY * nodesInZ + 20 - nodesInY * nodesInZ - 1;
                for (int i = 0; i <= 1; i++)
                {
                    connectivity[l] = new Dictionary<int, int>() { { 1, firstMasterNode - 1 }, { 2, firstMasterNode }, { 3, firstMasterNode + nodesInY * nodesInZ },
                                                             {4, firstMasterNode + nodesInY * nodesInZ - 1 },
                                                             { 5, firstSlaveNode - 1 }, { 6, firstSlaveNode }, { 7, firstSlaveNode + nodesInY * nodesInZ},
                                                             {8, firstSlaveNode + nodesInY * nodesInZ - 1 } };
                    l += 1;
                    firstSlaveNode += -1;
                    firstMasterNode = firstSlaveNode + nodesInX * nodesInY * nodesInZ + 20 - nodesInY * nodesInZ + 1;
                }
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

            string type = "Hex8";
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
                elementProperties[i] = new ElementProperties(E, A, type3, 5.0, 9, 1, 1);
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
            structuralSolution.NonLinearScheme.numberOfLoadSteps = 15;

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
                if (dof == 1862 || dof == 1874)
                {
                    externalForces3[dof - 1] = -externalStructuralLoad * xInterv2 * zInterv / 4.0;
                }
                else if (dof == 1865 || dof == 1868 || dof == 1871 ||
                    dof % 27 == 62 || dof % 75 == 74)
                {
                    externalForces3[dof - 1] = -externalStructuralLoad * xInterv2 * zInterv / 2.0;
                }
                else if (dof < 1800)
                {
                    externalForces3[dof - 1] = -externalStructuralLoad * xInterv2 * zInterv;
                }
                else if (dof == 1801 || dof == 1813 ||
                    dof == 1861 || dof == 1873)
                {
                    externalForces3[dof - 1] = 10 * externalStructuralLoad * yInterv * zInterv / 4.0;
                }
                else if (dof == 1804 || dof == 1807 ||
                    dof == 1810 || dof == 1864 || dof == 1867 || dof == 1870 || dof == 1816 ||
                    dof == 1831 || dof == 1846 || dof == 1828 || dof == 1843 || dof == 1858)
                {
                    externalForces3[dof - 1] = 10 * externalStructuralLoad * yInterv * zInterv / 2.0;
                }
                else
                {
                    externalForces3[dof - 1] = 10 * externalStructuralLoad * yInterv * zInterv;
                }
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
