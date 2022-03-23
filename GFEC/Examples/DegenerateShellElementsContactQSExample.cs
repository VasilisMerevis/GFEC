using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GFEC
{
    class DegenerateShellElementsContactQSExample
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
        const int ShellelementsNumber = 100;
        const int ContactElementsNumber = 116;
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
        const double externalStructuralLoad = -200000.0;

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
            for (int i = 0; i < nodesInY - 1; i++)
            {
                boundedDofs.Add(i * 5 + 3);
            }
            boundedDofs.Add(101);
            boundedDofs.Add(102);
            boundedDofs.Add(103);
            for (int i = 1; i < nodesInX - 1; i++)
            {
                boundedDofs.Add(i * nodesInY * 5 + 3);
                boundedDofs.Add((i * nodesInY + nodesInY - 1) * 5 + 3);
            }
            boundedDofs.Add(2101);
            boundedDofs.Add(2102);
            for (int i = 0; i < nodesInY; i++)
            {
                boundedDofs.Add((nodesInX - 1) * nodesInY * 5 + i * 5 + 3);
            }
            boundedDofs.Add(2201);
            boundedDofs.Add(2202);
            //boundedDofs.Add(2053);
            for(int i = 1; i <= 9; i++)
            {
                boundedDofs.Add(nodesNumberShellEmements * 5 / 2 + i);

            }
            for (int i = 1; i <= 9; i++)
            {
                boundedDofs.Add(nodesNumberShellEmements * 5 / 2 + 63 * 3 + i);

            }
            for (int i = 1; i <= 9; i++)
            {
                boundedDofs.Add(nodesNumberShellEmements * 5 / 2 + 126 * 3 + i);

            }
            structuralBoundaryConditions = boundedDofs.ToArray<int>();
        }
        private static void CreateStructuralLoadVector()
        {
            loadedStructuralDOFs = new List<int>();
            for (int i = 2; i <= 21; i++)
            {
                loadedStructuralDOFs.Add(nodesNumberShellEmements / 2 * 5 + i * 9);
            }
            for (int i = 2; i <= 21; i++)
            {
                loadedStructuralDOFs.Add(nodesNumberShellEmements / 2 * 5 + 189 + i * 9);
            }
            for (int i = 2; i <= 21; i++)
            {
                loadedStructuralDOFs.Add(nodesNumberShellEmements / 2 * 5 + 2 * 189 + i * 9);
            }
            externalForcesStructuralVector = new double[nodesNumberShellEmements / 2 * 5 +
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
            for (int i = 1; i <= nodesInX - 2; i += 2)
            {
                for (int j = 1; j <= nodesInY - 2; j += 2)
                {
                    int first = (i - 1) * nodesInY + j;
                    int first2 = first + nodesNumberShellEmements / 2;
                    connectivity[l] = new Dictionary<int, int>() { { 1, first }, { 2, first + nodesInY },
                            { 3,  first + nodesInY + nodesInY}, { 4,  first + nodesInY + nodesInY + 1},
                            { 5,  first + nodesInY + nodesInY + 2}, { 6,  first + nodesInY + 2},
                            { 7,  first + 2}, { 8,  first + 1},
                            { 9,  first + nodesInY + 1},
                            { 10, first2 }, { 11, first2 + nodesInY },
                            { 12,  first2 + nodesInY + nodesInY}, { 13,  first2 + nodesInY + nodesInY + 1},
                            { 14,  first2 + nodesInY + nodesInY + 2}, { 15,  first2 + nodesInY + 2},
                            { 16,  first2 + 2}, { 17,  first2 + 1},
                            { 18,  first2 + nodesInY + 1}};
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
            int masterNode1 = nodesNumberShellEmements/2 + 8 * 21 + 1;
            for (int i = nodesNumberShellEmements + 31; i <= nodesNumberShellEmements + 55; i += 6)
            {
                int slaveNode1 = i;
                connectivity[l] = new Dictionary<int, int>()
            {
                { 1, masterNode1 },
                { 2, masterNode1 + 21 },
                { 3, masterNode1 + 42 },
                { 4, masterNode1 + 42 + 1 },
                { 5, masterNode1 + 42 + 2},
                { 6, masterNode1 + 21 + 2},
                { 7, masterNode1 + 2 },
                { 8, masterNode1 + 1 },
                { 9, masterNode1 + 21 + 1 },
                { 10, slaveNode1 },
                { 11, slaveNode1 + 63 },
                { 12, slaveNode1 + 66 },
                { 13, slaveNode1 + 3}
            };
                l += 1;
                connectivity[l] = new Dictionary<int, int>()
            {
                { 1, masterNode1 + 42 },
                { 2, masterNode1 + 21 + 42 },
                { 3, masterNode1 + 42 + 42 },
                { 4, masterNode1 + 42 + 1 + 42 },
                { 5, masterNode1 + 42 + 2 + 42},
                { 6, masterNode1 + 21 + 2 + 42},
                { 7, masterNode1 + 2 + 42 },
                { 8, masterNode1 + 1 + 42 },
                { 9, masterNode1 + 21 + 1 + 42 },
                { 10, slaveNode1 },
                { 11, slaveNode1 + 63 },
                { 12, slaveNode1 + 66 },
                { 13, slaveNode1 + 3}
            };
                l += 1;
                //
                slaveNode1 += 63;
                connectivity[l] = new Dictionary<int, int>()
            {
                { 1, masterNode1 },
                { 2, masterNode1 + 21 },
                { 3, masterNode1 + 42 },
                { 4, masterNode1 + 42 + 1 },
                { 5, masterNode1 + 42 + 2},
                { 6, masterNode1 + 21 + 2},
                { 7, masterNode1 + 2 },
                { 8, masterNode1 + 1 },
                { 9, masterNode1 + 21 + 1 },
                { 10, slaveNode1 },
                { 11, slaveNode1 + 63 },
                { 12, slaveNode1 + 66 },
                { 13, slaveNode1 + 3}
            };
                l += 1;
                connectivity[l] = new Dictionary<int, int>()
            {
                { 1, masterNode1 + 42 },
                { 2, masterNode1 + 21 + 42 },
                { 3, masterNode1 + 42 + 42 },
                { 4, masterNode1 + 42 + 1 + 42 },
                { 5, masterNode1 + 42 + 2 + 42},
                { 6, masterNode1 + 21 + 2 + 42},
                { 7, masterNode1 + 2 + 42 },
                { 8, masterNode1 + 1 + 42 },
                { 9, masterNode1 + 21 + 1 + 42 },
                { 10, slaveNode1 },
                { 11, slaveNode1 + 63 },
                { 12, slaveNode1 + 66 },
                { 13, slaveNode1 + 3}
            };
                l += 1;
                //
                slaveNode1 = i + 3;
                connectivity[l] = new Dictionary<int, int>()
            {
                { 1, masterNode1 },
                { 2, masterNode1 + 21 },
                { 3, masterNode1 + 42 },
                { 4, masterNode1 + 42 + 1 },
                { 5, masterNode1 + 42 + 2},
                { 6, masterNode1 + 21 + 2},
                { 7, masterNode1 + 2 },
                { 8, masterNode1 + 1 },
                { 9, masterNode1 + 21 + 1 },
                { 10, slaveNode1 },
                { 11, slaveNode1 + 63 },
                { 12, slaveNode1 + 66 },
                { 13, slaveNode1 + 3}
            };
                l += 1;
                connectivity[l] = new Dictionary<int, int>()
            {
                { 1, masterNode1 + 42 },
                { 2, masterNode1 + 21 + 42 },
                { 3, masterNode1 + 42 + 42 },
                { 4, masterNode1 + 42 + 1 + 42 },
                { 5, masterNode1 + 42 + 2 + 42},
                { 6, masterNode1 + 21 + 2 + 42},
                { 7, masterNode1 + 2 + 42 },
                { 8, masterNode1 + 1 + 42 },
                { 9, masterNode1 + 21 + 1 + 42 },
                { 10, slaveNode1 },
                { 11, slaveNode1 + 63 },
                { 12, slaveNode1 + 66 },
                { 13, slaveNode1 + 3}
            };
                l += 1;
                //
                slaveNode1 += 63;
                connectivity[l] = new Dictionary<int, int>()
            {
                { 1, masterNode1 },
                { 2, masterNode1 + 21 },
                { 3, masterNode1 + 42 },
                { 4, masterNode1 + 42 + 1 },
                { 5, masterNode1 + 42 + 2},
                { 6, masterNode1 + 21 + 2},
                { 7, masterNode1 + 2 },
                { 8, masterNode1 + 1 },
                { 9, masterNode1 + 21 + 1 },
                { 10, slaveNode1 },
                { 11, slaveNode1 + 63 },
                { 12, slaveNode1 + 66 },
                { 13, slaveNode1 + 3}
            };
                l += 1;
                connectivity[l] = new Dictionary<int, int>()
            {
                { 1, masterNode1 + 42 },
                { 2, masterNode1 + 21 + 42 },
                { 3, masterNode1 + 42 + 42 },
                { 4, masterNode1 + 42 + 1 + 42 },
                { 5, masterNode1 + 42 + 2 + 42},
                { 6, masterNode1 + 21 + 2 + 42},
                { 7, masterNode1 + 2 + 42 },
                { 8, masterNode1 + 1 + 42 },
                { 9, masterNode1 + 21 + 1 + 42 },
                { 10, slaveNode1 },
                { 11, slaveNode1 + 63 },
                { 12, slaveNode1 + 66 },
                { 13, slaveNode1 + 3}
            };
                l += 1;
                masterNode1 += 2;
            }
            masterNode1 = nodesNumberShellEmements / 2 + 8 * 21 + 3;
            for (int i = nodesNumberShellEmements + 31; i <= nodesNumberShellEmements + 55; i += 6)
            {
                int slaveNode1 = i;
                connectivity[l] = new Dictionary<int, int>()
            {
                { 1, masterNode1 },
                { 2, masterNode1 + 21 },
                { 3, masterNode1 + 42 },
                { 4, masterNode1 + 42 + 1 },
                { 5, masterNode1 + 42 + 2},
                { 6, masterNode1 + 21 + 2},
                { 7, masterNode1 + 2 },
                { 8, masterNode1 + 1 },
                { 9, masterNode1 + 21 + 1 },
                { 10, slaveNode1 },
                { 11, slaveNode1 + 63 },
                { 12, slaveNode1 + 66 },
                { 13, slaveNode1 + 3}
            };
                l += 1;
                connectivity[l] = new Dictionary<int, int>()
            {
                { 1, masterNode1 + 42 },
                { 2, masterNode1 + 21 + 42 },
                { 3, masterNode1 + 42 + 42 },
                { 4, masterNode1 + 42 + 1 + 42 },
                { 5, masterNode1 + 42 + 2 + 42},
                { 6, masterNode1 + 21 + 2 + 42},
                { 7, masterNode1 + 2 + 42 },
                { 8, masterNode1 + 1 + 42 },
                { 9, masterNode1 + 21 + 1 + 42 },
                { 10, slaveNode1 },
                { 11, slaveNode1 + 63 },
                { 12, slaveNode1 + 66 },
                { 13, slaveNode1 + 3}
            };
                l += 1;
                //
                slaveNode1 += 63;
                connectivity[l] = new Dictionary<int, int>()
            {
                { 1, masterNode1 },
                { 2, masterNode1 + 21 },
                { 3, masterNode1 + 42 },
                { 4, masterNode1 + 42 + 1 },
                { 5, masterNode1 + 42 + 2},
                { 6, masterNode1 + 21 + 2},
                { 7, masterNode1 + 2 },
                { 8, masterNode1 + 1 },
                { 9, masterNode1 + 21 + 1 },
                { 10, slaveNode1 },
                { 11, slaveNode1 + 63 },
                { 12, slaveNode1 + 66 },
                { 13, slaveNode1 + 3}
            };
                l += 1;
                connectivity[l] = new Dictionary<int, int>()
            {
                { 1, masterNode1 + 42 },
                { 2, masterNode1 + 21 + 42 },
                { 3, masterNode1 + 42 + 42 },
                { 4, masterNode1 + 42 + 1 + 42 },
                { 5, masterNode1 + 42 + 2 + 42},
                { 6, masterNode1 + 21 + 2 + 42},
                { 7, masterNode1 + 2 + 42 },
                { 8, masterNode1 + 1 + 42 },
                { 9, masterNode1 + 21 + 1 + 42 },
                { 10, slaveNode1 },
                { 11, slaveNode1 + 63 },
                { 12, slaveNode1 + 66 },
                { 13, slaveNode1 + 3}
            };
                l += 1;
                //
                slaveNode1 = i + 3;
                connectivity[l] = new Dictionary<int, int>()
            {
                { 1, masterNode1 },
                { 2, masterNode1 + 21 },
                { 3, masterNode1 + 42 },
                { 4, masterNode1 + 42 + 1 },
                { 5, masterNode1 + 42 + 2},
                { 6, masterNode1 + 21 + 2},
                { 7, masterNode1 + 2 },
                { 8, masterNode1 + 1 },
                { 9, masterNode1 + 21 + 1 },
                { 10, slaveNode1 },
                { 11, slaveNode1 + 63 },
                { 12, slaveNode1 + 66 },
                { 13, slaveNode1 + 3}
            };
                l += 1;
                connectivity[l] = new Dictionary<int, int>()
            {
                { 1, masterNode1 + 42 },
                { 2, masterNode1 + 21 + 42 },
                { 3, masterNode1 + 42 + 42 },
                { 4, masterNode1 + 42 + 1 + 42 },
                { 5, masterNode1 + 42 + 2 + 42},
                { 6, masterNode1 + 21 + 2 + 42},
                { 7, masterNode1 + 2 + 42 },
                { 8, masterNode1 + 1 + 42 },
                { 9, masterNode1 + 21 + 1 + 42 },
                { 10, slaveNode1 },
                { 11, slaveNode1 + 63 },
                { 12, slaveNode1 + 66 },
                { 13, slaveNode1 + 3}
            };
                l += 1;
                //
                slaveNode1 += 63;
                connectivity[l] = new Dictionary<int, int>()
            {
                { 1, masterNode1 },
                { 2, masterNode1 + 21 },
                { 3, masterNode1 + 42 },
                { 4, masterNode1 + 42 + 1 },
                { 5, masterNode1 + 42 + 2},
                { 6, masterNode1 + 21 + 2},
                { 7, masterNode1 + 2 },
                { 8, masterNode1 + 1 },
                { 9, masterNode1 + 21 + 1 },
                { 10, slaveNode1 },
                { 11, slaveNode1 + 63 },
                { 12, slaveNode1 + 66 },
                { 13, slaveNode1 + 3}
            };
                l += 1;
                connectivity[l] = new Dictionary<int, int>()
            {
                { 1, masterNode1 + 42 },
                { 2, masterNode1 + 21 + 42 },
                { 3, masterNode1 + 42 + 42 },
                { 4, masterNode1 + 42 + 1 + 42 },
                { 5, masterNode1 + 42 + 2 + 42},
                { 6, masterNode1 + 21 + 2 + 42},
                { 7, masterNode1 + 2 + 42 },
                { 8, masterNode1 + 1 + 42 },
                { 9, masterNode1 + 21 + 1 + 42 },
                { 10, slaveNode1 },
                { 11, slaveNode1 + 63 },
                { 12, slaveNode1 + 66 },
                { 13, slaveNode1 + 3}
            };
                l += 1;
                masterNode1 += 2;
            }
            masterNode1 = nodesNumberShellEmements / 2 + 8 * 21 + 3;
            for (int i = nodesNumberShellEmements + 37; i <= nodesNumberShellEmements + 55; i += 6)
            {
                int slaveNode1 = i;
                connectivity[l] = new Dictionary<int, int>()
            {
                { 1, masterNode1 },
                { 2, masterNode1 + 21 },
                { 3, masterNode1 + 42 },
                { 4, masterNode1 + 42 + 1 },
                { 5, masterNode1 + 42 + 2},
                { 6, masterNode1 + 21 + 2},
                { 7, masterNode1 + 2 },
                { 8, masterNode1 + 1 },
                { 9, masterNode1 + 21 + 1 },
                { 10, slaveNode1 },
                { 11, slaveNode1 + 63 },
                { 12, slaveNode1 + 66 },
                { 13, slaveNode1 + 3}
            };
                l += 1;
                connectivity[l] = new Dictionary<int, int>()
            {
                { 1, masterNode1 + 42 },
                { 2, masterNode1 + 21 + 42 },
                { 3, masterNode1 + 42 + 42 },
                { 4, masterNode1 + 42 + 1 + 42 },
                { 5, masterNode1 + 42 + 2 + 42},
                { 6, masterNode1 + 21 + 2 + 42},
                { 7, masterNode1 + 2 + 42 },
                { 8, masterNode1 + 1 + 42 },
                { 9, masterNode1 + 21 + 1 + 42 },
                { 10, slaveNode1 },
                { 11, slaveNode1 + 63 },
                { 12, slaveNode1 + 66 },
                { 13, slaveNode1 + 3}
            };
                l += 1;
                //
                slaveNode1 += 63;
                connectivity[l] = new Dictionary<int, int>()
            {
                { 1, masterNode1 },
                { 2, masterNode1 + 21 },
                { 3, masterNode1 + 42 },
                { 4, masterNode1 + 42 + 1 },
                { 5, masterNode1 + 42 + 2},
                { 6, masterNode1 + 21 + 2},
                { 7, masterNode1 + 2 },
                { 8, masterNode1 + 1 },
                { 9, masterNode1 + 21 + 1 },
                { 10, slaveNode1 },
                { 11, slaveNode1 + 63 },
                { 12, slaveNode1 + 66 },
                { 13, slaveNode1 + 3}
            };
                l += 1;
                connectivity[l] = new Dictionary<int, int>()
            {
                { 1, masterNode1 + 42 },
                { 2, masterNode1 + 21 + 42 },
                { 3, masterNode1 + 42 + 42 },
                { 4, masterNode1 + 42 + 1 + 42 },
                { 5, masterNode1 + 42 + 2 + 42},
                { 6, masterNode1 + 21 + 2 + 42},
                { 7, masterNode1 + 2 + 42 },
                { 8, masterNode1 + 1 + 42 },
                { 9, masterNode1 + 21 + 1 + 42 },
                { 10, slaveNode1 },
                { 11, slaveNode1 + 63 },
                { 12, slaveNode1 + 66 },
                { 13, slaveNode1 + 3}
            };
                l += 1;
                //
                slaveNode1 = i + 3;
                connectivity[l] = new Dictionary<int, int>()
            {
                { 1, masterNode1 },
                { 2, masterNode1 + 21 },
                { 3, masterNode1 + 42 },
                { 4, masterNode1 + 42 + 1 },
                { 5, masterNode1 + 42 + 2},
                { 6, masterNode1 + 21 + 2},
                { 7, masterNode1 + 2 },
                { 8, masterNode1 + 1 },
                { 9, masterNode1 + 21 + 1 },
                { 10, slaveNode1 },
                { 11, slaveNode1 + 63 },
                { 12, slaveNode1 + 66 },
                { 13, slaveNode1 + 3}
            };
                l += 1;
                connectivity[l] = new Dictionary<int, int>()
            {
                { 1, masterNode1 + 42 },
                { 2, masterNode1 + 21 + 42 },
                { 3, masterNode1 + 42 + 42 },
                { 4, masterNode1 + 42 + 1 + 42 },
                { 5, masterNode1 + 42 + 2 + 42},
                { 6, masterNode1 + 21 + 2 + 42},
                { 7, masterNode1 + 2 + 42 },
                { 8, masterNode1 + 1 + 42 },
                { 9, masterNode1 + 21 + 1 + 42 },
                { 10, slaveNode1 },
                { 11, slaveNode1 + 63 },
                { 12, slaveNode1 + 66 },
                { 13, slaveNode1 + 3}
            };
                l += 1;
                //
                slaveNode1 += 63;
                connectivity[l] = new Dictionary<int, int>()
            {
                { 1, masterNode1 },
                { 2, masterNode1 + 21 },
                { 3, masterNode1 + 42 },
                { 4, masterNode1 + 42 + 1 },
                { 5, masterNode1 + 42 + 2},
                { 6, masterNode1 + 21 + 2},
                { 7, masterNode1 + 2 },
                { 8, masterNode1 + 1 },
                { 9, masterNode1 + 21 + 1 },
                { 10, slaveNode1 },
                { 11, slaveNode1 + 63 },
                { 12, slaveNode1 + 66 },
                { 13, slaveNode1 + 3}
            };
                l += 1;
                connectivity[l] = new Dictionary<int, int>()
            {
                { 1, masterNode1 + 42 },
                { 2, masterNode1 + 21 + 42 },
                { 3, masterNode1 + 42 + 42 },
                { 4, masterNode1 + 42 + 1 + 42 },
                { 5, masterNode1 + 42 + 2 + 42},
                { 6, masterNode1 + 21 + 2 + 42},
                { 7, masterNode1 + 2 + 42 },
                { 8, masterNode1 + 1 + 42 },
                { 9, masterNode1 + 21 + 1 + 42 },
                { 10, slaveNode1 },
                { 11, slaveNode1 + 63 },
                { 12, slaveNode1 + 66 },
                { 13, slaveNode1 + 3}
            };
                l += 1;
                masterNode1 += 2;
            }
            masterNode1 = nodesNumberShellEmements / 2 + 8 * 21 + 1;
            int sNode = nodesNumberShellEmements + 28;
            connectivity[l] = new Dictionary<int, int>()
            {
                { 1, masterNode1 },
                { 2, masterNode1 + 21 },
                { 3, masterNode1 + 42 },
                { 4, masterNode1 + 42 + 1 },
                { 5, masterNode1 + 42 + 2},
                { 6, masterNode1 + 21 + 2},
                { 7, masterNode1 + 2 },
                { 8, masterNode1 + 1 },
                { 9, masterNode1 + 21 + 1 },
                { 10, sNode },
                { 11, sNode + 63 },
                { 12, sNode + 66 },
                { 13, sNode + 3}
            };
            l += 1;
            sNode +=63;
            connectivity[l] = new Dictionary<int, int>()
            {
                { 1, masterNode1 },
                { 2, masterNode1 + 21 },
                { 3, masterNode1 + 42 },
                { 4, masterNode1 + 42 + 1 },
                { 5, masterNode1 + 42 + 2},
                { 6, masterNode1 + 21 + 2},
                { 7, masterNode1 + 2 },
                { 8, masterNode1 + 1 },
                { 9, masterNode1 + 21 + 1 },
                { 10, sNode },
                { 11, sNode + 63 },
                { 12, sNode + 66 },
                { 13, sNode + 3}
            };
            l += 1;
            //
            masterNode1 = nodesNumberShellEmements / 2 + 8 * 21 + 1 + 42;
            sNode = nodesNumberShellEmements + 28;
            connectivity[l] = new Dictionary<int, int>()
            {
                { 1, masterNode1 },
                { 2, masterNode1 + 21 },
                { 3, masterNode1 + 42 },
                { 4, masterNode1 + 42 + 1 },
                { 5, masterNode1 + 42 + 2},
                { 6, masterNode1 + 21 + 2},
                { 7, masterNode1 + 2 },
                { 8, masterNode1 + 1 },
                { 9, masterNode1 + 21 + 1 },
                { 10, sNode },
                { 11, sNode + 63 },
                { 12, sNode + 66 },
                { 13, sNode + 3}
            };
            l += 1;
            sNode += 63;
            connectivity[l] = new Dictionary<int, int>()
            {
                { 1, masterNode1 },
                { 2, masterNode1 + 21 },
                { 3, masterNode1 + 42 },
                { 4, masterNode1 + 42 + 1 },
                { 5, masterNode1 + 42 + 2},
                { 6, masterNode1 + 21 + 2},
                { 7, masterNode1 + 2 },
                { 8, masterNode1 + 1 },
                { 9, masterNode1 + 21 + 1 },
                { 10, sNode },
                { 11, sNode + 63 },
                { 12, sNode + 66 },
                { 13, sNode + 3}
            };
            return connectivity;
        }

        private static Dictionary<int, bool[]> CreateNodeFAT()
        {
            Dictionary<int, bool[]> nodeFAT = new Dictionary<int, bool[]>();
            for (int i = 1; i <= nodesNumberShellEmements / 2; i++)
            {
                nodeFAT[i] = new bool[] { true, true, true, true, true, false };
            }
            for (int i = nodesNumberShellEmements / 2 + 1; i <= nodesNumberShellEmements / 2 + nodesNumberSolidEmements; i++)
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
            string type = "IsoparamShell18";
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
                elementProperties[i] = new ElementProperties(E, A, type3, 10, 5, 1, 2);
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
            structuralSolution.NonLinearScheme.numberOfLoadSteps = 10;//

            double[] externalForces3 = externalForcesStructuralVector;
            int count = 1;
            foreach (var dof in loadedStructuralDOFs)
            {
                if (count < 20)
                {
                    externalForces3[dof - 1] = 2.0 * 0.01 * externalStructuralLoad;
                }
                else if(count == 20||
                    count == 60)
                {
                    externalForces3[dof - 1] = 0.01 * externalStructuralLoad;
                }
                else if (count == 40)
                {
                    externalForces3[dof - 1] = 2.0 * 0.01 * externalStructuralLoad;
                }
                else if (count > 20 && count <40)
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
                string name = "DegenElementsNLsolution" + i.ToString() + ".dat";
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