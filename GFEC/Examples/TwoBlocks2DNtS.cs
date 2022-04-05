﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GFEC
{
    public static class TwoBlocks2DNtS
    {

        public static ISolver structuralSolution;
        static int[] structuralBoundaryConditions;
        const double length = 1.0;
        const double gap = 0.01;
        const double thickness = 0.01;

        //External loads
        const double externalStructuralLoad = -3.0 * 1e5;

        static List<int> loadedStructuralDOFs;
        static double[] externalForcesStructuralVector;

        const double YoungMod = 1.0 * 1e9;
        const double YoungMod2 = 2.0 * 1e5;

        const double poissonRatio = 0.3;
        const double density = 8000.0;
        const double area = 1.0;
        const double contactArea = 0.005;



        private static void CreateStructuralBoundaryConditions()
        {
            List<int> boundedDofs = new List<int>();
            boundedDofs.Add(1);
            boundedDofs.Add(3);
            boundedDofs.Add(5);
            boundedDofs.Add(7);
            boundedDofs.Add(9);
            boundedDofs.Add(11);
            boundedDofs.Add(13);
            boundedDofs.Add(15);
            boundedDofs.Add(17);
            boundedDofs.Add(19);
            boundedDofs.Add(20);
            boundedDofs.Add(21);
            boundedDofs.Add(22);
            boundedDofs.Add(27);
            boundedDofs.Add(28);
            structuralBoundaryConditions = boundedDofs.ToArray<int>();
        }

        private static void CreateStructuralLoadVector()
        {
            loadedStructuralDOFs = new List<int>();
            loadedStructuralDOFs.Add(14);
            loadedStructuralDOFs.Add(18);
            externalForcesStructuralVector = new double[14 * 2];
        }

        private static Dictionary<int, INode> CreateNodes()
        {

            Dictionary<int, INode> nodes = new Dictionary<int, INode>();
            int k;
            k = 1;
            //First block
            nodes[k] = new Node(length, 0.0);
            k += 1;
            nodes[k] = new Node((3.0/2.0) * length, 0.0);
            k += 1;
            nodes[k] = new Node(2 * length, 0.0);
            k += 1;
            nodes[k] = new Node(length, (1.0/2.0) * length);
            k += 1;
            nodes[k] = new Node((3.0 / 2.0) * length, (1.0 / 2.0) * length);
            k += 1;
            nodes[k] = new Node(2 * length, (1.0 / 2.0) * length);
            k += 1;
            nodes[k] = new Node(length, length);
            k += 1;
            nodes[k] = new Node((3.0 / 2.0) * length, length);
            k += 1;
            nodes[k] = new Node(2.0 * length, length);
            k += 1;
            //Second block
            nodes[k] = new Node(0.5 * length, -2*length - gap);
            k += 1;
            nodes[k] = new Node(2.5 * length, -2*length - gap);
            k += 1;
            nodes[k] = new Node(0.5 * length, -gap);
            k += 1;
            nodes[k] = new Node(2.5 * length, -gap);
            k += 1;
            //Bars
            nodes[k] = new Node((3.0/2.0) * length, 2 * length);
            return nodes;
        }

        private static Dictionary<int, Dictionary<int, int>> CreateConnectivity()
        {

            Dictionary<int, Dictionary<int, int>> connectivity = new Dictionary<int, Dictionary<int, int>>();
            int k = 1;
            connectivity[k] = new Dictionary<int, int>() { { 1, 1 }, { 2, 2 }, { 3, 5 }, { 4, 4 } };
            k += 1;
            connectivity[k] = new Dictionary<int, int>() { { 1, 2 }, { 2, 3 }, { 3, 6 }, { 4, 5 } };
            k += 1;
            connectivity[k] = new Dictionary<int, int>() { { 1, 4 }, { 2, 5 }, { 3, 8 }, { 4, 7 } };
            k += 1;
            connectivity[k] = new Dictionary<int, int>() { { 1, 5 }, { 2, 6 }, { 3, 9 }, { 4, 8 } };
            k += 1;
            connectivity[k] = new Dictionary<int, int>() { { 1, 10 }, { 2, 11 }, { 3, 13 }, { 4, 12 } };
            k += 1;
            //bar elements
            connectivity[k] = new Dictionary<int, int>() { { 1, 8 }, { 2, 14 } };
            k += 1;
            //Contact elements
            connectivity[k] = new Dictionary<int, int>() { { 1, 12 }, { 2, 13 }, { 3, 1 } };
            k += 1;
            connectivity[k] = new Dictionary<int, int>() { { 1, 12 }, { 2, 13 }, { 3, 2 } };
            k += 1;
            connectivity[k] = new Dictionary<int, int>() { { 1, 12 }, { 2, 13 }, { 3, 3 } };
            return connectivity;
        }

        private static Dictionary<int, bool[]> CreateNodeFAT()
        {
            Dictionary<int, bool[]> nodeFAT = new Dictionary<int, bool[]>();
            for (int i = 1; i <= 14; i++)
            {
                nodeFAT[i] = new bool[] { true, true, false, false, false, false };
            }
            return nodeFAT;
        }
        private static Dictionary<int, IElementProperties> CreateElementProperties()
        {
            double E = YoungMod;
            double E2 = YoungMod2;

            double A = area;
            double CA = contactArea;

            string type = "Quad4";
            string type2 = "ContactNtS2D";
            string type3 = "Bar2D";

            Dictionary<int, IElementProperties> elementProperties = new Dictionary<int, IElementProperties>();
            for (int i = 1; i <= 5; i++)
            {
                elementProperties[i] = new ElementProperties(E, poissonRatio, A, thickness, density, type);

            }
            elementProperties[6] = new ElementProperties(E2, A, type3);
            elementProperties[6].Density = density;
            elementProperties[7] = new ElementProperties(E, CA, type2);
            elementProperties[7].Density = density;
            elementProperties[8] = new ElementProperties(E, CA, type2);
            elementProperties[8].Density = density;
            elementProperties[9] = new ElementProperties(E, CA, type2);
            elementProperties[9].Density = density;
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
            int countContactElements = elementsAssembly.CountElementsOfSameType(typeof(ContactNtN2D));
            ShowToGUI.PlotInitialGeometry(elementsAssembly);
            structuralSolution.LinearScheme = new LUFactorization();
            structuralSolution.NonLinearScheme.Tolerance = 1e-5;
            structuralSolution.ActivateNonLinearSolver = true;
            structuralSolution.NonLinearScheme.numberOfLoadSteps = 10;

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
            Dictionary<int, List<double[]>> gPointsStress = new Dictionary<int, List<double[]>>();
            Dictionary<int, List<double[]>> gPointsStrain = new Dictionary<int, List<double[]>>();
            Dictionary<int, List<double[]>> gPoints = new Dictionary<int, List<double[]>>();
            Dictionary<int, List<double[]>> nodalStress = new Dictionary<int, List<double[]>>();
            Dictionary<int, List<double[]>> nodalStrain = new Dictionary<int, List<double[]>>();
            for (int i = 1; i <= allStepsSolutions.Count; i++)
            {
                string name = "NodalCoordinates" + i.ToString() + ".dat";
                ExportToFile.ExportUpdatedNodalCoordinates(elementsAssembly, BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions.Single(m => m.Key == i).Value, elementsAssembly.BoundedDOFsVector), name);
                gPointsStress = elementsAssembly.GetElementsStresses(allStepsSolutions[i]);
                gPointsStrain = elementsAssembly.GetElementsStains(allStepsSolutions[i]);
                gPoints = elementsAssembly.GetElementsGaussPoints(allStepsSolutions[i]);
                nodalStress = elementsAssembly.GetElementsNodesStresses(allStepsSolutions[i]);
                nodalStrain = elementsAssembly.GetElementsNodesStains(allStepsSolutions[i]);
                string name1 = "GPointsStress" + i.ToString() + ".dat";
                string name2 = "GPointsStrain" + i.ToString() + ".dat";
                string name3 = "GPointsCoordinates" + i.ToString() + ".dat";
                string name4 = "NodalStress" + i.ToString() + ".dat";
                string name5 = "NodalStrain" + i.ToString() + ".dat";

                VectorOperations.PrintDictionaryofListsofVectorsToFile(gPointsStress, @"C:\Users\Public\Documents\" + name1);
                VectorOperations.PrintDictionaryofListsofVectorsToFile(gPointsStrain, @"C:\Users\Public\Documents\" + name2);
                VectorOperations.PrintDictionaryofListsofVectorsToFile(gPoints, @"C:\Users\Public\Documents\" + name3);
                VectorOperations.PrintDictionaryofListsofVectorsToFile(nodalStress, @"C:\Users\Public\Documents\" + name4);
                VectorOperations.PrintDictionaryofListsofVectorsToFile(nodalStrain, @"C:\Users\Public\Documents\" + name5);
            }
            elementsAssembly.UpdateDisplacements(solvector3);
            ShowToGUI.PlotFinalGeometry(elementsAssembly);
            double[] fullSolVector3 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(solvector3, elementsAssembly.BoundedDOFsVector);
            Dictionary<int, INode> finalNodes = Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullSolVector3);
            double[] xFinalNodalCoor = Assembly.NodalCoordinatesToVectors(finalNodes).Item1;
            double[] yFinalNodalCoor = Assembly.NodalCoordinatesToVectors(finalNodes).Item2;
            Dictionary<int, double[]> allStepsFullSolutions = new Dictionary<int, double[]>();
            Dictionary<int, Dictionary<int, double[]>> allStepsContactForces = new Dictionary<int, Dictionary<int, double[]>>();
            Dictionary<int, double[]> elementsInternalContactForcesVector;
            for (int i = 1; i <= allStepsSolutions.Count; i++)
            {
                elementsInternalContactForcesVector = new Dictionary<int, double[]>();
                elementsAssembly.UpdateDisplacements(allStepsSolutions[i]);
                elementsInternalContactForcesVector[7] = elementsAssembly.ElementsAssembly[7].CreateInternalGlobalForcesVector();
                elementsInternalContactForcesVector[8] = elementsAssembly.ElementsAssembly[8].CreateInternalGlobalForcesVector();
                elementsInternalContactForcesVector[9] = elementsAssembly.ElementsAssembly[9].CreateInternalGlobalForcesVector();
                allStepsContactForces[i] = elementsInternalContactForcesVector;
                string name = "ContactForces" + i.ToString() + ".dat";
                double[] Vector = new double[10];
                Vector[0] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == 7).Value[0] +
                    allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == 8).Value[0]+
                    allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == 9).Value[0];
                Vector[1] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == 7).Value[1] +
                    allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == 8).Value[1] +
                    allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == 9).Value[1];
                Vector[2] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == 7).Value[4];
                Vector[3] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == 7).Value[5];
                Vector[4] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == 8).Value[4];
                Vector[5] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == 8).Value[5];
                Vector[6] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == 9).Value[4];
                Vector[7] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == 9).Value[5];
                Vector[8] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == 7).Value[2] +
                    allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == 8).Value[2] +
                    allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == 9).Value[2];
                Vector[9] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == 7).Value[3] +
                    allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == 8).Value[3] +
                    allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == 9).Value[3];
                VectorOperations.PrintVectorToFile(Vector, @"C:\Users\Public\Documents\" + name);
            }

            for (int i = 0; i < allStepsSolutions.Count; i++)
            {
                allStepsFullSolutions.Add(i + 1, BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions.Single(m => m.Key == i + 1).Value, elementsAssembly.BoundedDOFsVector));
                int j = i + 1;
                string name = "solution" + j.ToString() + ".dat";
                VectorOperations.PrintVectorToFile(allStepsFullSolutions.Single(m => m.Key == i + 1).Value, @"C:\Users\Public\Documents\" + name);
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
