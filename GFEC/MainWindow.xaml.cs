using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using LiveCharts;
using LiveCharts.Wpf;
using Microsoft.Win32;

namespace GFEC
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        public SeriesCollection Graph { get; set; }
        public SeriesCollection Mesh { get; set; }
        private Results solverResults;
        private Dictionary<int, INode> nodes = new Dictionary<int, INode>();
        private Dictionary<int, Dictionary<int, int>> elementsConnectivity = new Dictionary<int, Dictionary<int, int>>();
        public string selectedExample;
        public SeriesCollection Something { get; set; }

        

        public MainWindow()
        {
            InitializeComponent();
            LoadComboBox();
            gnuplotImage.Source = null;


        }

        private void RunButton(object sender, RoutedEventArgs args)
        {
            //SolveSelectedExample();
            CoupledThermalStructural.diagramData = new ShowToGUI();
            CoupledThermalStructural.diagramData.TestEvent += TestEventMethod;
            CoupledThermalStructural.diagramData.TestEventMethod();
            selectedExample = ComboBox1.SelectedItem.ToString();
            Thread thread1 = new Thread(SolveSelectedExample);
            thread1.SetApartmentState(ApartmentState.STA);
            thread1.Start();
            
            
            //thread1.Join();
            //Graph = ShowToGUI.ShowResults(solverResults);


            //Dictionary<int, INode> nodes = new Dictionary<int, INode>();
            //nodes[1] = new Node(0.0, 0.01);
            //nodes[2] = new Node(0.3, 0.01);
            //nodes[3] = new Node(0.6, 0.01);
            //nodes[4] = new Node(0.6, 0.12);
            //nodes[5] = new Node(0.3, 0.12);
            //nodes[6] = new Node(0.0, 0.12);
            //nodes[7] = new Node(0.45, -0.11);
            //nodes[8] = new Node(0.75, -0.11);
            //nodes[9] = new Node(1.05, -0.11);
            //nodes[10] = new Node(0.45, 0.0);
            //nodes[11] = new Node(0.75, 0.0);
            //nodes[12] = new Node(1.05, 0.0);
            //Dictionary<int, Dictionary<int, int>> connectivity = new Dictionary<int, Dictionary<int, int>>();
            //connectivity[1] = new Dictionary<int, int>() { { 1, 1 }, { 2, 2 }, { 3, 5 }, { 4, 6 } };
            //connectivity[2] = new Dictionary<int, int>() { { 1, 2 }, { 2, 3 }, { 3, 4 }, { 4, 5 } };
            //connectivity[3] = new Dictionary<int, int>() { { 1, 7 }, { 2, 8 }, { 3, 11 }, { 4, 10 } };
            //connectivity[4] = new Dictionary<int, int>() { { 1, 8 }, { 2, 9 }, { 3, 12 }, { 4, 11 } };
            //connectivity[5] = new Dictionary<int, int>() { { 1, 10 }, { 2, 11 }, { 3, 3 } };
            //Mesh = ShowToGUI.DrawMesh(nodes, connectivity);

            DataContext = this;

            return;
        }

        private void LoadComboBox()
        {
            List<string> exampleList = new List<string>();
            exampleList.Add("LinearTrussExample");
            exampleList.Add("TwoQuadsExample");
            exampleList.Add("TwoBeamsInFrContactQuadsExample");
            exampleList.Add("ThermalExample");
            exampleList.Add("TwoThermalQuadsInContactExample");
            exampleList.Add("TwoThermalQuadsExample");
            exampleList.Add("TwoQuadsInContactNewExample");
            exampleList.Add("CoupledPhysicsExample");
            exampleList.Add("CoupledThermalStructural");
            exampleList.Add("CoupledThermalStructuralCNTs");
            exampleList.Add("CoupledThermalStructuralCNTs1b");
            exampleList.Add("CoupledThermalStructuralCNTs2");
            exampleList.Add("CoupledThermalStructuralCNTsInAngle");
            exampleList.Add("CoupledThermalStructuralCNTsInAngle2");
            exampleList.Add("CoupledThermalStructuralCNTsInAngle3");
            exampleList.Add("CoupledThermalStructuralCNTsInAngle4");
            exampleList.Add("CoupledThermalStructuralCNTsInAngle5");
            exampleList.Add("CoupledThermalStructural2");
            exampleList.Add("CNTExample");
            exampleList.Add("CNTsInParallelFinalExample");
            exampleList.Add("CNTsInAngleFinalExample");
            exampleList.Add("ElasticContrainsCNTsInAngleFinal");
            exampleList.Add("NewExampleContacts");
            exampleList.Add("DynamicExample");
            exampleList.Add("ImpactElasticAgainstRigid");
            exampleList.Add("NewDynamicExample");
            exampleList.Add("ImpactBetweenBars");
            exampleList.Add("ImpactElasticAgainstRigid2");
            exampleList.Add("ImpactCircle");
            exampleList.Add("ImpactCircle2");
            exampleList.Add("CantileverWithTriangElements");
            exampleList.Add("CantileverWithQuad8Elements");
            exampleList.Add("TwoBlocksInContact3D");
            exampleList.Add("Hxa8TestExample");
            exampleList.Add("ThreeTrusses");
            exampleList.Add("TwoBlocks2DNtN");
            exampleList.Add("TwoBlocks2DNtS");
            exampleList.Add("BendingBeamContact2d");
            exampleList.Add("BendingOveraRigidCylinder");
            exampleList.Add("TwoBlocksHigherOrderNTS");
            exampleList.Add("BendingBeamContact3d");
            exampleList.Add("BendingBeamContact3dWithFriction");
            exampleList.Add("BendingBeamContact3dWithFrictionRefinedMesh");
            exampleList.Add("BeamsInAngleContact3dWithFriction");
            exampleList.Add("Blocks3dContactSliding");
            exampleList.Add("BendingBeamContact3dWithFrictionRefinedMesh2");
            exampleList.Add("Blocks3dContactSlidingMeshRefined");
            exampleList.Add("BendingBeamContact3dWithFrictionQuadraticShapeFunctions");
            exampleList.Add("Blocks3dContactSlidingQuadratic");
            exampleList.Add("shell2DExample");
            exampleList.Add("Impactshell2DExample");
            exampleList.Add("Impact3dSolids");
            exampleList.Add("DegenerateShellElementsLinearExample");
            exampleList.Add("DegenerateShellElementsContactQSExample");
            exampleList.Add("Cantilever3dCheck");
            exampleList.Add("DegenerateShellElementsImpactExample");
            exampleList.Add("CNTsInParallelFinalSensitivityAnalysis");
            exampleList.Add("CNTsInParallelFinalSensitivityAnalysis2");
            exampleList.Add("SolidShellLinearExample");
            exampleList.Add("SolidShellElementsContactExample");
            ComboBox1.ItemsSource = exampleList;
        }

        private void SolveSelectedExample()
        {
            Results finalResults;
            Tuple<Dictionary<int, double[]>, Dictionary<int, double>> results;
            //string selectedExample = ComboBox1.SelectedItem.ToString();
            switch (selectedExample)
            {
                case "TwoQuadsExample":
                    finalResults = TwoQuadsExample.RunStaticExample();
                    break;
                case "LinearTrussExample":
                    finalResults = LinearTrussExample.RunExample();
                    break;
                case "TwoBeamsInFrContactQuadsExample":
                    finalResults = TwoBeamsInFrContactQuadsExample.RunDynamicExample();
                    break;
                case "ThermalExample":
                    finalResults = ThermalExample.RunStaticExample();
                    break;
                case "TwoThermalQuadsInContactExample":
                    finalResults = TwoThermalQuadsInContactExample.RunStaticExample();
                    break;
                case "TwoThermalQuadsExample":
                    finalResults = TwoThermalQuadsExample.RunStaticExample();
                    break;
                case "TwoQuadsInContactNewExample":
                    finalResults = TwoQuadsInContactNewExample.RunStaticExample();
                    break;
                case "CoupledPhysicsExample":
                    finalResults = CoupledPhysicsExample.RunStaticExample();
                    break;
                case "CoupledThermalStructural":
                    CoupledThermalStructural.structuralSolution = new StaticSolver();
                    CoupledThermalStructural.structuralSolution.NonLinearScheme = new LoadControlledNewtonRaphson();
                    CoupledThermalStructural.structuralSolution.NonLinearScheme.convergenceResult += NonLinearScheme_convergenceResult;
                    CoupledThermalStructural.diagramData = new ShowToGUI();
                    CoupledThermalStructural.diagramData.ShowDiagramInGUI += c_ShowDiagramInGUI;
                    finalResults = CoupledThermalStructural.RunStaticExample();
                    break;
                case "CoupledThermalStructural2":
                    CoupledThermalStructural2.structuralSolution = new StaticSolver();
                    CoupledThermalStructural2.structuralSolution.NonLinearScheme = new LoadControlledNewtonRaphson();
                    //CoupledThermalStructural2.structuralSolution.NonLinearScheme.convergenceResult += NonLinearScheme_convergenceResult;
                    CoupledThermalStructural2.diagramData = new ShowToGUI();
                    CoupledThermalStructural2.diagramData.ShowDiagramInGUI += c_ShowDiagramInGUI;
                    finalResults = CoupledThermalStructural2.RunStaticExample();
                    break;
                case "CoupledThermalStructuralCNTs":
                    CoupledThermalStructuralCNTs.structuralSolution = new StaticSolver(); 
                    CoupledThermalStructuralCNTs.structuralSolution.NonLinearScheme = new LoadControlledNewtonRaphson();
                    CoupledThermalStructuralCNTs.structuralSolution.NonLinearScheme.convergenceResult += NonLinearScheme_convergenceResult;
                    finalResults = CoupledThermalStructuralCNTs.RunStaticExample();
                    break;
                case "CoupledThermalStructuralCNTs1b":
                    CoupledThermalStructuralCNTs1b.structuralSolution = new StaticSolver();
                    CoupledThermalStructuralCNTs1b.structuralSolution.NonLinearScheme = new LoadControlledNewtonRaphson();
                    CoupledThermalStructuralCNTs1b.structuralSolution.NonLinearScheme.convergenceResult += NonLinearScheme_convergenceResult;
                    finalResults = CoupledThermalStructuralCNTs1b.RunStaticExample();
                    break;
                case "CoupledThermalStructuralCNTs2":
                    CoupledThermalStructuralCNTs2.structuralSolution = new StaticSolver();
                    CoupledThermalStructuralCNTs2.structuralSolution.NonLinearScheme = new LoadControlledNewtonRaphson();
                    CoupledThermalStructuralCNTs2.structuralSolution.NonLinearScheme.convergenceResult += NonLinearScheme_convergenceResult;
                    finalResults = CoupledThermalStructuralCNTs2.RunStaticExample();
                    break;
                case "CoupledThermalStructuralCNTsInAngle":
                    CoupledThermalStructuralCNTsInAngle.structuralSolution = new StaticSolver();
                    CoupledThermalStructuralCNTsInAngle.structuralSolution.NonLinearScheme = new LoadControlledNewtonRaphson();
                    CoupledThermalStructuralCNTsInAngle.structuralSolution.NonLinearScheme.convergenceResult += NonLinearScheme_convergenceResult;
                    finalResults = CoupledThermalStructuralCNTsInAngle.RunStaticExample();
                    break;
                case "CoupledThermalStructuralCNTsInAngle2":
                    CoupledThermalStructuralCNTsInAngle2.structuralSolution = new StaticSolver();
                    CoupledThermalStructuralCNTsInAngle2.structuralSolution.NonLinearScheme = new LoadControlledNewtonRaphson();
                    CoupledThermalStructuralCNTsInAngle2.structuralSolution.NonLinearScheme.convergenceResult += NonLinearScheme_convergenceResult;
                    CoupledThermalStructuralCNTsInAngle2.thermalSolution = new StaticSolver();
                    CoupledThermalStructuralCNTsInAngle2.thermalSolution.NonLinearScheme = new LoadControlledNewtonRaphson();
                    CoupledThermalStructuralCNTsInAngle2.thermalSolution.NonLinearScheme.convergenceResult += NonLinearScheme_convergenceResult;
                    finalResults = CoupledThermalStructuralCNTsInAngle2.RunStaticExample();
                    break;
                case "CoupledThermalStructuralCNTsInAngle3":
                    CoupledThermalStructuralCNTsInAngle3.structuralSolution = new StaticSolver();
                    CoupledThermalStructuralCNTsInAngle3.structuralSolution.NonLinearScheme = new LoadControlledNewtonRaphson();
                    CoupledThermalStructuralCNTsInAngle3.structuralSolution.NonLinearScheme.convergenceResult += NonLinearScheme_convergenceResult;
                    //finalResults = CoupledThermalStructuralCNTsInAngle3.RunStaticExample();
                    //CoupledThermalStructuralCNTsInAngle3.thermalSolution = new StaticSolver();
                    //CoupledThermalStructuralCNTsInAngle3.thermalSolution.NonLinearScheme = new LoadControlledNewtonRaphson();
                    //CoupledThermalStructuralCNTsInAngle3.thermalSolution.NonLinearScheme.convergenceResult += NonLinearScheme_convergenceResult;
                    finalResults = CoupledThermalStructuralCNTsInAngle3.RunStaticExample();
                    break;
                case "CoupledThermalStructuralCNTsInAngle4":
                    CoupledThermalStructuralCNTsInAngle4.structuralSolution = new StaticSolver();
                    CoupledThermalStructuralCNTsInAngle4.structuralSolution.NonLinearScheme = new LoadControlledNewtonRaphson();
                    CoupledThermalStructuralCNTsInAngle4.structuralSolution.NonLinearScheme.convergenceResult += NonLinearScheme_convergenceResult;
                    finalResults = CoupledThermalStructuralCNTsInAngle4.RunStaticExample();
                    break;
                case "CoupledThermalStructuralCNTsInAngle5":
                    CoupledThermalStructuralCNTsInAngle5.structuralSolution = new StaticSolver();
                    CoupledThermalStructuralCNTsInAngle5.structuralSolution.NonLinearScheme = new LoadControlledNewtonRaphson();
                    CoupledThermalStructuralCNTsInAngle5.structuralSolution.NonLinearScheme.convergenceResult += NonLinearScheme_convergenceResult;
                    CoupledThermalStructuralCNTsInAngle5.thermalSolution = new StaticSolver();
                    CoupledThermalStructuralCNTsInAngle5.thermalSolution.NonLinearScheme = new LoadControlledNewtonRaphson();
                    CoupledThermalStructuralCNTsInAngle5.thermalSolution.NonLinearScheme.convergenceResult += NonLinearScheme_convergenceResult;
                    finalResults = CoupledThermalStructuralCNTsInAngle5.RunStaticExample();
                    break;
                case "CNTExample":
                    CNTExample.structuralSolution = new StaticSolver();
                    CNTExample.structuralSolution.NonLinearScheme = new LoadControlledNewtonRaphson();
                    CNTExample.structuralSolution.NonLinearScheme.convergenceResult += NonLinearScheme_convergenceResult;
                    finalResults = CNTExample.RunStaticExample();
                    break;
                case "CNTsInParallelFinalExample":
                    CNTsInParallelFinalExample.structuralSolution = new StaticSolver();
                    CNTsInParallelFinalExample.structuralSolution.NonLinearScheme = new LoadControlledNewtonRaphson();
                    CNTsInParallelFinalExample.structuralSolution.NonLinearScheme.convergenceResult += NonLinearScheme_convergenceResult;
                    finalResults = CNTsInParallelFinalExample.RunStaticExample();
                    break;
                case "CNTsInAngleFinalExample":
                    CNTsInAngleFinalExample.structuralSolution = new StaticSolver();
                    CNTsInAngleFinalExample.structuralSolution.NonLinearScheme = new LoadControlledNewtonRaphson();
                    CNTsInAngleFinalExample.structuralSolution.NonLinearScheme.convergenceResult += NonLinearScheme_convergenceResult;
                    finalResults = CNTsInAngleFinalExample.RunStaticExample();
                    break;
                case "ElasticContrainsCNTsInAngleFinal":
                    ElasticContrainsCNTsInAngleFinal.structuralSolution = new StaticSolver();
                    ElasticContrainsCNTsInAngleFinal.structuralSolution.NonLinearScheme = new LoadControlledNewtonRaphson();
                    ElasticContrainsCNTsInAngleFinal.structuralSolution.NonLinearScheme.convergenceResult += NonLinearScheme_convergenceResult;
                    finalResults = ElasticContrainsCNTsInAngleFinal.RunStaticExample();
                    break;
                case "NewExampleContacts":
                    NewExampleContacts.structuralSolution = new StaticSolver();
                    NewExampleContacts.structuralSolution.NonLinearScheme = new LoadControlledNewtonRaphson();
                    NewExampleContacts.structuralSolution.NonLinearScheme.convergenceResult += NonLinearScheme_convergenceResult;
                    finalResults = NewExampleContacts.RunStaticExample();
                    break;
                case "DynamicExample":
                    finalResults = DynamicExample.RunExample();
                    break;
                case "ImpactElasticAgainstRigid":
                    finalResults = ImpactElasticAgainstRigid.RunExample();
                    break;
                case "NewDynamicExample":
                    finalResults = NewDynamicExample.RunExample();
                    break;
                case "ImpactBetweenBars":
                    finalResults = ImpactBetweenBars.RunExample();
                    break;
                case "ImpactElasticAgainstRigid2":
                    finalResults = ImpactElasticAgainstRigid2.RunExample();
                    break;
                case "ImpactCircle":
                    finalResults = ImpactCircle.RunExample();
                    break;
                case "CantileverWithTriangElements":
                    finalResults = CantileverWithTriangElements.RunExample();
                    break;
                case "CantileverWithQuad8Elements":
                    finalResults = CantileverWithQuad8Elements.RunExample();
                    break;
                case "ImpactCircle2":
                    finalResults = ImpactCircle2.RunExample();
                    break;
                case "TwoBlocksInContact3D":
                    TwoBlocksInContact3D.newSolu = new StaticSolver();
                    TwoBlocksInContact3D.newSolu.NonLinearScheme = new LoadControlledNewtonRaphson();
                    TwoBlocksInContact3D.newSolu.NonLinearScheme.convergenceResult += NonLinearScheme_convergenceResult;
                    finalResults = TwoBlocksInContact3D.RunStaticExample();
                    break;
                case "Hxa8TestExample":
                    Hxa8TestExample.newSolu = new StaticSolver();
                    Hxa8TestExample.newSolu.NonLinearScheme = new LoadControlledNewtonRaphson();
                    Hxa8TestExample.newSolu.NonLinearScheme.convergenceResult += NonLinearScheme_convergenceResult;
                    finalResults = Hxa8TestExample.RunStaticExample();
                    break;
                case "ThreeTrusses":
                    ThreeTrusses.structuralSolution = new StaticSolver();
                    ThreeTrusses.structuralSolution.NonLinearScheme = new LoadControlledNewtonRaphson();
                    ThreeTrusses.structuralSolution.NonLinearScheme.convergenceResult += NonLinearScheme_convergenceResult;
                    finalResults = ThreeTrusses.RunStaticExample();
                    break;
                case "TwoBlocks2DNtN":
                    TwoBlocks2DNtN.structuralSolution = new StaticSolver();
                    TwoBlocks2DNtN.structuralSolution.NonLinearScheme = new LoadControlledNewtonRaphson();
                    TwoBlocks2DNtN.structuralSolution.NonLinearScheme.convergenceResult += NonLinearScheme_convergenceResult;
                    finalResults = TwoBlocks2DNtN.RunStaticExample();
                    break;
                case "TwoBlocks2DNtS":
                    TwoBlocks2DNtS.structuralSolution = new StaticSolver();
                    TwoBlocks2DNtS.structuralSolution.NonLinearScheme = new LoadControlledNewtonRaphson();
                    TwoBlocks2DNtS.structuralSolution.NonLinearScheme.convergenceResult += NonLinearScheme_convergenceResult;
                    finalResults = TwoBlocks2DNtS.RunStaticExample();
                    break;
                case "BendingBeamContact2d":
                    BendingBeamContact2d.structuralSolution = new StaticSolver();
                    BendingBeamContact2d.structuralSolution.NonLinearScheme = new LoadControlledNewtonRaphson();
                    BendingBeamContact2d.structuralSolution.NonLinearScheme.convergenceResult += NonLinearScheme_convergenceResult;
                    finalResults = BendingBeamContact2d.RunStaticExample();
                    break;
                case "BendingOveraRigidCylinder":
                    BendingOveraRigidCylinder.structuralSolution = new StaticSolver();
                    BendingOveraRigidCylinder.structuralSolution.NonLinearScheme = new LoadControlledNewtonRaphson();
                    BendingOveraRigidCylinder.structuralSolution.NonLinearScheme.convergenceResult += NonLinearScheme_convergenceResult;
                    finalResults = BendingOveraRigidCylinder.RunStaticExample();
                    break;
                case "TwoBlocksHigherOrderNTS":
                    TwoBlocksHigherOrderNTS.structuralSolution = new StaticSolver();
                    TwoBlocksHigherOrderNTS.structuralSolution.NonLinearScheme = new LoadControlledNewtonRaphson();
                    TwoBlocksHigherOrderNTS.structuralSolution.NonLinearScheme.convergenceResult += NonLinearScheme_convergenceResult;
                    finalResults = TwoBlocksHigherOrderNTS.RunStaticExample();
                    break;
                case "BendingBeamContact3d":
                    BendingBeamContact3d.structuralSolution = new StaticSolver();
                    BendingBeamContact3d.structuralSolution.NonLinearScheme = new LoadControlledNewtonRaphson();
                    BendingBeamContact3d.structuralSolution.NonLinearScheme.convergenceResult += NonLinearScheme_convergenceResult;
                    finalResults = BendingBeamContact3d.RunStaticExample();
                    break;
                case "BendingBeamContact3dWithFriction":
                    BendingBeamContact3dWithFriction.structuralSolution = new StaticSolver();
                    BendingBeamContact3dWithFriction.structuralSolution.NonLinearScheme = new LoadControlledNewtonRaphson();
                    BendingBeamContact3dWithFriction.structuralSolution.NonLinearScheme.convergenceResult += NonLinearScheme_convergenceResult;
                    finalResults = BendingBeamContact3dWithFriction.RunStaticExample();
                    break;
                case "BendingBeamContact3dWithFrictionRefinedMesh":
                    BendingBeamContact3dWithFrictionRefinedMesh.structuralSolution = new StaticSolver();
                    BendingBeamContact3dWithFrictionRefinedMesh.structuralSolution.NonLinearScheme = new LoadControlledNewtonRaphson();
                    BendingBeamContact3dWithFrictionRefinedMesh.structuralSolution.NonLinearScheme.convergenceResult += NonLinearScheme_convergenceResult;
                    finalResults = BendingBeamContact3dWithFrictionRefinedMesh.RunStaticExample();
                    break;
                case "BeamsInAngleContact3dWithFriction":
                    BeamsInAngleContact3dWithFriction.structuralSolution = new StaticSolver();
                    //BeamsInAngleContact3dWithFriction.structuralSolution.NonLinearScheme = new ModifiedNewtonRaphson();
                    BeamsInAngleContact3dWithFriction.structuralSolution.NonLinearScheme = new LoadControlledNewtonRaphson();
                    BeamsInAngleContact3dWithFriction.structuralSolution.NonLinearScheme.convergenceResult += NonLinearScheme_convergenceResult;
                    finalResults = BeamsInAngleContact3dWithFriction.RunStaticExample();
                    break;
                case "Blocks3dContactSliding":
                    Blocks3dContactSliding.structuralSolution = new StaticSolver();
                    Blocks3dContactSliding.structuralSolution.NonLinearScheme = new LoadControlledNewtonRaphson();
                    Blocks3dContactSliding.structuralSolution.NonLinearScheme.convergenceResult += NonLinearScheme_convergenceResult;
                    finalResults = Blocks3dContactSliding.RunStaticExample();
                    break;
                case "BendingBeamContact3dWithFrictionRefinedMesh2":
                    BendingBeamContact3dWithFrictionRefinedMesh2.structuralSolution = new StaticSolver();
                    BendingBeamContact3dWithFrictionRefinedMesh2.structuralSolution.NonLinearScheme = new MMCPCGLoadControlledNewtonRaphson();
                    BendingBeamContact3dWithFrictionRefinedMesh2.structuralSolution.NonLinearScheme.convergenceResult += NonLinearScheme_convergenceResult;
                    finalResults = BendingBeamContact3dWithFrictionRefinedMesh2.RunStaticExample();
                    break;
                case "Blocks3dContactSlidingMeshRefined":
                    Blocks3dContactSlidingMeshRefined.structuralSolution = new StaticSolver();
                    Blocks3dContactSlidingMeshRefined.structuralSolution.NonLinearScheme = new LoadControlledNewtonRaphson();
                    Blocks3dContactSlidingMeshRefined.structuralSolution.NonLinearScheme.convergenceResult += NonLinearScheme_convergenceResult;
                    finalResults = Blocks3dContactSlidingMeshRefined.RunStaticExample();
                    break;
                case "BendingBeamContact3dWithFrictionQuadraticShapeFunctions":
                    BendingBeamContact3dWithFrictionQuadraticShapeFunctions.structuralSolution = new StaticSolver();
                    BendingBeamContact3dWithFrictionQuadraticShapeFunctions.structuralSolution.NonLinearScheme = new ModifiedNewtonRaphson();
                    BendingBeamContact3dWithFrictionQuadraticShapeFunctions.structuralSolution.NonLinearScheme.convergenceResult += NonLinearScheme_convergenceResult;
                    finalResults = BendingBeamContact3dWithFrictionQuadraticShapeFunctions.RunStaticExample();
                    break;
                case "Blocks3dContactSlidingQuadratic":
                    Blocks3dContactSlidingQuadratic.structuralSolution = new StaticSolver();
                    Blocks3dContactSlidingQuadratic.structuralSolution.NonLinearScheme = new LoadControlledNewtonRaphson();
                    Blocks3dContactSlidingQuadratic.structuralSolution.NonLinearScheme.convergenceResult += NonLinearScheme_convergenceResult;
                    finalResults = Blocks3dContactSlidingQuadratic.RunStaticExample();
                    break;
                case "shell2DExample":
                    shell2DExample.structuralSolution = new StaticSolver();
                    shell2DExample.structuralSolution.NonLinearScheme = new LoadControlledNewtonRaphson();
                    shell2DExample.structuralSolution.NonLinearScheme.convergenceResult += NonLinearScheme_convergenceResult;
                    finalResults = shell2DExample.RunStaticExample();
                    break;
                case "Impactshell2DExample":
                    finalResults = Impactshell2DExample.RunDynamicExample();
                    break;
                case "Impact3dSolids":
                    finalResults = Impact3dSolids.RunDynamicExample();
                    break;
                case "DegenerateShellElementsLinearExample":
                    DegenerateShellElementsLinearExample.structuralSolution = new StaticSolver();
                    DegenerateShellElementsLinearExample.structuralSolution.NonLinearScheme = new LoadControlledNewtonRaphson();
                    DegenerateShellElementsLinearExample.structuralSolution.NonLinearScheme.convergenceResult += NonLinearScheme_convergenceResult;
                    finalResults = DegenerateShellElementsLinearExample.RunStaticExample();
                    break;
                case "DegenerateShellElementsContactQSExample":
                    DegenerateShellElementsContactQSExample.structuralSolution = new StaticSolver();
                    DegenerateShellElementsContactQSExample.structuralSolution.NonLinearScheme = new LoadControlledNewtonRaphson();
                    DegenerateShellElementsContactQSExample.structuralSolution.NonLinearScheme.convergenceResult += NonLinearScheme_convergenceResult;
                    finalResults = DegenerateShellElementsContactQSExample.RunStaticExample();
                    break;
                case "Cantilever3dCheck":
                    Cantilever3dCheck.structuralSolution = new StaticSolver();
                    //Cantilever3dCheck.structuralSolution.NonLinearScheme = new LoadControlledNewtonRaphson();
                    //DegenerateShellElementsContactQSExample.structuralSolution.NonLinearScheme.convergenceResult += NonLinearScheme_convergenceResult;
                    finalResults = Cantilever3dCheck.RunStaticExample();
                    break;
                case "DegenerateShellElementsImpactExample":
                    finalResults = DegenerateShellElementsImpactExample.RunDynamicExample();
                    break;
                case "CNTsInParallelFinalSensitivityAnalysis":
                    CNTsInParallelFinalSensitivityAnalysis.structuralSolution = new StaticSolver();
                    CNTsInParallelFinalSensitivityAnalysis.structuralSolution.NonLinearScheme = new LoadControlledNewtonRaphson();
                    CNTsInParallelFinalSensitivityAnalysis.structuralSolution.NonLinearScheme.convergenceResult += NonLinearScheme_convergenceResult;
                    finalResults = CNTsInParallelFinalSensitivityAnalysis.RunStaticExample();
                    break;
                case "CNTsInParallelFinalSensitivityAnalysis2":
                    CNTsInParallelFinalSensitivityAnalysis2.structuralSolution = new StaticSolver();
                    CNTsInParallelFinalSensitivityAnalysis2.structuralSolution.NonLinearScheme = new LoadControlledNewtonRaphson();
                    CNTsInParallelFinalSensitivityAnalysis2.structuralSolution.NonLinearScheme.convergenceResult += NonLinearScheme_convergenceResult;
                    finalResults = CNTsInParallelFinalSensitivityAnalysis2.RunStaticExample();
                    break;
                case "SolidShellLinearExample":
                    SolidShellLinearExample.structuralSolution = new StaticSolver();
                    SolidShellLinearExample.structuralSolution.NonLinearScheme = new LoadControlledNewtonRaphson();
                    SolidShellLinearExample.structuralSolution.NonLinearScheme.convergenceResult += NonLinearScheme_convergenceResult;
                    finalResults = SolidShellLinearExample.RunStaticExample();
                    break;
                case "SolidShellElementsContactExample":
                    SolidShellElementsContactExample.structuralSolution = new StaticSolver();
                    SolidShellElementsContactExample.structuralSolution.NonLinearScheme = new LoadControlledNewtonRaphson();
                    SolidShellElementsContactExample.structuralSolution.NonLinearScheme.convergenceResult += NonLinearScheme_convergenceResult;
                    finalResults = SolidShellElementsContactExample.RunStaticExample();
                    break;
                default:
                    finalResults = TwoQuadsExample.RunStaticExample();
                    break;
            }
            //Results.Text = solution[0].ToString();

            solverResults = finalResults;
        }

        public void NonLinearScheme_convergenceResult(object sender, string e)
        {
            Application.Current.Dispatcher.Invoke(new Action(() =>
            {
                LogTool.Text = e;
            }));
        }

        private void PrintResultOnUI(object sender, string e)
        {
            Application.Current.Dispatcher.Invoke(new Action(() =>
            {
                LogTool.Text = e;
            }));
        }

        private void TestEventMethod(object sender, SeriesCollection e)
        {
            Application.Current.Dispatcher.Invoke(new Action(() =>
            {
                Something = e;
            }));
        }

        private void c_ShowDiagramInGUI(object sender, ShowDiagramInGUIArgs e)
        {
            Application.Current.Dispatcher.Invoke(new Action(() =>
            {
                Graph = e.DiagramData;
                Results.Text = "Test succeded";
            }));
            //Graph = e.DiagramData;
            //Results.Text = "Test succeded";
            //Environment.Exit(0);
        }

        private async void Import_Nodes_Button_Click(object sender, RoutedEventArgs args)
        {
            try
            {
                OpenFileDialog dialog1 = new OpenFileDialog();
                if (dialog1.ShowDialog() == true)
                {
                    StreamReader stream = new StreamReader(dialog1.FileName);
                    string file = await stream.ReadToEndAsync();
                    List<string> lines = new List<string>(file.Split(new string[] { "\r\n", "\r", "\n" }, StringSplitOptions.RemoveEmptyEntries));
                    lines.RemoveAt(0);

                    foreach (var line in lines)
                    {
                        // in case of first line ...
                        string[] fields = line.Split(new string[] { "\t" }, StringSplitOptions.None);
                        int nodeIndex = int.Parse(fields[0]);
                        var node = new Node(double.Parse(fields[1]), double.Parse(fields[2]));
                        nodes[nodeIndex] = node;
                    }
                }


            }
            catch (Exception ex)
            {
                throw ex;
            }
        }

        private async void Import_Connectivity_Button_Click(object sender, RoutedEventArgs e)
        {
            try
            {
                OpenFileDialog dialog1 = new OpenFileDialog();
                if (dialog1.ShowDialog() == true)
                {
                    StreamReader stream = new StreamReader(dialog1.FileName);
                    string file = await stream.ReadToEndAsync();
                    List<string> lines = new List<string>(file.Split(new string[] { "\r\n", "\r", "\n" }, StringSplitOptions.RemoveEmptyEntries));
                    lines.RemoveAt(0);

                    foreach (var line in lines)
                    {
                        // in case of first line ...
                        string[] fields = line.Split(new string[] { "\t" }, StringSplitOptions.None);
                        int elementIndex = int.Parse(fields[0]);

                        var element = new Dictionary<int, int>() {
                            { 1, int.Parse(fields[2]) },
                            { 2, int.Parse(fields[3]) },
                            { 3, int.Parse(fields[4]) },
                            { 4, int.Parse(fields[5]) }
                        };
                        elementsConnectivity[elementIndex] = element;
                    }
                }


            }
            catch (Exception ex)
            {
                throw ex;
            }
        }

        private void Button_Click(object sender, RoutedEventArgs e)
        {
            using (Game game = new Game(800, 600, "LearnOpenTK"))
            {
                game.Run(60.0);
            }



        }

        private void Button_Click_Gnuplot(object sender, RoutedEventArgs e)
        {
            GnuPlot.Set("terminal png size 500, 300");
            GnuPlot.Set("output 'gnuplot.png'");
            //GnuPlot.Plot("sin(x) + 2");
            //GnuPlot.Close();
            //gnuplotImage.Source = null;
            //gnuplotImage.Source = new BitmapImage(new Uri("pack://siteoforigin:,,/gnuplot.png"));



            //GnuPlot.Set("terminal png size 400, 300");
            //GnuPlot.Set("output 'gnuplot.png'");
            double[] X = new double[] { -15, -15, -15, -15, -15, -14, -14, -14, -14 };
            double[] Y = new double[] { 11, 12, 13, 14, 15, -15, -14, -13, -12, -11 };
            double[] Z = new double[] { 20394, 15745, 11885, 8771, 6330, 8771, 12155, 16469, 21818 };
            //double[,] Y = new double[,]
            //{
            //    { 0,0,0,1,2,2,1,0,0,0},
            //            { 0,0,2,3,3,3,3,2,0,0},
            //            { 0,2,3,4,4,4,4,3,2,0},
            //            { 2,3,4,5,5,5,5,4,3,2},
            //            { 3,4,5,6,7,7,6,5,4,3},
            //            { 3,4,5,6,7,7,6,5,4,3},
            //            { 2,3,4,5,5,5,5,4,3,2},
            //            { 0,2,3,4,4,4,4,3,2,0},
            //            { 0,0,2,3,3,3,3,2,0,0},
            //            { 0,0,0,1,2,2,1,0,0,0}
            //};
            //GnuPlot.Set("dgrid3d 50,50,2");
            //GnuPlot.Set("7,7,7");
            GnuPlot.Set("pm3d");
            GnuPlot.Set("dgrid3d");
            //GnuPlot.Set("contour");
            //GnuPlot.Set("map");
            //GnuPlot.Set("dgrid3d");
            //GnuPlot.Set("cntrparam levels 20", "isosamples 100");
            GnuPlot.Set("view map");
            //GnuPlot.Set("pm3d interpolate 10,10");

            GnuPlot.SPlot(X, Y, Z);
            GnuPlot.Set("output");

            GnuPlot.Close();

            while (true)
            {
                if (File.Exists(AppContext.BaseDirectory + "gnuplot.png") && new FileInfo(AppContext.BaseDirectory + "gnuplot.png").Length > 0)
                {
                    break;

                }
                Thread.Sleep(100);
            }
            GnuPlot.KillProcess();

            gnuplotImage.Source = null;
            gnuplotImage.Source = new BitmapImage(new Uri("file://" + AppContext.BaseDirectory + "gnuplot.png"));
        }

        private void Button_Test(object sender, RoutedEventArgs e)
        {
            double[] testVector = new double[75];
            for (int i = 0; i < 5; i++)
            {
                for (int j = 0; j < 15; j++)
                {
                    testVector[i * 15 + j] = i;
                }
            }
            string path = @"C:\Users\Public\Documents\ContourDataY.dat";
            ExportToFile.CreateContourDataForMatlab(testVector, testVector, testVector, 5, 15, path);
        }

        private void Button_ParallelTest(object sender, RoutedEventArgs e)
        {
            MultiThreadingExample test1 = new MultiThreadingExample();
            test1.timeElapsed += PrintResultOnUI;
            test1.RunExample();
        }
        //--------------------------------------------------------------------------
        //private double _axisMax;
        //private double _axisMin;


        //public Func<double, string> DateTimeFormatter { get; set; }
        //public double AxisStep { get; set; }
        //public double AxisUnit { get; set; }
        //public bool IsReading { get; set; }
        //public event PropertyChangedEventHandler PropertyChanged;
        //public void ConvergenceResults()
        //{
        //    var mapper = Mappers.Xy<ConvergenceValues>()
        //        .X(model => model.Iteration)   //use DateTime.Ticks as X
        //        .Y(model => model.ResidualNorm);           //use the value property as Y

        //    //lets save the mapper globally.
        //    Charting.For<ConvergenceValues>(mapper);

        //    //the values property will store our values array
        //    ChartValues = new ChartValues<ConvergenceValues>();

        //    //lets set how to display the X Labels
        //    //DateTimeFormatter = value => new DateTime((long)value).ToString("mm:ss");

        //    //AxisStep forces the distance between each separator in the X axis
        //    AxisStep = 1;// TimeSpan.FromSeconds(1).Ticks;
        //    //AxisUnit forces lets the axis know that we are plotting seconds
        //    //this is not always necessary, but it can prevent wrong labeling
        //    AxisUnit = 100;// TimeSpan.TicksPerSecond;
        //}

        //public double AxisMax
        //{
        //    get { return _axisMax; }
        //    set
        //    {
        //        _axisMax = value;
        //        OnPropertyChanged("AxisMax");
        //    }
        //}
        //public double AxisMin
        //{
        //    get { return _axisMin; }
        //    set
        //    {
        //        _axisMin = value;
        //        OnPropertyChanged("AxisMin");
        //    }
        //}

        //private void SetAxisLimits(int now)
        //{
        //    AxisMax = now + 1;
        //    AxisMin = now - 8;
        //}

        //protected virtual void OnPropertyChanged(string propertyName = null)
        //{
        //    if (PropertyChanged != null)
        //        PropertyChanged.Invoke(this, new PropertyChangedEventArgs(propertyName));
        //}

        //private void ImportOBJFile(object sender, RoutedEventArgs e)
        //{
        //    try
        //    {
        //        OpenFileDialog dialog1 = new OpenFileDialog();
        //        if (dialog1.ShowDialog() == true)
        //        {
        //            string selectedFilePath = dialog1.FileName;
        //            List<string> allLines = new List<string>(File.ReadAllLines(selectedFilePath));

        //            //List<string> lines = new List<string>(file.Split(new string[] { "\r\n", "\r", "\n" }, StringSplitOptions.RemoveEmptyEntries));
        //            allLines.RemoveRange(0, 4);
        //            int nodeIndex = 0;
        //            int connectivityIndex = 0;
        //            foreach (var line in allLines)
        //            {
        //                // in case of first line ...
        //                string separator = " ";
        //                string[] fields = line.Split(separator.ToCharArray()); //(new string[] { "\t" }, StringSplitOptions.RemoveEmptyEntries);
        //                if (fields[0] == "v")
        //                {
        //                    nodeIndex = nodeIndex + 1;
        //                    var node = new Node(double.Parse(fields[1]), double.Parse(fields[2]), double.Parse(fields[3]));
        //                    nodes[nodeIndex] = node;
        //                }
        //                else if (fields[0] == "f")
        //                {
        //                    connectivityIndex = connectivityIndex + 1;
        //                    string separatorForNode = "/";
        //                    int[] elementNodes = new int[4];
        //                    for (int i = 0; i < 4; i++)
        //                    {
        //                        string[] fieldsForNode = fields[i + 1].Split(separatorForNode.ToCharArray());
        //                        elementNodes[i] = Int16.Parse(fieldsForNode[0]);
        //                    }
        //                    elementsConnectivity[connectivityIndex] = new Dictionary<int, int>() { { 1, elementNodes[0] }, { 2, elementNodes[1] }, { 3, elementNodes[2] }, { 4, elementNodes[3] } };
        //                }

        //            }
        //        }


        //    }
        //    catch (Exception ex)
        //    {
        //        throw ex;
        //    }
        //}

        //private void Button_Click_1(object sender, RoutedEventArgs e)
        //{
        //    PlotOBJMesh plotMesh = new PlotOBJMesh();
        //    plotMesh.Window_Loaded();
        //    //ViewportGraphics = plotMesh.MainViewport;
        //    //ViewportGraphics.UpdateLayout();
        //    Window secondWindow = new Window();
        //    secondWindow.Show();
        //    secondWindow.KeyDown += plotMesh.Window_KeyDown;
        //    secondWindow.Content = plotMesh.MainViewport;
        //    //ViewportGraphics.InvalidateVisual();
        //    //ViewportGraphics.Children.Add(plotMesh.finalModel);
        //}
    }


}
