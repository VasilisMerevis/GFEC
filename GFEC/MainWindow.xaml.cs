using System;
using System.Collections.Generic;
using System.ComponentModel;
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
using System.Windows.Media.Media3D;
using System.Windows.Navigation;
using System.Windows.Shapes;
using LiveCharts;
using LiveCharts.Configurations;
using LiveCharts.Wpf;
using Microsoft.Win32;

namespace GFEC
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window, INotifyPropertyChanged
    {
        public SeriesCollection Graph { get; set; }
        public SeriesCollection Mesh { get; set; }
        private Results solverResults;
        private Dictionary<int, INode> nodes = new Dictionary<int, INode>();
        private Dictionary<int, Dictionary<int, int>> elementsConnectivity = new Dictionary<int, Dictionary<int, int>>();
        public string selectedExample;
        public SeriesCollection Something { get; set; }
        public ChartValues<ConvergenceValues> ChartValues { get; set; }
        private int LoadStepNumber { get; set; }



        public MainWindow()
        {
            InitializeComponent();
            LoadComboBox();
            gnuplotImage.Source = null;
            ConvergenceResults();

            
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
                    CNTsInAngleFinalExample.thermalSolution = new StaticSolver();
                    CNTsInAngleFinalExample.thermalSolution.NonLinearScheme = new LoadControlledNewtonRaphson();
                    CNTsInAngleFinalExample.thermalSolution.NonLinearScheme.convergenceResult += NonLinearScheme_convergenceResult;
                    
                    finalResults = CNTsInAngleFinalExample.RunStaticExample();

                    break;
                default:
                    finalResults = TwoQuadsExample.RunStaticExample();
                    break;
            }
            //Results.Text = solution[0].ToString();

            solverResults = finalResults;
        }

        public void NonLinearScheme_convergenceResult(object sender, ConvergenceValues e)
        {
            Application.Current.Dispatcher.Invoke(new Action(() =>
            { 
                LogTool.Text = "Load Step " + e.LoadStep + "-Iteration " + e.Iteration + " : Convergence State: " + e.ConvergenceResult + " with residual " + e.ResidualNorm;
                if (e.LoadStep > LoadStepNumber)
                {
                    ChartValues.Clear();
                }
                //ChartValues.Add(new ConvergenceValues
                //{
                //    Iteration = e.Iteration,
                //    ResidualNorm = e.ResidualNorm
                //});
                
                LoadStepNumber = e.LoadStep;
            }));
            SetAxisLimits(e.Iteration);
            
            if (ChartValues.Count > 50) ChartValues.RemoveAt(0);
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
            /// Declare scene objects.
            Viewport3D myViewport3D = new Viewport3D();
            Model3DGroup myModel3DGroup = new Model3DGroup();
            GeometryModel3D myGeometryModel = new GeometryModel3D();
            ModelVisual3D myModelVisual3D = new ModelVisual3D();
            // Defines the camera used to view the 3D object. In order to view the 3D object,
            // the camera must be positioned and pointed such that the object is within view
            // of the camera.
            OrthographicCamera myPCamera = new OrthographicCamera();

            // Specify where in the 3D scene the camera is.
            myPCamera.Position = new Point3D(0, 0, 2);

            // Specify the direction that the camera is pointing.
            myPCamera.LookDirection = new Vector3D(0, 0, -1);

            // Define camera's horizontal field of view in degrees.
            //myPCamera.FieldOfView = 60;

            // Asign the camera to the viewport
            myViewport3D.Camera = myPCamera;
            // Define the lights cast in the scene. Without light, the 3D object cannot
            // be seen. Note: to illuminate an object from additional directions, create
            // additional lights.
            DirectionalLight myDirectionalLight = new DirectionalLight();
            myDirectionalLight.Color = Colors.White;
            myDirectionalLight.Direction = new Vector3D(-0.61, -0.5, -0.61);

            myModel3DGroup.Children.Add(myDirectionalLight);

            // The geometry specifes the shape of the 3D plane. In this sample, a flat sheet
            // is created.
            MeshGeometry3D myMeshGeometry3D = new MeshGeometry3D();

            // Create a collection of normal vectors for the MeshGeometry3D.
            Vector3DCollection myNormalCollection = new Vector3DCollection();
            myNormalCollection.Add(new Vector3D(0, 0, 1));
            myNormalCollection.Add(new Vector3D(0, 0, 1));
            myNormalCollection.Add(new Vector3D(0, 0, 1));
            myNormalCollection.Add(new Vector3D(0, 0, 1));
            myNormalCollection.Add(new Vector3D(0, 0, 1));
            myNormalCollection.Add(new Vector3D(0, 0, 1));
            myMeshGeometry3D.Normals = myNormalCollection;

            // Create a collection of vertex positions for the MeshGeometry3D.
            Point3DCollection myPositionCollection = new Point3DCollection();
            myPositionCollection.Add(new Point3D(-0.5, -0.5, 0.5));
            myPositionCollection.Add(new Point3D(0.5, -0.5, 0.5));
            myPositionCollection.Add(new Point3D(0.5, 0.5, 0.5));
            myPositionCollection.Add(new Point3D(0.5, 0.5, 0.5));
            myPositionCollection.Add(new Point3D(-0.5, 0.5, 0.5));
            myPositionCollection.Add(new Point3D(-0.5, -0.5, 0.5));
            myMeshGeometry3D.Positions = myPositionCollection;

            // Create a collection of texture coordinates for the MeshGeometry3D.
            PointCollection myTextureCoordinatesCollection = new PointCollection();
            myTextureCoordinatesCollection.Add(new Point(0, 0));
            myTextureCoordinatesCollection.Add(new Point(1, 0));
            myTextureCoordinatesCollection.Add(new Point(1, 1));
            myTextureCoordinatesCollection.Add(new Point(1, 1));
            myTextureCoordinatesCollection.Add(new Point(0, 1));
            myTextureCoordinatesCollection.Add(new Point(0, 0));
            myMeshGeometry3D.TextureCoordinates = myTextureCoordinatesCollection;

            // Create a collection of triangle indices for the MeshGeometry3D.
            Int32Collection myTriangleIndicesCollection = new Int32Collection();
            myTriangleIndicesCollection.Add(0);
            myTriangleIndicesCollection.Add(1);
            myTriangleIndicesCollection.Add(2);
            myTriangleIndicesCollection.Add(3);
            myTriangleIndicesCollection.Add(4);
            myTriangleIndicesCollection.Add(5);
            myMeshGeometry3D.TriangleIndices = myTriangleIndicesCollection;

            // Apply the mesh to the geometry model.
            myGeometryModel.Geometry = myMeshGeometry3D;

            // The material specifies the material applied to the 3D object. In this sample a
            // linear gradient covers the surface of the 3D object.

            // Create a horizontal linear gradient with four stops.
            LinearGradientBrush myHorizontalGradient = new LinearGradientBrush();
            myHorizontalGradient.StartPoint = new Point(0, 0.5);
            myHorizontalGradient.EndPoint = new Point(1, 0.5);
            myHorizontalGradient.GradientStops.Add(new GradientStop(Colors.Yellow, 0.0));
            myHorizontalGradient.GradientStops.Add(new GradientStop(Colors.Red, 0.25));
            myHorizontalGradient.GradientStops.Add(new GradientStop(Colors.Blue, 0.75));
            myHorizontalGradient.GradientStops.Add(new GradientStop(Colors.LimeGreen, 1.0));

            // Define material and apply to the mesh geometries.
            DiffuseMaterial myMaterial = new DiffuseMaterial(myHorizontalGradient);
            myGeometryModel.Material = myMaterial;

            // Apply a transform to the object. In this sample, a rotation transform is applied,
            // rendering the 3D object rotated.
            RotateTransform3D myRotateTransform3D = new RotateTransform3D();
            AxisAngleRotation3D myAxisAngleRotation3d = new AxisAngleRotation3D();
            myAxisAngleRotation3d.Axis = new Vector3D(0, 3, 0);
            myAxisAngleRotation3d.Angle = 40;
            myRotateTransform3D.Rotation = myAxisAngleRotation3d;
            myGeometryModel.Transform = myRotateTransform3D;

            // Add the geometry model to the model group.
            myModel3DGroup.Children.Add(myGeometryModel);

            // Add the group of models to the ModelVisual3d.
            myModelVisual3D.Content = myModel3DGroup;

            //
            myViewport3D.Children.Add(myModelVisual3D);

            // Apply the viewport to the page so it will be rendered.
            //this.Content = myViewport3D;
            Window testWindow = new Window();
            testWindow.Content = myViewport3D;
            testWindow.Show();
            //myViewport3D2.DataContext = myViewport3D;

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

        
        private void Button_ParallelTest(object sender, RoutedEventArgs e)
        {
            MultiThreadingExample test1 = new MultiThreadingExample();
            test1.timeElapsed += PrintResultOnUI;
            test1.RunExample();
        }

        private double _axisMax;
        private double _axisMin;

        
        public Func<double, string> DateTimeFormatter { get; set; }
        public double AxisStep { get; set; }
        public double AxisUnit { get; set; }
        public bool IsReading { get; set; }
        public event PropertyChangedEventHandler PropertyChanged;
        public void ConvergenceResults()
        {
            var mapper = Mappers.Xy<ConvergenceValues>()
                .X(model => model.Iteration)   //use DateTime.Ticks as X
                .Y(model => model.ResidualNorm);           //use the value property as Y

            //lets save the mapper globally.
            Charting.For<ConvergenceValues>(mapper);

            //the values property will store our values array
            ChartValues = new ChartValues<ConvergenceValues>();

            //lets set how to display the X Labels
            //DateTimeFormatter = value => new DateTime((long)value).ToString("mm:ss");

            //AxisStep forces the distance between each separator in the X axis
            AxisStep = 1;// TimeSpan.FromSeconds(1).Ticks;
            //AxisUnit forces lets the axis know that we are plotting seconds
            //this is not always necessary, but it can prevent wrong labeling
            AxisUnit = 100;// TimeSpan.TicksPerSecond;
        }

        public double AxisMax
        {
            get { return _axisMax; }
            set
            {
                _axisMax = value;
                OnPropertyChanged("AxisMax");
            }
        }
        public double AxisMin
        {
            get { return _axisMin; }
            set
            {
                _axisMin = value;
                OnPropertyChanged("AxisMin");
            }
        }

        private void SetAxisLimits(int now)
        {
            AxisMax = now + 1;
            AxisMin = now - 8;
        }

        protected virtual void OnPropertyChanged(string propertyName = null)
        {
            if (PropertyChanged != null)
                PropertyChanged.Invoke(this, new PropertyChangedEventArgs(propertyName));
        }

    }


}
