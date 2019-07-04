﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
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

namespace GFEC
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        public MainWindow()
        {
            InitializeComponent();
            LoadComboBox();
        }

        private void RunButton(object sender, RoutedEventArgs args)
        {
            SolveSelectedExample();
        }

        private void LoadComboBox()
        {
            List<string> exampleList = new List<string>();
            exampleList.Add("LinearTrussExample");
            exampleList.Add("TwoQuadsExample");
            
            ComboBox1.ItemsSource = exampleList;
        }

        private void SolveSelectedExample()
        {
            double[] solution;
            string selectedExample = ComboBox1.SelectedItem.ToString();
            switch (selectedExample)
            {
                case "TwoQuadsExample":
                    solution = TwoQuadsExample.RunStaticExample();
                    break;
                case "LinearTrussExample":
                    solution = LinearTrussExample.RunExample();
                    break;
                default:
                    solution = TwoQuadsExample.RunStaticExample();
                    break;
            }
            Results.Text = solution[0].ToString();

        }
    }

    
}
