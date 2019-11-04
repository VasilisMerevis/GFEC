﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GFEC
{
    class ContactNtN2DTh : IElement
    {
        public Dictionary<int, INode> Nodes { get; }
        public IElementProperties Properties { get; set; }
        public Dictionary<int, bool[]> ElementFreedomSignature { get; } = new Dictionary<int, bool[]>();
        public List<int> ElementFreedomList { get; set; }
        public double[] DisplacementVector { get; set; }
        public double[] AccelerationVector { get; set; }
        private double PenaltyFactor { get; set; }
        private double ContactArea { get; set; }

        public ContactNtN2DTh(IElementProperties properties, Dictionary<int, INode> nodes)
        {
            Properties = properties;
            this.Nodes = nodes;
            ElementFreedomSignature[1] = new bool[] { true, false, false, false, false, false };
            ElementFreedomSignature[2] = new bool[] { true, false, false, false, false, false };
            DisplacementVector = new double[2];
            properties.SectionArea = ContactArea;
        }

        private double CalculateConductivity()
        {
            double cc = 1.0;
            double cH = cc * ContactArea;
            return cH;
        }

        private double CalculateTemperatureJump()
        {
            double theta1 = 100.0;
            double theta2 = 0.0;
            double gH = (theta2 + DisplacementVector[1]) - (theta1 + DisplacementVector[0]);
            return gH;
        }

        public double[,] CreateGlobalStiffnessMatrix()
        {
            double cH = 19.2;
            double[,] stiffMatrix = new double[,]
            {
                {1.0 * cH, -1.0 * cH },
                {-1.0 * cH, 1.0 * cH }
            };
            return stiffMatrix;
        }

        public double[] CreateInternalGlobalForcesVector()
        {
            double cH = CalculateConductivity();
            double gH = CalculateTemperatureJump();
            double[] heatFlux = new double[]
            {
                -1.0 * cH * gH ,
                 1.0 * cH * gH
            };
            return heatFlux;
        }

        public double[,] CreateMassMatrix()
        {
            return new double[4, 4];
        }

        public double[,] CreateDampingMatrix()
        {
            return new double[4, 4];
        }
    }
}
