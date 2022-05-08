using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GFEC
{
    public class Object3D
    {
        PlotOBJMesh Mesh3D { get; set; }
        Dictionary<int, INode> Nodes { get; set; }
        Dictionary<int, Dictionary<int, int>> ElementsList { get; set; }

        public Object3D(Dictionary<int, INode> nodes, Dictionary<int, Dictionary<int, int>> elementsList)
        {
            Mesh3D = new PlotOBJMesh();
            Nodes = nodes;
            ElementsList = elementsList;
        }

        public void Transfrom()
        {
            Dictionary<int, Dictionary<int, int>> transformedList = new Dictionary<int, Dictionary<int, int>>();
            int k = 0;
            foreach (var element in ElementsList)
            {
                Dictionary<int, int> face1NodesList = new Dictionary<int, int>();
                Dictionary<int, int> face2NodesList = new Dictionary<int, int>();
                Dictionary<int, int> face3NodesList = new Dictionary<int, int>();
                Dictionary<int, int> face4NodesList = new Dictionary<int, int>();
                Dictionary<int, int> face5NodesList = new Dictionary<int, int>();
                Dictionary<int, int> face6NodesList = new Dictionary<int, int>();

                face1NodesList.Add(1, element.Value[1]);
                face1NodesList.Add(2, element.Value[2]);
                face1NodesList.Add(3, element.Value[3]);
                face1NodesList.Add(4, element.Value[4]);

                face2NodesList.Add(1, element.Value[5]);
                face2NodesList.Add(2, element.Value[6]);
                face2NodesList.Add(3, element.Value[7]);
                face2NodesList.Add(4, element.Value[8]);

                face3NodesList.Add(1, element.Value[1]);
                face3NodesList.Add(2, element.Value[2]);
                face3NodesList.Add(3, element.Value[5]);
                face3NodesList.Add(4, element.Value[6]);

                face4NodesList.Add(1, element.Value[4]);
                face4NodesList.Add(2, element.Value[3]);
                face4NodesList.Add(3, element.Value[8]);
                face4NodesList.Add(4, element.Value[7]);

                face5NodesList.Add(1, element.Value[4]);
                face5NodesList.Add(2, element.Value[1]);
                face5NodesList.Add(3, element.Value[8]);
                face5NodesList.Add(4, element.Value[5]);

                face6NodesList.Add(1, element.Value[3]);
                face6NodesList.Add(2, element.Value[2]);
                face6NodesList.Add(3, element.Value[7]);
                face6NodesList.Add(4, element.Value[6]);

                k = k + 1;
                transformedList.Add(k, face1NodesList);
                k = k + 1;
                transformedList.Add(k, face2NodesList);
                k = k + 1;
                transformedList.Add(k, face3NodesList);
                k = k + 1;
                transformedList.Add(k, face4NodesList);
                k = k + 1;
                transformedList.Add(k, face5NodesList);
                k = k + 1;
                transformedList.Add(k, face6NodesList);
            }
        }
    }
}
