using System;
using System.IO;
using System.Collections.Generic;
using System.Text;
using CF_Slander;

namespace Send_and_collect_data_sf18_v3
{
    static public class Common
    {
        static public List<Vector4> PrepareDataForCF18(string path)
        {
            List<Vector4> data = new List<Vector4>();

            using (var reader = new StreamReader(path))
            {
                while (!reader.EndOfStream)
                {
                    var line = reader.ReadLine();
                    var values = line.Split(';');

                    if (values.Length == 4)
                    {
                        var vec = new Vector4()
                        {
                            X = float.Parse(values[0]),
                            Y = float.Parse(values[1]),
                            Z = float.Parse(values[2]),
                            W = float.Parse(values[3])

                        };
                        data.Add(vec);
                    }
                }
            }
            return data;
        }
    }
    

}
