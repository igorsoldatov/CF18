using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Utils;

namespace CF_Slander
{
    /// <summary>
    /// Еще один класс для создания сообщений нужного 
    /// Здесь метод более тупой: просто по частям собираем отправной массив
    /// </summary>
    public class CFMessage
    {
        private const int Header = 57;
        private const int Spacer = 0;

        public float[] NVector { get; private set; }
        public byte[] OutMessage { get; private set; }
        public CFMessage(Vec3 nVector)
        {
            NVector = new float[] { (float)nVector.x, (float)nVector.y, (float)nVector.z };
            OutMessage = GenerateMessage(NVector);
        }

        public CFMessage(float x, float y, float z)
        {
            NVector = new float[] { x, y, z };
            OutMessage = GenerateMessage(NVector);
        }

        private byte[] GenerateMessage(float[] vec)
        {
            int iw = sizeof(int);
            int fw = sizeof(float);
            int i = 0; 

            byte[] data = new byte[16 * iw + 3 * fw + 10 * iw];

            i = 0;
            byte[] converted = BitConverter.GetBytes(Header);
            if (BitConverter.IsLittleEndian)
            {
                Array.Reverse(converted);
            }

            for (int j = 0; j < iw; ++j)
            { 
                data[i * iw + j] = converted[j];
            }
            
            for (i = 1; i < 16; i ++)
            {
                converted = BitConverter.GetBytes(Spacer);
                if (BitConverter.IsLittleEndian)
                {
                    Array.Reverse(converted);
                }

                for (int j = 0; j < iw; ++j)
                {
                    data[i * iw + j] = converted[j];
                }
            }

            for (i = 16; i < 19; i++)
            {
                converted = BitConverter.GetBytes(vec[i-16]);
                if (BitConverter.IsLittleEndian)
                {
                    Array.Reverse(converted);
                }

                for (int j = 0; j < iw; ++j)
                {
                    data[i * fw + j] = converted[j];
                }
            }

            for (i = 19; i< 29; i++)
            {
                converted = BitConverter.GetBytes(Spacer);
                if (BitConverter.IsLittleEndian)
                {
                    Array.Reverse(converted);
                }

                for (int j = 0; j < iw; ++j)
                {
                    data[i * iw + j] = converted[j];
                }
            }

            return data;
        }
    }
}
