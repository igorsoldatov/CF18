using System;
using System.Net.Sockets;
using System.Text;
using System.Threading;
using System.ComponentModel;
using System.Runtime.InteropServices;
using Utils;
using System.Numerics;

namespace CF_Slander
{
    /// <summary>
    /// Класс для формирования и отправки сообщений на рабочую станцию
    /// </summary>
    public class CF_sender
    {
        public UdpClient SendingUdpClient { get; private set; }
        public string Address { get; private set; }
        public int Port { get; private set; }
        

        /// <summary>
        /// Конструктор
        /// </summary>
        /// <param name="address">адрес рабочей станции</param>
        /// <param name="port">порт на который слать</param>
        public CF_sender(string address, int port)
        {
            Address = address;
            Port = port;
            SendingUdpClient = new UdpClient();
        }

        /// <summary>
        /// Создание сообщения из координат вектора перегрузки и его отправка
        /// </summary>
        /// <param name="nx">nx</param>
        /// <param name="ny">ny</param>
        /// <param name="nz">nz</param>
        public void Send(float nx, float ny, float nz)
        {
            CF_18 data_info = new CF_18();

            data_info.fPX = nx;
            data_info.fPY = ny;
            data_info.fPZ = nz;

            byte[] toBytes = StructToBytes(data_info);

            SendingUdpClient.Send(toBytes, toBytes.Length, Address, Port);
        }

        /// <summary>
        /// Создание сообщения из координат вектора перегрузки в формате Vec3 и его отправка
        /// </summary>
        /// <param name="nVec">вектор перегрузки</param>
        public void SendNVec(Vec3 nVec)
        {
            CF_18 data_info = new CF_18();

            data_info.fPX = (float)nVec.x;
            data_info.fPY = (float)nVec.y;
            data_info.fPZ = (float)nVec.z;

            byte[] toBytes = StructToBytes(data_info);
            SendingUdpClient.Send(toBytes, toBytes.Length, Address, Port);
        }

        public void SendGXYZ(Vector4 xyzVec)
        {
            CF_18 data_info = new CF_18();

            data_info.fFPX = xyzVec.W;
            data_info.fUX = xyzVec.X;
            data_info.fUY = xyzVec.Y;
            data_info.fUZ = xyzVec.Z;

            byte[] toBytes = StructToBytes(data_info);
            SendingUdpClient.Send(toBytes, toBytes.Length, Address, Port);

        }

        /// <summary>
        /// Еще один вариант создания сообщения нужного формата
        /// из координат вектора перегрузки  и его отправка
        /// </summary>
        /// <param name="nx">nx</param>
        /// <param name="ny">ny</param>
        /// <param name="nz">nz</param>
        public void TrySendNVec(float nx, float ny, float nz)
        {
            CFMessage cfmes = new CFMessage(nx, ny, nz);
            var data = cfmes.OutMessage;
            SendingUdpClient.Send(data, data.Length, Address, Port);
        }
        /// <summary>
        /// Еще один вариант создания сообщения нужного формата
        /// из вектора Vec3 и его отправка
        /// </summary>
        /// <param name="nVec">вектор перегрузки</param>
        public void TrySendNVec(Vec3 nVec)
        {
            CFMessage cfmes = new CFMessage(nVec);
            var data = cfmes.OutMessage;
            SendingUdpClient.Send(data, data.Length, Address, Port);
        }

        


        /// <summary>
        /// Получает представление класса в виде массива байтов
        /// </summary>
        /// <param name="str">объект CF_18</param>
        /// <returns>массив байтов в формате Big-Endian, готовый для дальнейшей отправки</returns>
        private byte[] StructToBytes(object str)
        {
            int width = sizeof(int);
            int size = Marshal.SizeOf(str.GetType());
            //Console.WriteLine(size + "\n");
            byte[] arr = new byte[size];
            IntPtr ptr = Marshal.AllocHGlobal(size);
            Marshal.StructureToPtr(str, ptr, true);
            Marshal.Copy(ptr, arr, 0, size);
            Marshal.FreeHGlobal(ptr);

            if (BitConverter.IsLittleEndian)
            {
                //Console.WriteLine("little-endian, trying to convert to big-endian");
                arr = reverse(arr, size, width);
            }
            else
            {
                //Console.WriteLine("most likely big-endian, let's try to send as it is");
            }

            return arr;
        }

        /// <summary>
        /// Конвертатирует представление строки 
        /// в виде массива байтов из Little-Endian в Big-Endian
        /// </summary>
        /// <param name="_arr">строка в виде массива байтов</param>
        /// <param name="n">длина массива</param>
        /// <param name="k">размер элемента</param>
        /// <returns>Сконвертированный массив байтов</returns>
        private byte[] reverse(byte[] _arr, int n, int k)
        {
            for (int i = 0; i < n; i += k)
            {
                int left = i;

                int right = Math.Min(i + k - 1, n - 1);
                byte temp;

                while (left < right)
                {
                    temp = _arr[left];
                    _arr[left] = _arr[right];
                    _arr[right] = temp;
                    left += 1;
                    right -= 1;
                }
            }

            return _arr;
        }

        /// <summary>
        /// Представление отправляемой строки в виде класса
        /// </summary>
        [StructLayout(LayoutKind.Sequential)]
        private class CF_18
        {
            [Description("управляющая команда")] public int command = 57;                   //header        0

            [Description("общая длина пакета в байтах")] public int total_lengh = 0;        //spacer #1     1
            [Description("идентификатор отправителя")] public int sender_id = 0;            //spacer #2     2
            [Description("идентификатор получателя")] public int reseiver_id = 0;           //spacer #3     3
            [Description("цикл модели в мс")] public uint time_stamp = 0;                   //spacer #4     4
            [Description("дата")] public int jdate = 0;                                     //spacer #5     5
            [Description("время")] public int jtime = 0;                                    //spacer #6     6
            [Description("текущий номер сообщения")] public uint message_no = 0;            //spacer #7     7
            [Description("общее количество сообщений")] public uint total_messages = 0;     //spacer #8     8
            [Description("текущее количество записей")] public uint no_of_records = 0;      //spacer #9     9
            [Description("общее количество  записей")] public uint total_no_records = 0;    //spacer #10    10
            [Description("зарезервировано")] public uint reserved1 = 0;                     //spacer #11    11
            [Description("зарезервировано")] public uint reserved2 = 0;                     //spacer #12    12
            [Description("зарезервировано")] public uint reserved3 = 0;                     //spacer #13    13
            [Description("зарезервировано")] public uint reserved4 = 0;                     //spacer #14    14
            [Description("зарезервировано")] public uint reserved5 = 0;                     //spacer #15    15

            [Description("Перегрузка по оси X")] public float fPX = 0f;                      //fx            16
            [Description("Перегрузка по оси Y")] public float fPY = 0f;                      //fy            17
            [Description("Перегрузка по оси Z")] public float fPZ = 0f;                      //fz            18

            [Description("угол крена")] public float fUX = 0f;                               //spacer #1     19
            [Description("угол курса")] public float fUY = 0f;                               //spacer #2     20
            [Description("угол тангажа")] public float fUZ = 0f;                             //spacer #3     21
            [Description("первая производная перегрузки по оси X")] public float fFPX = 0f;  //spacer #4     22
            [Description("первая производная перегрузки по оси Y")] public float fFPY = 0f;  //spacer #5     23
            [Description("первая производная перегрузки по оси Z")] public float fFPZ = 0f;  //spacer #6     24
            [Description("вторая производная перегрузки по оси X")] public float fSPX = 0f;  //spacer #7     25
            [Description("вторая производная перегрузки по оси Y")] public float fSPY = 0f;  //spacer #8     26
            [Description("вторая производная перегрузки по оси Z")] public float fSPZ = 0f;  //spacer #9     27
            [Description(" переменная старта")] public int STAT = 0;                        //spacer #10    28
        }
    }
}
