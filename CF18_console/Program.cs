using System;
using System.IO;
using System.IO.Ports;
using System.Text;
using System.Threading;
using System.Collections.Generic;
using CF_Slander;
using System.Numerics;
using System.Timers;
using System.Diagnostics;

namespace CommunicationProtocol
{
    class Program
    {
        static Object locker = new Object();
        static int CursorReceiveData = 7;
        static int CursorCalibration = 10;
        static int CursorSendData = 14;

        static byte imu_total = 3;
        static int init_time = 20;
        static int baud_rate = 500000;
        static string com_port = "COM3";

        static int gyro_dataset_size = 10000;
        static int magnet_dataset_size = 10000;
        static string time_stamp = DateTime.Now.ToString("dd.MM.yyyy.HH.mm.ss");
        static string data_folder = System.IO.Directory.GetCurrentDirectory() + "\\raw_cf18_data-" + time_stamp;

        static int MyPort = 5917;
        static string MyAdress = "";
        static string WorkStationAdress = "192.168.62.84"; //для проверки -- "127.0.0.1"
        static int WorkStationPort = 5918;
        static CF_sender MesSender = new CF_sender(WorkStationAdress, WorkStationPort);
        static int time_delay = 15;

        // Calibration
        static private System.Timers.Timer time_calib;
        static bool calibration_finished = false;
        static int prev_ms = 0;
        static int calibration_itr = 0;
        static List<Vector4> calibration_data = new List<Vector4>();

        static void Collect_IMU_raw_data(ArduinoController controller, int dataset_size, int side, string title)
        {
            var csv_files = new Dictionary<int, StringBuilder>();
            for (byte i = 0; i < imu_total; ++i)
            {
                StringBuilder csv_file = new StringBuilder();
                csv_files.Add(i, csv_file);
            }

            int data_counter = 0;
            while (true)
            {
                if (controller.ReadRawData(out Dictionary<int, RawData> data, out UInt32 time, out ulong package) > 0)
                {
                    data_counter++;
                    lock (locker)
                    {
                        Console.SetCursorPosition(0, CursorReceiveData);
                        Console.Write(new string(' ', Console.WindowWidth));
                        Console.Write("\rReceived data, package {0}: {1}", side, data_counter);
                    }

                    foreach (var csv_pair in csv_files)
                    {
                        var imu_id = csv_pair.Key;
                        var csv_file = csv_pair.Value;
                        //AddRawDataInFile(imu_id, data[imu_id], csv_file);
                        AddRawDataInFile_DifferentialNet(imu_id, data[imu_id], csv_file);
                    }

                    if (data_counter >= dataset_size)
                    {
                        lock (locker)
                        {
                            //Console.WriteLine();
                            Console.SetCursorPosition(0, CursorReceiveData);
                            Console.Write(new string(' ', Console.WindowWidth));
                        }

                        foreach (var csv_pair in csv_files)
                        {
                            var imu_id = csv_pair.Key;
                            var csv_file = csv_pair.Value;

                            string path = string.Format("{0}\\{1}-{2}-{3}-.csv", data_folder, title, imu_id, side);
                            System.IO.File.WriteAllText(path, csv_file.ToString());
                        }
                        break;
                    }
                }
            }
        }

        static void Collect_Magnetometr_Calibration_Data(ArduinoController controller)
        {
            Collect_IMU_raw_data(controller, magnet_dataset_size, 0, "magnet");
        }

        static void Collect_Gyro_Accel_Calibration_Data(ArduinoController controller)
        {
            int sides = 1000;

            //controller.Wait(timeout, string.Format("Зафиксируйте датчик на стороне {0}/{1}", 1, sides));
            Console.Beep();

            for (int side_itr = 1; side_itr <= sides; ++side_itr)
            {
                Collect_IMU_raw_data(controller, gyro_dataset_size, side_itr, "cf");

                /*
                if (side_itr != sides)
                {
                    controller.Wait(timeout, string.Format("Переверните датчик на строну {0}/{1}", side_itr + 1, sides));
                    Console.Beep();
                }
                */
            }
        }

        static void AddRawDataInFile(int imu, RawData row, StringBuilder file)
        {
            if (file.Length == 0)
            {
                file.AppendLine("imu,server_time,arduino_time,ax,ay,az,gx,gy,gz,mx,my,mz");
            }
            var server_time = DateTime.Now.ToString("yyyy-MM-ddTHH:mm:ss.ffffff");
            var newLine = string.Format("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11}",
                                imu, server_time, row.time,
                                row.ax, row.ay, row.az,
                                row.gx, row.gy, row.gz,
                                row.mx, row.my, row.mz);
            file.AppendLine(newLine);
        }

        static void AddRawDataInFile_DifferentialNet(int imu, RawData row, StringBuilder file)
        {
            if (file.Length == 0)
            {
                file.AppendLine("server_time;(ax,ay,az);(gx,gy,gz)");
            }
            var server_time = DateTime.Now.ToString("yyyy-MM-ddTHH:mm:ss.ffffff");
            var newLine = string.Format("{0};({1},{2},{3});({4},{5},{6})",
                                server_time,
                                row.ax, row.ay, row.az,
                                row.gx, row.gy, row.gz);
            file.AppendLine(newLine);
        }

        static void PrintAttentionMessage(string message)
        {
            var default_color = Console.ForegroundColor;
            Console.ForegroundColor = ConsoleColor.Green;
            Console.WriteLine(message);
            Console.ForegroundColor = default_color;
        }

        static void PrintErrorMessage(string message)
        {
            var default_color = Console.ForegroundColor;
            Console.ForegroundColor = ConsoleColor.Red;
            Console.WriteLine(message);
            Console.ForegroundColor = default_color;
        }

        static void Collect_calibration_data()
        {
            var default_color = Console.ForegroundColor;

            var controller = new ArduinoController(imu_total);
            if (controller.Connect(com_port, baud_rate) < 0)
            {
                Console.WriteLine("No connection");
            }

            Console.WriteLine("Controller initialization, please wait {0} second...", init_time);
            controller.Wait(init_time);

            if (controller.SetParameterToArduino(Parameters.Mode, (byte)Modes.RawData) <= 0)
            {
                Console.WriteLine("Can't set arduino parameters. Calibration process is stoped.");
                return;
            }

            var motionProc = new MotionProcessor(AccelRange.ACCEL_RANGE_16G, GyroRange.GYRO_RANGE_2000DPS);

            Console.WriteLine("\nCollection gyroscope and accelerometr data:");
            PrintAttentionMessage("Attention: IMU sensor must be hard fixed on each of 6 side!!!");

            Collect_Gyro_Accel_Calibration_Data(controller);

            Console.WriteLine("\nCalibration data was written to path:");
            Console.WriteLine(data_folder);
        }

        static string ChosePort()
        {
            string port = "";
            var ports = SerialPort.GetPortNames();
            if (ports.Length == 0)
            {
                PrintErrorMessage("COM ports wasn't found.");
            }
            else if (ports.Length == 1)
            {
                port = ports[0];
                PrintAttentionMessage(string.Format("Choosen port: {0}", port));
            }
            else
            {
                string com_ports_str = "";
                foreach (var p in ports)
                {
                    if (com_ports_str.Length != 0)
                    {
                        com_ports_str += ", ";
                    }
                    com_ports_str += p;
                }
                PrintAttentionMessage(string.Format("Available COM ports: {0}", com_ports_str));
                Console.WriteLine("Chose COM port: ");
                port = Console.ReadLine();
                if (port.Length == 0)
                {
                    port = ports[1];
                }
            }
            return port;
        }

        static List<Vector4> PrepareDataForCF18(string path)
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

        static void AddCFControlDataInFile(Vector4 row, StringBuilder file)
        {
            if (file.Length == 0)
            {
                file.AppendLine("Timestamp;X;Y;Z;GF");
            }
            var time_stamp = DateTime.Now.ToString("yyyy-MM-ddTHH:mm:ss.ffffff");
            var newLine = string.Format("{0};{1};{2};{3};{4}", time_stamp, row.X, row.Y, row.Z, row.W);
            file.AppendLine(newLine);
        }


        static void SendDataToCF18(int delay)
        {
            while (true)
            {
                if (calibration_finished)
                {
                    break;
                }
                Thread.Sleep(1000);
            }

            var data = PrepareDataForCF18(@"./wxyz.csv");

            lock (locker)
            {
                Console.SetCursorPosition(0, CursorSendData - 1);
                Console.WriteLine("Press any key to start sending data to ЦФ-18");
            }
            Console.ReadKey();

            StringBuilder csv_file = new StringBuilder();

            prev_ms = DateTime.Now.Millisecond;
            var curr_ms = 0;
            var pack_itr = 0;
            var e = Stopwatch.StartNew();

            foreach (var vec in data)
            {
                //curr_ms = DateTime.Now.Millisecond;
                MesSender.SendGXYZ(vec);
                AddCFControlDataInFile(vec, csv_file);
                lock (locker)
                {
                    Console.SetCursorPosition(0, CursorSendData);
                    Console.Write(new string(' ', Console.WindowWidth));
                    Console.Write("Sended package: {0}, dt: {1}, x: {2:0.000}, y: {3:0.000}, z: {4:0.000}, gf: {5:0.000}",
                        pack_itr, e.ElapsedMilliseconds, vec.X, vec.Y, vec.Z, vec.W);
                    //Console.WriteLine("Package sended: {0}", pack_itr);
                    //prev_ms = curr_ms;
                    e = Stopwatch.StartNew();
                    ++pack_itr;
                }
                Thread.Sleep(50);
            }
            string time_stamp = DateTime.Now.ToString("yyyy-MM-ddTHH.mm.ss.ffffff");
            string path = string.Format("{0}\\cf_control_data-{1}.csv", data_folder, time_stamp);
            System.IO.File.WriteAllText(path, csv_file.ToString());
            e.Stop();
        }

        private static void SendCalibrationComands()
        {
            if (calibration_itr < calibration_data.Count)
            {
                var vec = calibration_data[calibration_itr];
                MesSender.SendGXYZ(vec);
                int curr_ms = DateTime.Now.Millisecond;
                lock (locker)
                {
                    Console.SetCursorPosition(0, CursorCalibration);
                    Console.Write(new string(' ', Console.WindowWidth));
                    Console.Write("Package sended: {0}, dt: {1}, x: {2:0.000}, y: {3:0.000}, z: {4:0.000}, gf: {5:0.000}",
                        calibration_itr, curr_ms - prev_ms, vec.X, vec.Y, vec.Z, vec.W);
                    prev_ms = curr_ms;
                }
                calibration_itr++;
            }
            else
            {
                lock (locker)
                {
                    Console.SetCursorPosition(0, CursorCalibration + 1);
                    Console.WriteLine("Calibration finished");
                }
                calibration_finished = true;
                time_calib.Stop();
            }
        }

        private static void OnTimedEvent_Calibration(Object source, ElapsedEventArgs e)
        {
            SendCalibrationComands();
        }

        private static void SendDataForCalibrationToCF18()
        {
            lock (locker)
            {
                Console.SetCursorPosition(0, CursorCalibration - 1);
                Console.Write("Press any key to start calibration");
            }
            Console.ReadKey();

            prev_ms = DateTime.Now.Millisecond;

            var calibration_data_tmp = PrepareDataForCF18(@"./calib.csv");
            calibration_data.Clear();
            foreach (var vec in calibration_data_tmp)
            {
                for (int i = 0; i < 1000; ++i)
                {
                    calibration_data.Add(vec);
                }
            }
            calibration_itr = 0;

            SendCalibrationComands();
            time_calib = new System.Timers.Timer(10);
            time_calib.Elapsed += OnTimedEvent_Calibration;
            time_calib.Start();
        }

        static void Main(string[] args)
        {
            Console.Write("Chose time delay: ");
            var delay = Console.Read();
            com_port = ChosePort();
            System.IO.DirectoryInfo di = System.IO.Directory.CreateDirectory(data_folder);

            new Thread(() =>
            {
                Thread.CurrentThread.IsBackground = true;
                Collect_calibration_data();
            }).Start();

            Thread.Sleep(25000);

            SendDataForCalibrationToCF18();
            SendDataToCF18(delay);

            Console.ReadKey();
        }
    }
}
