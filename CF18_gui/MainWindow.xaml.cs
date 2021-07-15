using System;
using System.IO;
using System.IO.Ports;
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
using System.Threading;
using System.Timers;
using CommunicationProtocol;
using CF_Slander;
using System.Windows.Threading;

namespace Send_and_collect_data_sf18_v3
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        private Object locker = new object();

        private string WorkStationAdress = "192.168.62.84";
        private int WorkStationPort = 5918;
        private CF_sender MesSender;
        private ArduinoController controller;
        private MotionProcessor motionProc;

        private byte imu_total = 3;
        private int init_time = 20;
        private int baud_rate = 500000;
        private int dataset_size = 10000;
        private string data_folder = Directory.GetCurrentDirectory() + "\\raw_cf18_data-" + DateTime.Now.ToString("dd.MM.yyyy.HH.mm.ss");


        private System.Timers.Timer time_calib;
        private bool collecting_data = false;
        private int prev_ms = 0;
        private int command_timeout;
        private int command_itr = 0;
        private string data_name = "cf_data";
        private List<Vector4> commands_list = new List<Vector4>();
        private Dictionary<int, StringBuilder> imu_data_files = new Dictionary<int, StringBuilder>();
        private StringBuilder command_file;
        private ProgressBar current_progress;

        private bool CollectionData { 
            get { return collecting_data; } 
            set {
                collecting_data = value;
                Dispatcher.BeginInvoke(new Action(() =>
                {
                    btnCalibration.IsEnabled = !collecting_data;
                    btnBoltanka.IsEnabled = !collecting_data;
                }));
            }
        }

        private void ChoosePort()
        {
            string[] ports = SerialPort.GetPortNames();
            if (ports.Length > 0)
            {
                foreach (string port in ports)
                {
                    comboBoxCOMPorts.Items.Add(port);
                }
                comboBoxCOMPorts.SelectedItem = ports[0];
            }
        }

        public MainWindow()
        {
            InitializeComponent();

            btnCalibration.IsEnabled = false;
            btnBoltanka.IsEnabled = false;

            DirectoryInfo di = Directory.CreateDirectory(data_folder);
            MesSender = new CF_sender(WorkStationAdress, WorkStationPort);
            ChoosePort();
        }

        private void ConnectArduino(string com_port)
        {
            controller = new ArduinoController(imu_total);
            if (controller.Connect(com_port, baud_rate) < 0)
            {
                Dispatcher.BeginInvoke(new Action(() =>
                {
                    btnConnect.Content = "Not connected";
                    btnConnect.Background = Brushes.Red;
                }));
            }

            for (int i = init_time; i > 0; --i)
            {
                btnConnect.Dispatcher.Invoke(DispatcherPriority.Normal, new Action(() =>
                {
                    btnConnect.Content = string.Format("Conneting... {0}", i);
                }));
                Thread.Sleep(1000);
            }

            if (controller.SetParameterToArduino(Parameters.Mode, (byte)Modes.RawData) <= 0)
            {
                Dispatcher.BeginInvoke(new Action(() =>
                {
                    btnConnect.Content = "Not connected";
                    btnConnect.Background = Brushes.Red;
                }));
                return;
            }

            motionProc = new MotionProcessor(AccelRange.ACCEL_RANGE_16G, GyroRange.GYRO_RANGE_2000DPS);

            Dispatcher.BeginInvoke(new Action(() =>
            {
                btnConnect.Content = "Connected";
                btnConnect.Background = Brushes.Green;
                btnCalibration.IsEnabled = true;
                btnBoltanka.IsEnabled = true;
            }));

            ThreadPool.QueueUserWorkItem(o =>
            {
                Collect_IMU_raw_data(controller);
            });
        }

        private void InitCollectionData()
        {
            imu_data_files.Clear();
            for (byte i = 0; i < imu_total; ++i)
            {
                StringBuilder csv_file = new StringBuilder();
                imu_data_files.Add(i, csv_file);
            }

            command_file = new StringBuilder();
        }

        private void FlushData()
        {
            string server_time = DateTime.Now.ToString("yyyy.MM.dd.HH.mm.ss.ffffff");
            string path = "";
            foreach (KeyValuePair<int, StringBuilder> csv_pair in imu_data_files)
            {
                int imu_id = csv_pair.Key;
                StringBuilder csv_file = csv_pair.Value;

                path = string.Format("{0}\\{1}-{2}-{3}-.csv", data_folder, data_name, imu_id, server_time);
                File.WriteAllText(path, csv_file.ToString());
            }
            imu_data_files.Clear();


            if (command_file != null)
            {
                path = string.Format("{0}\\commands-{1}-{2}-.csv", data_folder, data_name, server_time);
                File.WriteAllText(path, command_file.ToString());
                command_file = null;
            }
        }

        private void Collect_IMU_raw_data(ArduinoController controller)
        {
            int data_counter = 0;
            while (true)
            {
                lock (locker)
                {
                    if (CollectionData)
                    {
                        if (controller.ReadRawData(out Dictionary<int, RawData> data, out UInt32 time, out ulong package) > 0)
                        {
                            data_counter++;
                            lblReceiveData.Dispatcher.Invoke(new Action(() =>
                            {
                                lblReceiveData.Content = string.Format("Received data: {0}", data_counter);
                            }));

                            string server_time = DateTime.Now.ToString("yyyy-MM-ddTHH:mm:ss.ffffff");
                            foreach (KeyValuePair<int, StringBuilder> csv_pair in imu_data_files)
                            {
                                int imu_id = csv_pair.Key;
                                StringBuilder csv_file = csv_pair.Value;
                                AddRawDataInFile(imu_id, data[imu_id], server_time, csv_file);
                            }
                        }
                    } else
                    {
                        if (imu_data_files.Count > 0)
                        {
                            FlushData();
                            data_counter = 0;
                        }
                    }
                }
            }
        }

        private void AddRawDataInFile(int imu, RawData row, string server_time, StringBuilder file)
        {
            if (file.Length == 0)
            {
                _ = file.AppendLine("imu,server_time,arduino_time,ax,ay,az,gx,gy,gz,mx,my,mz");
            }
            string newLine = string.Format("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11}",
                                imu, server_time, row.time,
                                row.ax, row.ay, row.az,
                                row.gx, row.gy, row.gz,
                                row.mx, row.my, row.mz);
            _ = file.AppendLine(newLine);
        }



        private void AddCFControlDataInFile(Vector4 row, StringBuilder file)
        {
            if (file.Length == 0)
            {
                file.AppendLine("Timestamp;X;Y;Z;GF");
            }
            var time_stamp = DateTime.Now.ToString("yyyy-MM-ddTHH:mm:ss.ffffff");
            var newLine = string.Format("{0};{1};{2};{3};{4}", time_stamp, row.X, row.Y, row.Z, row.W);
            file.AppendLine(newLine);
        }

        private void SendComands()
        {
            if (command_itr < commands_list.Count)
            {
                var vec = commands_list[command_itr];
                MesSender.SendGXYZ(vec);
                AddCFControlDataInFile(vec, command_file);
                int curr_ms = DateTime.Now.Millisecond;
                prev_ms = curr_ms;
                command_itr++;
                Dispatcher.BeginInvoke(new Action(() =>
                {
                    current_progress.Value = command_itr;
                }));
            }
            else
            {
                lock (locker)
                {
                    CollectionData = false;
                }                    
                time_calib.Stop();
            }
        }

        private void OnTimedEvent_SendCommand(Object source, ElapsedEventArgs e)
        {
            SendComands();
        }

        private void SendDataToCF18(string name)
        {
            lock (locker)
            {
                InitCollectionData();
                data_name = name;
                CollectionData = true;
            }
            
            prev_ms = DateTime.Now.Millisecond;
            command_itr = 0;
            current_progress.Minimum = command_itr;
            current_progress.Maximum = commands_list.Count;

            SendComands();
            time_calib = new System.Timers.Timer(10);
            time_calib.Elapsed += OnTimedEvent_SendCommand;
            time_calib.Start();
        }



        private void btnConnect_Click(object sender, RoutedEventArgs e)
        {
            if (controller == null)
            {
                string com_port = comboBoxCOMPorts.SelectedItem.ToString();
                ThreadPool.QueueUserWorkItem(o =>
                {
                    ConnectArduino(com_port);
                });
                comboBoxCOMPorts.IsEnabled = false;
            }
        }

        private void btnCalibration_Click(object sender, RoutedEventArgs e)
        {
            commands_list.Clear();
            var data_tmp = Common.PrepareDataForCF18(@"./calib.csv");
            foreach (Vector4 vec in data_tmp)
            {
                for (int i = 0; i < 2000; ++i)
                {
                    commands_list.Add(vec);
                }
            }

            current_progress = progressCalibration;
            command_timeout = 10;
                        
            SendDataToCF18("calibration");
        }

        private void btnBoltanka_Click(object sender, RoutedEventArgs e)
        {
            commands_list.Clear();
            commands_list = Common.PrepareDataForCF18(@"./wxyz.csv");

            current_progress = progressBoltanka;
            command_timeout = 55;

            SendDataToCF18("boltanca");
        }
    }
}
