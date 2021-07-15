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
        private string WorkStationAdress = "192.168.62.84";
        private int WorkStationPort = 5918;
        private CF_sender MesSender;
        private ArduinoController controller;
        private MotionProcessor motionProc;

        private byte imu_total = 3;
        private int init_time = 20;
        private int baud_rate = 500000;
        private string data_folder = Directory.GetCurrentDirectory() + "\\raw_cf18_data-" + DateTime.Now.ToString("dd.MM.yyyy.HH.mm.ss");


        private System.Timers.Timer time_calib;
        private bool controling_finished = false;
        private int prev_ms = 0;
        private int command_timeout;
        private int command_itr = 0;
        private List<Vector4> commands_list = new List<Vector4>();
        private ProgressBar current_progress;

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
            var controller = new ArduinoController(imu_total);
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
        }

        private void btnConnect_Click(object sender, RoutedEventArgs e)
        {
            string com_port = comboBoxCOMPorts.SelectedItem.ToString();
            ThreadPool.QueueUserWorkItem(o =>
            {
                ConnectArduino(com_port);
            });
            comboBoxCOMPorts.IsEnabled = false;
        }

        private void SendComands()
        {
            if (command_itr < commands_list.Count)
            {
                var vec = commands_list[command_itr];
                MesSender.SendGXYZ(vec);
                int curr_ms = DateTime.Now.Millisecond;
                prev_ms = curr_ms;
                command_itr++;
                Dispatcher.BeginInvoke(new Action(() =>
                {
                    current_progress.Value = command_itr;
                    //progressCalibration.
                }));
            }
            else
            {
                controling_finished = true;
                time_calib.Stop();
            }
        }

        private void OnTimedEvent_SendCommand(Object source, ElapsedEventArgs e)
        {
            SendComands();
        }

        private void SendDataToCF18()
        {
            prev_ms = DateTime.Now.Millisecond;

            command_itr = 0;

            current_progress.Minimum = command_itr;
            current_progress.Maximum = commands_list.Count;

            SendComands();
            time_calib = new System.Timers.Timer(10);
            time_calib.Elapsed += OnTimedEvent_SendCommand;
            time_calib.Start();
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

            SendDataToCF18();
        }

        private void btnBoltanka_Click(object sender, RoutedEventArgs e)
        {
            commands_list.Clear();
            commands_list = Common.PrepareDataForCF18(@"./wxyz.csv");

            current_progress = progressBoltanka;
            command_timeout = 55;

            SendDataToCF18();
        }
    }
}
